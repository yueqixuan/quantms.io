import gzip
import logging
import os
import re
import time
from contextlib import contextmanager
from pathlib import Path
from typing import Any, Dict, Generator, Optional, Union, List

import duckdb
import pandas as pd

from quantmsio.core.duckdb import MzTabIndexer
from quantmsio.core.project import create_uuid_filename
from quantmsio.operate.tools import get_modification_details
from quantmsio.utils.pride_utils import (
    get_quantmsio_modifications,
    fetch_modifications_from_mztab_line,
)


def generate_modification_list(modification_str: str, modifications):
    if pd.isna(modification_str):
        return None
    modifications = get_quantmsio_modifications(modification_str, modifications)
    modifications_string = ""
    for _, value in modifications.items():
        modifications_string += "|".join(map(str, value["position"]))
        modifications_string = (
            modifications_string + "-" + value["unimod_accession"] + ","
        )
    modifications_string = modifications_string[:-1]  # Remove last comma
    modification_list = modifications_string.split(",")
    return modification_list


def fetch_modifications_from_mztab_line(line: str, _modifications: dict) -> dict:
    """
    get the modifications from a mztab line. An mzTab modification could be a fixed or variable modification.
    The structure of a fixed is the following:
      MTD	fixed_mod[1]	[UNIMOD, UNIMOD:4, Carbamidomethyl, ]
      MTD	fixed_mod[1]-site	C
      MTD	fixed_mod[1]-position	Anywhere
    while the structure of a variable modification is the following:
      MTD	var_mod[1]	[UNIMOD, UNIMOD:21, Phospho, ]
      MTD	var_mod[1]-site	S
      MTD   var_mod[1]-position	Anywhere

    :param line: mztab line
    :param _modifications: modifications dictionary
    :return: modification dictionary
    """
    line = line.strip()
    line_parts = line.split("\t")
    if line_parts[0] == "MTD" and "_mod[" in line_parts[1]:
        if "site" not in line_parts[1] and "position" not in line_parts[1]:
            values = line_parts[2].replace("[", "").replace("]", "").split(",")
            accession = values[1].strip()
            name = values[2].strip()
            index = line_parts[1].split("[")[1].split("]")[0]
            _modifications[accession] = [name, index, None, None]
        elif "site" in line_parts[1]:
            index = line_parts[1].split("[")[1].split("]")[0]
            accession = None
            for (
                key,
                value,
            ) in (
                _modifications.items()
            ):  # for name, age in dictionary.iteritems():  (for Python 2.x)
                if value[1] == index:
                    accession = key
            if accession is None:
                raise Exception("The accession for the modification is None")
            _modifications[accession][2] = line_parts[2]
        elif "position" in line_parts[1]:
            index = line_parts[1].split("[")[1].split("]")[0]
            accession = None
            for key, value in _modifications.items():
                if value[1] == index:
                    accession = key
            if accession is None:
                raise Exception("The accession for the modification is None")
            _modifications[accession][3] = line_parts[2]
    return _modifications


class MzTab(MzTabIndexer):
    def __init__(
        self,
        mztab_path: Union[Path, str],
        msstats_path: Union[Path, str] = None,
        duckdb_max_memory: str = "16GB",
        duckdb_threads: int = 4,
    ) -> None:
        super().__init__(
            mztab_path=mztab_path,
            msstats_path=msstats_path,
            max_memory=duckdb_max_memory,
            worker_threads=duckdb_threads,
        )
        self.logger = logging.getLogger(f"{__name__}.{self.__class__.__name__}")

        # Column caching
        self._psms_columns = None
        self._pep_columns = None

        # Metadata caching
        self._metadata_cache = {}
        self._metadata_parsed = False

    def _parse_metadata_once(self) -> None:
        """Parse metadata section once and cache results using DuckDB."""
        if self._metadata_parsed:
            return

        self._metadata_cache = {
            "ms_runs": {},
            "score_names": {},
            "modifications": {},
            "mods_map": {},
        }

        try:
            metadata_df = self.get_metadata()
            if metadata_df.empty:
                self.logger.warning("Metadata table is empty in DuckDB.")
                self._metadata_parsed = True
                return

            # Using vectorized operations for filtering with escaped patterns
            is_location = metadata_df["key"].str.endswith("-location", na=False)
            is_score = metadata_df["key"].str.contains(
                r"psm_search_engine_score", regex=True, na=False
            ) & ~metadata_df["key"].str.contains(r"\-", regex=True, na=False)
            is_mod = metadata_df["key"].str.contains(r"_mod\[", regex=True, na=False)

            # Parse MS runs
            for _, row in metadata_df[is_location].iterrows():
                ms_run_id = row["key"].split("-")[0]
                if pd.notna(row["value"]) and isinstance(row["value"], str):
                    filename = row["value"].split("//")[-1].split(".")[0]
                    self._metadata_cache["ms_runs"][ms_run_id] = filename

            # Parse score names
            for _, row in metadata_df[is_score].iterrows():
                if pd.notna(row["value"]) and isinstance(row["value"], str):
                    score_values = (
                        row["value"].replace("[", "").replace("]", "").split(",")
                    )
                    if len(score_values) >= 3:
                        score_name = score_values[2].strip().split(":")[0]
                        self._metadata_cache["score_names"][score_name] = row[
                            "key"
                        ].replace("psm_", "")

            # Parse modifications
            _modifications = self._metadata_cache.get("modifications", {})
            for _, row in metadata_df[is_mod].iterrows():
                if pd.notna(row["value"]) and isinstance(row["value"], str):
                    line = f"MTD\t{row['key']}\t{row['value']}"
                    _modifications = fetch_modifications_from_mztab_line(
                        line, _modifications
                    )
            self._metadata_cache["modifications"] = _modifications

        except (duckdb.Error, AttributeError, Exception) as e:
            self.logger.error(f"Error parsing metadata from DuckDB: {e}")

        self._metadata_parsed = True

    def extract_ms_runs(self):
        """Extract MS runs with caching."""
        self._parse_metadata_once()
        return self._metadata_cache["ms_runs"].copy()

    def get_score_names(self):
        """Extract score names with caching."""
        self._parse_metadata_once()
        return self._metadata_cache["score_names"].copy()

    @staticmethod
    def generate_positions(start, end) -> list:
        start = start.split(",")
        end = end.split(",")
        return [start + ":" + end for start, end in zip(start, end)]

    def get_modifications(self):
        """Get modifications with caching."""
        self._parse_metadata_once()
        return self._metadata_cache["modifications"].copy()

    def get_mods_map(self):
        """Get modifications map with caching."""
        self._parse_metadata_once()

        # If mods_map is empty, parse it from cached modifications
        if not self._metadata_cache.get("mods_map") and self._metadata_cache.get(
            "modifications"
        ):
            mods_map = {}
            for accession, value in self._metadata_cache["modifications"].items():
                name, _, site, _ = value
                if name and accession and site:
                    mods_map[name] = [accession.upper(), site]
                    mods_map[accession.upper()] = [name, site]
            self._metadata_cache["mods_map"] = mods_map

        return self._metadata_cache.get("mods_map", {}).copy()

    @staticmethod
    def generate_modifications_details(seq, mods_map, automaton, select_mods):
        seq = seq.replace(".", "")
        peptidoform, modification_details = get_modification_details(
            seq, mods_map, automaton, select_mods
        )
        if len(modification_details) == 0:
            return [peptidoform, None]
        else:
            return [peptidoform, modification_details]

    def stream_section(
        self, section: str, chunk_size: int = 10000
    ) -> Generator[pd.DataFrame, None, None]:
        """Stream a section in chunks for memory-efficient processing."""
        section_map = {
            "PSM": MzTabIndexer._MZTAB_INDEXER_TABLE_PSMS,
            "PRT": MzTabIndexer._MZTAB_INDEXER_TABLE_PROTEINS,
        }
        table_name = section_map.get(section.upper())

        if not table_name:
            raise ValueError(f"Invalid section: {section}")

        try:
            # Ensure we have a valid connection
            if not self._duckdb:
                self.logger.error("No database connection available")
                return

            # Get total count for logging
            total_count = self._duckdb.execute(
                f"SELECT COUNT(*) FROM {table_name}"
            ).fetchone()[0]
            self.logger.debug(
                f"Starting to stream {section} section with {total_count:,} total rows"
            )

            # Use DuckDB for streaming
            offset = 0
            rows_processed = 0

            while True:
                self.logger.debug(f"Fetching chunk at offset {offset:,}")
                query = f"SELECT * FROM {table_name} LIMIT {chunk_size} OFFSET {offset}"
                chunk_df = self._duckdb.execute(query).df()

                if chunk_df.empty:
                    self.logger.debug("Reached end of data")
                    break

                rows_processed += len(chunk_df)
                self.logger.debug(
                    f"Processed {rows_processed:,}/{total_count:,} rows ({(rows_processed/total_count)*100:.1f}%)"
                )

                yield chunk_df
                offset += chunk_size

            self.logger.debug(
                f"Finished streaming {section} section. Total rows processed: {rows_processed:,}"
            )

        except Exception as e:
            self.logger.error(f"Failed to stream section {section}: {e}", exc_info=True)
            raise  # Re-raise the exception after logging it
