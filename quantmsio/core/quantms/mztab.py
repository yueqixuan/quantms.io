import gzip
import logging
import os
import time
from contextlib import contextmanager
from pathlib import Path
from typing import Any, Dict, Generator, Optional, Union

import duckdb
import pandas as pd

from quantmsio.core.project import create_uuid_filename
from quantmsio.operate.tools import get_modification_details
from quantmsio.utils.pride_utils import get_quantmsio_modifications


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


class MzTab:
    def __init__(
        self,
        mztab_path: Union[Path, str],
        enable_duckdb: bool = False,
        duckdb_max_memory: str = "16GB",
        duckdb_threads: int = 4,
    ) -> None:
        self.mztab_path = Path(mztab_path)
        self.logger = logging.getLogger("quantmsio.core.mztab")
        self.enable_duckdb = enable_duckdb

        # Detect if file is gzipped
        self.is_gzipped = str(self.mztab_path).endswith(".gz")

        # For gzipped files, we'll need to use a different strategy since seeking is not supported
        # We'll decompress to a temporary file when needed for complex operations
        self._temp_decompressed_file = None

        # Position tracking for efficient re-reading
        self._psm_pos = None
        self._psm_len = None
        self._psm_end_pos = None
        self._pep_pos = None
        self._pep_len = None
        self._pep_end_pos = None
        self._prt_pos = None
        self._prt_len = None
        self._prt_end_pos = None

        # Column caching
        self._psms_columns = None
        self._pep_columns = None

        # Metadata caching
        self._metadata_cache = {}
        self._metadata_parsed = False

        # DuckDB integration
        self._duckdb: Optional[duckdb.DuckDBPyConnection] = None
        self._duckdb_name: Optional[str] = None
        self._duckdb_max_memory = duckdb_max_memory
        self._duckdb_threads = duckdb_threads

        # Initialize DuckDB if enabled and file is large enough to benefit
        if self.enable_duckdb and self._should_use_duckdb():
            try:
                self._setup_duckdb()
            except Exception as e:
                self.logger.warning(
                    f"Failed to setup DuckDB, falling back to pandas: {e}"
                )
                self.enable_duckdb = False

    def _should_use_duckdb(self) -> bool:
        """Determine if DuckDB should be used based on file size and complexity."""
        try:
            file_size = self.mztab_path.stat().st_size
            # Use DuckDB for files larger than 100MB or if they have complex structure
            return file_size > 100 * 1024 * 1024  # 100MB threshold
        except Exception:
            return False

    def _setup_duckdb(self) -> None:
        """Initialize DuckDB with optimized settings for mzTab processing."""
        if not self.mztab_path.exists():
            raise FileNotFoundError(f"mzTab file not found: {self.mztab_path}")

        self._duckdb_name = create_uuid_filename("mztab-duckdb", ".db")
        start_time = time.time()

        self._duckdb = duckdb.connect(self._duckdb_name)
        self._duckdb.execute("SET enable_progress_bar=true")
        self._duckdb.execute(f"SET max_memory='{self._duckdb_max_memory}'")
        self._duckdb.execute(f"SET worker_threads='{self._duckdb_threads}'")

        # Load mzTab sections into separate tables if they exist
        self._load_mztab_sections_to_duckdb()

        elapsed = time.time() - start_time
        self.logger.info(f"DuckDB setup completed in {elapsed:.2f} seconds")

    def _load_mztab_sections_to_duckdb(self) -> None:
        """Load mzTab data sections (PSM, PEP, PRT) into DuckDB tables."""
        sections = ["PSM", "PEP", "PRT"]
        section_headers = {"PSM": "PSH", "PEP": "PEH", "PRT": "PRH"}

        for section in sections:
            try:
                header = section_headers[section]
                if self._has_section(header):
                    # Create temporary CSV file for the section
                    temp_csv = self._extract_section_to_temp_csv(section, header)
                    if temp_csv and temp_csv.exists():
                        # Load into DuckDB
                        table_name = section.lower()
                        self._duckdb.execute(
                            f"""
                            CREATE TABLE IF NOT EXISTS {table_name} AS 
                            SELECT * FROM read_csv_auto('{temp_csv}', header=true, sep='\t')
                        """
                        )
                        # Clean up temporary file
                        temp_csv.unlink()
                        self.logger.debug(
                            f"Loaded {section} section into DuckDB table '{table_name}'"
                        )
            except Exception as e:
                self.logger.warning(
                    f"Failed to load {section} section into DuckDB: {e}"
                )

    def _has_section(self, header: str) -> bool:
        """Check if a section exists in the mzTab file."""
        try:
            with self._safe_file_open() as f:
                for line in f:
                    if line.startswith(header):
                        return True
                    if not line.startswith("MTD") and not line.startswith(header[:2]):
                        break
            return False
        except Exception:
            return False

    def _extract_section_to_temp_csv(self, section: str, header: str) -> Optional[Path]:
        """Extract a section to a temporary CSV file for DuckDB loading."""
        try:
            temp_file = Path(create_uuid_filename(f"mztab_{section.lower()}", ".csv"))

            # Use seekable file path for consistent handling
            seekable_path = self._get_seekable_file_path()
            with open(seekable_path, "r", encoding="utf-8") as infile, open(
                temp_file, "w", encoding="utf-8"
            ) as outfile:

                # Find the header line
                header_found = False
                for line in infile:
                    if line.startswith(header):
                        # Write the header (remove the PSH/PEH/PRH prefix)
                        outfile.write("\t".join(line.strip().split("\t")[1:]) + "\n")
                        header_found = True
                        break

                if not header_found:
                    temp_file.unlink()
                    return None

                # Write data lines
                for line in infile:
                    if line.startswith(section):
                        # Write the data (remove the PSM/PEP/PRT prefix)
                        outfile.write("\t".join(line.strip().split("\t")[1:]) + "\n")
                    elif not line.startswith(section):
                        # End of section
                        break

            return temp_file
        except Exception as e:
            self.logger.warning(f"Failed to extract {section} section: {e}")
            return None

    @contextmanager
    def _safe_file_open(self, mode="r", encoding="utf-8"):
        """Context manager for safe file operations with gzip support."""
        try:
            self._validate_file()
            if self.is_gzipped:
                # For gzipped files, use gzip.open
                if "b" not in mode:
                    mode = mode + "t"  # Ensure text mode for gzipped files
                f = gzip.open(self.mztab_path, mode, encoding=encoding)
            else:
                # For regular files, use standard open
                f = open(self.mztab_path, mode, encoding=encoding)
            try:
                yield f
            finally:
                f.close()
        except Exception as e:
            self.logger.error(f"File operation failed: {e}")
            raise

    def _validate_file(self) -> None:
        """Validate that the mzTab file exists and is not empty."""
        if not self.mztab_path.exists():
            raise FileNotFoundError(f"mzTab file not found: {self.mztab_path}")
        if self.mztab_path.stat().st_size == 0:
            raise ValueError("mzTab file is empty")

    def _get_seekable_file_path(self) -> Path:
        """Get a seekable file path. For gzipped files, create temporary decompressed file."""
        if not self.is_gzipped:
            return self.mztab_path

        # For gzipped files, create a temporary decompressed file
        if (
            self._temp_decompressed_file is None
            or not self._temp_decompressed_file.exists()
        ):
            self._temp_decompressed_file = Path(
                create_uuid_filename("mztab_temp", ".mztab")
            )
            self.logger.info(
                f"Decompressing gzipped file to temporary file: {self._temp_decompressed_file}"
            )

            try:
                with gzip.open(
                    self.mztab_path, "rt", encoding="utf-8"
                ) as gz_file, open(
                    self._temp_decompressed_file, "w", encoding="utf-8"
                ) as temp_file:
                    # Copy content from gzipped file to temporary file
                    for line in gz_file:
                        temp_file.write(line)
                self.logger.info("Decompression completed successfully")
            except Exception as e:
                self.logger.error(f"Failed to decompress gzipped file: {e}")
                if (
                    self._temp_decompressed_file
                    and self._temp_decompressed_file.exists()
                ):
                    self._temp_decompressed_file.unlink()
                    self._temp_decompressed_file = None
                raise

        return self._temp_decompressed_file

    def _parse_metadata_once(self) -> None:
        """Parse metadata section once and cache results."""
        if self._metadata_parsed:
            return

        self._metadata_cache = {
            "ms_runs": {},
            "score_names": {},
            "modifications": {},
            "mods_map": {},
        }

        try:
            with self._safe_file_open() as f:
                for line in f:
                    line = line.strip()
                    if not line.startswith("MTD"):
                        break

                    parts = line.split("\t")
                    if len(parts) < 3:
                        continue

                    # Parse MS runs
                    if parts[1].endswith("-location"):
                        ms_run_id = parts[1].split("-")[0]
                        filename = parts[2].split("//")[-1].split(".")[0]
                        self._metadata_cache["ms_runs"][ms_run_id] = filename

                    # Parse score names
                    elif "psm_search_engine_score" in parts[1]:
                        score_values = (
                            parts[2].replace("[", "").replace("]", "").split(",")
                        )
                        if len(score_values) >= 3:
                            score_name = score_values[2].strip()
                            if ":" in score_name:
                                score_name = score_name.split(":")[0]
                            self._metadata_cache["score_names"][score_name] = parts[
                                1
                            ].replace("psm_", "")

                            # Parse modifications
                    elif "_mod[" in parts[1]:
                        self._metadata_cache["modifications"] = (
                            fetch_modifications_from_mztab_line(
                                line, self._metadata_cache["modifications"]
                            )
                        )

        except Exception as e:
            self.logger.warning(f"Failed to parse metadata: {e}")

        self._metadata_parsed = True

    def _get_pos(self, header):
        if header == "PSH" and self._pep_pos is not None:
            return self._pep_end_pos
        elif header == "PEH" and self._prt_pos is not None:
            return self._prt_end_pos
        else:
            return 0

    def __extract_len(self, header):
        map_tag = {"PSH": "PSM", "PEH": "PEP", "PRH": "PRT"}

        # For operations requiring seeking, use seekable file path
        seekable_path = self._get_seekable_file_path()

        with open(seekable_path, "r", encoding="utf-8") as f:
            pos = self._get_pos(header)
            f.seek(pos)
            line = f.readline()
            while line and not line.startswith(header):
                pos = f.tell()
                line = f.readline()

            if not line:
                return 0, pos, pos

            if header == "PSH":
                self._psms_columns = line.split("\n")[0].split("\t")
            if header == "PEH":
                self._pep_columns = line.split("\n")[0].split("\t")

            line = f.readline()
            fle_len = 0
            while line and line.startswith(map_tag[header]):
                fle_len += 1
                line = f.readline()
            end_pos = f.tell()

        return fle_len, pos, end_pos

    def __load_second(self, header, **kwargs):
        # For operations requiring seeking, use seekable file path
        seekable_path = self._get_seekable_file_path()

        # Check if this is chunked reading
        if "chunksize" in kwargs and kwargs["chunksize"] is not None:
            # For chunked reading, we need to keep the file handle open
            f = open(seekable_path, "r", encoding="utf-8")
            if header == "PSH":
                f.seek(self._psm_pos)
                return pd.read_csv(
                    f, sep="\t", nrows=self._psm_len, low_memory=False, **kwargs
                )
            elif header == "PEH":
                f.seek(self._pep_pos)
                return pd.read_csv(
                    f, sep="\t", nrows=self._pep_len, low_memory=False, **kwargs
                )
            else:
                f.seek(self._prt_pos)
                return pd.read_csv(
                    f, sep="\t", nrows=self._prt_len, low_memory=False, **kwargs
                )
        else:
            # Non-chunked reading - use context manager
            with open(seekable_path, "r", encoding="utf-8") as f:
                if header == "PSH":
                    f.seek(self._psm_pos)
                    return pd.read_csv(
                        f, sep="\t", nrows=self._psm_len, low_memory=False, **kwargs
                    )
                elif header == "PEH":
                    f.seek(self._pep_pos)
                    return pd.read_csv(
                        f, sep="\t", nrows=self._pep_len, low_memory=False, **kwargs
                    )
                else:
                    f.seek(self._prt_pos)
                    return pd.read_csv(
                        f, sep="\t", nrows=self._prt_len, low_memory=False, **kwargs
                    )

    def __set_table_config(self, header, length, pos, end_pos):
        if header == "PSH":
            self._psm_pos = pos
            self._psm_len = length
            self._psm_end_pos = end_pos
        elif header == "PEH":
            self._pep_pos = pos
            self._pep_len = length
            self._pep_end_pos = end_pos
        else:
            self._prt_pos = pos
            self._prt_len = length
            self._prt_end_pos = end_pos

    def skip_and_load_csv(self, header, **kwargs):
        """Load CSV data with improved caching and error handling."""
        if self._psm_pos is not None and header == "PSH":
            return self.__load_second(header, **kwargs)
        if self._pep_pos is not None and header == "PEH":
            return self.__load_second(header, **kwargs)
        if self._prt_pos is not None and header == "PRH":
            return self.__load_second(header, **kwargs)

        # Check if this is a chunked read (iterator)
        if "chunksize" in kwargs and kwargs["chunksize"] is not None:
            # For chunked reading, we need to return a pandas iterator
            # Cannot use DuckDB for this, fall back to pandas
            fle_len, pos, end_pos = self.__extract_len(header)
            seekable_path = self._get_seekable_file_path()
            f = open(seekable_path, "r", encoding="utf-8")
            f.seek(pos)
            self.__set_table_config(header, fle_len, pos, end_pos)
            return pd.read_csv(f, nrows=fle_len, sep="\t", low_memory=False, **kwargs)

        # Use DuckDB if available and appropriate for non-chunked reads
        if self.enable_duckdb and self._duckdb:
            table_map = {"PSH": "psm", "PEH": "pep", "PRH": "prt"}
            table_name = table_map.get(header)
            if table_name:
                try:
                    # Check if table exists
                    result = self._duckdb.execute(
                        f"SELECT name FROM pragma_table_info('{table_name}')"
                    ).fetchall()
                    if result:
                        # Table exists, use DuckDB query
                        columns = kwargs.get("usecols", None)
                        if columns:
                            cols_str = ", ".join(f'"{col}"' for col in columns)
                            df = self._duckdb.execute(
                                f"SELECT {cols_str} FROM {table_name}"
                            ).df()
                        else:
                            df = self._duckdb.execute(
                                f"SELECT * FROM {table_name}"
                            ).df()
                        return df
                except Exception as e:
                    self.logger.warning(
                        f"DuckDB query failed, falling back to pandas: {e}"
                    )

        # Fallback to original method
        fle_len, pos, end_pos = self.__extract_len(header)
        seekable_path = self._get_seekable_file_path()
        with open(seekable_path, "r", encoding="utf-8") as f:
            f.seek(pos)
            self.__set_table_config(header, fle_len, pos, end_pos)
            return pd.read_csv(f, nrows=fle_len, sep="\t", low_memory=False, **kwargs)

    def extract_ms_runs(self):
        """Extract MS runs with caching."""
        self._parse_metadata_once()
        return self._metadata_cache["ms_runs"].copy()

    def get_protein_map(self, protein_str=None):
        """
        return: a dict about protein score with improved performance
        """
        try:
            if self.enable_duckdb and self._duckdb:
                # Try DuckDB approach first
                try:
                    if protein_str:
                        query = f"""
                        SELECT ambiguity_members, MIN("best_search_engine_score[1]") as min_score
                        FROM prt 
                        WHERE ambiguity_members LIKE '%{protein_str}%'
                        GROUP BY ambiguity_members
                        """
                    else:
                        query = """
                        SELECT ambiguity_members, MIN("best_search_engine_score[1]") as min_score
                        FROM prt 
                        GROUP BY ambiguity_members
                        """

                    result_df = self._duckdb.execute(query).df()
                    return dict(
                        zip(result_df["ambiguity_members"], result_df["min_score"])
                    )
                except Exception as e:
                    self.logger.warning(f"DuckDB protein map query failed: {e}")

            # Fallback to original pandas approach
            prt = self.skip_and_load_csv(
                "PRH",
                usecols=["ambiguity_members", "best_search_engine_score[1]"],
            )
            if protein_str:
                prt = prt[
                    prt["ambiguity_members"].str.contains(f"{protein_str}", na=False)
                ]
            prt_score = prt.groupby("ambiguity_members").min()
            protein_map = prt_score.to_dict()["best_search_engine_score[1]"]
            return protein_map

        except Exception as e:
            self.logger.error(f"Failed to get protein map: {e}")
            return {}

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

        # If mods_map is empty, parse it separately
        if not self._metadata_cache["mods_map"]:
            try:
                with self._safe_file_open() as f:
                    mods_map = {}
                    lines = []

                    # Read all MTD lines first
                    for line in f:
                        if line.startswith("MTD"):
                            lines.append(line.strip())
                        else:
                            break

                    # Parse modifications in two passes
                    i = 0
                    while i < len(lines):
                        line = lines[i]
                        if "_mod[" in line:
                            line_parts = line.split("\t")
                            if (
                                len(line_parts) >= 3
                                and "site" not in line_parts[1]
                                and "position" not in line_parts[1]
                            ):
                                values = (
                                    line_parts[2]
                                    .replace("[", "")
                                    .replace("]", "")
                                    .split(",")
                                )
                                if len(values) >= 3:
                                    accession = values[1].strip()
                                    name = values[2].strip()

                                    # Look for corresponding site line
                                    for j in range(i + 1, min(i + 3, len(lines))):
                                        if j < len(lines) and "site" in lines[j]:
                                            site_parts = lines[j].split("\t")
                                            if len(site_parts) >= 3:
                                                site = site_parts[2]
                                                mods_map[name] = [
                                                    accession.upper(),
                                                    site,
                                                ]
                                                mods_map[accession.upper()] = [
                                                    name,
                                                    site,
                                                ]
                                                break
                        i += 1

                    self._metadata_cache["mods_map"] = mods_map
            except Exception as e:
                self.logger.warning(f"Failed to parse mods_map: {e}")

        return self._metadata_cache["mods_map"].copy()

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

    # New methods for enhanced functionality

    def get_file_statistics(self) -> Dict[str, Any]:
        """Get comprehensive file statistics."""
        stats = {
            "file_size_mb": round(self.mztab_path.stat().st_size / (1024 * 1024), 2),
            "sections": {},
        }

        try:
            if self.enable_duckdb and self._duckdb:
                # Get statistics using DuckDB
                for table in ["psm", "pep", "prt"]:
                    try:
                        count_result = self._duckdb.execute(
                            f"SELECT COUNT(*) as count FROM {table}"
                        ).fetchone()
                        if count_result:
                            stats["sections"][table.upper()] = {
                                "row_count": count_result[0]
                            }
                    except Exception:
                        stats["sections"][table.upper()] = {"row_count": 0}
            else:
                # Fallback to position-based counting
                for header in ["PSH", "PEH", "PRH"]:
                    try:
                        length, _, _ = self.__extract_len(header)
                        section_name = {"PSH": "PSM", "PEH": "PEP", "PRH": "PRT"}[
                            header
                        ]
                        stats["sections"][section_name] = {"row_count": length}
                    except (ValueError, KeyError, IndexError) as e:
                        # Log specific errors for debugging but continue processing other sections
                        self.logger.debug(
                            f"Failed to extract length for section {header}: {e}"
                        )
                        continue
                    except Exception as e:
                        # Log unexpected errors but continue processing other sections
                        self.logger.warning(
                            f"Unexpected error extracting length for section {header}: {e}"
                        )
                        continue

            # Add metadata statistics
            self._parse_metadata_once()
            stats["metadata"] = {
                "ms_runs_count": len(self._metadata_cache["ms_runs"]),
                "modifications_count": len(self._metadata_cache["modifications"]),
                "score_types_count": len(self._metadata_cache["score_names"]),
            }

        except Exception as e:
            self.logger.warning(f"Failed to get file statistics: {e}")

        return stats

    def query_proteins_efficiently(
        self, protein_list: list, columns: list = None
    ) -> pd.DataFrame:
        """Efficiently query multiple proteins at once."""
        if not protein_list:
            return pd.DataFrame()

        try:
            if self.enable_duckdb and self._duckdb:
                # Use DuckDB for efficient querying
                protein_conditions = []
                for protein in protein_list:
                    protein_conditions.append(f"ambiguity_members LIKE '%{protein}%'")

                where_clause = " OR ".join(protein_conditions)
                cols = ", ".join(f'"{col}"' for col in columns) if columns else "*"

                query = f"SELECT {cols} FROM prt WHERE {where_clause}"
                return self._duckdb.execute(query).df()
            else:
                # Fallback to pandas
                prt = self.skip_and_load_csv("PRH", usecols=columns)
                mask = prt["ambiguity_members"].str.contains(
                    "|".join(protein_list), na=False
                )
                return prt[mask]

        except Exception as e:
            self.logger.error(f"Failed to query proteins efficiently: {e}")
            return pd.DataFrame()

    def stream_section(
        self, section: str, chunk_size: int = 10000
    ) -> Generator[pd.DataFrame, None, None]:
        """Stream a section in chunks for memory-efficient processing."""
        section_map = {"PSM": "PSH", "PEP": "PEH", "PRT": "PRH"}
        header = section_map.get(section.upper())

        if not header:
            raise ValueError(f"Invalid section: {section}")

        try:
            if self.enable_duckdb and self._duckdb:
                # Use DuckDB for streaming
                table_name = section.lower()
                offset = 0

                while True:
                    query = (
                        f"SELECT * FROM {table_name} LIMIT {chunk_size} OFFSET {offset}"
                    )
                    chunk_df = self._duckdb.execute(query).df()

                    if chunk_df.empty:
                        break

                    yield chunk_df
                    offset += chunk_size
            else:
                # Fallback to pandas chunking
                full_df = self.skip_and_load_csv(header)
                for i in range(0, len(full_df), chunk_size):
                    yield full_df.iloc[i : i + chunk_size]

        except Exception as e:
            self.logger.error(f"Failed to stream section {section}: {e}")

    def __enter__(self):
        """Context manager entry."""
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Context manager exit with cleanup."""
        self.cleanup()

    def cleanup(self):
        """Clean up resources."""
        if self._duckdb:
            try:
                self._duckdb.close()
            except Exception as e:
                self.logger.warning(f"Failed to close DuckDB connection: {e}")
            self._duckdb = None

        if self._duckdb_name and Path(self._duckdb_name).exists():
            try:
                os.remove(self._duckdb_name)
            except Exception as e:
                self.logger.warning(
                    f"Failed to remove DuckDB file {self._duckdb_name}: {e}"
                )
            self._duckdb_name = None

        # Clean up temporary decompressed file
        if self._temp_decompressed_file and self._temp_decompressed_file.exists():
            try:
                self._temp_decompressed_file.unlink()
            except Exception as e:
                self.logger.warning(
                    f"Failed to remove temporary file {self._temp_decompressed_file}: {e}"
                )
            self._temp_decompressed_file = None

    def __del__(self):
        """Destructor to ensure cleanup."""
        self.cleanup()
