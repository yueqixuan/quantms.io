import codecs
import os
from pathlib import Path
from typing import Union, Optional
import logging

import pandas as pd
from quantmsio.utils.pride_utils import get_quantmsio_modifications
from quantmsio.operate.tools import get_modification_details


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
    def __init__(self, mztab_path: Union[Path, str]) -> None:
        self.mztab_path = mztab_path
        self.logger = logging.getLogger("quantmsio.core.mztab")
        # psm pos
        self._psm_pos = None
        # psm len
        self._psm_len = None
        self._psm_end_pos = None
        # pep pos
        self._pep_pos = None
        # pep len
        self._pep_len = None
        self._pep_end_pos = None
        # prt pos
        self._prt_pos = None
        # prt len
        self._prt_len = None
        self._prt_end_pos = None
        # load psms columns
        self._psms_columns = None
        # load pep columns
        self._pep_columns = None

    def _get_pos(self, pattern: str) -> int:
        """Get position of pattern in file."""
        self.logger.debug(f"🔍 Searching for pattern: {pattern}")
        if os.stat(self.mztab_path).st_size == 0:
            raise ValueError("File is empty")
        f = open(self.mztab_path)
        pos = 0
        line = f.readline()
        while not line.startswith(pattern):
            pos = f.tell()
            line = f.readline()
            if not line:
                break
        f.close()
        return pos

    def __extract_len(self, header):
        map_tag = {"PSH": "PSM", "PEH": "PEP", "PRH": "PRT"}
        if os.stat(self.mztab_path).st_size == 0:
            raise ValueError("File is empty")
        f = open(self.mztab_path)
        pos = self._get_pos(header)
        f.seek(pos)
        line = f.readline()
        while not line.startswith(header):
            pos = f.tell()
            line = f.readline()

        if header == "PSH":
            self._psms_columns = line.split("\n")[0].split("\t")
        if header == "PEH":
            self._pep_columns = line.split("\n")[0].split("\t")

        line = f.readline()
        fle_len = 0
        while line.startswith(map_tag[header]):
            fle_len += 1
            line = f.readline()
        end_pos = f.tell()
        f.close()
        return fle_len, pos, end_pos

    def __load_second(self, header, **kwargs):
        f = open(self.mztab_path)
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

    def skip_and_load_csv(
        self, pattern: str, chunksize: Optional[int] = None, usecols=None
    ):
        """Skip to pattern and load as CSV."""
        import pandas as pd

        self.logger.debug(f"📖 Loading CSV data after pattern: {pattern}")
        pos = self._get_pos(pattern)
        chunks_processed = 0
        total_rows = 0

        for chunk in pd.read_csv(
            self.mztab_path,
            sep="\t",
            skiprows=lambda x: x < pos,
            chunksize=chunksize if chunksize else None,
            usecols=usecols,
        ):
            chunks_processed += 1
            total_rows += len(chunk)
            if chunks_processed % 5 == 0:  # Log every 5 chunks
                self.logger.debug(
                    f"⏳ Loaded {chunks_processed} chunks, {total_rows:,} rows so far..."
                )
            yield chunk

    def extract_ms_runs(self):
        """Extract MS run information."""
        self.logger.debug("🔍 Extracting MS run information...")
        ms_runs = {}
        if os.stat(self.mztab_path).st_size == 0:
            raise ValueError("File is empty")
        f = open(self.mztab_path)
        for line in f:
            if line.startswith("MTD"):
                if "ms_run[" in line:
                    if "location" in line:
                        key = line.split("\t")[1].split("-")[0].strip()
                        value = line.split("\t")[2].split("/")[-1].strip()
                        ms_runs[key] = value
            elif line.startswith("PRH"):
                break
        f.close()
        self.logger.debug(f"✓ Found {len(ms_runs)} MS runs")
        return ms_runs

    def get_protein_map(self, protein_str=None):
        """
        return: a dict about protein score
        """
        self.logger.debug("🔍 Extracting protein global q-value map...")
        protein_map = {}
        rows_processed = 0

        for df in self.skip_and_load_csv("PRT", chunksize=100000):
            rows_processed += len(df)
            if "accession" in df.columns and "global_qvalue" in df.columns:
                temp_dict = df.set_index("accession")["global_qvalue"].to_dict()
                protein_map.update(temp_dict)
            if rows_processed % 100000 == 0:
                self.logger.debug(f"⏳ Processed {rows_processed:,} protein rows...")

        self.logger.debug(f"✓ Extracted {len(protein_map)} protein q-values")
        if protein_str:
            protein_map = {k: v for k, v in protein_map.items() if protein_str in k}
        return protein_map

    def get_score_names(self):
        """Extract score names."""
        self.logger.debug("🔍 Extracting score names...")
        score_names = {}
        if os.stat(self.mztab_path).st_size == 0:
            raise ValueError("File is empty")
        f = open(self.mztab_path)
        for line in f:
            if line.startswith("MTD"):
                if "psm_search_engine_score" in line:
                    key = line.split("[")[1].split("]")[0]
                    value = f"search_engine_score[{key}]"
                    score_names[line.split("\t")[2].strip()] = value
            elif line.startswith("PSH"):
                break
        f.close()
        self.logger.debug(f"✓ Found {len(score_names)} score names")
        return score_names

    @staticmethod
    def generate_positions(start, end) -> list:
        start = start.split(",")
        end = end.split(",")
        return [start + ":" + end for start, end in zip(start, end)]

    def get_modifications(self):
        """Extract modifications."""
        self.logger.debug("🔍 Extracting modifications...")
        modifications = {}
        if os.stat(self.mztab_path).st_size == 0:
            raise ValueError("File is empty")
        f = open(self.mztab_path)
        for line in f:
            if line.startswith("MTD"):
                if "fixed_mod" in line or "variable_mod" in line:
                    key = line.split("[")[1].split("]")[0]
                    value = line.split("\t")[2].strip()
                    modifications[key] = value
            elif line.startswith("PRT"):
                break
        f.close()
        self.logger.debug(f"✓ Found {len(modifications)} modifications")
        return modifications

    def get_mods_map(self):
        """Extract modifications map."""
        self.logger.debug("🔍 Extracting modifications map...")
        mods_map = {}
        modifications = self.get_modifications()
        for key, value in modifications.items():
            if "[" in value and "]" in value:
                name = value.split("[")[0].strip()
                mass = value.split("[")[1].split("]")[0].strip()
                mods_map[name] = float(mass)
        self.logger.debug(f"✓ Created map for {len(mods_map)} modifications")
        return mods_map

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
