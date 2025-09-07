"""
MaxQuant data processing module for proteomics data conversion.
Combines core functionality with advanced features in a single, clean interface.
"""

import logging
import re
import zipfile
from pathlib import Path
from typing import List, Union

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
from pyopenms import AASequence

from quantmsio.core.format import PG_SCHEMA, FEATURE_SCHEMA, PSM_SCHEMA
from quantmsio.core.common import (
    MAXQUANT_FEATURE_MAP,
    MAXQUANT_FEATURE_USECOLS,
    MAXQUANT_PSM_MAP,
    MAXQUANT_PSM_USECOLS,
)
from quantmsio.core.quantms.feature import Feature
from quantmsio.core.sdrf import SDRFHandler
from quantmsio.operate.tools import (
    get_ahocorasick,
    get_protein_accession,
)
from quantmsio.utils.file_utils import (
    extract_protein_list,
    ParquetBatchWriter,
)

logging.basicConfig(format="%(asctime)s - %(message)s", level=logging.INFO)
intensity_normalize_pattern = r"Reporter intensity corrected \d+"
intensity_pattern = r"Reporter intensity \d+"


def clean_peptidoform(peptidoform):
    """Clean MaxQuant peptidoform string for pyopenms compatibility."""
    if not isinstance(peptidoform, str):
        return ""

    peptidoform = peptidoform.strip("_")

    modification_mapping = {
        "ac": "Acetyl",
        "ox": "Oxidation",
        "me": "Methyl",
        "ph": "Phospho",
        "de": "Deamidated",
        "cam": "Carbamidomethyl",
        "dim": "Dimethyl",
        "tri": "Trimethyl",
        "ub": "GlyGly",
        "su": "Sumo",
    }

    for short_name, full_name in modification_mapping.items():
        pattern = f"\\({re.escape(short_name)}\\)"
        replacement = f"({full_name})"
        peptidoform = re.sub(pattern, replacement, peptidoform)

    return peptidoform


def calculate_theoretical_mz(row):
    """Calculate theoretical m/z using pyopenms."""
    try:
        cleaned_peptidoform = clean_peptidoform(row["peptidoform"])
        if not cleaned_peptidoform:
            return None
        sequence = AASequence.fromString(cleaned_peptidoform)
        charge = row["precursor_charge"]
        if charge > 0:
            return sequence.getMZ(charge)
        return sequence.getMonoWeight()
    except Exception:
        return None


def convert_maxquant_flag(value):
    """Convert MaxQuant '+' flag to binary integer."""
    return 1 if value == "+" else 0


def parse_modifications_from_peptidoform(peptidoform: str) -> list:
    """Parse modifications from peptidoform string using pyopenms."""
    if not isinstance(peptidoform, str):
        return None

    try:
        cleaned_peptidoform = clean_peptidoform(peptidoform)
        if not cleaned_peptidoform:
            return None

        sequence = AASequence.fromString(cleaned_peptidoform)
        modifications = {}

        if sequence.hasNTerminalModification():
            mod = sequence.getNTerminalModification()
            accession = mod.getUniModAccession()
            short_name = mod.getId()
            if not short_name:
                short_name = (
                    mod.getFullName()
                    if hasattr(mod, "getFullName")
                    else f"UniMod:{accession}"
                )
            if short_name not in modifications:
                modifications[short_name] = {
                    "name": short_name,
                    "accession": accession,
                    "positions": [],
                }
            modifications[short_name]["positions"].append(
                {"position": "N-term.0", "scores": []}
            )

        if sequence.hasCTerminalModification():
            mod = sequence.getCTerminalModification()
            accession = mod.getUniModAccession()
            short_name = mod.getId()
            if not short_name:
                short_name = (
                    mod.getFullName()
                    if hasattr(mod, "getFullName")
                    else f"UniMod:{accession}"
                )
            if short_name not in modifications:
                modifications[short_name] = {
                    "name": short_name,
                    "accession": accession,
                    "positions": [],
                }
            modifications[short_name]["positions"].append(
                {"position": f"C-term.{sequence.size()+1}", "scores": []}
            )

        for i in range(sequence.size()):
            residue = sequence.getResidue(i)
            if residue.isModified():
                mod = residue.getModification()
                accession = mod.getUniModAccession()
                short_name = mod.getId()
                if not short_name:
                    short_name = (
                        mod.getFullName()
                        if hasattr(mod, "getFullName")
                        else f"UniMod:{accession}"
                    )
                if short_name not in modifications:
                    modifications[short_name] = {
                        "name": short_name,
                        "accession": accession,
                        "positions": [],
                    }
                modifications[short_name]["positions"].append(
                    {"position": f"{residue.getOneLetterCode()}.{i+1}", "scores": []}
                )

        return list(modifications.values()) if modifications else None

    except Exception:
        return None


class MaxQuant:
    """MaxQuant data processor for proteomics data conversion."""

    def __init__(self):
        self.mods_map = {}
        self._automaton = None
        self._intensity_normalize_names = []
        self._intensity_names = []
        self._sample_map = {}
        self.experiment_type = "LFQ"

    def get_mods_map(self, line: str):
        """Extract modification mapping from column headers."""
        mods_map = {}
        for col in line.split("\t"):
            if " Probabilities" in col:
                mod_name = col.replace(" Probabilities", "")
                mods_map[mod_name] = mod_name
        return mods_map

    # ===== CORE FILE READING FUNCTIONS =====

    def read_evidence(self, file_path):
        """Read evidence.txt file."""
        return pd.read_csv(file_path, sep="\t", low_memory=False)

    def read_protein_groups(self, file_path):
        """Read proteinGroups.txt file."""
        return pd.read_csv(file_path, sep="\t", low_memory=False)

    def read_msms(self, file_path):
        """Read msms.txt file."""
        return pd.read_csv(file_path, sep="\t", low_memory=False)

    # ===== ADVANCED BATCH PROCESSING =====

    def _parse_peptide_safely(self, peptidoform):
        """Safe peptide parsing with fallback strategies."""
        cleaned_peptidoform = clean_peptidoform(peptidoform)
        if not cleaned_peptidoform:
            raise ValueError(f"Empty peptidoform after cleaning: {peptidoform}")

        strategies = [
            ("fromString", lambda x: AASequence.fromString(x)),
            ("constructor", lambda x: AASequence(x)),
        ]

        for strategy_name, strategy_func in strategies:
            try:
                seq = strategy_func(cleaned_peptidoform)
                logging.debug(
                    f"Successfully parsed {peptidoform} using {strategy_name} strategy"
                )
                return seq
            except Exception as e:
                logging.debug(f"{strategy_name} strategy failed for {peptidoform}: {e}")
                continue

        logging.warning(
            f"All pyopenms parsing strategies failed for: {peptidoform} (cleaned: {cleaned_peptidoform})"
        )
        raise ValueError(f"Cannot parse peptide sequence: {peptidoform}")

    def generate_calculated_mz_batch(self, df):
        """Calculate m/z with error handling."""

        def safe_get_sequence(peptidoform):
            """Safely get AASequence object."""
            try:
                cleaned_peptidoform = clean_peptidoform(peptidoform)
                if not cleaned_peptidoform:
                    return None
                return AASequence.fromString(cleaned_peptidoform)
            except Exception as e1:
                try:
                    return AASequence(cleaned_peptidoform)
                except Exception as e2:
                    logging.warning(
                        f"Cannot parse peptide {peptidoform} (cleaned: {cleaned_peptidoform}): fromString={e1}, constructor={e2}"
                    )
                    return None

        uniq_p = df["peptidoform"].unique()
        sequence_map = {}
        failed_peptides = []

        for peptide in uniq_p:
            sequence = safe_get_sequence(peptide)
            if sequence is not None:
                sequence_map[peptide] = sequence
            else:
                failed_peptides.append(peptide)
                sequence_map[peptide] = None

        if failed_peptides:
            logging.warning(
                f"Failed to parse {len(failed_peptides)} peptides: {failed_peptides[:5]}{'...' if len(failed_peptides) > 5 else ''}"
            )

        def calculate_mz_for_row(row):
            sequence = sequence_map.get(row["peptidoform"])
            if sequence is None:
                return 0.0
            charge = row["precursor_charge"]
            if charge > 0:
                return sequence.getMZ(charge)
            return sequence.getMonoWeight()

        df.loc[:, "calculated_mz"] = df.apply(calculate_mz_for_row, axis=1)

    def iter_batch(
        self,
        file_path: Union[Path, str],
        label: str = "feature",
        chunksize: int = 100000,
        protein_str: str = None,
    ):
        """Iterate over file in batches."""
        for df in pd.read_csv(
            file_path, sep="\t", low_memory=False, chunksize=chunksize
        ):
            if protein_str and "Proteins" in df.columns:
                df = df[df["Proteins"].str.contains(protein_str, na=False)]
            yield df

    @staticmethod
    def open_from_zip_archive(zip_file, file_name, **kwargs):
        """Open file from zip archive."""
        with zipfile.ZipFile(zip_file) as z:
            with z.open(file_name) as f:
                df = pd.read_csv(f, sep="\t", low_memory=False, **kwargs)
        return df

    def read_zip_file(self, zip_path: str, **kwargs):
        """Read evidence.txt from zip file."""
        filepath = Path(zip_path)
        df = self.open_from_zip_archive(
            zip_path, f"{filepath.stem}/evidence.txt", **kwargs
        )
        return df

    # ===== CORE PROCESSING FUNCTIONS =====

    def process_evidence_to_feature_table(self, df: pd.DataFrame) -> pa.Table:
        """Process evidence DataFrame to feature table."""

        column_mapping = {
            "Sequence": "sequence",
            "Modified sequence": "peptidoform",
            "Charge": "precursor_charge",
            "PEP": "posterior_error_probability",
            "m/z": "observed_mz",
            "Raw file": "reference_file_name",
            "MS/MS scan number": "scan",
            "Retention time": "rt",
            "Leading razor protein": "anchor_protein",
            "Score": "andromeda_score",
            "Delta score": "andromeda_delta_score",
            "PIF": "parent_ion_fraction",
            "Intensity": "intensity",
        }

        processed_df = df.rename(columns=column_mapping)

        processed_df["calculated_mz"] = processed_df.apply(
            calculate_theoretical_mz, axis=1
        )

        processed_df["precursor_charge"] = processed_df["precursor_charge"].astype(
            "int32"
        )
        processed_df["posterior_error_probability"] = processed_df[
            "posterior_error_probability"
        ].astype("float32")
        processed_df["observed_mz"] = processed_df["observed_mz"].astype("float32")
        processed_df["calculated_mz"] = processed_df["calculated_mz"].astype("float32")
        processed_df["scan"] = processed_df["scan"].astype(str)

        processed_df["rt"] = (processed_df["rt"] * 60).astype("float32")

        processed_df["is_decoy"] = (
            df["Reverse"].apply(convert_maxquant_flag).astype("int32")
        )

        processed_df["pg_accessions"] = (
            df["Protein group IDs"].fillna("").astype(str).str.split(";")
        )
        processed_df["gg_names"] = (
            df["Gene names"].fillna("").astype(str).str.split(";")
        )

        def create_intensity_struct(row):
            if (
                "intensity" in row
                and pd.notna(row["intensity"])
                and row["intensity"] > 0
            ):
                return [
                    {
                        "sample_accession": row["reference_file_name"],
                        "channel": "MS",
                        "intensity": float(row["intensity"]),
                    }
                ]
            else:
                return None

        processed_df["intensities"] = processed_df.apply(
            create_intensity_struct, axis=1
        )

        processed_df["unique"] = (
            df["Proteins"].apply(self._is_unique_protein).astype("int32")
        )

        processed_df["modifications"] = processed_df["peptidoform"].apply(
            parse_modifications_from_peptidoform
        )

        processed_df["additional_scores"] = processed_df.apply(
            self._create_processed_additional_scores, axis=1
        )

        processed_df["cv_params"] = processed_df["parent_ion_fraction"].apply(
            self._create_processed_cv_params
        )

        processed_df["ion_mobility"] = None
        processed_df["pg_global_qvalue"] = None
        processed_df["start_ion_mobility"] = None
        processed_df["stop_ion_mobility"] = None
        processed_df["gg_accessions"] = None
        processed_df["scan_reference_file_name"] = None
        processed_df["rt_start"] = None
        processed_df["rt_stop"] = None
        processed_df["additional_intensities"] = None

        schema_fields = [field.name for field in FEATURE_SCHEMA]

        final_df = pd.DataFrame()
        for field in schema_fields:
            if field in processed_df.columns:
                final_df[field] = processed_df[field]
            else:
                final_df[field] = None

        final_df = final_df[schema_fields]

        return pa.Table.from_pandas(
            final_df, schema=FEATURE_SCHEMA, preserve_index=False
        )

    def process_protein_groups_to_pg_table(
        self, df: pd.DataFrame, sdrf_path: str = None
    ) -> pa.Table:
        """Process protein groups DataFrame to PG table."""
        df = self._apply_pg_column_mapping(df)
        df = self._process_pg_list_columns(df)
        df = self._process_pg_flag_columns(df)
        df = self._add_anchor_protein(df)
        df = self._process_pg_intensities(df, sdrf_path)
        df = self._process_pg_additional_scores(df)
        df["reference_file_name"] = "proteinGroups.txt"
        return self._create_pg_final_table(df)

    def _apply_pg_column_mapping(self, df: pd.DataFrame) -> pd.DataFrame:
        """Apply column mapping for protein groups."""
        pg_column_mapping = {
            "Protein IDs": "pg_accessions",
            "Fasta headers": "pg_names",
            "Gene names": "gg_accessions",
            "Q-value": "global_qvalue",
            "Reverse": "is_decoy",
            "Potential contaminant": "contaminant",
            "Score": "andromeda_score",
            "Sequence coverage [%]": "sequence_coverage",
            "Mol. weight [kDa]": "molecular_weight",
            "MS/MS count": "msms_count",
        }
        df.rename(columns=pg_column_mapping, inplace=True)
        return df

    def _process_pg_list_columns(self, df: pd.DataFrame) -> pd.DataFrame:
        """Process columns that should be split into lists."""
        list_columns = ["pg_accessions", "pg_names", "gg_accessions"]
        for col in list_columns:
            if col in df.columns:
                df[col] = df[col].fillna("").astype(str).str.split(";")
        return df

    def _process_pg_flag_columns(self, df: pd.DataFrame) -> pd.DataFrame:
        """Process flag columns (is_decoy, contaminant)."""
        flag_columns = ["is_decoy", "contaminant"]
        for col in flag_columns:
            if col in df.columns:
                df[col] = df[col].apply(convert_maxquant_flag).astype("int32")
        return df

    def _add_anchor_protein(self, df: pd.DataFrame) -> pd.DataFrame:
        """Add anchor protein (first protein in the group)."""
        df["anchor_protein"] = df["pg_accessions"].apply(self._get_anchor_protein)
        return df

    def _get_anchor_protein(self, pg_accessions):
        """Get the first protein from a list of protein accessions."""
        return (
            pg_accessions[0]
            if isinstance(pg_accessions, list) and len(pg_accessions) > 0
            else None
        )

    def _process_pg_intensities(
        self, df: pd.DataFrame, sdrf_path: str = None
    ) -> pd.DataFrame:
        """Process intensity data for protein groups."""
        tmt_channels = self.get_tmt_channels_from_sdrf(sdrf_path)
        intensity_cols = [col for col in df.columns if col.startswith("Intensity ")]
        lfq_cols = [col for col in df.columns if col.startswith("LFQ intensity ")]

        intensity_results = df.apply(
            lambda row: self._create_pg_intensity_struct(
                row, tmt_channels, intensity_cols, lfq_cols
            ),
            axis=1,
        )

        df["intensities"] = intensity_results.str[0]
        df["additional_intensities"] = intensity_results.str[1]

        return df

    def _create_pg_intensity_struct(self, row, tmt_channels, intensity_cols, lfq_cols):
        """Create intensity structure for a single protein group row."""
        intensities = []
        additional_intensities = []

        sample_accession = self._get_sample_accession(row)
        intensities.extend(self._process_basic_intensities(row, intensity_cols))
        additional_intensities.extend(self._process_lfq_intensities(row, lfq_cols))

        if tmt_channels:
            tmt_intensities, tmt_additional = self._process_tmt_intensities(
                row, tmt_channels, sample_accession
            )
            intensities.extend(tmt_intensities)
            additional_intensities.extend(tmt_additional)

        return (
            intensities if intensities else None,
            additional_intensities if additional_intensities else None,
        )

    def _get_sample_accession(self, row) -> str:
        """Get sample accession from row data."""
        for col_name in ["reference_file_name", "Raw file"]:
            if col_name in row.index and pd.notna(row[col_name]):
                return row[col_name]
        return "Unknown"

    def _process_basic_intensities(self, row, intensity_cols) -> list:
        """Process basic intensity columns."""
        intensities = []
        for col in intensity_cols:
            if col in row.index and pd.notna(row[col]) and row[col] > 0:
                sample = col.replace("Intensity ", "")
                intensities.append(
                    {
                        "sample_accession": sample,
                        "channel": "LFQ",
                        "intensity": float(row[col]),
                    }
                )
        return intensities

    def _process_lfq_intensities(self, row, lfq_cols) -> list:
        """Process LFQ intensity columns."""
        additional_intensities = []
        for col in lfq_cols:
            if col in row.index and pd.notna(row[col]) and row[col] > 0:
                sample = col.replace("LFQ intensity ", "")
                additional_intensities.append(
                    {
                        "sample_accession": sample,
                        "channel": "LFQ",
                        "intensities": [
                            {
                                "intensity_name": "lfq_intensity",
                                "intensity_value": float(row[col]),
                            }
                        ],
                    }
                )
        return additional_intensities

    def _process_tmt_intensities(self, row, tmt_channels, sample_accession) -> tuple:
        """Process TMT intensity columns."""
        intensities = []
        additional_intensities = []

        for i in range(min(8, len(tmt_channels))):
            reporter_col = f"Reporter intensity {i}"
            corrected_col = f"Reporter intensity corrected {i}"

            if (
                reporter_col in row.index
                and pd.notna(row[reporter_col])
                and row[reporter_col] > 0
            ):

                channel_name = tmt_channels[i]

                intensities.append(
                    {
                        "sample_accession": sample_accession,
                        "channel": channel_name,
                        "intensity": float(row[reporter_col]),
                    }
                )

                if corrected_col in row.index and pd.notna(row[corrected_col]):
                    additional_intensities.append(
                        {
                            "sample_accession": sample_accession,
                            "channel": channel_name,
                            "intensities": [
                                {
                                    "intensity_name": "corrected_intensity",
                                    "intensity_value": float(row[corrected_col]),
                                }
                            ],
                        }
                    )

        return intensities, additional_intensities

    def _process_pg_additional_scores(self, df: pd.DataFrame) -> pd.DataFrame:
        """Process additional scores for protein groups."""
        df["additional_scores"] = df.apply(self._create_pg_additional_scores, axis=1)
        return df

    def _create_pg_additional_scores(self, row) -> list:
        """Create additional scores structure for a single row."""
        if not pd.notna(row["andromeda_score"]):
            return None

        return [
            {
                "score_name": "andromeda_score",
                "score_value": float(row["andromeda_score"]),
            },
            {
                "score_name": "sequence_coverage",
                "score_value": float(row["sequence_coverage"]),
            },
            {
                "score_name": "molecular_weight",
                "score_value": float(row["molecular_weight"]),
            },
            {
                "score_name": "msms_count",
                "score_value": float(row["msms_count"]),
            },
        ]

    def _create_pg_final_table(self, df: pd.DataFrame) -> pa.Table:
        """Create final PyArrow table with proper schema."""
        schema_fields = [field.name for field in PG_SCHEMA]
        final_df = pd.DataFrame()

        for field in schema_fields:
            if field in df.columns:
                final_df[field] = df[field]

        for field in schema_fields:
            if field not in final_df.columns:
                final_df[field] = None

        final_df = final_df[schema_fields]
        return pa.Table.from_pandas(final_df, schema=PG_SCHEMA, preserve_index=False)

    def process_msms_to_psm_table(self, df: pd.DataFrame) -> pa.Table:
        """Process MS/MS DataFrame to PSM table."""

        psm_column_mapping = {
            "Sequence": "sequence",
            "Modified sequence": "peptidoform",
            "Charge": "precursor_charge",
            "PEP": "posterior_error_probability",
            "m/z": "observed_mz",
            "Raw file": "reference_file_name",
            "Scan number": "scan",
            "Retention time": "rt",
            "Proteins": "protein_accessions",
            "Score": "andromeda_score",
            "Delta score": "andromeda_delta_score",
            "PIF": "parent_ion_fraction",
            "Reverse": "is_decoy",
        }

        df.rename(columns=psm_column_mapping, inplace=True)

        df["calculated_mz"] = df.apply(calculate_theoretical_mz, axis=1)

        df["precursor_charge"] = df["precursor_charge"].astype("int32")
        df["posterior_error_probability"] = df["posterior_error_probability"].astype(
            "float32"
        )
        df["observed_mz"] = df["observed_mz"].astype("float32")
        df["calculated_mz"] = df["calculated_mz"].astype("float32")
        df["scan"] = df["scan"].astype(str)

        df["rt"] = (df["rt"] * 60).astype("float32")
        df["is_decoy"] = df["is_decoy"].apply(convert_maxquant_flag).astype("int32")

        if "protein_accessions" in df.columns:
            df["protein_accessions"] = (
                df["protein_accessions"].fillna("").astype(str).str.split(";")
            )

        df["modifications"] = df["peptidoform"].apply(
            parse_modifications_from_peptidoform
        )

        df["additional_scores"] = df.apply(
            lambda row: (
                [
                    {
                        "score_name": "andromeda_score",
                        "score_value": float(row["andromeda_score"]),
                    },
                    {
                        "score_name": "andromeda_delta_score",
                        "score_value": float(row["andromeda_delta_score"]),
                    },
                ]
                if pd.notna(row["andromeda_score"])
                else None
            ),
            axis=1,
        )

        if "parent_ion_fraction" in df.columns:
            df["cv_params"] = df["parent_ion_fraction"].apply(
                lambda pif: (
                    [{"cv_name": "parent_ion_fraction", "cv_value": str(pif)}]
                    if pd.notna(pif)
                    else None
                )
            )
        else:
            df["cv_params"] = None

        df["ion_mobility"] = None
        df["predicted_rt"] = None
        df["number_peaks"] = None
        df["mz_array"] = None
        df["intensity_array"] = None

        schema_fields = [field.name for field in PSM_SCHEMA]
        final_df = df[[col for col in schema_fields if col in df.columns]]
        for field in schema_fields:
            if field not in final_df.columns:
                final_df[field] = None
        final_df = final_df[schema_fields]

        return pa.Table.from_pandas(final_df, schema=PSM_SCHEMA, preserve_index=False)

    # ===== HIGH-LEVEL FILE WRITING METHODS =====

    def extract_col_msg(self, col_df, label: str = "feature"):
        """Extract column mapping for MaxQuant files."""
        available_columns = set(col_df.columns)

        if label == "feature":
            intensity_normalize_names = []
            intensity_names = []
            use_cols = MAXQUANT_FEATURE_USECOLS.copy()
            use_map = MAXQUANT_FEATURE_MAP.copy()

            use_cols = [col for col in use_cols if col in available_columns]
            use_map = {k: v for k, v in use_map.items() if k in available_columns}

            for col in col_df.columns:
                if re.search(intensity_normalize_pattern, col, re.IGNORECASE):
                    use_cols.append(col)
                    intensity_normalize_names.append(col)
                elif re.search(intensity_pattern, col, re.IGNORECASE):
                    use_cols.append(col)
                    intensity_names.append(col)

            self._intensity_normalize_names = intensity_normalize_names
            self._intensity_names = intensity_names

            if "Intensity" in available_columns:
                use_cols.append("Intensity")
        else:
            use_cols = MAXQUANT_PSM_USECOLS.copy()
            use_map = MAXQUANT_PSM_MAP.copy()
            use_cols = [col for col in use_cols if col in available_columns]
            use_map = {k: v for k, v in use_map.items() if k in available_columns}

        line = "\t".join(col_df.columns)
        self.mods_map = self.get_mods_map(line)
        self._automaton = get_ahocorasick(self.mods_map)

        for key in self.mods_map.keys():
            col = f"{key} Probabilities"
            if col in col_df.columns:
                use_cols.append(col)

        return use_map, use_cols

    def iter_batch_with_mapping(
        self,
        file_path: Union[Path, str],
        label: str = "feature",
        chunksize: int = 100000,
        protein_str: str = None,
    ):
        """Iterate through file batches with column mapping."""
        col_df = pd.read_csv(file_path, sep="\t", nrows=0)
        use_map, use_cols = self.extract_col_msg(col_df, label=label)

        for df in pd.read_csv(
            file_path, sep="\t", usecols=use_cols, low_memory=False, chunksize=chunksize
        ):
            df.rename(columns=use_map, inplace=True)
            if protein_str:
                df = df[df["mp_accessions"].str.contains(f"{protein_str}", na=False)]
            df = self.main_operate(df)
            yield df

    def main_operate(self, df: pd.DataFrame):
        """Main data processing operations."""
        self.generate_peptidoform(df)
        self.generate_calculated_mz(df)
        self.generate_modification_details(df)

        df["posterior_error_probability"] = pd.to_numeric(
            df["posterior_error_probability"], errors="coerce"
        )
        df = df[
            (df["posterior_error_probability"] < 0.05)
            | (df["posterior_error_probability"].isna())
        ].copy()

        df["is_decoy"] = (
            df["is_decoy"].map({None: 0, np.nan: 0, "+": 1}).astype("int32")
        )
        df["andromeda_score"] = pd.to_numeric(df["andromeda_score"], errors="coerce")
        df["andromeda_delta_score"] = pd.to_numeric(
            df["andromeda_delta_score"], errors="coerce"
        )

        df["additional_scores"] = df[
            ["andromeda_score", "andromeda_delta_score"]
        ].apply(self._create_additional_scores, axis=1)
        df.loc[:, "cv_params"] = df["parent_ion_fraction"].apply(self._create_cv_params)

        df["predicted_rt"] = pd.Series([np.nan] * len(df), dtype="float32")
        df["ion_mobility"] = pd.Series([np.nan] * len(df), dtype="float32")
        return df

    def _create_additional_scores(self, row):
        """Create additional scores structure for a row."""
        return [
            {
                "score_name": "andromeda_score",
                "score_value": row["andromeda_score"],
            },
            {
                "score_name": "andromeda_delta_score",
                "score_value": row["andromeda_delta_score"],
            },
        ]

    def _create_cv_params(self, score):
        """Create CV params structure for a score."""
        return [{"cv_name": "parent_ion_fraction", "cv_value": str(score)}]

    def _split_protein_accessions(self, protein_accessions):
        """Split protein accessions string into list."""
        return (
            protein_accessions.split(";")
            if isinstance(protein_accessions, str) and protein_accessions
            else []
        )

    def _is_unique_protein(self, proteins):
        """Check if protein is unique (no semicolon)."""
        return 1 if isinstance(proteins, str) and ";" not in proteins else 0

    def _create_processed_additional_scores(self, row):
        """Create additional scores for processed data."""
        if pd.notna(row["andromeda_score"]) and pd.notna(row["andromeda_delta_score"]):
            return [
                {
                    "score_name": "andromeda_score",
                    "score_value": float(row["andromeda_score"]),
                },
                {
                    "score_name": "andromeda_delta_score",
                    "score_value": float(row["andromeda_delta_score"]),
                },
            ]
        return None

    def _create_processed_cv_params(self, pif):
        """Create CV params for processed data."""
        return (
            [{"cv_name": "parent_ion_fraction", "cv_value": str(pif)}]
            if pd.notna(pif)
            else None
        )

    def generate_peptidoform(self, df):
        """Generate cleaned peptidoform."""
        df["peptidoform"] = df["peptidoform"].apply(clean_peptidoform)

    def generate_calculated_mz(self, df):
        """Calculate m/z values."""
        df["calculated_mz"] = df.apply(calculate_theoretical_mz, axis=1)

    def generate_modification_details(self, df):
        """Generate modification details."""
        df["modifications"] = df["peptidoform"].apply(
            parse_modifications_from_peptidoform
        )

    def transform_psm(self, df: pd.DataFrame):
        """Transform PSM data for schema compliance."""
        if "mp_accessions" in df.columns:
            df.rename(columns={"mp_accessions": "protein_accessions"}, inplace=True)

        if "protein_accessions" in df.columns:
            df["protein_accessions"] = df["protein_accessions"].apply(
                self._split_protein_accessions
            )

        if "scan" in df.columns:
            df["scan"] = df["scan"].astype(str)

        df.loc[:, "intensity_array"] = None
        df.loc[:, "mz_array"] = None
        df.loc[:, "number_peaks"] = None

    def transform_feature(self, df: pd.DataFrame, sdrf_path: str = None):
        """Transform feature data for schema compliance."""
        df = self._apply_feature_column_mapping(df)
        df = self._process_feature_protein_groups(df)
        df = self._process_feature_unique_peptides(df)
        df = self._process_feature_gene_names(df)
        df = self._add_feature_anchor_protein(df)
        df = self._convert_feature_scan_to_string(df)
        df = self._initialize_feature_schema_fields(df)
        df = self._process_feature_modifications(df)
        df = self._process_feature_intensities(df, sdrf_path)

        return df

    def _apply_feature_column_mapping(self, df: pd.DataFrame) -> pd.DataFrame:
        """Apply column mapping for features."""
        from ..common import MAXQUANT_FEATURE_MAP

        columns_to_map = {}
        for old_col, new_col in MAXQUANT_FEATURE_MAP.items():
            if old_col in df.columns:
                columns_to_map[old_col] = new_col

        if columns_to_map:
            df.rename(columns=columns_to_map, inplace=True)

        return df

    def _process_feature_protein_groups(self, df: pd.DataFrame) -> pd.DataFrame:
        """Process protein group columns for features."""
        if "pg_accessions" not in df.columns:
            if "Protein group IDs" in df.columns:
                df["pg_accessions"] = (
                    df["Protein group IDs"].fillna("").astype(str).str.split(";")
                )
            else:
                df["pg_accessions"] = [[]] * len(df)

        return df

    def _process_feature_unique_peptides(self, df: pd.DataFrame) -> pd.DataFrame:
        """Process unique peptide detection for features."""

        def is_unique_peptide(pg_accessions):
            """Check if peptide is unique."""
            if not isinstance(pg_accessions, list):
                return 1
            if len(pg_accessions) > 1:
                return 0
            if len(pg_accessions) == 1 and (
                ";" in str(pg_accessions[0]) or "," in str(pg_accessions[0])
            ):
                return 0
            return 1

        df.loc[:, "unique"] = df["pg_accessions"].apply(is_unique_peptide)
        df["pg_accessions"] = df["pg_accessions"].apply(get_protein_accession)

        return df

    def _process_feature_gene_names(self, df: pd.DataFrame) -> pd.DataFrame:
        """Process gene names for features."""
        if "gg_names" in df.columns:
            df["gg_names"] = df["gg_names"].str.split(";")
        else:
            df["gg_names"] = [[]] * len(df)

        return df

    def _add_feature_anchor_protein(self, df: pd.DataFrame) -> pd.DataFrame:
        """Add anchor protein for features."""
        df.loc[:, "anchor_protein"] = df["pg_accessions"].str[0]
        return df

    def _convert_feature_scan_to_string(self, df: pd.DataFrame) -> pd.DataFrame:
        """Convert scan column to string for features."""
        if "scan" in df.columns:
            df["scan"] = df["scan"].astype(str)
        return df

    def _initialize_feature_schema_fields(self, df: pd.DataFrame) -> pd.DataFrame:
        """Initialize nullable FEATURE schema fields."""
        df.loc[:, "gg_accessions"] = [[]] * len(df)
        df["pg_global_qvalue"] = pd.Series([np.nan] * len(df), dtype="float32")
        df.loc[:, "scan_reference_file_name"] = None
        df["start_ion_mobility"] = pd.Series([np.nan] * len(df), dtype="float32")
        df["stop_ion_mobility"] = pd.Series([np.nan] * len(df), dtype="float32")
        df["rt_start"] = pd.Series([np.nan] * len(df), dtype="float32")
        df["rt_stop"] = pd.Series([np.nan] * len(df), dtype="float32")

        return df

    def _process_feature_modifications(self, df: pd.DataFrame) -> pd.DataFrame:
        """Process modifications field for features."""
        if "modifications" in df.columns:
            df["modifications"] = df["modifications"].apply(
                self._ensure_modifications_list
            )
        return df

    def _ensure_modifications_list(self, modifications):
        """Ensure modifications is a list."""
        return modifications if isinstance(modifications, list) else []

    def _process_feature_intensities(
        self, df: pd.DataFrame, sdrf_path: str = None
    ) -> pd.DataFrame:
        """Process intensity data for features."""
        if "intensity" in df.columns:
            intensity_results = df.apply(
                lambda row: self._create_feature_intensity_struct(row, sdrf_path),
                axis=1,
            )

            df["intensities"] = intensity_results.str[0]
            df["additional_intensities"] = intensity_results.str[1]
        else:
            df.loc[:, "intensities"] = [[]] * len(df)
            df.loc[:, "additional_intensities"] = [[]] * len(df)

        return df

    def _create_feature_intensity_struct(self, row, sdrf_path: str = None):
        """Create intensity structure for a single feature row."""
        intensities = []
        additional_intensities = []

        sample_accession = self._get_feature_sample_accession(row)
        intensity_value = self._get_feature_primary_intensity(row)
        if intensity_value is not None and intensity_value > 0:
            intensities.append(
                {
                    "sample_accession": sample_accession,
                    "channel": "LFQ",
                    "intensity": float(intensity_value),
                }
            )

        tmt_channels = self.get_tmt_channels_from_sdrf(sdrf_path)
        if tmt_channels:
            tmt_intensities, tmt_additional = self._process_feature_tmt_intensities(
                row, tmt_channels, sample_accession
            )
            intensities.extend(tmt_intensities)
            additional_intensities.extend(tmt_additional)

        return (
            intensities if intensities else None,
            additional_intensities if additional_intensities else None,
        )

    def _get_feature_sample_accession(self, row) -> str:
        """Get sample accession for feature row."""
        if "reference_file_name" in row.index:
            return row["reference_file_name"]
        elif "Raw file" in row.index:
            return row["Raw file"]
        else:
            return "Unknown"

    def _get_feature_primary_intensity(self, row):
        """Get primary intensity value for feature row."""
        if "intensity" in row.index and pd.notna(row["intensity"]):
            return row["intensity"]
        elif "Intensity" in row.index and pd.notna(row["Intensity"]):
            return row["Intensity"]
        return None

    def _process_feature_tmt_intensities(self, row, tmt_channels, sample_accession):
        """Process TMT intensities for feature row."""
        intensities = []
        additional_intensities = []

        for i in range(min(8, len(tmt_channels))):
            reporter_col = f"Reporter intensity {i}"
            corrected_col = f"Reporter intensity corrected {i}"

            if (
                reporter_col in row.index
                and pd.notna(row[reporter_col])
                and row[reporter_col] > 0
            ):

                channel_name = tmt_channels[i]

                intensities.append(
                    {
                        "sample_accession": sample_accession,
                        "channel": channel_name,
                        "intensity": float(row[reporter_col]),
                    }
                )

                if corrected_col in row.index and pd.notna(row[corrected_col]):
                    additional_intensities.append(
                        {
                            "sample_accession": sample_accession,
                            "channel": channel_name,
                            "intensities": [
                                {
                                    "intensity_name": "corrected_intensity",
                                    "intensity_value": float(row[corrected_col]),
                                }
                            ],
                        }
                    )

        return intensities, additional_intensities

    def _init_sdrf(self, sdrf_path: Union[Path, str]):
        """Initialize SDRF handler."""
        sdrf = SDRFHandler(sdrf_path)
        self.experiment_type = sdrf.get_experiment_type_from_sdrf()
        self._sample_map = sdrf.get_sample_map_run()

    def get_tmt_channels_from_sdrf(self, sdrf_path: str = None):
        """Get TMT channel mapping from SDRF."""
        try:
            sdrf_path_to_use = None

            if sdrf_path:
                sdrf_path_to_use = sdrf_path
            elif (
                hasattr(self, "sdrf_handler")
                and self.sdrf_handler
                and hasattr(self.sdrf_handler, "sdrf_path")
                and self.sdrf_handler.sdrf_path
            ):
                sdrf_path_to_use = self.sdrf_handler.sdrf_path
            elif hasattr(self, "_current_sdrf_path") and self._current_sdrf_path:
                sdrf_path_to_use = self._current_sdrf_path

            if sdrf_path_to_use:
                sdrf_df = pd.read_csv(sdrf_path_to_use, sep="\t")
                if "comment[label]" in sdrf_df.columns:
                    all_labels = sdrf_df["comment[label]"].dropna().tolist()
                    unique_labels = sorted(set(all_labels))
                    tmt_labels = [
                        label
                        for label in unique_labels
                        if any(
                            keyword in label.upper()
                            for keyword in ["TMT", "ITRAQ", "PLEX"]
                        )
                    ]
                    if tmt_labels:
                        return tmt_labels[:8]

            return []

        except Exception as e:
            print(f"Warning: Could not read TMT channels from SDRF: {e}")
            return []

    def write_psm_to_file(
        self, msms_path: str, output_path: str, chunksize: int = 1000000
    ):
        """Write PSM data to Parquet file."""
        batch_writer = ParquetBatchWriter(output_path, PSM_SCHEMA)

        try:
            for df in self.iter_batch_with_mapping(
                msms_path, label="psm", chunksize=chunksize
            ):
                if df is not None and not df.empty:
                    self.transform_psm(df)
                    records = df.to_dict("records")

                    if records and "intensities" in records[0]:
                        first_intensities = records[0]["intensities"]
                        print(
                            f"DEBUG: First record intensities type: {type(first_intensities)}"
                        )
                        print(f"DEBUG: First record intensities: {first_intensities}")
                        print(
                            f"DEBUG: Has intensities data: {len(first_intensities) > 0}"
                        )

                    batch_writer.write_batch(records)
        finally:
            if batch_writer:
                batch_writer.close()
            logging.info(f"PSM file written to {output_path}")

    def write_feature_to_file(
        self,
        evidence_path: str,
        sdrf_path: str,
        output_path: str,
        chunksize: int = 1000000,
        protein_file=None,
    ):
        """Write feature data to Parquet file."""
        self._init_sdrf(sdrf_path)
        protein_list = extract_protein_list(protein_file) if protein_file else None
        protein_str = "|".join(protein_list) if protein_list else None

        batch_writer = ParquetBatchWriter(output_path, FEATURE_SCHEMA)

        try:
            for df in self.iter_batch_with_mapping(
                evidence_path, "feature", chunksize, protein_str=protein_str
            ):
                if df is not None and not df.empty:
                    self.transform_feature(df, sdrf_path)

                    Feature.convert_to_parquet_format(df)

                    records = df.to_dict("records")
                    batch_writer.write_batch(records)
        finally:
            if batch_writer:
                batch_writer.close()
            logging.info(f"Feature file written to {output_path}")

    def write_protein_groups_to_file(
        self, protein_groups_path: str, sdrf_path: str, output_path: str
    ):
        """Write protein groups data to Parquet file."""
        self._init_sdrf(sdrf_path)

        df = self.read_protein_groups(protein_groups_path)
        table = self.process_protein_groups_to_pg_table(df, sdrf_path)

        pq.write_table(table, output_path)
        logging.info(f"Protein groups file written to {output_path}")


# ===== CONVENIENCE FUNCTIONS FOR BACKWARD COMPATIBILITY =====


# Backward compatibility functions - simplified using class methods directly
def read_evidence(file_path):
    """Backward compatibility function."""
    return MaxQuant().read_evidence(file_path)


def read_protein_groups(file_path):
    """Backward compatibility function."""
    return MaxQuant().read_protein_groups(file_path)


def read_msms(file_path):
    """Backward compatibility function."""
    return MaxQuant().read_msms(file_path)


def process_evidence_to_feature_table(df: pd.DataFrame) -> pa.Table:
    """Backward compatibility function."""
    return MaxQuant().process_evidence_to_feature_table(df)


def process_protein_groups_to_pg_table(
    df: pd.DataFrame, sdrf_path: str = None
) -> pa.Table:
    """Backward compatibility function."""
    return MaxQuant().process_protein_groups_to_pg_table(df, sdrf_path)


def process_msms_to_psm_table(df: pd.DataFrame) -> pa.Table:
    """Backward compatibility function."""
    return MaxQuant().process_msms_to_psm_table(df)
