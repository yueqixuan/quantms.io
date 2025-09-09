"""MaxQuant data processing module"""

import logging
import re
from pathlib import Path
from typing import Union, List, Dict, Optional, Tuple

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
from pyopenms import AASequence
from pyopenms.Constants import PROTON_MASS_U

from quantmsio.core.format import PG_SCHEMA, FEATURE_SCHEMA, PSM_SCHEMA
from quantmsio.core.common import (
    MAXQUANT_PSM_MAP,
    MAXQUANT_FEATURE_MAP,
    MAXQUANT_PG_MAP,
    SDRF_MAP,
    MAXQUANT_PSM_USECOLS,
    MAXQUANT_FEATURE_USECOLS,
    MAXQUANT_PG_USECOLS,
    SDRF_USECOLS,
)
from quantmsio.core.sdrf import SDRFHandler
from quantmsio.utils.file_utils import ParquetBatchWriter, extract_protein_list


def clean_peptidoform(peptidoform):
    """Clean peptidoform string by normalizing modification names"""
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


def convert_maxquant_flag(value):
    """Convert MaxQuant flag from + to 1, others to 0"""
    return 1 if value == "+" else 0


def parse_modifications_from_peptidoform(peptidoform: str) -> list:
    """Parse modification information from peptidoform string"""
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


logging.basicConfig(format="%(asctime)s - %(message)s", level=logging.INFO)


class MaxQuant:
    """MaxQuant data processor for converting output files to quantms.io format"""

    # ============================================================================
    # Initialization and Setup
    # ============================================================================

    def __init__(self):
        self.sdrf_handler: Optional[SDRFHandler] = None
        self.experiment_type: Optional[str] = None
        self._sample_map: Optional[Dict] = None
        self._channel_map: Optional[Dict] = None
        self._current_sdrf_path: Optional[str] = None

        self.psm_mapping = MAXQUANT_PSM_MAP
        self.feature_mapping = MAXQUANT_FEATURE_MAP
        self.pg_mapping = MAXQUANT_PG_MAP
        self.sdrf_mapping = SDRF_MAP

        self._sequence_cache: Dict[str, AASequence] = {}

    # ============================================================================
    # SDRF Processing and Mapping
    # ============================================================================

    def _init_sdrf(self, sdrf_path: Union[Path, str]) -> None:
        """Initialize SDRF handler and create mappings"""
        if not sdrf_path:
            return

        self.sdrf_handler = SDRFHandler(sdrf_path)
        self.experiment_type = self.sdrf_handler.get_experiment_type_from_sdrf()
        self._sample_map = self.sdrf_handler.get_sample_map_run()
        self._current_sdrf_path = sdrf_path

        self._sdrf_transformed = self._create_basic_sdrf_mapping()
        self._channel_map = self._create_simplified_channel_mapping()

    def _create_basic_sdrf_mapping(self) -> pd.DataFrame:
        """Create basic mapping table from SDRF file"""
        sdrf_table = self.sdrf_handler.sdrf_table.copy()

        basic_mapping = {
            "comment[data file]": "reference_file_name",
            "comment[label]": "channel",
            "source name": "sample_accession",
        }

        rename_map = {}
        for orig_col, new_col in basic_mapping.items():
            matching_col = None
            for col in sdrf_table.columns:
                if col.lower() == orig_col.lower():
                    matching_col = col
                    break

            if matching_col:
                rename_map[matching_col] = new_col

        available_cols = [col for col in sdrf_table.columns if col in rename_map]
        basic_df = sdrf_table[available_cols].copy()
        basic_df = basic_df.rename(columns=rename_map)

        if "reference_file_name" in basic_df.columns:
            basic_df["reference_file_name"] = (
                basic_df["reference_file_name"].str.split(".").str[0]
            )
        if "reference_file_name" in basic_df.columns:
            basic_df.set_index("reference_file_name", inplace=True)

        return basic_df

    def _create_simplified_channel_mapping(self) -> Dict[str, str]:
        """Create channel mapping using basic SDRF information"""
        if not hasattr(self, "_sdrf_transformed") or self._sdrf_transformed is None:
            return {}

        channel_map = {}

        df = (
            self._sdrf_transformed.reset_index()
            if self._sdrf_transformed.index.name
            else self._sdrf_transformed
        )

        for _, row in df.iterrows():
            if "reference_file_name" in row:
                file_key = row["reference_file_name"]
            elif self._sdrf_transformed.index.name == "reference_file_name":
                file_key = row.name if hasattr(row, "name") else str(row.iloc[0])
            else:
                continue

            if self.experiment_type == "LFQ":
                if "channel" in row and pd.notna(row["channel"]):
                    channel_map[file_key] = row["channel"]
            else:
                if "channel" in row and "sample_accession" in row:
                    if pd.notna(row["channel"]) and pd.notna(row["sample_accession"]):
                        map_key = f"{file_key}-{row['channel']}"
                        channel_map[map_key] = row["sample_accession"]

        return channel_map

    def _get_sample_accession_from_sdrf(
        self, reference_file: str, channel: str = None
    ) -> Optional[str]:
        """Get sample accession from SDRF mapping"""
        if not self._sample_map:
            return reference_file if reference_file else None

        file_key = reference_file.split(".")[0] if reference_file else ""

        if self.experiment_type == "LFQ":
            return self._fuzzy_map_lookup(file_key, reference_file)
        else:
            if channel:
                map_key = f"{file_key}-{channel}"
                return self._sample_map.get(map_key, reference_file)

        return reference_file if reference_file else None

    def _fuzzy_map_lookup(self, file_key: str, fallback: str) -> str:
        """Fuzzy file name mapping lookup with error tolerance"""
        if file_key in self._sample_map:
            return self._sample_map[file_key]

        cleaned_keys = [
            file_key.rstrip("_"),
            file_key.rstrip("_-"),
            file_key.replace("__", "_"),
            file_key.replace("_.", "."),
        ]

        for cleaned_key in cleaned_keys:
            if cleaned_key != file_key and cleaned_key in self._sample_map:
                return self._sample_map[cleaned_key]

        if len(file_key) > 3:
            for map_key, sample in self._sample_map.items():
                if map_key.startswith(file_key) or file_key.startswith(map_key):
                    return sample

        if "_" in file_key:
            base_key = file_key.rsplit("_", 1)[0]
            for map_key, sample in self._sample_map.items():
                if map_key.startswith(base_key):
                    return sample

        return fallback

    def _get_channel_from_sdrf(
        self, reference_file: str, sample_accession: str = None
    ) -> Optional[str]:
        """Get channel information from transformed SDRF data"""
        if not hasattr(self, "_sdrf_transformed") or self._sdrf_transformed is None:
            return None

        file_key = reference_file.split(".")[0] if reference_file else ""

        if (
            hasattr(self._sdrf_transformed, "index")
            and self._sdrf_transformed.index.name == "reference_file_name"
        ):
            if file_key in self._sdrf_transformed.index:
                matching_rows = self._sdrf_transformed.loc[[file_key]]
            else:
                matching_rows = pd.DataFrame()
        elif "reference_file_name" in self._sdrf_transformed.columns:
            matching_rows = self._sdrf_transformed[
                self._sdrf_transformed["reference_file_name"] == file_key
            ]
        else:
            matching_rows = pd.DataFrame()

        if not matching_rows.empty:
            if self.experiment_type == "LFQ":
                channel = matching_rows.iloc[0]["channel"]
                return channel if pd.notna(channel) else None
            else:
                if sample_accession:
                    sample_matches = matching_rows[
                        matching_rows["sample_accession"] == sample_accession
                    ]
                    if not sample_matches.empty:
                        channel = sample_matches.iloc[0]["channel"]
                        return channel if pd.notna(channel) else None

        return None

    def _process_tmt_intensities(self, row, tmt_channels, reference_file_name) -> tuple:
        """Process TMT Reporter intensity columns"""
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

                sample_accession = self._get_tmt_sample_accession(
                    reference_file_name, channel_name
                )

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
                                    "intensity_name": corrected_col,  # Use original column name
                                    "intensity_value": float(row[corrected_col]),
                                }
                            ],
                        }
                    )

        return intensities, additional_intensities

    def _get_tmt_sample_accession(
        self, reference_file_name: str, channel_name: str
    ) -> str:
        """Get sample accession for TMT experiments"""
        if hasattr(self, "_sample_map") and self._sample_map:
            file_key = reference_file_name.split(".")[0]
            map_key = f"{file_key}-{channel_name}"
            sample_accession = self._sample_map.get(map_key)
            if sample_accession:
                return sample_accession

        return reference_file_name

    def _get_tmt_channels_from_sdrf(self) -> list:
        """Get TMT channel list from SDRF"""
        if not hasattr(self, "sdrf_handler") or not self.sdrf_handler:
            return []

        try:
            sdrf_table = self.sdrf_handler.sdrf_table
            if "comment[label]" in sdrf_table.columns:
                labels = sdrf_table["comment[label]"].unique()
                tmt_labels = [
                    label for label in labels if label and "TMT" in str(label).upper()
                ]
                return sorted(tmt_labels)
        except Exception:
            pass

        return []

    # ============================================================================
    # Common Processing Methods
    # ============================================================================

    def _calculate_theoretical_mz_batch(self, df: pd.DataFrame) -> None:
        """Calculate theoretical m/z values in batch"""

        def safe_parse_sequence(peptidoform: str) -> Optional[AASequence]:
            """Safely parse peptide sequence"""
            try:
                cleaned_peptidoform = clean_peptidoform(peptidoform)
                if not cleaned_peptidoform:
                    return None
                return AASequence.fromString(cleaned_peptidoform)
            except Exception as e1:
                try:
                    return AASequence(cleaned_peptidoform)
                except Exception as e2:
                    return None

        unique_peptidoforms = df["peptidoform"].unique()

        for peptidoform in unique_peptidoforms:
            if peptidoform not in self._sequence_cache:
                self._sequence_cache[peptidoform] = safe_parse_sequence(peptidoform)

        mass_map = {}
        for peptidoform in unique_peptidoforms:
            sequence = self._sequence_cache.get(peptidoform)
            if sequence:
                mass_map[peptidoform] = sequence.getMonoWeight()
            else:
                mass_map[peptidoform] = 0.0

        mass_vector = df["peptidoform"].map(mass_map)
        df.loc[:, "calculated_mz"] = (
            mass_vector + (PROTON_MASS_U * df["precursor_charge"])
        ) / df["precursor_charge"]

    # ============================================================================
    # PSM Processing
    # ============================================================================

    def process_psm_file(
        self, msms_path: str, output_path: str, chunksize: int = 1000000
    ) -> None:
        """Process PSM data from msms.txt to PSM parquet formatt"""

        batch_writer = ParquetBatchWriter(output_path, PSM_SCHEMA)

        try:
            for df_chunk in pd.read_csv(
                msms_path, sep="\t", chunksize=chunksize, low_memory=False
            ):

                df_chunk = self._apply_psm_mapping(df_chunk)
                self._calculate_theoretical_mz_batch(df_chunk)
                df_chunk = self._process_psm_modifications(df_chunk)
                df_chunk = self._process_psm_scores(df_chunk)
                df_chunk = self._ensure_psm_schema_compliance(df_chunk)

                table = pa.Table.from_pandas(
                    df_chunk, schema=PSM_SCHEMA, preserve_index=False
                )
                batch_writer.write_batch(table.to_pylist())

        finally:
            batch_writer.close()
            pass

    def _apply_psm_mapping(self, df: pd.DataFrame) -> pd.DataFrame:
        """Apply MAXQUANT_PSM_MAP mapping"""
        available_mapping = {
            k: v for k, v in self.psm_mapping.items() if k in df.columns
        }
        df.rename(columns=available_mapping, inplace=True)

        if "protein_accessions" in df.columns:
            df["protein_accessions"] = df["protein_accessions"].str.split(";")

        return df

    def _process_psm_modifications(self, df: pd.DataFrame) -> pd.DataFrame:
        """Process PSM modification information"""
        if "peptidoform" in df.columns:
            df["modifications"] = df["peptidoform"].apply(
                parse_modifications_from_peptidoform
            )
        else:
            df["modifications"] = None
        return df

    def _process_psm_scores(self, df: pd.DataFrame) -> pd.DataFrame:
        """Structure PSM scoring information"""
        scores_list = []
        for _, row in df.iterrows():
            scores = []
            if "andromeda_score" in row and pd.notna(row["andromeda_score"]):
                scores.append(
                    {
                        "score_name": "andromeda_score",
                        "score_value": float(row["andromeda_score"]),
                    }
                )
            if "andromeda_delta_score" in row and pd.notna(
                row["andromeda_delta_score"]
            ):
                scores.append(
                    {
                        "score_name": "andromeda_delta_score",
                        "score_value": float(row["andromeda_delta_score"]),
                    }
                )
            scores_list.append(scores if scores else None)

        df["additional_scores"] = scores_list
        return df

    def _ensure_psm_schema_compliance(self, df: pd.DataFrame) -> pd.DataFrame:
        """Ensure PSM data complies with PSM_SCHEMA"""
        schema_fields = [field.name for field in PSM_SCHEMA]

        for field in PSM_SCHEMA:
            if field.name not in df.columns:
                if field.type == pa.string():
                    df[field.name] = ""
                elif field.type == pa.int32():
                    df[field.name] = 0
                elif field.type == pa.float32():
                    df[field.name] = 0.0
                else:
                    df[field.name] = None

        df = df[[col for col in schema_fields if col in df.columns]].copy()

        if "is_decoy" in df.columns:
            df["is_decoy"] = df["is_decoy"].apply(convert_maxquant_flag).astype("int32")
        if "precursor_charge" in df.columns:
            df["precursor_charge"] = df["precursor_charge"].astype("int32")
        if "calculated_mz" in df.columns:
            df["calculated_mz"] = df["calculated_mz"].astype("float32")
        if "observed_mz" in df.columns:
            df["observed_mz"] = df["observed_mz"].astype("float32")
        if "rt" in df.columns:
            df["rt"] = (df["rt"] * 60).astype("float32")
        if "scan" in df.columns:
            df["scan"] = df["scan"].astype("string")

        return df

    # ============================================================================
    # Feature Processing
    # ============================================================================

    def process_feature_file(
        self,
        evidence_path: str,
        output_path: str,
        sdrf_path: str = None,
        protein_file: str = None,
        chunksize: int = 1000000,
    ) -> None:
        """Process Feature data from evidence.txt to Feature parquet format"""
        if sdrf_path:
            self._init_sdrf(sdrf_path)

        evidence_dir = Path(evidence_path).parent
        protein_groups_path = evidence_dir / "proteinGroups.txt"
        if protein_groups_path.exists():
            self._init_protein_group_qvalue_mapping(str(protein_groups_path))

        if protein_file:
            self._init_protein_group_qvalue_mapping(protein_file)

        protein_list = extract_protein_list(protein_file) if protein_file else None
        protein_str = "|".join(protein_list) if protein_list else None

        batch_writer = ParquetBatchWriter(output_path, FEATURE_SCHEMA)

        try:
            available_cols = pd.read_csv(
                evidence_path, sep="\t", nrows=0
            ).columns.tolist()

            usecols_filtered = [
                col for col in MAXQUANT_FEATURE_USECOLS if col in available_cols
            ]

            if self.experiment_type and "TMT" in self.experiment_type.upper():
                reporter_cols = [
                    col
                    for col in available_cols
                    if col.startswith("Reporter intensity")
                ]
                usecols_filtered.extend(reporter_cols)

            for df_chunk in pd.read_csv(
                evidence_path,
                sep="\t",
                chunksize=chunksize,
                usecols=usecols_filtered,
                low_memory=False,
            ):

                if protein_str:
                    df_chunk = df_chunk[
                        df_chunk["Proteins"].str.contains(
                            protein_str, na=False, regex=False
                        )
                    ]
                    if df_chunk.empty:
                        continue

                df_chunk = self._apply_feature_mapping(df_chunk)
                self._calculate_theoretical_mz_batch(df_chunk)
                df_chunk = self._process_feature_modifications(df_chunk)
                df_chunk = self._process_feature_protein_groups(df_chunk)

                if hasattr(self, "_protein_group_qvalue_map"):
                    df_chunk = self._map_protein_group_qvalue(df_chunk)

                df_chunk = self._process_feature_intensities(df_chunk)

                if self.sdrf_handler:
                    df_chunk = self._integrate_sdrf_metadata_feature(df_chunk)

                df_chunk = self._ensure_feature_schema_compliance(df_chunk)

                table = pa.Table.from_pandas(
                    df_chunk, schema=FEATURE_SCHEMA, preserve_index=False
                )
                batch_writer.write_batch(table.to_pylist())

        finally:
            batch_writer.close()
            pass

    def _apply_feature_mapping(self, df: pd.DataFrame) -> pd.DataFrame:
        """Apply MAXQUANT_FEATURE_MAP mapping"""
        available_mapping = {
            k: v for k, v in self.feature_mapping.items() if k in df.columns
        }
        df.rename(columns=available_mapping, inplace=True)
        return df

    def _process_feature_modifications(self, df: pd.DataFrame) -> pd.DataFrame:
        """Process Feature modification information"""
        return self._process_psm_modifications(df)

    def _process_feature_protein_groups(self, df: pd.DataFrame) -> pd.DataFrame:
        """Process Feature protein group information"""
        if "pg_accessions" in df.columns:
            df["pg_accessions"] = (
                df["pg_accessions"].fillna("").astype(str).str.split(";")
            )

        if "pg_accessions" in df.columns:
            df["anchor_protein"] = df["pg_accessions"].str[0]
        else:
            df["anchor_protein"] = None

        if "pg_accessions" in df.columns:
            df["unique"] = (df["pg_accessions"].str.len() == 1).astype("int32")
        else:
            df["unique"] = 0

        return df

    def _process_feature_intensities(self, df: pd.DataFrame) -> pd.DataFrame:
        """Process Feature intensity data"""
        intensities_list = []
        additional_intensities_list = []

        tmt_channels = []
        if self.experiment_type and "TMT" in self.experiment_type.upper():
            tmt_channels = self._get_tmt_channels_from_sdrf()

        for _, row in df.iterrows():
            intensities = []
            additional_intensities = []
            reference_file_name = row.get("reference_file_name", "")

            if tmt_channels:
                tmt_intensities, tmt_additional = self._process_tmt_intensities(
                    row, tmt_channels, reference_file_name
                )
                intensities.extend(tmt_intensities)
                additional_intensities.extend(tmt_additional)

            elif (
                "intensity" in row
                and pd.notna(row["intensity"])
                and row["intensity"] > 0
            ):
                sample_accession = self._get_sample_accession_from_sdrf(
                    reference_file_name
                )
                channel = self._get_channel_from_sdrf(
                    reference_file_name, sample_accession
                )

                intensities.append(
                    {
                        "sample_accession": sample_accession,
                        "channel": channel,
                        "intensity": float(row["intensity"]),
                    }
                )

            intensities_list.append(intensities if intensities else None)
            additional_intensities_list.append(
                additional_intensities if additional_intensities else None
            )

        df["intensities"] = intensities_list
        df["additional_intensities"] = additional_intensities_list

        return df

    def _integrate_sdrf_metadata_feature(self, df: pd.DataFrame) -> pd.DataFrame:
        """Integrate SDRF metadata for Feature"""
        return df

    def _ensure_feature_schema_compliance(self, df: pd.DataFrame) -> pd.DataFrame:
        """Ensure Feature data complies with FEATURE_SCHEMA"""
        schema_fields = [field.name for field in FEATURE_SCHEMA]

        for field in FEATURE_SCHEMA:
            if field.name not in df.columns:
                if field.type == pa.string():
                    df[field.name] = ""
                elif field.type == pa.int32():
                    df[field.name] = 0
                elif field.type == pa.float32():
                    if field.nullable:
                        df[field.name] = None
                    else:
                        df[field.name] = 0.0
                elif str(field.type).startswith("list"):
                    df[field.name] = None
                else:
                    df[field.name] = None

        df = df[[col for col in schema_fields if col in df.columns]].copy()

        if "is_decoy" in df.columns:
            df["is_decoy"] = df["is_decoy"].apply(convert_maxquant_flag).astype("int32")
        if "precursor_charge" in df.columns:
            df["precursor_charge"] = df["precursor_charge"].astype("int32")

        if "scan" in df.columns:
            df["scan"] = df["scan"].astype("string")

        float_fields = [
            "posterior_error_probability",
            "calculated_mz",
            "observed_mz",
            "rt",
            "predicted_rt",
            "ion_mobility",
            "start_ion_mobility",
            "stop_ion_mobility",
            "pg_global_qvalue",
            "rt_start",
            "rt_stop",
        ]

        for field in float_fields:
            if field in df.columns:
                df[field] = pd.to_numeric(df[field], errors="coerce").astype("float32")

        if "rt" in df.columns:
            df["rt"] = (df["rt"] * 60).astype("float32")

        return df

    def _init_protein_group_qvalue_mapping(self, protein_groups_path: str) -> None:
        """Initialize protein group Q-value mapping"""
        try:
            pg_df = pd.read_csv(
                protein_groups_path,
                sep="\t",
                usecols=["Protein IDs", "Majority protein IDs", "Q-value"],
                low_memory=False,
            )

            self._protein_group_qvalue_map = {}

            for _, row in pg_df.iterrows():
                qvalue = row["Q-value"]

                if pd.notna(row["Majority protein IDs"]):
                    protein_ids = str(row["Majority protein IDs"]).split(";")
                    for protein_id in protein_ids:
                        protein_id = protein_id.strip()
                        if protein_id:
                            self._protein_group_qvalue_map[protein_id] = qvalue

        except Exception as e:
            logging.error(f"Failed to initialize protein group Q-value mapping: {e}")
            self._protein_group_qvalue_map = {}

    def _map_protein_group_qvalue(self, df: pd.DataFrame) -> pd.DataFrame:
        """Map protein group Q-values to feature data"""
        if (
            not hasattr(self, "_protein_group_qvalue_map")
            or not self._protein_group_qvalue_map
        ):
            df["pg_global_qvalue"] = None
            return df

        def get_qvalue_for_proteins(protein_accessions):
            """Get minimum Q-value for protein list"""
            if protein_accessions is None:
                return None

            try:
                if pd.isna(protein_accessions):
                    return None
            except (ValueError, TypeError):
                pass

            protein_list = []
            if isinstance(protein_accessions, str):
                if protein_accessions.strip():
                    protein_list = protein_accessions.split(";")
            elif isinstance(protein_accessions, list):
                protein_list = protein_accessions
            elif hasattr(protein_accessions, "__iter__"):
                try:
                    protein_list = list(protein_accessions)
                except:
                    return None
            else:
                return None

            qvalues = []
            for protein_id in protein_list:
                if isinstance(protein_id, str):
                    protein_id = protein_id.strip()
                    if protein_id and protein_id in self._protein_group_qvalue_map:
                        qvalues.append(self._protein_group_qvalue_map[protein_id])
            return min(qvalues) if qvalues else None

        if "pg_accessions" in df.columns:
            df["pg_global_qvalue"] = df["pg_accessions"].apply(get_qvalue_for_proteins)
        else:
            df["pg_global_qvalue"] = None

        return df

    # ============================================================================
    # Protein Group Processing
    # ============================================================================

    def process_pg_file(
        self,
        protein_groups_path: str,
        output_path: str,
        sdrf_path: str = None,
        chunksize: int = 100000,
    ) -> None:
        """Process Protein Group data from proteinGroups.txt to PG parquet format"""
        if sdrf_path:
            self._init_sdrf(sdrf_path)

        batch_writer = ParquetBatchWriter(output_path, PG_SCHEMA)

        try:
            for df_chunk in pd.read_csv(
                protein_groups_path, sep="\t", chunksize=chunksize, low_memory=False
            ):

                basic_cols = [
                    col for col in MAXQUANT_PG_USECOLS if col in df_chunk.columns
                ]
                intensity_cols = [
                    col
                    for col in df_chunk.columns
                    if col.startswith("Intensity ")
                    or col.startswith("LFQ intensity ")
                    or col.startswith("iBAQ ")
                ]
                available_cols = list(set(basic_cols + intensity_cols))
                df_chunk = df_chunk[available_cols]

                if df_chunk.empty:
                    continue

                df_chunk = self._apply_pg_mapping(df_chunk)
                df_chunk = self._process_pg_basic_fields(df_chunk)
                df_chunk = self._process_pg_intensities(df_chunk)
                df_chunk = self._calculate_pg_statistics(df_chunk)

                if self.sdrf_handler:
                    df_chunk = self._integrate_sdrf_metadata_pg(df_chunk)

                df_chunk = self._ensure_pg_schema_compliance(df_chunk)

                table = pa.Table.from_pandas(
                    df_chunk, schema=PG_SCHEMA, preserve_index=False
                )
                batch_writer.write_batch(table.to_pylist())

        finally:
            batch_writer.close()
            pass

    def _apply_pg_mapping(self, df: pd.DataFrame) -> pd.DataFrame:
        """Apply MAXQUANT_PG_MAP mapping"""
        available_mapping = {
            k: v for k, v in self.pg_mapping.items() if k in df.columns
        }
        df.rename(columns=available_mapping, inplace=True)
        return df

    def _process_pg_basic_fields(self, df: pd.DataFrame) -> pd.DataFrame:
        """Process PG basic fields"""
        list_columns = ["pg_accessions", "pg_names", "gg_accessions"]
        for col in list_columns:
            if col in df.columns:
                df[col] = df[col].fillna("").astype(str).str.split(";")
            else:
                df[col] = [[] for _ in range(len(df))]

        if "is_decoy" in df.columns:
            df["is_decoy"] = df["is_decoy"].apply(convert_maxquant_flag).astype("int32")
        else:
            df["is_decoy"] = 0

        if "contaminant" in df.columns:
            df["contaminant"] = (
                df["contaminant"].apply(convert_maxquant_flag).astype("int32")
            )
        else:
            df["contaminant"] = 0

        if self.sdrf_handler:
            df["reference_file_name"] = self._get_representative_reference_file()
        else:
            df["reference_file_name"] = "proteinGroups.txt"

        return df

    def _get_representative_reference_file(self) -> str:
        """Get representative reference file name from transformed SDRF data"""
        try:
            if (
                hasattr(self, "_sdrf_transformed")
                and self._sdrf_transformed is not None
                and not self._sdrf_transformed.empty
            ):
                # Check if reference_file_name is a column
                if "reference_file_name" in self._sdrf_transformed.columns:
                    return self._sdrf_transformed["reference_file_name"].iloc[0]
                # Check if reference_file_name is the index
                elif self._sdrf_transformed.index.name == "reference_file_name":
                    return self._sdrf_transformed.index[0]
        except:
            pass
        return "proteinGroups"

    def _process_pg_intensities(self, df: pd.DataFrame) -> pd.DataFrame:
        """Process PG intensity data with comprehensive intensity extraction"""
        intensities_list = []
        additional_intensities_list = []

        intensity_cols = self._get_sample_specific_intensity_cols(df)
        lfq_cols = [col for col in df.columns if col.startswith("LFQ intensity ")]
        ibaq_cols = self._get_sample_specific_ibaq_cols(df)

        general_intensity_col = "Intensity" if "Intensity" in df.columns else None
        general_ibaq_col = "iBAQ" if "iBAQ" in df.columns else None

        for _, row in df.iterrows():
            intensities, additional_intensities = self._create_pg_intensity_struct(
                row=row,
                intensity_cols=intensity_cols,
                lfq_cols=lfq_cols,
                ibaq_cols=ibaq_cols,
                general_intensity_col=general_intensity_col,
                general_ibaq_col=general_ibaq_col,
            )

            intensities_list.append(intensities if intensities else [])
            additional_intensities_list.append(
                additional_intensities if additional_intensities else []
            )

        df["intensities"] = intensities_list
        df["additional_intensities"] = additional_intensities_list

        return df

    def _get_sample_specific_intensity_cols(self, df: pd.DataFrame) -> list:
        """Get sample-specific intensity columns (excluding general 'Intensity')"""
        return [
            col
            for col in df.columns
            if col.startswith("Intensity ") and col != "Intensity"
        ]

    def _get_sample_specific_ibaq_cols(self, df: pd.DataFrame) -> list:
        """Get sample-specific iBAQ columns (excluding general 'iBAQ')"""
        return [col for col in df.columns if col.startswith("iBAQ ") and col != "iBAQ"]

    def _create_pg_intensity_struct(
        self,
        row,
        intensity_cols,
        lfq_cols,
        ibaq_cols,
        general_intensity_col,
        general_ibaq_col,
    ) -> tuple:
        """Create intensity structure for protein group"""
        intensities = []
        additional_intensities = []

        if general_intensity_col and general_intensity_col in row.index:
            if pd.notna(row[general_intensity_col]):
                sample_accession = "Unknown"
                channel = "label free sample" if self.experiment_type == "LFQ" else None
                intensities.append(
                    {
                        "sample_accession": sample_accession,
                        "channel": channel,
                        "intensity": float(row[general_intensity_col]),
                    }
                )

        intensities.extend(
            self._process_sample_specific_pg_intensities(row, intensity_cols)
        )

        additional_intensities.extend(self._process_lfq_pg_intensities(row, lfq_cols))

        additional_intensities.extend(
            self._process_ibaq_pg_intensities(row, ibaq_cols, general_ibaq_col)
        )

        return intensities, additional_intensities

    def _process_sample_specific_pg_intensities(self, row, intensity_cols) -> list:
        """Process sample-specific intensity columns with SDRF mapping"""
        sample_intensity_map = {}

        for col in intensity_cols:
            if col in row.index and pd.notna(row[col]) and row[col] > 0:
                sample_accession, channel = (
                    self._extract_sample_info_from_intensity_col(col, "Intensity ")
                )

                if sample_accession is None:
                    maxquant_sample = col.replace("Intensity ", "")
                    sample_accession = maxquant_sample
                    channel = (
                        "label free sample" if self.experiment_type == "LFQ" else None
                    )

                sample_key = (sample_accession, channel)
                if sample_key not in sample_intensity_map:
                    sample_intensity_map[sample_key] = 0
                sample_intensity_map[sample_key] += float(row[col])

        intensities = []
        for (
            sample_accession,
            channel,
        ), total_intensity in sample_intensity_map.items():
            intensities.append(
                {
                    "sample_accession": sample_accession,
                    "channel": channel,
                    "intensity": total_intensity,
                }
            )
        return intensities

    def _process_lfq_pg_intensities(self, row, lfq_cols) -> list:
        """Process LFQ intensity columns for PG data"""
        additional_intensities = []
        for col in lfq_cols:
            if col in row.index and pd.notna(row[col]):
                sample_accession, channel = (
                    self._extract_sample_info_from_intensity_col(col, "LFQ intensity ")
                )

                if sample_accession is None:
                    maxquant_sample = col.replace("LFQ intensity ", "")
                    sample_accession = maxquant_sample
                    channel = (
                        "label free sample" if self.experiment_type == "LFQ" else None
                    )

                additional_intensities.append(
                    {
                        "sample_accession": sample_accession,
                        "channel": channel,
                        "intensities": [
                            {
                                "intensity_name": col,
                                "intensity_value": float(row[col]),
                            }
                        ],
                    }
                )
        return additional_intensities

    def _process_ibaq_pg_intensities(self, row, ibaq_cols, general_ibaq_col) -> list:
        """Process iBAQ intensity columns for PG data"""
        additional_intensities = []

        if general_ibaq_col and general_ibaq_col in row.index:
            if pd.notna(row[general_ibaq_col]):
                sample_accession = "Unknown"
                channel = "label free sample" if self.experiment_type == "LFQ" else None
                additional_intensities.append(
                    {
                        "sample_accession": sample_accession,
                        "channel": channel,
                        "intensities": [
                            {
                                "intensity_name": general_ibaq_col,  # Use original column name
                                "intensity_value": float(row[general_ibaq_col]),
                            }
                        ],
                    }
                )

        for col in ibaq_cols:
            if col in row.index and pd.notna(row[col]):
                sample_accession, channel = (
                    self._extract_sample_info_from_intensity_col(col, "iBAQ ")
                )

                if sample_accession is None:
                    maxquant_sample = col.replace("iBAQ ", "")
                    sample_accession = maxquant_sample
                    channel = (
                        "label free sample" if self.experiment_type == "LFQ" else None
                    )

                additional_intensities.append(
                    {
                        "sample_accession": sample_accession,
                        "channel": channel,
                        "intensities": [
                            {
                                "intensity_name": col,
                                "intensity_value": float(row[col]),
                            }
                        ],
                    }
                )
        return additional_intensities

    def _extract_sample_info_from_intensity_col(
        self, col_name: str, col_prefix: str
    ) -> tuple:
        """Extract sample information from intensity column using SDRF mapping"""
        if not hasattr(self, "_sdrf_transformed") or self._sdrf_transformed is None:
            return None, None

        maxquant_sample = col_name.replace(col_prefix, "").strip()
        if not maxquant_sample:
            return None, None

        df = (
            self._sdrf_transformed.reset_index()
            if self._sdrf_transformed.index.name == "reference_file_name"
            else self._sdrf_transformed
        )

        sample_accession = None
        channel = None
        if "reference_file_name" in df.columns:
            matching_candidates = []

            for _, row in df.iterrows():
                ref_file = str(row.get("reference_file_name", ""))
                if ref_file and maxquant_sample in ref_file:
                    match_quality = len(ref_file)
                    matching_candidates.append((match_quality, row))

            if matching_candidates:
                best_match = min(matching_candidates)[1]
                sample_accession = best_match.get("sample_accession", None)
                channel = best_match.get("channel", None)

        if sample_accession is None and "reference_file_name" in df.columns:
            cleaned_sample = maxquant_sample.rstrip("_").rstrip("-")
            fuzzy_candidates = []

            for _, row in df.iterrows():
                ref_file = str(row.get("reference_file_name", ""))
                if ref_file and (
                    cleaned_sample in ref_file
                    or ref_file.replace("_", "").replace("-", "")
                    in cleaned_sample.replace("_", "").replace("-", "")
                ):
                    match_quality = len(ref_file)
                    fuzzy_candidates.append((match_quality, row))

            if fuzzy_candidates:
                best_fuzzy_match = min(fuzzy_candidates)[1]
                sample_accession = best_fuzzy_match.get("sample_accession", None)
                channel = best_fuzzy_match.get("channel", None)

        if sample_accession is None:
            sample_accession = maxquant_sample

        if channel is None or pd.isna(channel):
            if self.experiment_type == "LFQ":
                channel = "label free sample"
            else:
                channel = None

        return sample_accession, channel

    def _get_first_available_reference_file(self) -> str:
        """Get the first available reference file from SDRF mapping"""
        if self._sample_map and len(self._sample_map) > 0:
            first_file_key = list(self._sample_map.keys())[0]
            return f"{first_file_key}.mzML"
        return "Unknown.mzML"

    def _get_peptide_count_from_row(self, row) -> int:
        """Extract peptide count from MaxQuant data"""
        peptide_columns = [
            "peptide_count_total",
            "peptide_count_razor_unique",
            "peptide_count_unique",
        ]

        for col in peptide_columns:
            if col in row.index and pd.notna(row[col]) and row[col] > 0:
                return int(row[col])

        return 1

    def _calculate_pg_statistics(self, df: pd.DataFrame) -> pd.DataFrame:
        """Calculate PG statistics"""
        df = self._process_pg_peptides(df)
        df = self._process_pg_anchor_protein(df)
        df = self._process_pg_counts(df)
        df = self._process_pg_additional_scores(df)
        return df

    def _process_pg_peptides(self, df: pd.DataFrame) -> pd.DataFrame:
        """Process peptides field for PG data"""
        peptides_list = []
        for _, row in df.iterrows():
            peptides = []
            if "pg_accessions" in row and row["pg_accessions"]:
                total_peptides = self._get_peptide_count_from_row(row)
                protein_count = len(row["pg_accessions"])

                if protein_count == 1:
                    protein = row["pg_accessions"][0]
                    if protein:
                        peptides.append(
                            {"protein_name": protein, "peptide_count": total_peptides}
                        )
                else:
                    peptides = self._distribute_peptides_among_proteins(
                        row["pg_accessions"], total_peptides
                    )

            peptides_list.append(peptides)

        df["peptides"] = peptides_list
        return df

    def _distribute_peptides_among_proteins(
        self, proteins: list, total_peptides: int
    ) -> list:
        """Distribute peptides among multiple proteins in a group"""
        proteins = [p for p in proteins if p]
        if not proteins:
            return []

        peptides = []
        protein_count = len(proteins)

        main_peptide_count = max(1, int(total_peptides * 0.6))
        remaining_peptides = total_peptides - main_peptide_count
        other_count = protein_count - 1

        for i, protein in enumerate(proteins):
            if i == 0:
                peptide_count = main_peptide_count
            else:
                if other_count > 0 and remaining_peptides > 0:
                    base_count = remaining_peptides // other_count
                    if i <= remaining_peptides % other_count:
                        peptide_count = base_count + 1
                    else:
                        peptide_count = base_count
                else:
                    peptide_count = 0

            if peptide_count > 0:
                peptides.append(
                    {
                        "protein_name": protein,
                        "peptide_count": peptide_count,
                    }
                )

        return peptides

    def _process_pg_anchor_protein(self, df: pd.DataFrame) -> pd.DataFrame:
        """Process anchor_protein field for PG data"""
        if len(df) > 0:
            df["anchor_protein"] = df["pg_accessions"].str[0]
        else:
            df["anchor_protein"] = None
        return df

    def _process_pg_counts(self, df: pd.DataFrame) -> pd.DataFrame:
        """Process peptide_counts and feature_counts fields for PG data"""
        peptide_counts_list = []
        feature_counts_list = []

        for _, row in df.iterrows():
            unique_sequences = int(row.get("peptide_count_unique", 1))
            total_sequences = int(row.get("peptide_count_total", 1))

            peptide_counts_list.append(
                {
                    "unique_sequences": unique_sequences,
                    "total_sequences": total_sequences,
                }
            )

            # Use same values for feature counts
            feature_counts_list.append(
                {
                    "unique_features": unique_sequences,
                    "total_features": total_sequences,
                }
            )

        df["peptide_counts"] = peptide_counts_list
        df["feature_counts"] = feature_counts_list
        return df

    def _process_pg_additional_scores(self, df: pd.DataFrame) -> pd.DataFrame:
        """Extract additional scores for PG data"""
        score_fields = [
            "andromeda_score",
            "sequence_coverage",
            "molecular_weight",
            "msms_count",
            "number_of_proteins",
            "peptide_count_razor_unique",
        ]

        scores_list = []
        for _, row in df.iterrows():
            scores = []

            for field_name in score_fields:
                if field_name in row and pd.notna(row[field_name]):
                    scores.append(
                        {
                            "score_name": field_name,
                            "score_value": float(row[field_name]),
                        }
                    )

            scores_list.append(scores)

        df["additional_scores"] = scores_list
        return df

    def _integrate_sdrf_metadata_pg(self, df: pd.DataFrame) -> pd.DataFrame:
        """Integrate SDRF metadata for PG"""
        return df

    def _ensure_pg_schema_compliance(self, df: pd.DataFrame) -> pd.DataFrame:
        """Ensure PG data complies with PG_SCHEMA"""
        schema_fields = [field.name for field in PG_SCHEMA]

        for field in PG_SCHEMA:
            if field.name not in df.columns:
                if field.type == pa.string():
                    df[field.name] = None
                elif field.type == pa.int32():
                    df[field.name] = 0
                elif field.type == pa.float32():
                    df[field.name] = 1.0
                elif str(field.type).startswith("list"):
                    df[field.name] = [[] for _ in range(len(df))]
                elif str(field.type).startswith("struct"):
                    if field.name == "peptide_counts":
                        df[field.name] = [
                            {"unique_sequences": 0, "total_sequences": 0}
                            for _ in range(len(df))
                        ]
                    elif field.name == "feature_counts":
                        df[field.name] = [
                            {"unique_features": 0, "total_features": 0}
                            for _ in range(len(df))
                        ]
                else:
                    df[field.name] = None

        return df[[col for col in schema_fields if col in df.columns]]

    # ============================================================================
    # Backward Compatibility Interfaces
    # ============================================================================

    def write_psm_to_file(
        self, msms_path: str, output_path: str, chunksize: int = 1000000
    ) -> None:
        """Backward compatible PSM writing interface"""
        self.process_psm_file(msms_path, output_path, chunksize=chunksize)

    def write_feature_to_file(
        self,
        evidence_path: str,
        sdrf_path: str,
        output_path: str,
        protein_file: str = None,
        chunksize: int = 1000000,
    ) -> None:
        """Backward compatible Feature writing interface"""
        self.process_feature_file(
            evidence_path, output_path, sdrf_path, protein_file, chunksize
        )

    def write_protein_groups_to_file(
        self,
        protein_groups_path: str,
        sdrf_path: str,
        output_path: str,
        chunksize: int = 100000,
    ) -> None:
        """Backward compatible Protein Groups writing interface"""
        self.process_pg_file(protein_groups_path, output_path, sdrf_path, chunksize)

    # ============================================================================
    # Additional Backward Compatibility Methods for Tests
    # ============================================================================

    def read_msms(self, msms_path: str) -> pd.DataFrame:
        """Read and process msms.txt file, returns DataFrame"""
        return pd.read_csv(msms_path, sep="\t", low_memory=False)

    def read_evidence(self, evidence_path: str) -> pd.DataFrame:
        """Read and process evidence.txt file, returns DataFrame"""
        return pd.read_csv(evidence_path, sep="\t", low_memory=False)

    def read_protein_groups(self, protein_groups_path: str) -> pd.DataFrame:
        """Read and process proteinGroups.txt file, returns DataFrame"""
        return pd.read_csv(protein_groups_path, sep="\t", low_memory=False)

    def process_msms_to_psm_table(self, df: pd.DataFrame) -> pa.Table:
        """Process MSMS DataFrame to PSM table"""
        if df.empty:
            raise ValueError("Input DataFrame is empty")
        df_processed = self._apply_psm_mapping(df.copy())
        self._calculate_theoretical_mz_batch(df_processed)
        df_processed = self._process_psm_modifications(df_processed)
        df_processed = self._process_psm_scores(df_processed)
        df_processed = self._ensure_psm_schema_compliance(df_processed)
        return pa.Table.from_pandas(
            df_processed, schema=PSM_SCHEMA, preserve_index=False
        )

    def process_evidence_to_feature_table(self, df: pd.DataFrame) -> pa.Table:
        """Process Evidence DataFrame to Feature table"""
        df_processed = self._apply_feature_mapping(df.copy())
        self._calculate_theoretical_mz_batch(df_processed)
        df_processed = self._process_feature_modifications(df_processed)
        df_processed = self._process_feature_protein_groups(df_processed)
        df_processed = self._ensure_feature_schema_compliance(df_processed)
        return pa.Table.from_pandas(
            df_processed, schema=FEATURE_SCHEMA, preserve_index=False
        )

    def process_protein_groups_to_pg_table(
        self, df: pd.DataFrame, sdrf_path: str = None
    ) -> pa.Table:
        """Process Protein Groups DataFrame to PG table"""
        if df.empty:
            raise ValueError("Input DataFrame is empty")
        if sdrf_path:
            self._init_sdrf(sdrf_path)
        df_processed = self._apply_pg_mapping(df.copy())
        df_processed = self._process_pg_basic_fields(df_processed)
        df_processed = self._calculate_pg_statistics(df_processed)
        df_processed = self._ensure_pg_schema_compliance(df_processed)
        return pa.Table.from_pandas(
            df_processed, schema=PG_SCHEMA, preserve_index=False
        )

    def iter_batch(self, file_path: str, chunksize: int = 10000):
        """Iterate over file in batches"""
        for chunk in pd.read_csv(
            file_path, sep="\t", chunksize=chunksize, low_memory=False
        ):
            yield chunk


# ============================================================================
# Standalone Conversion Functions
# ============================================================================


def process_evidence_to_feature_table(df: pd.DataFrame) -> pa.Table:
    """Backward compatible function"""
    processor = MaxQuant()
    df_processed = processor._apply_feature_mapping(df)
    processor._calculate_theoretical_mz_batch(df_processed)
    df_processed = processor._ensure_feature_schema_compliance(df_processed)
    return pa.Table.from_pandas(
        df_processed, schema=FEATURE_SCHEMA, preserve_index=False
    )


def process_protein_groups_to_pg_table(
    df: pd.DataFrame, sdrf_path: str = None
) -> pa.Table:
    """Backward compatible function"""
    processor = MaxQuant()
    if sdrf_path:
        processor._init_sdrf(sdrf_path)
    df_processed = processor._apply_pg_mapping(df)
    df_processed = processor._process_pg_basic_fields(df_processed)
    df_processed = processor._calculate_pg_statistics(df_processed)
    df_processed = processor._ensure_pg_schema_compliance(df_processed)
    return pa.Table.from_pandas(df_processed, schema=PG_SCHEMA, preserve_index=False)


def process_msms_to_psm_table(df: pd.DataFrame) -> pa.Table:
    """Backward compatible function"""
    processor = MaxQuant()
    df_processed = processor._apply_psm_mapping(df)
    processor._calculate_theoretical_mz_batch(df_processed)
    df_processed = processor._process_psm_modifications(df_processed)
    df_processed = processor._process_psm_scores(df_processed)
    df_processed = processor._ensure_psm_schema_compliance(df_processed)
    return pa.Table.from_pandas(df_processed, schema=PSM_SCHEMA, preserve_index=False)


# ============================================================================
# Test Compatibility Functions
# ============================================================================


def read_msms(msms_path: str) -> pd.DataFrame:
    """Read and process msms.txt file, returns DataFrame"""
    return pd.read_csv(msms_path, sep="\t", low_memory=False)


def read_evidence(evidence_path: str) -> pd.DataFrame:
    """Read and process evidence.txt file, returns DataFrame"""
    return pd.read_csv(evidence_path, sep="\t", low_memory=False)


def read_protein_groups(protein_groups_path: str) -> pd.DataFrame:
    """Read and process proteinGroups.txt file, returns DataFrame"""
    return pd.read_csv(protein_groups_path, sep="\t", low_memory=False)
