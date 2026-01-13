import logging
import re
from pathlib import Path
import pandas as pd
import pyarrow as pa

from qpx.core.common import FEATURE_SCHEMA

# MsstatsIN functionality now available in MzTabIndexer
from qpx.core.quantms.mztab import MzTabIndexer
from qpx.operate.tools import get_protein_accession
from qpx.utils.file_utils import (
    extract_protein_list,
    ParquetBatchWriter,
)

# Pre-compiled regular expression patterns for better performance
MOD_PATTERN = re.compile(r"_mod\[")
BRACKET_MOD_PATTERN = re.compile(r"\[.*?\]")
PAREN_MOD_PATTERN = re.compile(r"\(.*?\)")
NON_AMINO_ACID_PATTERN = re.compile(r"[^A-Z]")


class Feature:
    """Feature processor using composition pattern.

    This class processes feature data from an MzTabIndexer instance, providing
    feature-specific functionality without inheriting from the indexer.
    """

    def __init__(self, mztab_indexer):
        """Initialize Feature processor with file paths.

        Args:
            mztab_file_path: Path to mzTab file
            sdrf_path: Path to SDRF file
            msstats_in_path: Path to MSstats input file
        """
        self.logger = logging.getLogger(f"{__name__}.{self.__class__.__name__}")

        self._indexer: MzTabIndexer = mztab_indexer
        self._ms_runs = self._extract_ms_runs()
        self._protein_global_qvalue_map = self._get_protein_map()
        self._score_names = self._get_score_names()
        self._mods_map = self._get_mods_map()

    def cleanup(self):
        """Clean up resources."""
        if hasattr(self._indexer, "close"):
            self._indexer.close()

    def extract_psm_msg(self, protein_str=None):
        """Extract PSM messages for merging with features.

        Returns:
            tuple: (map_dict, pep_dict) containing PSM data mappings
        """
        # Get PSMs from the indexer
        try:
            self.logger.info("Get PSMs from the indexer.")
            psms_df = self._indexer.get_psms()
            if psms_df.empty:
                return {}, {}

            self.logger.info("Get metadata from the indexer.")
            metadata_df = self._indexer.get_metadata()

            self.logger.info(
                "Calling get_ms_run_reference_map to retrieve MS run reference map."
            )
            ms_reference_map = self.get_ms_run_reference_map(psms_df, metadata_df)
            psms_df["reference_file_name"] = psms_df["spectra_ref_file"].map(
                ms_reference_map
            )

            # Apply protein filter if specified
            if protein_str:
                psms_df = psms_df[
                    psms_df["accession"].str.contains(protein_str, na=False)
                ]

            map_dict = {}
            pep_dict = {}

            # Create mapping dictionaries
            self.logger.info("Start creating mapping dictionaries.")
            for _, row in psms_df.iterrows():
                # Map key: (reference_file_name, peptidoform, precursor_charge)
                reference_file_name = row.get("reference_file_name", "")
                peptidoform = row.get(
                    "opt_global_cv_MS:1000889_peptidoform_sequence", ""
                )
                charge = str(row.get("charge", "0"))

                map_key = (reference_file_name, peptidoform, charge)

                # Store PSM features for merging
                map_dict[map_key] = [
                    row.get("opt_global_Posterior_Error_Probability_score"),
                    row.get("calc_mass_to_charge"),
                    row.get("exp_mass_to_charge"),
                    row.get("accession"),
                    (
                        1
                        if row.get("opt_global_cv_MS:1002217_decoy_peptide") == "1"
                        else 0
                    ),
                    [{}],  # additional_scores placeholder
                    [{}],  # cv_params placeholder
                ]

                # Store peptide-charge mapping for best scan
                pep_key = (peptidoform, charge)
                if pep_key not in pep_dict:
                    scan_info = row.get("spectra_ref", "").split(":")
                    pep_dict[pep_key] = [
                        peptidoform,
                        scan_info[0] if len(scan_info) > 0 else "",
                        scan_info[1] if len(scan_info) > 1 else "",
                    ]

            self.logger.info("Finished creating mapping dictionaries.")

            return map_dict, pep_dict

        except Exception as e:
            self.logger.warning(f"Could not extract PSM messages: {e}")
            return {}, {}

    def get_ms_run_reference_map(self, psms, metadata):
        ms_run_list = psms["spectra_ref_file"].drop_duplicates().to_list()
        mapping = metadata.set_index("key")["value"].to_dict()

        reference_map = dict()
        for ms_run in ms_run_list:
            key = ms_run + "-location"
            if key in mapping:
                reference_map[ms_run] = Path(mapping[key]).stem

        return reference_map

    def _parse_peptidoform(self, peptidoform):
        """
        Parse peptidoform string to extract sequence and modifications.

        Args:
            peptidoform (str): The peptidoform string to parse

        Returns:
            tuple: (sequence, mod_list) where mod_list contains dicts with mod_name and position
        """
        sequence = ""
        mod_list = []
        state = "SEQ"
        mod_name_buffer = ""
        pos = 1

        for char in peptidoform:
            if state == "SEQ":
                if char == "(":
                    state = "MOD"
                    mod_name_buffer = ""
                elif char == ".":
                    continue
                else:
                    sequence += char
                    pos += 1
            elif state == "MOD":
                if char == ")":
                    mod_pos = 0 if pos == 1 else pos - 1
                    mod_list.append({"mod_name": mod_name_buffer, "position": mod_pos})
                    state = "SEQ"
                else:
                    mod_name_buffer += char

        return sequence, mod_list

    def _match_modification_position(self, mod, rule, sequence):
        """
        Match a modification with position rules and return the JSON position if valid.

        Args:
            mod (dict): Modification with mod_name and position
            rule (dict): Rule with site and position information
            sequence (str): Peptide sequence

        Returns:
            str or None: JSON position string if match is valid, None otherwise
        """
        site = rule["site"]
        position_rule = rule["position"]
        aa_pos = mod["position"]

        if aa_pos == 0 and position_rule.lower() in ["any n-term", "n-term.0"]:
            return "N-term.0"

        elif site != "X" and aa_pos > 0:
            aa = sequence[aa_pos - 1]
            if aa == site:
                return f"{aa}.{aa_pos}"

        return None

    def _build_modification_entry(self, mod_name, accession, mod_positions):
        """
        Build a modification entry for the results.

        Args:
            mod_name (str): Name of the modification
            accession (str): Modification accession
            mod_positions (list): List of position dictionaries

        Returns:
            dict: Modification entry for results
        """
        return {
            "name": mod_name,
            "accession": accession,
            "positions": mod_positions,
        }

    def _create_position_entry(self, json_pos, mod_name, select_mods):
        """
        Create a position entry with scores.

        Args:
            json_pos (str): Position in JSON format
            mod_name (str): Name of the modification
            select_mods (list): List of selected modifications

        Returns:
            dict: Position entry with scores
        """
        return {
            "position": json_pos,
            "scores": (
                # Note: localization_probability scores are set to None as no data source
                # is currently available. Future implementation should integrate with
                # modification localization scoring algorithms or external databases.
                [
                    {
                        "score_name": "localization_probability",
                        "score_value": None,
                    }
                ]
                if mod_name in select_mods
                else None
            ),
        }

    def _process_modifications_from_dict(
        self, mod_list, sequence, modifications_dict, select_mods
    ):
        """
        Process modifications from the modifications dictionary.

        Args:
            mod_list (list): List of modifications from peptidoform
            sequence (str): Peptide sequence
            modifications_dict (dict): Dictionary of modification information
            select_mods (list): List of selected modifications

        Returns:
            list: List of processed modification results
        """
        mod_results = []

        for mod_name, mod_info in modifications_dict.items():
            accession = mod_info["accession"]
            site_positions = mod_info["site_positions"]
            mod_positions = []

            for rule in site_positions:
                for mod in mod_list:
                    if mod["mod_name"] != mod_name:
                        continue

                    json_pos = self._match_modification_position(mod, rule, sequence)
                    if json_pos:
                        position_entry = self._create_position_entry(
                            json_pos, mod_name, select_mods
                        )
                        mod_positions.append(position_entry)

            if mod_positions:
                mod_entry = self._build_modification_entry(
                    mod_name, accession, mod_positions
                )
                mod_results.append(mod_entry)

        return mod_results

    def generate_modifications_details(self, peptidoform, modifications_dict):
        """
        Generate detailed modification information from peptidoform and modifications dictionary.

        Args:
            peptidoform (str): The peptidoform string containing sequence and modifications
            modifications_dict (dict): Dictionary containing modification information

        Returns:
            list or None: List of modification details or None if no modifications found
        """
        select_mods = list(modifications_dict.keys())

        # Parse peptidoform to extract sequence and modifications
        sequence, mod_list = self._parse_peptidoform(peptidoform)

        # Process modifications from dictionary
        mod_results = self._process_modifications_from_dict(
            mod_list, sequence, modifications_dict, select_mods
        )

        return mod_results if mod_results else None

    def _extract_ms_runs(self) -> dict:
        """Extract MS runs from metadata."""
        try:
            metadata_df = self._indexer.get_metadata()
            ms_runs = {}
            for _, row in metadata_df.iterrows():
                if row["key"].startswith("ms_run[") and row["key"].endswith(
                    "]-location"
                ):
                    ms_run_id = row["key"].split("[")[1].split("]")[0]
                    ms_runs[ms_run_id] = row["value"]
            return ms_runs
        except Exception as e:
            self.logger.warning(f"Could not extract MS runs: {e}")
            return {}

    def _get_protein_map(self) -> dict:
        """Get protein global q-value map."""
        try:
            proteins_df = self._indexer.get_proteins()
            protein_map = {}
            # Check if required columns exist in proteins table
            # Note: opt_global_q-value presence varies by data source
            if (
                "accession" in proteins_df.columns
                and "opt_global_q-value" in proteins_df.columns
            ):
                for _, row in proteins_df.iterrows():
                    if pd.notna(row["accession"]) and pd.notna(
                        row["opt_global_q-value"]
                        # row["opt_global_qvalue"]
                    ):
                        protein_map[row["accession"]] = row["opt_global_q-value"]
                        # protein_map[row["accession"]] = row["opt_global_qvalue"]
            return protein_map
        except Exception as e:
            self.logger.warning(f"Could not get protein map: {e}")
            return {}

    def _get_score_names(self) -> dict:
        """Get score names from metadata."""
        try:
            metadata_df = self._indexer.get_metadata()
            score_names = {}
            for _, row in metadata_df.iterrows():
                if "search_engine_score" in row["key"] and "[1]" in row["key"]:
                    score_names[row["value"]] = row["key"]
            return score_names
        except Exception as e:
            self.logger.warning(f"Could not get score names: {e}")
            return {}

    def _get_mods_map(self) -> dict:
        """Get modifications map."""
        try:
            metadata_df = self._indexer.get_metadata()

            mod_df = metadata_df[
                metadata_df["key"].str.contains(MOD_PATTERN, regex=True, na=False)
            ]
            mod_dict = mod_df.set_index("key")["value"].to_dict()

            mod_groups = list(set(i.split("-", 1)[0] for i in mod_df["key"].to_list()))

            modifications = dict()
            for mod_group in mod_groups:

                suffixes = ["", "-site", "-position"]

                accession = None
                name = None
                site = None
                position = None

                for suffix in suffixes:
                    key = mod_group + suffix
                    if key in mod_dict:

                        if suffix == "":
                            row_values = (
                                mod_dict[key]
                                .replace("[", "")
                                .replace("]", "")
                                .split(",")
                            )
                            if len(row_values) >= 3:
                                accession = row_values[1].strip()
                                name = row_values[2].strip()

                        if suffix == "-site":
                            site = mod_dict[key]

                        if suffix == "-position":
                            position = mod_dict[key]

                if name is not None:
                    if name not in modifications:
                        modifications[name] = {
                            "accession": accession,
                            "site_positions": [],
                        }
                    modifications[name]["site_positions"].append(
                        {"site": site, "position": position}
                    )

            return modifications
        except Exception as e:
            self.logger.warning(f"Could not get mods map: {e}")
            return {}

    def _create_file_metadata(self):
        """Create file metadata structure according to feature.avsc schema"""
        import uuid
        from datetime import datetime
        from qpx import __version__

        return {
            "qpx_version": __version__,
            "creator": "QPX",
            "file_type": "feature_file",
            "creation_date": datetime.now().isoformat(),
            "uuid": str(uuid.uuid4()),
            "scan_format": "scan",  # Default scan format
            "software_provider": "QPX",
        }

    def transform_msstats_in(self, file_num=10, protein_str=None):

        # Determine experiment type (LFQ vs TMT)
        experiment_type = self._indexer.get_msstats_experiment_type()

        # Use the enhanced MSstats analysis methods from MzTabIndexer
        for batch in self._indexer.iter_msstats_files(file_batch_size=file_num):
            if batch is not None and not batch.empty:
                # Apply protein filter if specified
                if protein_str:
                    batch = batch[
                        batch["ProteinName"].str.contains(protein_str, na=False)
                    ]

                if not batch.empty:

                    # Unique peptide indicator from PSM Table
                    unique_peptide_df = self._indexer.get_unique_from_psm_table()

                    batch = pd.merge(
                        batch,
                        unique_peptide_df,
                        on=["pg_accessions", "peptidoform"],
                        how="left",
                    )

                    # Aggregate data to create feature-level records with intensities array
                    aggregated_features = self._aggregate_msstats_to_features(
                        batch, experiment_type
                    )
                    if not aggregated_features.empty:
                        yield aggregated_features

    def _aggregate_msstats_to_features(self, msstats_batch, experiment_type):
        """
        Aggregate MSstats data into feature-level records with proper intensities structure.
        Groups by (PeptideSequence, ProteinName, Charge, reference_file_name) and creates intensities array.
        """

        # Group by feature identifier (peptidoform + charge + reference file + protein)
        grouping_cols = ["peptidoform", "pg_accessions", "reference_file_name"]

        # Add charge column if available, otherwise use default
        if "charge" in msstats_batch.columns:
            grouping_cols.append("charge")
        else:
            # Add a default charge if not available
            msstats_batch["charge"] = 3
            grouping_cols.append("charge")

        features_list = []

        for group_key, group_data in msstats_batch.groupby(grouping_cols):
            if len(grouping_cols) == 4:
                peptidoform, protein_name, reference_file_name, precursor_charge = (
                    group_key
                )
            else:
                peptidoform, protein_name, reference_file_name = group_key
                precursor_charge = 3  # default charge

            # Create intensities array from group data
            intensities = []
            for _, row in group_data.iterrows():
                intensity_entry = {
                    "sample_accession": row.get("sample_accession", ""),
                    "channel": str(
                        row["channel"]
                        if "channel" in row and row["channel"] is not None
                        else ("LFQ" if experiment_type == "LFQ" else "Unknown")
                    ),
                    "intensity": float(row.get("intensity", 0.0)),
                }
                intensities.append(intensity_entry)

            # Extract other metadata from the first row of the group
            first_row = group_data.iloc[0]

            # Create feature record
            feature_record = {
                "peptidoform": peptidoform,
                "precursor_charge": int(precursor_charge),
                "reference_file_name": reference_file_name,
                "intensities": intensities,
                "pg_accessions": [protein_name] if protein_name else [],
                "anchor_protein": protein_name or "",
                "rt": first_row.get("rt", None),
                "unique": float(first_row.get("unique", 1)),
                # Will add more fields in subsequent processing steps
            }

            features_list.append(feature_record)

        # Convert to DataFrame
        if features_list:
            return pd.DataFrame(features_list)
        else:
            return pd.DataFrame()

    @staticmethod
    def merge_msstats_and_psm(msstats, map_dict):
        map_features = [
            "posterior_error_probability",
            "calculated_mz",
            "observed_mz",
            "mp_accessions",
            "is_decoy",
            "additional_scores",
            "cv_params",
        ]

        def merge_psm(rows, index):
            key = (
                rows["reference_file_name"],
                rows["peptidoform"],
                rows["precursor_charge"],
            )
            if key in map_dict:
                return map_dict[key][index]
            else:
                return None

        for i, feature in enumerate(map_features):
            msstats.loc[:, feature] = msstats[
                ["reference_file_name", "peptidoform", "precursor_charge"]
            ].apply(
                lambda rows: merge_psm(rows, i),
                axis=1,
            )

    def generate_feature(self, file_num=10, protein_str=None):

        feature_count = 0
        for msstats in self.generate_feature_report(file_num, protein_str):
            feature = self.transform_feature(msstats)

            feature_count += len(feature)
            self.logger.info(
                f"Generated {len(feature)} features, the total to {feature_count}."
            )

            yield feature

    def generate_feature_report(self, file_num=10, protein_str=None):
        map_dict, pep_dict = self.extract_psm_msg(protein_str)
        for msstats in self.transform_msstats_in(file_num, protein_str):
            if not msstats.empty:
                # Merge PSM data with MSstats aggregated features
                self.merge_msstats_and_psm_for_features(msstats, map_dict)
                # Add additional metadata fields
                self.add_additional_msg(msstats, pep_dict)
                # Convert data types for parquet format
                self.convert_to_parquet_format(msstats)
                yield msstats

    def merge_msstats_and_psm_for_features(self, msstats, map_dict):
        """Merge PSM data with aggregated feature data"""
        if msstats.empty:
            return

        map_features = [
            "posterior_error_probability",
            "calculated_mz",
            "observed_mz",
            "mp_accessions",
            "is_decoy",
            "additional_scores",
            "cv_params",
        ]

        def merge_psm(rows, index):
            # Use the anchor protein as the main protein identifier for PSM lookup
            # protein_key = rows.get("anchor_protein", "")
            key = (
                rows["reference_file_name"],
                rows["peptidoform"],
                str(rows["precursor_charge"]),
            )
            if key in map_dict:
                return map_dict[key][index]
            else:
                return None

        # Apply PSM data merging
        for i, feature in enumerate(map_features):
            if feature not in msstats.columns:
                msstats.loc[:, feature] = msstats[
                    ["reference_file_name", "peptidoform", "precursor_charge"]
                ].apply(
                    lambda rows: merge_psm(rows, i),
                    axis=1,
                )

    @staticmethod
    def transform_feature(df):
        return pa.Table.from_pandas(df, schema=FEATURE_SCHEMA)

    def write_feature_to_file(
        self,
        output_path,
        file_num=10,
        protein_file=None,
    ):
        protein_list = extract_protein_list(protein_file) if protein_file else None
        protein_str = "|".join(protein_list) if protein_list else None

        # Create file metadata for parquet file
        file_metadata = self._create_file_metadata()

        # Use the new generic batch writer with file metadata
        batch_writer = ParquetBatchWriter(
            output_path, FEATURE_SCHEMA, file_metadata=file_metadata
        )

        try:
            for feature_df in self.generate_feature(file_num, protein_str):
                if feature_df.num_rows > 0:
                    batch_writer.write_batch(feature_df.to_pylist())
        finally:
            batch_writer.close()

            if Path(output_path).exists():
                self.logger.info(
                    f"[Writer] Successfully wrote Feature to: {output_path}"
                )

            # Clean up the temporary MzTabIndexer
            self._indexer.cleanup_duckdb()

    @staticmethod
    def generate_best_scan(rows, pep_dict):
        key = (rows["peptidoform"], rows["precursor_charge"])
        if key in pep_dict:
            return [pep_dict[key][1], pep_dict[key][2]]
        else:
            return [None, None]

    def add_additional_msg(self, msstats, pep_dict):
        """Add additional metadata fields to the feature records"""

        # Add protein global qvalue (note: field name is pg_global_qvalue)
        if "anchor_protein" in msstats.columns:
            msstats.loc[:, "pg_global_qvalue"] = msstats["anchor_protein"].map(
                self._protein_global_qvalue_map
            )

        # Add best scan information
        if "peptidoform" in msstats.columns and "precursor_charge" in msstats.columns:
            msstats[["scan_reference_file_name", "scan"]] = msstats[
                ["peptidoform", "precursor_charge"]
            ].apply(
                lambda rows: self.generate_best_scan(rows, pep_dict),
                axis=1,
                result_type="expand",
            )

        # Process modifications
        if hasattr(self, "_mods_map") and self._mods_map is not None:
            msstats.loc[:, "modifications"] = msstats["peptidoform"].apply(
                lambda x: self.generate_modifications_details(x, self._mods_map)
            )
        else:
            # Add empty modifications
            msstats.loc[:, "modifications"] = None

        # Extract sequence from peptidoform (remove modifications)
        msstats.loc[:, "sequence"] = msstats["peptidoform"].apply(
            lambda x: self._extract_sequence_from_peptidoform(x) if pd.notna(x) else x
        )

        # Extract protein accession from anchor_protein
        if "anchor_protein" in msstats.columns:
            msstats["anchor_protein"] = msstats["anchor_protein"].apply(
                get_protein_accession
            )

        # Add additional fields with default values
        msstats.loc[:, "additional_intensities"] = None
        msstats.loc[:, "predicted_rt"] = None
        msstats.loc[:, "gg_accessions"] = None
        msstats.loc[:, "gg_names"] = None
        msstats.loc[:, "rt_start"] = None
        msstats.loc[:, "rt_stop"] = None
        msstats.loc[:, "ion_mobility"] = None
        msstats.loc[:, "start_ion_mobility"] = None
        msstats.loc[:, "stop_ion_mobility"] = None

    def _extract_sequence_from_peptidoform(self, peptidoform):
        """Extract plain sequence from peptidoform by removing modifications"""
        if not peptidoform:
            return peptidoform

        # Remove modifications in brackets [modification]
        sequence = BRACKET_MOD_PATTERN.sub("", peptidoform)
        # Remove terminal modifications like (modification)
        sequence = PAREN_MOD_PATTERN.sub("", sequence)
        # Remove any remaining special characters commonly used in modifications
        sequence = NON_AMINO_ACID_PATTERN.sub("", sequence.upper())

        return sequence if sequence else peptidoform

    @staticmethod
    def _filter_null_values(value):
        """
        Filter null values for nullable complex fields.

        Args:
            value: The value to check

        Returns:
            The original value if it's valid, None otherwise
        """
        if value is None:
            return None

        # For string or float types, check if they are pandas NaN
        if isinstance(value, (str, float)) and pd.isna(value):
            return None

        return value

    @staticmethod
    def _ensure_list_type(value):
        """
        Ensure value is a list type.

        Args:
            value: The value to check

        Returns:
            The original value if it's a list, empty list otherwise
        """
        return value if isinstance(value, list) else []

    @staticmethod
    def _ensure_dict_type(value):
        """
        Ensure value is a dict type.

        Args:
            value: The value to check

        Returns:
            The original value if it's a dict, empty dict otherwise
        """
        return value if isinstance(value, dict) else {}

    @staticmethod
    def _ensure_proper_list(value):
        """
        Ensure value is a proper list, converting single values to lists when appropriate.

        Args:
            value: The value to process

        Returns:
            - Original value if it's already a list
            - Single-item list [value] if value is not NaN and not empty string
            - Empty list otherwise
        """
        if isinstance(value, list):
            return value
        elif pd.notna(value) and value != "":
            return [value]
        else:
            return []

    @staticmethod
    def _convert_float_columns(res):
        """Convert float columns with proper NaN handling."""
        import pandas as pd

        float_columns = [
            "pg_global_qvalue",
            "calculated_mz",
            "observed_mz",
            "posterior_error_probability",
            "predicted_rt",
            "rt",
            "ion_mobility",
            "start_ion_mobility",
            "stop_ion_mobility",
            "rt_start",
            "rt_stop",
        ]

        for col in float_columns:
            if col in res.columns:
                res[col] = pd.to_numeric(res[col], errors="coerce").astype("float32")

    @staticmethod
    def _convert_integer_columns(res):
        """Convert integer columns with proper handling of NaN values."""
        import pandas as pd

        int_columns = ["precursor_charge", "unique", "is_decoy"]
        for col in int_columns:
            if col in res.columns:
                res[col] = pd.to_numeric(res[col], errors="coerce").astype("Int32")

    @staticmethod
    def _convert_string_columns(res):
        """Convert string columns."""
        string_columns = [
            "sequence",
            "peptidoform",
            "reference_file_name",
            "anchor_protein",
            "scan",
            "scan_reference_file_name",
        ]
        for col in string_columns:
            if col in res.columns:
                res[col] = res[col].astype(str)

    @staticmethod
    def _convert_list_columns(res):
        """Handle list columns (protein accessions, gene accessions, etc.)."""
        list_columns = ["pg_accessions", "gg_accessions", "gg_names"]
        for col in list_columns:
            if col in res.columns:
                # Ensure these are proper lists
                res[col] = res[col].apply(Feature._ensure_proper_list)

    @staticmethod
    def _convert_complex_columns(res):
        """Handle complex structured columns."""
        complex_columns = [
            "intensities",
            "additional_intensities",
            "modifications",
            "additional_scores",
            "cv_params",
            "file_metadata",
        ]

        for col in complex_columns:
            if col in res.columns:
                # Ensure proper structure for complex fields
                if (
                    col == "intensities"
                    or col == "modifications"
                    or col == "additional_scores"
                ):
                    res[col] = res[col].apply(Feature._ensure_list_type)
                elif col == "file_metadata":
                    # file_metadata should be a dict for each record
                    res[col] = res[col].apply(Feature._ensure_dict_type)
                else:
                    # For nullable complex fields, set None where appropriate
                    res[col] = res[col].apply(Feature._filter_null_values)

    @staticmethod
    def _ensure_required_fields(res):
        """Ensure all required fields exist with default values if missing."""
        required_fields = {
            "sequence": "",
            "peptidoform": "",
            "precursor_charge": 0,
            "is_decoy": 0,
            "calculated_mz": 0.0,
            "observed_mz": 0.0,
            "reference_file_name": "",
            "scan": "",
            "anchor_protein": "",
            "intensities": [],
            "file_metadata": {},
        }

        for field, default_value in required_fields.items():
            if field not in res.columns:
                res[field] = default_value

    @staticmethod
    def convert_to_parquet_format(res):
        """
        Convert DataFrame columns to appropriate types for Parquet format according to feature.avsc schema.
        This handles all field types including the new intensities structure and file_metadata.

        Parameters:
        -----------
        res : pandas.DataFrame
            DataFrame to convert

        Returns:
        --------
        None (modifies DataFrame in-place)
        """
        if res.empty:
            return

        # Convert different column types using helper methods
        Feature._convert_float_columns(res)
        Feature._convert_integer_columns(res)
        Feature._convert_string_columns(res)
        Feature._convert_list_columns(res)
        Feature._convert_complex_columns(res)
        Feature._ensure_required_fields(res)
