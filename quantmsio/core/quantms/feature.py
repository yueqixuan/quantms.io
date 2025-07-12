import logging
import tempfile
from pathlib import Path
from typing import Union

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

from quantmsio.core.common import FEATURE_SCHEMA

# MsstatsIN functionality now available in MzTabIndexer
from quantmsio.core.quantms.mztab import MzTabIndexer
from quantmsio.core.quantms.psm import Psm
from quantmsio.core.sdrf import SDRFHandler
from quantmsio.operate.tools import get_ahocorasick, get_protein_accession
from quantmsio.utils.file_utils import (
    close_file,
    extract_protein_list,
    save_slice_file,
    ParquetBatchWriter,
)


class Feature:
    """Feature processor using composition pattern.

    This class processes feature data from an MzTabIndexer instance, providing
    feature-specific functionality without inheriting from the indexer.
    """

    def __init__(self, mztab_indexer: MzTabIndexer, sdrf_path, msstats_in_path):
        """Initialize Feature processor with an MzTabIndexer instance.

        Args:
            mztab_indexer: An initialized MzTabIndexer instance
            sdrf_path: Path to SDRF file
            msstats_in_path: Path to MSstats input file
        """
        self._indexer = mztab_indexer
        self._msstats_in = msstats_in_path
        self._sdrf_path = sdrf_path
        self._ms_runs = self._extract_ms_runs()
        self._protein_global_qvalue_map = self._get_protein_map()
        self._score_names = self._get_score_names()
        self.experiment_type = SDRFHandler(sdrf_path).get_experiment_type_from_sdrf()
        self._mods_map = self._get_mods_map()
        self._automaton = get_ahocorasick(self._mods_map)
        self.logger = logging.getLogger(f"{__name__}.{self.__class__.__name__}")

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
            if (
                "accession" in proteins_df.columns
                and "opt_global_qvalue" in proteins_df.columns
            ):
                for _, row in proteins_df.iterrows():
                    if pd.notna(row["accession"]) and pd.notna(
                        row["opt_global_qvalue"]
                    ):
                        protein_map[row["accession"]] = row["opt_global_qvalue"]
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
            modifications = {}
            for _, row in metadata_df.iterrows():
                if "fixed_mod[" in row["key"] or "var_mod[" in row["key"]:
                    if "site" not in row["key"] and "position" not in row["key"]:
                        values = (
                            row["value"].replace("[", "").replace("]", "").split(",")
                        )
                        if len(values) >= 3:
                            accession = values[1].strip()
                            name = values[2].strip()
                            modifications[accession] = [
                                name,
                                row["key"].split("[")[1].split("]")[0],
                                None,
                                None,
                            ]

            mods_map = {}
            for accession, (name, index, site, position) in modifications.items():
                mods_map[name] = [accession, site or "X"]
            return mods_map
        except Exception as e:
            self.logger.warning(f"Could not get mods map: {e}")
            return {}

    def transform_msstats_in(
        self, file_num=10, protein_str=None, duckdb_max_memory="16GB", duckdb_threads=4
    ):
        # Check if msstats data is already loaded in the indexer
        if not self._indexer._msstats_path:
            # Add msstats data to the existing indexer
            self._indexer.add_msstats_table(self._msstats_in)

        # Use the enhanced MSstats analysis methods from MzTabIndexer
        for batch in self._indexer.iter_msstats_files(file_batch_size=file_num):
            if batch is not None and not batch.empty:
                # Apply protein filter if specified
                if protein_str:
                    batch = batch[
                        batch["ProteinName"].str.contains(protein_str, na=False)
                    ]

                if not batch.empty:
                    # Transform to the expected format for backward compatibility
                    batch_transformed = batch.rename(
                        columns={
                            "ProteinName": "ProteinName",
                            "Reference": "Reference",
                            "Intensity": "Intensity",
                            "PeptideSequence": "PeptideSequence",
                            "Channel": "Channel",
                        }
                    )
                    yield batch_transformed

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

    def generate_feature(
        self, file_num=10, protein_str=None, duckdb_max_memory="16GB", duckdb_threads=4
    ):
        for msstats in self.generate_feature_report(
            file_num, protein_str, duckdb_max_memory, duckdb_threads
        ):
            feature = self.transform_feature(msstats)
            yield feature

    def generate_feature_report(
        self, file_num=10, protein_str=None, duckdb_max_memory="16GB", duckdb_threads=4
    ):
        map_dict, pep_dict = self.extract_psm_msg(2000000, protein_str)
        for msstats in self.transform_msstats_in(
            file_num, protein_str, duckdb_max_memory, duckdb_threads
        ):
            self.merge_msstats_and_psm(msstats, map_dict)
            self.add_additional_msg(msstats, pep_dict)
            self.convert_to_parquet_format(msstats)
            yield msstats

    @staticmethod
    def slice(df, partitions):
        cols = df.columns
        if not isinstance(partitions, list):
            raise Exception(f"{partitions} is not a list")
        if len(partitions) == 0:
            raise Exception(f"{partitions} is empty")
        for partion in partitions:
            if partion not in cols:
                raise Exception(f"{partion} does not exist")
        for key, df in df.groupby(partitions):
            yield key, df

    def generate_slice_feature(
        self,
        partitions,
        file_num=10,
        protein_str=None,
        duckdb_max_memory="16GB",
        duckdb_threads=4,
    ):
        for msstats in self.generate_feature_report(
            file_num, protein_str, duckdb_max_memory, duckdb_threads
        ):
            for key, df in self.slice(msstats, partitions):
                feature = self.transform_feature(df)
                yield key, feature

    @staticmethod
    def transform_feature(df):
        return pa.Table.from_pandas(df, schema=FEATURE_SCHEMA)

    def write_feature_to_file(
        self,
        output_path,
        file_num=10,
        protein_file=None,
        duckdb_max_memory="16GB",
        duckdb_threads=4,
    ):
        protein_list = extract_protein_list(protein_file) if protein_file else None
        protein_str = "|".join(protein_list) if protein_list else None

        # Use the new generic batch writer
        batch_writer = ParquetBatchWriter(output_path, FEATURE_SCHEMA)

        try:
            for feature_df in self.generate_feature(
                file_num, protein_str, duckdb_max_memory, duckdb_threads
            ):
                if not feature_df.empty:
                    # The schema is applied when creating the table
                    records = feature_df.to_dict("records")
                    batch_writer.write_batch(records)
        finally:
            batch_writer.close()
            self.logger.info(f"Feature file written to {output_path}")

    def write_features_to_file(
        self,
        output_folder,
        filename,
        partitions,
        file_num=10,
        protein_file=None,
        duckdb_max_memory="16GB",
        duckdb_threads=4,
    ):
        logger = logging.getLogger("quantmsio.core.feature")

        # Log input and output paths
        logger.info(f"Input mzTab file: {self._indexer._mztab_path}")
        logger.info(f"Output folder: {output_folder}")
        logger.info(f"Base filename: {filename}")
        if protein_file:
            logger.info(f"Protein filter file: {protein_file}")

        pqwriters = {}
        protein_list = extract_protein_list(protein_file) if protein_file else None
        protein_str = "|".join(protein_list) if protein_list else None
        for key, feature in self.generate_slice_feature(
            partitions, file_num, protein_str, duckdb_max_memory, duckdb_threads
        ):
            pqwriters = save_slice_file(
                feature, pqwriters, output_folder, key, filename
            )
        close_file(pqwriters)

    @staticmethod
    def generate_best_scan(rows, pep_dict):
        key = (rows["peptidoform"], rows["precursor_charge"])
        if key in pep_dict:
            return [pep_dict[key][1], pep_dict[key][2]]
        else:
            return [None, None]

    def add_additional_msg(self, msstats, pep_dict):
        select_mods = list(self._mods_map.keys())
        msstats.loc[:, "pg_global_qvalue"] = msstats["mp_accessions"].map(
            self._protein_global_qvalue_map
        )
        msstats[["scan_reference_file_name", "scan"]] = msstats[
            ["peptidoform", "precursor_charge"]
        ].apply(
            lambda rows: self.generate_best_scan(rows, pep_dict),
            axis=1,
            result_type="expand",
        )
        msstats[["peptidoform", "modifications"]] = msstats[["peptidoform"]].apply(
            lambda row: self.generate_modifications_details(
                row["peptidoform"], self._mods_map, self._automaton, select_mods
            ),
            axis=1,
            result_type="expand",
        )
        msstats["mp_accessions"] = msstats["mp_accessions"].apply(get_protein_accession)
        msstats.loc[:, "additional_intensities"] = None
        msstats.loc[:, "predicted_rt"] = None
        msstats.loc[:, "gg_accessions"] = None
        msstats.loc[:, "gg_names"] = None
        msstats.loc[:, "rt_start"] = None
        msstats.loc[:, "rt_stop"] = None
        msstats.loc[:, "ion_mobility"] = None
        msstats.loc[:, "start_ion_mobility"] = None
        msstats.loc[:, "stop_ion_mobility"] = None

    @staticmethod
    def convert_to_parquet_format(res):
        """
        Convert DataFrame columns to appropriate types for Parquet format.
        This is optimized to handle NaN values and use vectorized operations.

        Parameters:
        -----------
        res : pandas.DataFrame
            DataFrame to convert

        Returns:
        --------
        None (modifies DataFrame in-place)
        """
        # Convert float columns in a single pass
        float_columns = [
            "pg_global_qvalue",
            "calculated_mz",
            "observed_mz",
            "posterior_error_probability",
        ]
        for col in float_columns:
            if col in res.columns:
                res[col] = pd.to_numeric(res[col], errors="coerce")

        # Convert integer columns with proper handling of NaN values
        res["unique"] = pd.to_numeric(res["unique"], errors="coerce").astype("Int32")

        # Use numpy for faster conversion of precursor_charge
        if "precursor_charge" in res.columns:
            res["precursor_charge"] = pd.to_numeric(
                res["precursor_charge"], errors="coerce"
            ).astype("Int32")

        # Convert is_decoy more efficiently
        if "is_decoy" in res.columns:
            res["is_decoy"] = pd.to_numeric(res["is_decoy"], errors="coerce").astype(
                "Int32"
            )

        # Convert string columns
        res["scan"] = res["scan"].astype(str)
        res["scan_reference_file_name"] = res["scan_reference_file_name"].astype(str)

        # Handle rt column
        if "rt" in res.columns:
            res["rt"] = pd.to_numeric(res["rt"], errors="coerce")
        else:
            res.loc[:, "rt"] = None
