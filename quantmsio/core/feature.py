from pathlib import Path
from typing import Union, Optional, Callable
import logging
import re

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import duckdb
from quantmsio.operate.tools import get_ahocorasick, get_protein_accession
from quantmsio.utils.file_utils import extract_protein_list, save_slice_file, close_file
from quantmsio.core.mztab import MzTab
from quantmsio.core.psm import Psm
from quantmsio.core.sdrf import SDRFHandler
from quantmsio.core.msstats_in import MsstatsIN
from quantmsio.core.common import FEATURE_SCHEMA
from quantmsio.utils.logger import get_logger


class Feature(MzTab):
    def __init__(self, mztab_path: Union[Path, str], sdrf_path, msstats_in_path):
        super(Feature, self).__init__(mztab_path)
        self._msstats_in = msstats_in_path
        self._sdrf_path = sdrf_path
        self._ms_runs = self.extract_ms_runs()
        self._protein_global_qvalue_map = self.get_protein_map()
        self._score_names = self.get_score_names()
        self.experiment_type = SDRFHandler(sdrf_path).get_experiment_type_from_sdrf()
        self._mods_map = self.get_mods_map()
        self._automaton = get_ahocorasick(self._mods_map)
        self.logger = get_logger("quantmsio.core.feature")

    def extract_psm_msg(self, chunksize=2000000, protein_str=None):
        psm = Psm(self.mztab_path)
        pep_dict = psm.extract_from_pep(chunksize=100000)
        map_dict = {}
        for psm_chunk in psm.iter_psm_table(chunksize, protein_str):
            for key, df in psm_chunk.groupby(
                ["reference_file_name", "peptidoform", "precursor_charge"]
            ):
                df.reset_index(drop=True, inplace=True)
                temp_df = df.iloc[df["posterior_error_probability"].idxmin()]
                if key not in map_dict:
                    map_dict[key] = [None for _ in range(7)]
                pep_value = temp_df["posterior_error_probability"]
                if map_dict[key][0] is None or float(map_dict[key][0]) > float(
                    pep_value
                ):
                    map_dict[key][0] = pep_value
                    map_dict[key][1] = temp_df["calculated_mz"]
                    map_dict[key][2] = temp_df["observed_mz"]
                    map_dict[key][3] = temp_df["mp_accessions"]
                    map_dict[key][4] = temp_df["is_decoy"]
                    map_dict[key][5] = temp_df["additional_scores"]
                    map_dict[key][6] = temp_df["cv_params"]
        return map_dict, pep_dict

    def transform_msstats_in(
        self, file_num=10, protein_str=None, duckdb_max_memory="16GB", duckdb_threads=4
    ):
        msstats_in = MsstatsIN(
            self._msstats_in, self._sdrf_path, duckdb_max_memory, duckdb_threads
        )
        for msstats in msstats_in.generate_msstats_in(file_num, protein_str):
            yield msstats
        msstats_in.destroy_duckdb_database()

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
        output_path: str,
        file_num: int = 50,
        protein_file: Optional[str] = None,
        duckdb_max_memory: Optional[str] = None,
        duckdb_threads: Optional[int] = None,
        progress_callback: Optional[Callable[[int, int, str], None]] = None,
    ) -> None:
        """
        Write features to a single parquet file.
        
        Args:
            output_path: Path to write the output file
            file_num: Number of files to process at once
            protein_file: Optional protein file path
            duckdb_max_memory: Maximum memory for DuckDB
            duckdb_threads: Number of threads for DuckDB
            progress_callback: Optional callback for progress updates
        """
        try:
            if progress_callback:
                progress_callback(0, 100, "Starting feature conversion")
            
            # Configure DuckDB
            if duckdb_max_memory or duckdb_threads:
                self.logger.debug("Configuring DuckDB settings...")
                if duckdb_max_memory:
                    self.logger.debug(f"Setting DuckDB max memory to {duckdb_max_memory}")
                    duckdb.config.set_memory_limit(duckdb_max_memory)
                if duckdb_threads:
                    self.logger.debug(f"Setting DuckDB threads to {duckdb_threads}")
                    duckdb.config.set_threads(duckdb_threads)
            
            # Extract protein information
            if progress_callback:
                progress_callback(10, 100, "Processing protein information")
            self.logger.debug("Processing protein information...")
            protein_list = None
            protein_str = None
            if protein_file:
                try:
                    protein_list = extract_protein_list(protein_file)
                    if protein_list:
                        # Join with pipe character for regex OR matching
                        protein_str = "|".join(re.escape(p) for p in protein_list)
                        self.logger.debug(f"Found {len(protein_list)} proteins to filter")
                    else:
                        self.logger.warning("No proteins found in protein file")
                except Exception as e:
                    self.logger.warning(f"Error processing protein file: {str(e)}")
            
            # Extract PSM messages
            if progress_callback:
                progress_callback(20, 100, "Extracting PSM data")
            self.logger.debug("Extracting PSM data...")
            map_dict, pep_dict = self.extract_psm_msg(2000000, protein_str)
            
            # Process features
            pqwriter = None
            processed_count = 0
            total_count = 0
            
            for msstats in self.transform_msstats_in(
                file_num, protein_str, duckdb_max_memory, duckdb_threads
            ):
                if total_count == 0:
                    total_count = len(msstats)
                
                if progress_callback:
                    progress = min(30 + (processed_count / total_count * 50), 80)
                    progress_callback(int(progress), 100, "Processing features")
                
                self.logger.debug(f"Processing batch of {len(msstats)} features...")
                self.merge_msstats_and_psm(msstats, map_dict)
                self.add_additional_msg(msstats, pep_dict)
                self.convert_to_parquet_format(msstats)
                
                feature = self.transform_feature(msstats)
                if not pqwriter:
                    self.logger.debug("Creating Parquet writer...")
                    pqwriter = pq.ParquetWriter(output_path, feature.schema)
                
                pqwriter.write_table(feature)
                processed_count += len(msstats)
            
            if progress_callback:
                progress_callback(90, 100, "Finalizing output file")
            
            self.logger.debug("Closing Parquet writer...")
            close_file(pqwriter=pqwriter)
            
            if progress_callback:
                progress_callback(100, 100, "Feature conversion completed")
            self.logger.debug("Feature conversion completed successfully")
            
        except Exception as e:
            self.logger.error(f"Error during feature conversion: {str(e)}")
            raise

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
