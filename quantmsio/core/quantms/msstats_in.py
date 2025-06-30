from pathlib import Path
from typing import Generator, List, Optional, Union

import pandas as pd

from quantmsio.core.common import MSSTATS_MAP, MSSTATS_USECOLS
from quantmsio.core.duckdb import DuckDB
from quantmsio.core.sdrf import SDRFHandler
from quantmsio.operate.tools import get_protein_accession
from quantmsio.utils.constants import ITRAQ_CHANNEL, TMT_CHANNELS
from quantmsio.utils.pride_utils import clean_peptidoform_sequence


class MsstatsIN(DuckDB):
    def __init__(
        self,
        report_path: Union[Path, str],
        sdrf_path: Union[Path, str],
        duckdb_max_memory="16GB",
        duckdb_threads=4,
    ):
        super(MsstatsIN, self).__init__(report_path, duckdb_max_memory, duckdb_threads)
        self._sdrf = SDRFHandler(sdrf_path)
        self.experiment_type = self._sdrf.get_experiment_type_from_sdrf()
        self._sample_map = self._sdrf.get_sample_map_run()
        self._optimized_setup_done = False

    def __enter__(self):
        """Context manager entry."""
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Context manager exit with cleanup."""
        self.destroy_duckdb_database()

    def _setup_optimized_processing(self):
        """Create optimized database views and tables for processing."""
        if self._optimized_setup_done:
            return

        try:
            # First, detect available columns in the report table
            self._detect_available_columns()

            # Create channel mapping table for TMT/iTRAQ
            if self.experiment_type != "LFQ":
                self._create_channel_mapping_table()

            # Create sample mapping table
            self._create_sample_mapping_table()

            # Create optimized views
            self._create_processed_view()

            self._optimized_setup_done = True
        except Exception as e:
            print(f"Warning: Could not setup optimized processing: {e}")
            # Fall back to original processing

    def _detect_available_columns(self):
        """Detect what columns are available in the report table."""
        try:
            # Try different methods to get column information
            try:
                columns_query = "PRAGMA table_info('report')"
                columns_info = self._duckdb.execute(columns_query).df()
                available_columns = set(columns_info["name"].tolist())
            except Exception as e:
                # Fallback: get columns from a sample query
                self.logger.warning(
                    f"Failed to get column info via PRAGMA: {e}, trying fallback method"
                )
                try:
                    sample_query = "SELECT * FROM report LIMIT 1"
                    sample_data = self._duckdb.execute(sample_query).df()
                    available_columns = set(sample_data.columns.tolist())
                except Exception as e2:
                    self.logger.warning(f"Fallback column detection also failed: {e2}")
                    # Default column set based on common msstats format
                    available_columns = {
                        "Reference",
                        "ProteinName",
                        "PeptideSequence",
                        "Intensity",
                    }

            self._available_columns = available_columns

        except Exception as e:
            # Default column set based on common msstats format
            self._available_columns = {
                "Reference",
                "ProteinName",
                "PeptideSequence",
                "Intensity",
            }

    def _create_channel_mapping_table(self):
        """Create a channel mapping table in DuckDB for efficient lookups."""
        if "TMT" in self.experiment_type:
            channels = TMT_CHANNELS[self.experiment_type]
        else:
            channels = ITRAQ_CHANNEL[self.experiment_type]

        # Create mapping data
        mapping_data = [
            {"channel_num": i + 1, "channel_name": channel}
            for i, channel in enumerate(channels)
        ]

        # Drop existing table if it exists
        self._duckdb.execute("DROP TABLE IF EXISTS channel_mapping")

        # Insert into DuckDB using pandas DataFrame
        mapping_df = pd.DataFrame(mapping_data)
        self._duckdb.execute("CREATE TABLE channel_mapping AS SELECT * FROM mapping_df")

    def _create_sample_mapping_table(self):
        """Create sample mapping table for efficient joins."""
        sample_data = [
            {"file_channel": k, "sample_accession": v}
            for k, v in self._sample_map.items()
        ]

        # Drop existing table if it exists
        self._duckdb.execute("DROP TABLE IF EXISTS sample_mapping")

        # Insert into DuckDB using pandas DataFrame
        sample_df = pd.DataFrame(sample_data)
        self._duckdb.execute("CREATE TABLE sample_mapping AS SELECT * FROM sample_df")

    def _create_processed_view(self):
        """Create a processed view with all transformations applied."""

        # Drop existing view if it exists
        self._duckdb.execute("DROP VIEW IF EXISTS processed_msstats")

        # Check what columns are available
        has_channel = "Channel" in self._available_columns
        has_precursor_charge = "PrecursorCharge" in self._available_columns
        has_charge = "Charge" in self._available_columns
        has_retention_time = "RetentionTime" in self._available_columns

        # Build charge column expression
        if has_precursor_charge and has_charge:
            charge_expr = "COALESCE(PrecursorCharge, Charge)"
        elif has_precursor_charge:
            charge_expr = "PrecursorCharge"
        elif has_charge:
            charge_expr = "Charge"
        else:
            charge_expr = "3"  # Default charge

        # Build retention time expression
        rt_expr = "RetentionTime" if has_retention_time else "NULL"

        if self.experiment_type == "LFQ" or not has_channel:
            # LFQ view or when Channel column is missing
            sql = f"""
            CREATE VIEW processed_msstats AS
            SELECT 
                ProteinName as pg_accessions_raw,
                split_part(Reference, '.', 1) as reference_file_name,
                Intensity as intensity,
                PeptideSequence as peptidoform,
                {rt_expr} as rt,
                'LFQ' as channel,
                {charge_expr} as precursor_charge,
                -- Extract protein accession
                CASE 
                    WHEN ProteinName LIKE '%|%|%' 
                    THEN split_part(ProteinName, '|', 2)
                    ELSE ProteinName
                END as anchor_protein,
                -- Check if unique protein
                CASE 
                    WHEN ProteinName LIKE '%;%' OR ProteinName LIKE '%,%' 
                    THEN 0 ELSE 1 
                END as unique_protein,
                -- Create peptide map key
                split_part(Reference, '.', 1) || PeptideSequence || CAST({charge_expr} AS VARCHAR) as peptide_map_key,
                -- Sample mapping
                COALESCE(sm.sample_accession, split_part(Reference, '.', 1) || '-LFQ') as sample_accession
            FROM report r
            LEFT JOIN sample_mapping sm ON (split_part(r.Reference, '.', 1) || '-LFQ') = sm.file_channel
            """
        else:
            # TMT/iTRAQ view
            sql = f"""
            CREATE VIEW processed_msstats AS
            SELECT 
                ProteinName as pg_accessions_raw,
                split_part(Reference, '.', 1) as reference_file_name,
                Intensity as intensity,
                PeptideSequence as peptidoform,
                {rt_expr} as rt,
                COALESCE(cm.channel_name, CAST(r.Channel AS VARCHAR)) as channel,
                {charge_expr} as precursor_charge,
                -- Extract protein accession
                CASE 
                    WHEN ProteinName LIKE '%|%|%' 
                    THEN split_part(ProteinName, '|', 2)
                    ELSE ProteinName
                END as anchor_protein,
                -- Check if unique protein
                CASE 
                    WHEN ProteinName LIKE '%;%' OR ProteinName LIKE '%,%' 
                    THEN 0 ELSE 1 
                END as unique_protein,
                -- Create peptide map key
                split_part(Reference, '.', 1) || PeptideSequence || CAST({charge_expr} AS VARCHAR) as peptide_map_key,
                -- Sample mapping
                COALESCE(sm.sample_accession, split_part(Reference, '.', 1) || '-' || COALESCE(cm.channel_name, CAST(r.Channel AS VARCHAR))) as sample_accession
            FROM report r
            LEFT JOIN channel_mapping cm ON r.Channel = cm.channel_num
            LEFT JOIN sample_mapping sm ON (split_part(r.Reference, '.', 1) || '-' || COALESCE(cm.channel_name, CAST(r.Channel AS VARCHAR))) = sm.file_channel
            """

        self._duckdb.execute(sql)

    def get_runs(self):
        references = self._duckdb.sql("SELECT Reference FROM report").df()
        references = references["Reference"].str.split(".").str[0]
        return list(set(references))

    def iter_runs(self, file_num=10, columns: list = None):
        references = self.get_runs()
        ref_list = [
            references[i : i + file_num] for i in range(0, len(references), file_num)
        ]
        for refs in ref_list:
            batch_df = self.query_field("Reference", refs, columns)
            yield batch_df

    def get_file_statistics(self) -> pd.DataFrame:
        """Get comprehensive file-level statistics using SQL."""
        self._setup_optimized_processing()

        if not self._optimized_setup_done:
            # Fallback to basic statistics
            return self._get_basic_file_statistics()

        sql = """
        SELECT 
            reference_file_name,
            channel,
            COUNT(*) as feature_count,
            COUNT(DISTINCT pg_accessions_raw) as protein_count,
            COUNT(DISTINCT peptidoform) as peptide_count,
            AVG(intensity) as avg_intensity,
            MIN(intensity) as min_intensity,
            MAX(intensity) as max_intensity,
            SUM(CASE WHEN unique_protein = 1 THEN 1 ELSE 0 END) as unique_features,
            COUNT(DISTINCT sample_accession) as sample_count
        FROM processed_msstats 
        GROUP BY reference_file_name, channel
        ORDER BY reference_file_name, channel
        """

        return self._duckdb.execute(sql).df()

    def _get_basic_file_statistics(self) -> pd.DataFrame:
        """Fallback basic file statistics without optimized views."""

        # First check what columns are available
        try:
            columns_query = "PRAGMA table_info('report')"
            columns_info = self._duckdb.execute(columns_query).df()
            available_columns = set(columns_info["name"].tolist())
        except Exception as e:
            # If pragma doesn't work, try a different approach
            self.logger.warning(
                f"Failed to get column info via PRAGMA: {e}, trying fallback method"
            )
            try:
                sample_query = "SELECT * FROM report LIMIT 1"
                sample_data = self._duckdb.execute(sample_query).df()
                available_columns = set(sample_data.columns.tolist())
            except Exception as e2:
                # Fallback to basic assumptions
                self.logger.warning(f"Fallback column detection also failed: {e2}")
                available_columns = {
                    "Reference",
                    "ProteinName",
                    "PeptideSequence",
                    "Intensity",
                }

        # Build the query based on available columns and experiment type
        if self.experiment_type == "LFQ" or "Channel" not in available_columns:
            sql = """
            SELECT 
                split_part(Reference, '.', 1) as reference_file_name,
                'LFQ' as channel,
                COUNT(*) as feature_count,
                COUNT(DISTINCT ProteinName) as protein_count,
                COUNT(DISTINCT PeptideSequence) as peptide_count,
                AVG(Intensity) as avg_intensity,
                MIN(Intensity) as min_intensity,
                MAX(Intensity) as max_intensity
            FROM report 
            GROUP BY reference_file_name, channel
            ORDER BY reference_file_name, channel
            """
        else:
            sql = """
            SELECT 
                split_part(Reference, '.', 1) as reference_file_name,
                CAST(Channel AS VARCHAR) as channel,
                COUNT(*) as feature_count,
                COUNT(DISTINCT ProteinName) as protein_count,
                COUNT(DISTINCT PeptideSequence) as peptide_count,
                AVG(Intensity) as avg_intensity,
                MIN(Intensity) as min_intensity,
                MAX(Intensity) as max_intensity
            FROM report 
            GROUP BY reference_file_name, channel
            ORDER BY reference_file_name, channel
            """

        return self._duckdb.execute(sql).df()

    def get_protein_file_matrix(
        self, protein_filter: Optional[str] = None
    ) -> pd.DataFrame:
        """Get protein presence matrix across files using SQL."""
        self._setup_optimized_processing()

        where_clause = ""
        params = []
        if protein_filter:
            where_clause = "WHERE pg_accessions_raw LIKE ?"
            params.append(f"%{protein_filter}%")

        if self._optimized_setup_done:
            sql = f"""
            SELECT 
                anchor_protein,
                COUNT(DISTINCT reference_file_name) as file_count,
                string_agg(DISTINCT reference_file_name, ';') as files_present,
                COUNT(*) as total_features,
                AVG(intensity) as avg_intensity
            FROM processed_msstats 
            {where_clause}
            GROUP BY anchor_protein
            HAVING COUNT(DISTINCT reference_file_name) > 1
            ORDER BY file_count DESC, total_features DESC
            """
        else:
            # Fallback query
            sql = f"""
            SELECT 
                ProteinName as anchor_protein,
                COUNT(DISTINCT split_part(Reference, '.', 1)) as file_count,
                string_agg(DISTINCT split_part(Reference, '.', 1), ';') as files_present,
                COUNT(*) as total_features,
                AVG(Intensity) as avg_intensity
            FROM report 
            {where_clause}
            GROUP BY ProteinName
            HAVING COUNT(DISTINCT split_part(Reference, '.', 1)) > 1
            ORDER BY file_count DESC, total_features DESC
            """

        return self._duckdb.execute(sql, params).df()

    def iter_files_optimized(
        self, file_batch_size: int = 10
    ) -> Generator[pd.DataFrame, None, None]:
        """
        Optimized file iteration using SQL window functions for better memory management.
        """
        self._setup_optimized_processing()

        if not self._optimized_setup_done:
            # Fallback to original method
            yield from self.iter_runs(file_batch_size)
            return

        # Get file list with row numbers for batching
        files_sql = """
        SELECT DISTINCT reference_file_name,
               ROW_NUMBER() OVER (ORDER BY reference_file_name) as file_rank
        FROM processed_msstats
        """
        files_df = self._duckdb.execute(files_sql).df()

        total_files = len(files_df)

        for batch_start in range(0, total_files, file_batch_size):
            batch_end = min(batch_start + file_batch_size, total_files)

            # Get files for this batch
            batch_files = files_df.iloc[batch_start:batch_end][
                "reference_file_name"
            ].tolist()
            file_list_str = "', '".join(batch_files)

            # Use SQL to get batch data efficiently
            batch_sql = f"""
            SELECT *
            FROM processed_msstats 
            WHERE reference_file_name IN ('{file_list_str}')
            ORDER BY reference_file_name, peptidoform, channel
            """

            yield self._duckdb.execute(batch_sql).df()

    def aggregate_channels_sql(
        self, file_batch_size: int = 10
    ) -> Generator[pd.DataFrame, None, None]:
        """
        Aggregate channels per peptide using pure SQL for TMT/iTRAQ experiments.
        """
        self._setup_optimized_processing()

        if not self._optimized_setup_done:
            # Fallback to original method
            for batch in self.generate_msstats_in(file_batch_size):
                yield batch
            return

        if self.experiment_type == "LFQ":
            # For LFQ, just return the processed data
            for batch in self.iter_files_optimized(file_batch_size):
                batch["intensities"] = batch.apply(
                    lambda row: [
                        {
                            "sample_accession": row["sample_accession"],
                            "channel": row["channel"],
                            "intensity": row["intensity"],
                        }
                    ],
                    axis=1,
                )
                yield batch[
                    [
                        "reference_file_name",
                        "pg_accessions_raw",
                        "peptidoform",
                        "precursor_charge",
                        "anchor_protein",
                        "unique_protein",
                        "intensities",
                    ]
                ]
        else:
            # For TMT/iTRAQ, aggregate channels using SQL
            for batch_files in self._get_file_batches(file_batch_size):
                file_list_str = "', '".join(batch_files)

                sql = f"""
                WITH channel_aggregated AS (
                    SELECT 
                        reference_file_name,
                        pg_accessions_raw,
                        peptidoform,
                        precursor_charge,
                        anchor_protein,
                        unique_protein,
                        -- Aggregate all channels for each peptide using array aggregation
                        array_agg(STRUCT_PACK(
                            sample_accession := sample_accession,
                            channel := channel,
                            intensity := intensity
                        )) as intensities,
                        -- Keep one representative channel for deduplication
                        MIN(channel) as representative_channel
                    FROM processed_msstats 
                    WHERE reference_file_name IN ('{file_list_str}')
                    GROUP BY reference_file_name, pg_accessions_raw, peptidoform, precursor_charge, 
                             anchor_protein, unique_protein
                )
                SELECT 
                    reference_file_name,
                    pg_accessions_raw,
                    peptidoform, 
                    precursor_charge,
                    anchor_protein,
                    unique_protein,
                    representative_channel as channel,
                    intensities
                FROM channel_aggregated
                ORDER BY reference_file_name, peptidoform
                """

                try:
                    result = self._duckdb.execute(sql).df()
                    if len(result) > 0:
                        yield result
                except Exception as e:
                    print(
                        f"Warning: SQL aggregation failed ({e}), falling back to pandas"
                    )
                    # Fallback to pandas aggregation for this batch
                    for batch in self.iter_files_optimized(len(batch_files)):
                        if len(batch) > 0:
                            # Convert to expected format using pandas
                            batch = self._pandas_channel_aggregation(batch)
                            yield batch

    def _pandas_channel_aggregation(self, batch: pd.DataFrame) -> pd.DataFrame:
        """Fallback pandas-based channel aggregation."""
        intensities_map = {}

        def get_intensities_map(row):
            key = row["peptide_map_key"]
            if key not in intensities_map:
                intensities_map[key] = []
            intensities_map[key].append(
                {
                    "sample_accession": row["sample_accession"],
                    "channel": row["channel"],
                    "intensity": row["intensity"],
                }
            )

        batch.apply(get_intensities_map, axis=1)

        # Remove duplicates and map intensities
        result = batch.drop_duplicates(
            subset=["reference_file_name", "peptidoform", "precursor_charge"]
        ).copy()
        result["intensities"] = result["peptide_map_key"].map(intensities_map)

        return result[
            [
                "reference_file_name",
                "pg_accessions_raw",
                "peptidoform",
                "precursor_charge",
                "anchor_protein",
                "unique_protein",
                "channel",
                "intensities",
            ]
        ]

    def _get_file_batches(self, batch_size: int) -> Generator[List[str], None, None]:
        """Get file names in batches."""
        files = self._duckdb.execute(
            "SELECT DISTINCT split_part(Reference, '.', 1) as ref FROM report ORDER BY ref"
        ).df()
        file_list = files["ref"].tolist()

        for i in range(0, len(file_list), batch_size):
            yield file_list[i : i + batch_size]

    def generate_msstats_in(self, file_num=10, protein_str=None):
        """Original method maintained for backward compatibility."""
        msstats_map = MSSTATS_MAP.copy()
        usecols = list(MSSTATS_USECOLS)
        if self.experiment_type == "LFQ":
            usecols.remove("Channel")
            usecols.remove("RetentionTime")
            usecols += ["PrecursorCharge"]
            msstats_map["PrecursorCharge"] = "precursor_charge"
        else:
            usecols += ["Charge"]
            msstats_map["Charge"] = "precursor_charge"
        for msstats in self.iter_runs(file_num=file_num, columns=usecols):
            if self.experiment_type == "LFQ":
                msstats.loc[:, "Channel"] = "LFQ"
                msstats.loc[:, "RetentionTime"] = None
            if protein_str:
                msstats = msstats[
                    msstats["ProteinName"].str.contains(f"{protein_str}", na=False)
                ]
            msstats.rename(columns=msstats_map, inplace=True)
            self.transform_msstats_in(msstats)
            self.transform_experiment(msstats)
            yield msstats

    def generate_msstats_in_optimized(
        self, file_num: int = 10, protein_str: Optional[str] = None
    ) -> Generator[pd.DataFrame, None, None]:
        """
        Optimized version of generate_msstats_in using SQL aggregation.
        """
        # Add protein filter if specified
        if protein_str:
            print(f"Filtering for protein: {protein_str}")

        for batch in self.aggregate_channels_sql(file_num):
            if len(batch) > 0:
                # Apply final transformations
                if "pg_accessions_raw" in batch.columns:
                    batch["pg_accessions"] = batch["pg_accessions_raw"].apply(
                        get_protein_accession
                    )
                    batch["sequence"] = batch["peptidoform"].apply(
                        clean_peptidoform_sequence
                    )

                    # Filter by protein if specified
                    if protein_str:
                        batch = batch[
                            batch["pg_accessions_raw"].str.contains(
                                protein_str, na=False
                            )
                        ]

                    # Rename columns to match expected output
                    batch.rename(columns={"unique_protein": "unique"}, inplace=True)

                if len(batch) > 0:
                    yield batch

    def transform_msstats_in(self, msstats):
        msstats["reference_file_name"] = (
            msstats["reference_file_name"].str.split(".").str[0]
        )
        msstats.loc[:, "sequence"] = msstats["peptidoform"].apply(
            clean_peptidoform_sequence
        )
        if self.experiment_type != "LFQ":
            if "TMT" in self.experiment_type:
                msstats["channel"] = msstats["channel"].apply(
                    lambda row: TMT_CHANNELS[self.experiment_type][row - 1]
                )
            else:
                msstats["channel"] = msstats["channel"].apply(
                    lambda row: ITRAQ_CHANNEL[self.experiment_type][row - 1]
                )
        msstats.loc[:, "unique"] = msstats["pg_accessions"].apply(
            lambda x: 0 if ";" in x or "," in x else 1
        )
        msstats["pg_accessions"] = msstats["pg_accessions"].apply(get_protein_accession)
        msstats.loc[:, "anchor_protein"] = msstats["pg_accessions"].str[0]

    def transform_experiment(self, msstats):
        intensities_map = {}
        select_cols = [
            "map",
            "reference_file_name",
            "peptidoform",
            "precursor_charge",
            "channel",
            "intensity",
        ]
        if self.experiment_type != "LFQ":
            msstats.loc[:, "map"] = (
                msstats["reference_file_name"]
                + msstats["peptidoform"]
                + msstats["precursor_charge"].astype(str)
            )

            def get_intensities_map(rows):
                key = rows["map"]
                sample_key = rows["reference_file_name"] + "-" + rows["channel"]
                if key not in intensities_map:
                    intensities_map[key] = []
                intensities_map[key].append(
                    {
                        "sample_accession": self._sample_map[sample_key],
                        "channel": rows["channel"],
                        "intensity": rows["intensity"],
                    }
                )

            msstats[select_cols].apply(get_intensities_map, axis=1)
            msstats.drop_duplicates(
                subset=["reference_file_name", "peptidoform", "precursor_charge"],
                inplace=True,
            )
            msstats.reset_index(inplace=True, drop=True)
            msstats.loc[:, "intensities"] = msstats["map"].map(intensities_map)
            msstats.drop(["map"], inplace=True, axis=1)
        else:
            msstats.loc[:, "intensities"] = msstats[
                ["reference_file_name", "channel", "intensity"]
            ].apply(
                lambda rows: [
                    {
                        "sample_accession": self._sample_map[
                            rows["reference_file_name"] + "-" + rows["channel"]
                        ],
                        "channel": rows["channel"],
                        "intensity": rows["intensity"],
                    }
                ],
                axis=1,
            )

    def __del__(self):
        """Cleanup database views and tables."""
        try:
            if hasattr(self, "_duckdb") and self._duckdb and self._optimized_setup_done:
                self._duckdb.execute("DROP VIEW IF EXISTS processed_msstats")
                self._duckdb.execute("DROP TABLE IF EXISTS channel_mapping")
                self._duckdb.execute("DROP TABLE IF EXISTS sample_mapping")
                self._duckdb.execute("DROP TABLE IF EXISTS protein_groups")
                self._duckdb.execute("DROP VIEW IF EXISTS processed_msstats_with_pg")
            # Always call parent cleanup to close connection and remove database file
            if hasattr(self, "_duckdb") and self._duckdb:
                self.destroy_duckdb_database()
        except Exception as e:
            import logging

            logging.getLogger("quantmsio.core.msstats_in").warning(
                f"Exception during __del__ cleanup: {e}"
            )
