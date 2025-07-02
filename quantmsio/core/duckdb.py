import logging
import os
import shutil
import tempfile
import time
import gzip
from pathlib import Path
from typing import Generator, Optional, Union, Dict, List
from functools import lru_cache

import duckdb
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

from quantmsio.core.project import create_uuid_filename
from quantmsio.utils.constants import (
    PROTEIN_GROUP,
    RUN,
    MODIFIED_SEQUENCE,
    Q_VALUE,
    PG_Q_VALUE,
    PRECURSOR_QUANTITY,
    OPT_GLOBAL_RESULT_TYPE,
    INDISTINGUISHABLE_GROUP,
    SINGLE_PROTEIN_MZTAB,
    PROTEIN_DETAILS_MZTAB,
    ACCESSION,
    AMBIGUITY_MEMBERS,
    ANCHOR_PROTEIN,
    SPECTRA_REF,
    SPECTRA_REF_FILE,
    SPECTRA_REF_SCAN,
    MSSTATS_PROTEIN_NAME,
    MSSTATS_PEPTIDE_SEQUENCE,
    MSSTATS_REFERENCE,
    MSSTATS_REFERENCE_NAME,
    MSSTATS_INTENSITY,
)


class DuckDB:
    """Base DuckDB class for database operations."""

    def __init__(self, database_name: str):
        """Initialize DuckDB connection.

        Args:
            database_name: Name/path of the database file
        """
        self._duckdb_name = str(database_name)
        self._duckdb = None
        self.logger = logging.getLogger(f"{__name__}.{self.__class__.__name__}")

    def initialize_database(
        self, max_memory: Optional[str] = None, worker_threads: Optional[int] = None
    ) -> None:
        """Initialize the database with configuration.

        Args:
            max_memory: Maximum memory to use (e.g. '4GB')
            worker_threads: Number of worker threads to use
        """
        s = time.time()

        self._duckdb = duckdb.connect(self._duckdb_name)
        self._duckdb.execute("SET enable_progress_bar=true")

        if max_memory is not None:
            self._duckdb.execute(f"SET max_memory='{max_memory}'")
        if worker_threads is not None:
            self._duckdb.execute(f"SET worker_threads='{worker_threads}'")

        msg = self._duckdb.execute(
            "SELECT * FROM duckdb_settings() where name in ('worker_threads', 'max_memory')"
        ).df()
        self.logger.info(f"duckdb uses {str(msg['value'][0])} threads.")
        self.logger.info(f"duckdb uses {str(msg['value'][1])} of memory.")

        et = time.time() - s
        self.logger.info(f"Time to initialize duckdb {et} seconds")

    def create_table_from_file(
        self, table_name: str, file_path: str, indices: Optional[list[str]] = None
    ) -> None:
        """Create a table from a file and optionally create indices.

        Args:
            table_name: Name of the table to create
            file_path: Path to the source file
            indices: Optional list of column names to create indices on
        """
        if not self._duckdb:
            raise RuntimeError(
                "Database not initialized. Call initialize_database first."
            )

        self._duckdb.execute(
            f"CREATE TABLE {table_name} AS SELECT * FROM '{file_path}'"
        )

        if indices:
            for column in indices:
                # Replace periods with underscores in the index name
                safe_column = column.replace(".", "_")
                self._duckdb.execute(
                    f'CREATE INDEX IF NOT EXISTS idx_{table_name}_{safe_column} ON {table_name} ("{column}")'
                )

    def get_unique_values(self, table_name: str, column: str) -> list:
        """Get unique values from a column in a table.

        Args:
            table_name: Name of the table
            column: Column name to get unique values from

        Returns:
            List of unique values
        """
        if not self._duckdb:
            raise RuntimeError("Database not initialized")

        result = self._duckdb.execute(
            f'SELECT DISTINCT "{column}" FROM {table_name}'
        ).fetchall()
        return [x[0] for x in result if x[0] is not None]

    def query_to_df(self, query: str) -> "pd.DataFrame":
        """Execute a query and return results as a DataFrame.

        Args:
            query: SQL query to execute

        Returns:
            DataFrame with query results
        """
        if not self._duckdb:
            raise RuntimeError("Database not initialized")

        return self._duckdb.execute(query).df()

    def destroy_database(self) -> None:
        """Close the database connection and clean up."""
        if self._duckdb:
            self._duckdb.close()
            self._duckdb = None
        if self._duckdb_name and os.path.exists(self._duckdb_name):
            os.remove(self._duckdb_name)


class DiannDuckDB(DuckDB):
    """DuckDB implementation specific to DIA-NN data."""

    def __init__(
        self,
        diann_report_path: Union[Path, str],
        max_memory: str = "16GB",
        worker_threads: int = 4,
        pg_matrix_path: Optional[Union[Path, str]] = None,
        cache_size: int = 128,
    ):
        """Initialize DiannDuckDB.

        Args:
            diann_report_path: Path to DIA-NN report file
            max_memory: Maximum memory to use
            worker_threads: Number of worker threads
            pg_matrix_path: Optional path to protein groups matrix file
            cache_size: Size of the LRU cache for frequently accessed data (default: 128)
        """
        self._report_path = str(diann_report_path)
        self._pg_matrix_path = str(pg_matrix_path) if pg_matrix_path else None
        self._cache_size = cache_size
        database_name = create_uuid_filename("diann-report", ".db")
        super().__init__(database_name)

        # Initialize database and create report table
        self.initialize_database(max_memory, worker_threads)
        self.create_table_from_file("report", self._report_path, [PROTEIN_GROUP, RUN])

        # Load protein groups matrix if provided
        if self._pg_matrix_path:
            self.create_table_from_file(
                "pg_matrix", self._pg_matrix_path, [PROTEIN_GROUP]
            )

        # Initialize cache for common statistics
        self._init_cache()

    def _init_cache(self):
        """Initialize cache decorators for frequently accessed methods."""
        # Cache for unique values
        self.get_unique_values = lru_cache(maxsize=self._cache_size)(
            super().get_unique_values
        )

        # Cache for protein group matrix
        self._cached_get_protein_group_matrix = lru_cache(maxsize=self._cache_size)(
            self._get_protein_group_matrix_impl
        )

        # Cache for common statistics
        self._cached_get_statistics = lru_cache(maxsize=self._cache_size)(
            self._get_statistics_impl
        )

    def _get_protein_group_matrix_impl(
        self, protein_groups_tuple: Optional[tuple[str, ...]] = None
    ) -> pd.DataFrame:
        """Implementation of protein group matrix retrieval.

        Args:
            protein_groups_tuple: Optional tuple of protein groups to filter by

        Returns:
            DataFrame containing protein group matrix data
        """
        if not self._pg_matrix_path:
            raise RuntimeError(
                "Protein groups matrix was not loaded. Provide pg_matrix_path during initialization."
            )

        query = "SELECT * FROM pg_matrix"
        if protein_groups_tuple:
            query += f' WHERE "{PROTEIN_GROUP}" = ANY({list(protein_groups_tuple)})'

        return self.query_to_df(query)

    def get_protein_group_matrix(
        self, protein_groups: Optional[List[str]] = None
    ) -> pd.DataFrame:
        """Get protein group matrix data with caching.

        Args:
            protein_groups: Optional list of protein groups to filter by

        Returns:
            DataFrame containing protein group matrix data
        """
        # Convert list to tuple for caching (lists are not hashable)
        protein_groups_tuple = tuple(protein_groups) if protein_groups else None
        return self._cached_get_protein_group_matrix(protein_groups_tuple)

    def _get_statistics_impl(self, q_value_threshold: float = 0.01) -> Dict[str, float]:
        """Implementation of common statistics calculation.

        Args:
            q_value_threshold: Q-value threshold for filtering

        Returns:
            Dictionary containing common statistics
        """
        stats = {}

        # Basic counts
        stats.update(
            self.query_to_df(
                f"""
            SELECT 
                COUNT(*) as total_psms,
                COUNT(DISTINCT "{PROTEIN_GROUP}") as total_proteins,
                COUNT(DISTINCT "{MODIFIED_SEQUENCE}") as unique_peptides,
                COUNT(DISTINCT "{RUN}") as total_runs
            FROM report
        """
            )
            .iloc[0]
            .to_dict()
        )

        # Q-value based statistics
        stats.update(
            self.query_to_df(
                f"""
            SELECT 
                COUNT(CASE WHEN CAST("{Q_VALUE}" AS FLOAT) <= {q_value_threshold} THEN 1 END) as psms_passing_qvalue,
                COUNT(DISTINCT CASE WHEN CAST("{PG_Q_VALUE}" AS FLOAT) <= {q_value_threshold} THEN "{PROTEIN_GROUP}" END) as proteins_passing_qvalue
            FROM report
        """
            )
            .iloc[0]
            .to_dict()
        )

        # Intensity statistics
        stats.update(
            self.query_to_df(
                f"""
            SELECT 
                MIN(CAST("{PRECURSOR_QUANTITY}" AS FLOAT)) as min_intensity,
                MAX(CAST("{PRECURSOR_QUANTITY}" AS FLOAT)) as max_intensity,
                AVG(CAST("{PRECURSOR_QUANTITY}" AS FLOAT)) as avg_intensity
            FROM report
            WHERE "{PRECURSOR_QUANTITY}" IS NOT NULL
        """
            )
            .iloc[0]
            .to_dict()
        )

        return stats

    def get_statistics(self, q_value_threshold: float = 0.01) -> Dict[str, float]:
        """Get common statistics with caching.

        Args:
            q_value_threshold: Q-value threshold for filtering

        Returns:
            Dictionary containing common statistics
        """
        return self._cached_get_statistics(q_value_threshold)

    def clear_cache(self):
        """Clear all cached data."""
        self.get_unique_values.cache_clear()
        self._cached_get_protein_group_matrix.cache_clear()
        self._cached_get_statistics.cache_clear()

    def destroy_database(self) -> None:
        """Close the database connection and clean up."""
        self.clear_cache()
        super().destroy_database()

    def get_unique_references(self, column: str) -> list:
        """Get unique values from a column in the report.

        Args:
            column: Column name to get unique values from

        Returns:
            List of unique values
        """
        return self.get_unique_values("report", column)

    def get_report_from_database(self, refs: list, columns: str) -> pd.DataFrame:
        """Get report data for specific references and columns.

        Args:
            refs: List of reference values to filter by
            columns: Comma-separated string of column names to select

        Returns:
            DataFrame containing the requested data
        """
        query = f"""
        SELECT {columns}
        FROM report
        WHERE "{RUN}" = ANY({refs})
        """
        return self.query_to_df(query)

    def iter_file(
        self, field: str, file_num: int = 10, columns: Optional[list[str]] = None
    ) -> "Generator[tuple[list, pd.DataFrame], None, None]":
        """Iterate over report files in batches.

        Args:
            field: Field to group by
            file_num: Number of files per batch
            columns: Optional list of columns to select

        Yields:
            Tuple of (list of references, DataFrame with data)
        """
        references = self.get_unique_references(field)
        ref_list = [
            references[i : i + file_num] for i in range(0, len(references), file_num)
        ]
        for refs in ref_list:
            cols = "*"
            if columns:
                cols = ", ".join(f'"{c}"' for c in columns)
                cols = cols.replace('"unique"', "unique")

            query = f"""
            SELECT {cols}
            FROM report
            WHERE "{field}" = ANY('{refs}')
            """
            yield refs, self.query_to_df(query)

    def query_field(
        self, field: str, queries: list, columns: Optional[list[str]] = None
    ) -> pd.DataFrame:
        """Query report by field values.

        Args:
            field: Field to query on
            queries: List of values to match
            columns: Optional list of columns to select

        Returns:
            DataFrame with matching rows
        """
        field_conditions = [f"{field} LIKE '%{q}%'" for q in queries]
        where_clause = " OR ".join(field_conditions)
        cols = ", ".join(columns) if columns else "*"
        query = f"SELECT {cols} FROM report WHERE {where_clause}"
        return self.query_to_df(query)


class MzTabIndexer(DuckDB):
    """DuckDB implementation specific to mzTab data."""

    _MZTAB_INDEXER_FOLDER = "mztab_indexer_"

    _MZTAB_INDEXER_TABLE_METADATA = "metadata"
    _MZTAB_INDEXER_TABLE_PROTEINS = "proteins"
    _MZTAB_INDEXER_TABLE_PROTEIN_DETAILS = "protein_details"
    _MZTAB_INDEXER_TABLE_PSMS = "psms"
    _MZTAB_INDEXER_TABLE_MSSTATS = "msstats"

    def __init__(
        self,
        mztab_path: Optional[Union[Path, str]] = None,
        msstats_path: Optional[Union[Path, str]] = None,
        database_path: Optional[Union[Path, str]] = None,
        max_memory: str = "16GB",
        worker_threads: int = 8,
        batch_size: int = 100000,
        protein_columns: Optional[List[str]] = None,
        psm_columns: Optional[List[str]] = None,
    ):
        """Initialize MzTabIndexer.

        Args:
            mztab_path: Path to mzTab file. Required if database_path is not provided.
            msstats_path: Optional path to MSstats file.
            database_path: Optional path to an existing DuckDB database file.
            max_memory: Maximum memory to use.
            worker_threads: Number of worker threads.
            batch_size: Number of rows to process in each batch.
            protein_columns: Optional list of protein columns to include (if None, include all).
            psm_columns: Optional list of PSM columns to include (if None, include all).
        """
        if not database_path and (
            not mztab_path or not MzTabIndexer.validate_file(mztab_path)
        ):
            raise ValueError(
                "A valid mztab_path must be provided if database_path is not set."
            )

        db_name = (
            database_path
            if database_path
            else create_uuid_filename("mztab-quantms", ".db")
        )
        super().__init__(db_name)

        self._mztab_path = str(mztab_path) if mztab_path else None
        self._msstats_path = str(msstats_path) if msstats_path else None
        self._batch_size = batch_size
        self._temp_dir = None
        self._protein_columns = protein_columns
        self._psm_columns = psm_columns

        self.initialize_database(max_memory, worker_threads)

        if not database_path:
            self._create_temp_dir()
            self.create_tables()

    @staticmethod
    def validate_file(mztab_path: Union[str, Path]) -> bool:
        """Validate that the mzTab file exists and is not empty."""
        if not Path(mztab_path).exists():
            raise FileNotFoundError(f"mzTab file not found: {mztab_path}")
        if Path(mztab_path).stat().st_size == 0:
            raise ValueError("mzTab file is empty")
        return mztab_path.stat().st_size > 0

    def _create_temp_dir(self):
        """Create temporary directory for parquet files."""
        self._temp_dir = tempfile.mkdtemp(prefix=self._MZTAB_INDEXER_FOLDER)

    def _cleanup_temp_dir(self):
        """Clean up temporary directory."""
        if self._temp_dir and os.path.exists(self._temp_dir):
            shutil.rmtree(self._temp_dir)

    def __enter__(self):
        """Context manager entry."""
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Context manager exit with cleanup."""
        self._cleanup_temp_dir()
        self.destroy_database()

    def _process_mztab_to_parquet(self):
        """Process mzTab file to parquet files in batches."""
        # Determine if file is gzipped
        is_gzipped = self._mztab_path.endswith(".gz")

        metadata_parquet = Path(self._temp_dir) / "metadata.parquet"
        proteins_parquet = Path(self._temp_dir) / "proteins.parquet"
        protein_details_parquet = Path(self._temp_dir) / "protein_details.parquet"
        psms_parquet = Path(self._temp_dir) / "psms.parquet"

        # Initialize parquet writers
        proteins_writer = None
        protein_details_writer = None
        psms_writer = None

        # Open file (gzipped or not)
        open_func = gzip.open if is_gzipped else open
        with open_func(self._mztab_path, "rt") as f:
            # Process metadata section
            metadata = []
            proteins = []
            protein_details = []
            psms = []

            current_section = None

            for line in f:
                line = line.strip()
                if not line:
                    continue

                if line.startswith("MTD"):
                    current_section = "metadata"
                    metadata.append(self._parse_metadata_line(line))
                elif line.startswith("PRH"):
                    current_section = "proteins"
                    # Get protein header
                    self._protein_header = self._parse_header_line(line)
                elif line.startswith("PRT"):
                    protein_dict = self._parse_protein_line(line)
                    if OPT_GLOBAL_RESULT_TYPE in protein_dict:
                        if protein_dict[OPT_GLOBAL_RESULT_TYPE] in [
                            INDISTINGUISHABLE_GROUP,
                            SINGLE_PROTEIN_MZTAB,
                        ]:
                            proteins.append(protein_dict)
                        elif (
                            protein_dict[OPT_GLOBAL_RESULT_TYPE]
                            == PROTEIN_DETAILS_MZTAB
                        ):
                            protein_details.append(protein_dict)
                    if len(proteins) >= self._batch_size:
                        proteins_writer = self._write_protein_batch(
                            proteins, proteins_parquet, proteins_writer
                        )
                        proteins = []
                    if len(protein_details) >= self._batch_size:
                        protein_details_writer = self._write_protein_details_batch(
                            protein_details,
                            protein_details_parquet,
                            protein_details_writer,
                        )
                        protein_details = []
                elif line.startswith("PSH"):
                    current_section = "psms"
                    # Get PSM header
                    self._psm_header = self._parse_header_line(line)
                elif line.startswith("PSM"):
                    psms.append(self._parse_psm_line(line))
                    if len(psms) >= self._batch_size:
                        psms_writer = self._write_psm_batch(
                            psms, psms_parquet, psms_writer
                        )
                        psms = []

            # Write any remaining batches
            if metadata:
                self._write_metadata_batch(metadata, metadata_parquet)
            if proteins:
                self._write_protein_batch(proteins, proteins_parquet, proteins_writer)
            if protein_details:
                self._write_protein_details_batch(
                    protein_details, protein_details_parquet, protein_details_writer
                )
            if psms:
                self._write_psm_batch(psms, psms_parquet, psms_writer)

        return metadata_parquet, proteins_parquet, protein_details_parquet, psms_parquet

    def _parse_metadata_line(self, line: str) -> Dict:
        """Parse metadata line into dictionary."""
        parts = line.split("\t")
        return {"key": parts[1], "value": parts[2] if len(parts) > 2 else None}

    def _parse_header_line(self, line: str) -> List[str]:
        """Parse header line into column names."""
        return line.split("\t")[1:]

    def _parse_protein_line(self, line: str) -> Dict:
        """Parse protein line into dictionary."""
        values = line.split("\t")
        protein_dict = dict(zip(self._protein_header, values[1:]))

        if OPT_GLOBAL_RESULT_TYPE in protein_dict:
            result_type = protein_dict[OPT_GLOBAL_RESULT_TYPE]
            if result_type == INDISTINGUISHABLE_GROUP:
                anchor = protein_dict.get(ACCESSION)
                ambiguity_members = protein_dict.get(AMBIGUITY_MEMBERS)
                if ambiguity_members and ambiguity_members != "null":
                    members = sorted(ambiguity_members.split(","))
                    protein_dict[ACCESSION] = ",".join(members)
                protein_dict[ANCHOR_PROTEIN] = anchor
            elif result_type == SINGLE_PROTEIN_MZTAB:
                protein_dict[ANCHOR_PROTEIN] = protein_dict.get(ACCESSION)

        return protein_dict

    def _parse_psm_line(self, line: str) -> Dict:
        """Parse PSM line into dictionary."""
        values = line.split("\t")
        psm_dict = dict(zip(self._psm_header, values[1:]))
        if (
            ACCESSION in psm_dict
            and psm_dict[ACCESSION]
            and psm_dict[ACCESSION] != "null"
        ):
            accessions = sorted(psm_dict[ACCESSION].split(","))
            psm_dict[ACCESSION] = ",".join(accessions)

        # Split spectra_ref into file and scan info for plex experiments
        if SPECTRA_REF in psm_dict and psm_dict[SPECTRA_REF]:
            spectra_ref_parts = psm_dict[SPECTRA_REF].split("_controllerType")
            psm_dict[SPECTRA_REF_FILE] = spectra_ref_parts[0]
            if len(spectra_ref_parts) > 1:
                psm_dict[SPECTRA_REF_SCAN] = "_controllerType" + spectra_ref_parts[1]
            else:
                psm_dict[SPECTRA_REF_SCAN] = None
            del psm_dict[SPECTRA_REF]

        return psm_dict

    def _write_metadata_batch(self, metadata: List[Dict], output_path: Path):
        """Write metadata batch to parquet file."""
        df = pd.DataFrame(metadata)
        df.to_parquet(output_path, index=False)

    def _write_protein_batch(
        self, proteins: List[Dict], output_path: Path, writer=None
    ):
        """Write protein batch to parquet file."""
        df = pd.DataFrame(proteins)

        # Filter columns if specified
        if self._protein_columns:
            df = df.reindex(columns=self._protein_columns)

        table = pa.Table.from_pandas(df)

        if writer is None:
            # First batch or new writer needed
            writer = pq.ParquetWriter(output_path, table.schema)

        writer.write_table(table)
        return writer

    def _write_protein_details_batch(
        self, protein_details: List[Dict], output_path: Path, writer=None
    ):
        """Write protein details batch to parquet file."""
        df = pd.DataFrame(protein_details)
        table = pa.Table.from_pandas(df)

        if writer is None:
            # First batch or new writer needed
            writer = pq.ParquetWriter(output_path, table.schema)

        writer.write_table(table)
        return writer

    def _write_psm_batch(self, psms: List[Dict], output_path: Path, writer=None):
        """Write PSM batch to parquet file."""
        df = pd.DataFrame(psms)

        # Filter columns if specified
        if self._psm_columns:
            df = df.reindex(columns=self._psm_columns)

        table = pa.Table.from_pandas(df)

        if writer is None:
            # First batch or new writer needed
            writer = pq.ParquetWriter(output_path, table.schema)

        writer.write_table(table)
        return writer

    def create_tables(self):
        """Create DuckDB tables from parquet files."""
        # Ensure temp directory exists
        if not self._temp_dir:
            self._create_temp_dir()

        # Process mzTab to parquet files
        self.logger.info("Processing mzTab file to parquet format...")
        metadata_parquet, proteins_parquet, protein_details_parquet, psms_parquet = (
            self._process_mztab_to_parquet()
        )
        self.logger.info("Finished processing mzTab file to parquet format.")

        # Create tables from parquet files if they don't exist
        tables = self._duckdb.execute("SHOW TABLES").fetchall()
        table_names = [table[0] for table in tables]

        if self._MZTAB_INDEXER_TABLE_METADATA not in table_names:
            self.logger.info(f"Creating table: {self._MZTAB_INDEXER_TABLE_METADATA}")
            self._duckdb.execute(
                f"CREATE TABLE {self._MZTAB_INDEXER_TABLE_METADATA} AS SELECT * FROM read_parquet('{metadata_parquet}')"
            )

        if (
            self._MZTAB_INDEXER_TABLE_PROTEINS not in table_names
            and os.path.exists(proteins_parquet)
            and os.path.getsize(proteins_parquet) > 0
        ):
            self.logger.info(f"Creating table: {self._MZTAB_INDEXER_TABLE_PROTEINS}")
            self._duckdb.execute(
                f"CREATE TABLE {self._MZTAB_INDEXER_TABLE_PROTEINS} AS SELECT * FROM read_parquet('{proteins_parquet}')"
            )

        if (
            self._MZTAB_INDEXER_TABLE_PROTEIN_DETAILS not in table_names
            and os.path.exists(protein_details_parquet)
            and os.path.getsize(protein_details_parquet) > 0
        ):
            self.logger.info(
                f"Creating table: {self._MZTAB_INDEXER_TABLE_PROTEIN_DETAILS}"
            )
            self._duckdb.execute(
                f"CREATE TABLE {self._MZTAB_INDEXER_TABLE_PROTEIN_DETAILS} AS SELECT * FROM read_parquet('{protein_details_parquet}')"
            )

        if (
            self._MZTAB_INDEXER_TABLE_PSMS not in table_names
            and os.path.exists(psms_parquet)
            and os.path.getsize(psms_parquet) > 0
        ):
            self.logger.info(f"Creating table: {self._MZTAB_INDEXER_TABLE_PSMS}")
            self._duckdb.execute(
                f"CREATE TABLE {self._MZTAB_INDEXER_TABLE_PSMS} AS SELECT * FROM read_parquet('{psms_parquet}')"
            )

        # Create MSstats table if path provided and table doesn't exist
        if self._msstats_path and self._MZTAB_INDEXER_TABLE_MSSTATS not in table_names:
            self.logger.info(f"Creating table: {self._MZTAB_INDEXER_TABLE_MSSTATS}")
            self.add_msstats_table(self._msstats_path)

    def get_metadata(self) -> pd.DataFrame:
        """Get metadata as DataFrame."""
        return self.query_to_df(f"SELECT * FROM {self._MZTAB_INDEXER_TABLE_METADATA}")

    def get_proteins(self) -> pd.DataFrame:
        """Get proteins as DataFrame."""
        return self.query_to_df(f"SELECT * FROM {self._MZTAB_INDEXER_TABLE_PROTEINS}")

    def get_protein_details(self) -> pd.DataFrame:
        """Get protein details as DataFrame."""
        return self.query_to_df(
            f"SELECT * FROM {self._MZTAB_INDEXER_TABLE_PROTEIN_DETAILS}"
        )

    def get_psms(self) -> pd.DataFrame:
        """Get PSMs as DataFrame."""
        return self.query_to_df(f"SELECT * FROM {self._MZTAB_INDEXER_TABLE_PSMS}")

    def get_msstats(self) -> Optional[pd.DataFrame]:
        """Get MSstats as DataFrame if available."""
        if self._msstats_path:
            return self.query_to_df(
                f"SELECT * FROM {self._MZTAB_INDEXER_TABLE_MSSTATS}"
            )
        return None

    def get_protein_by_accession(self, accession: str) -> pd.DataFrame:
        """Get protein data by accession.

        Args:
            accession: Protein accession

        Returns:
            DataFrame with protein data
        """
        query = f"""
        SELECT *
        FROM {self._MZTAB_INDEXER_TABLE_PROTEINS}
        WHERE {ACCESSION} = '{accession}'
        """
        return self.query_to_df(query)

    def get_psms_by_protein(self, accession: str) -> pd.DataFrame:
        """Get PSMs for a specific protein.

        Args:
            accession: Protein accession

        Returns:
            DataFrame with PSM data
        """
        query = f"""
        SELECT *
        FROM {self._MZTAB_INDEXER_TABLE_PSMS}
        WHERE {ACCESSION} LIKE '%{accession}%'
        """
        return self.query_to_df(query)

    def get_metadata_value(self, key: str) -> Optional[str]:
        """Get metadata value by key.

        Args:
            key: Metadata key

        Returns:
            Metadata value if found, None otherwise
        """
        query = f"""
        SELECT value
        FROM {self._MZTAB_INDEXER_TABLE_METADATA}
        WHERE key = '{key}'
        LIMIT 1
        """
        result = self.query_to_df(query)
        return result["value"].iloc[0] if not result.empty else None

    def get_protein_count(self) -> int:
        """Get total number of proteins."""
        query = f"SELECT COUNT(*) as count FROM {self._MZTAB_INDEXER_TABLE_PROTEINS}"
        result = self.query_to_df(query)
        return result["count"].iloc[0]

    def get_protein_details_count(self) -> int:
        """Get total number of protein details."""
        query = (
            f"SELECT COUNT(*) as count FROM {self._MZTAB_INDEXER_TABLE_PROTEIN_DETAILS}"
        )
        result = self.query_to_df(query)
        return result["count"].iloc[0]

    def get_psm_count(self) -> int:
        """Get total number of PSMs."""
        query = f"SELECT COUNT(*) as count FROM {self._MZTAB_INDEXER_TABLE_PSMS}"
        result = self.query_to_df(query)
        return result["count"].iloc[0]

    def get_msstats_by_peptide(self, peptide: str) -> Optional[pd.DataFrame]:
        """Get MSstats data for a specific peptide.

        Args:
            peptide: Peptide sequence

        Returns:
            DataFrame with MSstats data if available
        """
        if not self._msstats_path:
            return None

        query = f"""
        SELECT *
        FROM {self._MZTAB_INDEXER_TABLE_MSSTATS}
        WHERE {MSSTATS_PEPTIDE_SEQUENCE} = '{peptide}'
        """
        return self.query_to_df(query)

    def get_msstats_runs(self) -> Optional[List[str]]:
        """Get list of runs from MSstats data.

        Returns:
            List of run names if MSstats data is available
        """
        if not self._msstats_path:
            return None

        query = f"""
        SELECT DISTINCT {MSSTATS_REFERENCE}
        FROM {self._MZTAB_INDEXER_TABLE_MSSTATS}
        ORDER BY {MSSTATS_REFERENCE}
        """
        result = self.query_to_df(query)
        return result[MSSTATS_REFERENCE].tolist() if not result.empty else []

    def get_msstats_protein_intensities(self, protein: str) -> Optional[pd.DataFrame]:
        """Get protein intensities across runs from MSstats data.

        Args:
            protein: Protein accession or name

        Returns:
            DataFrame with protein intensities if available
        """
        if not self._msstats_path:
            return None

        query = f"""
        SELECT {MSSTATS_REFERENCE}, {MSSTATS_INTENSITY}
        FROM {self._MZTAB_INDEXER_TABLE_MSSTATS}
        WHERE {MSSTATS_PROTEIN_NAME} LIKE '%{protein}%'
        ORDER BY {MSSTATS_REFERENCE}
        """
        return self.query_to_df(query)

    def add_msstats_table(self, msstats_path: str):
        """
        Adds or replaces an msstats table in an existing database.
        If the table already exists, it will be replaced with the new data.
        """
        if not self._duckdb:
            raise RuntimeError("Database not initialized.")

        self.logger.info("Attempting to add/replace msstats table.")
        self._msstats_path = msstats_path

        # Check if the msstats table already exists
        self.logger.debug("Checking for existing msstats table")
        tables = self._duckdb.execute("SHOW TABLES").fetchall()
        table_names = [table[0] for table in tables]
        if self._MZTAB_INDEXER_TABLE_MSSTATS in table_names:
            self.logger.debug("Dropping existing msstats table")
            self._duckdb.execute(f"DROP TABLE {self._MZTAB_INDEXER_TABLE_MSSTATS}")

        # Create temp directory if it doesn't exist
        if not self._temp_dir:
            self.logger.debug("Creating temporary directory")
            self._create_temp_dir()

        try:
            # Process MSstats to parquet with transformations
            self.logger.debug("Starting MSstats to parquet processing")
            msstats_parquet = self._process_msstats_to_parquet()

            # Create table from parquet
            self.logger.debug("Creating DuckDB table from parquet file")
            self._duckdb.execute(
                f"CREATE TABLE {self._MZTAB_INDEXER_TABLE_MSSTATS} AS SELECT * FROM read_parquet('{msstats_parquet}')"
            )

            # Verify table was created and has data
            row_count = self._duckdb.execute(
                f"SELECT COUNT(*) FROM {self._MZTAB_INDEXER_TABLE_MSSTATS}"
            ).fetchone()[0]
            self.logger.debug(f"Created MSstats table with {row_count} rows")

            # Get the columns of the newly created table
            self.logger.debug("Getting table columns for indexing")
            table_columns = [
                col[0]
                for col in self._duckdb.execute(
                    f"SELECT name FROM pragma_table_info('{self._MZTAB_INDEXER_TABLE_MSSTATS}')"
                ).fetchall()
            ]

            # Define columns to be indexed
            columns_to_index = [MSSTATS_PROTEIN_NAME]

            # Create indices on the msstats table for existing columns
            self.logger.debug("Creating indices on columns")
            for column in columns_to_index:
                if column in table_columns:
                    self.logger.debug(f"Creating index for column: {column}")
                    safe_column_name = column.replace(".", "_")
                    self._duckdb.execute(
                        f'CREATE INDEX IF NOT EXISTS idx_msstats_{safe_column_name} ON {self._MZTAB_INDEXER_TABLE_MSSTATS} ("{column}")'
                    )
                else:
                    self.logger.warning(
                        f"Column '{column}' not found in msstats table, skipping index creation."
                    )

            self.logger.info(
                f"Successfully added/replaced msstats table from {msstats_path}"
            )
        except Exception as e:
            self.logger.error(f"Error adding MSstats table: {str(e)}")
            raise
        finally:
            # Clean up temp directory if we created it
            if self._temp_dir:
                self.logger.debug("Cleaning up temporary directory")
                self._cleanup_temp_dir()
                self._temp_dir = None

    def _process_msstats_to_parquet(self) -> Path:
        """Process MSstats file to parquet format with transformations.

        Returns:
            Path to the created parquet file
        """
        self.logger.debug(
            f"Starting MSstats processing to parquet: {self._msstats_path}"
        )
        msstats_parquet = Path(self._temp_dir) / "msstats.parquet"
        writer = None
        total_rows = 0

        try:
            # Process file in chunks to avoid memory issues
            self.logger.debug(
                f"Reading MSstats file in chunks of size {self._batch_size}"
            )
            for chunk_idx, df in enumerate(
                pd.read_csv(self._msstats_path, chunksize=self._batch_size)
            ):
                chunk_size = len(df)
                total_rows += chunk_size

                # Sort and join protein names
                if MSSTATS_PROTEIN_NAME in df.columns:
                    df[MSSTATS_PROTEIN_NAME] = df[MSSTATS_PROTEIN_NAME].apply(
                        lambda x: (
                            ";".join(sorted(str(x).split(";"))) if pd.notna(x) else x
                        )
                    )

                # Create Reference_Name by removing file extensions and controller info
                if MSSTATS_REFERENCE in df.columns:
                    df[MSSTATS_REFERENCE_NAME] = df[MSSTATS_REFERENCE].apply(
                        lambda x: (
                            x.split("_controllerType")[0].rsplit(".", 1)[0]
                            if isinstance(x, str)
                            else x
                        )
                    )

                # Convert to pyarrow table
                table = pa.Table.from_pandas(df)

                if writer is None:
                    # First chunk - create writer with schema
                    writer = pq.ParquetWriter(
                        msstats_parquet,
                        table.schema,
                        compression="ZSTD",  # Use ZSTD compression for better compression ratio
                        compression_level=3,  # Moderate compression level for balance
                        use_dictionary=True,  # Enable dictionary encoding for strings
                        dictionary_pagesize_limit=1048576,  # 1MB dictionary page size
                    )

                # Write chunk to parquet
                writer.write_table(table)

                # Log progress
                self.logger.debug(f"Processed chunk {chunk_idx + 1} ({len(df)} rows)")
        except Exception as e:
            self.logger.error(f"Error processing MSstats file: {str(e)}")
            if writer:
                writer.close()
            raise
        finally:
            if writer:
                self.logger.debug("Closing parquet writer")
                writer.close()

        self.logger.info(
            f"Completed MSstats processing. Total rows processed: {total_rows}"
        )
        return msstats_parquet
