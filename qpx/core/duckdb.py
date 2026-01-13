import logging
import time
from pathlib import Path
from typing import Generator, Optional, Union, Dict, List, Any
from functools import lru_cache

import duckdb
import pandas as pd

from qpx.core.project import create_uuid_filename
from qpx.utils.constants import (
    PROTEIN_GROUP,
    RUN,
    MODIFIED_SEQUENCE,
    Q_VALUE,
    PG_Q_VALUE,
    PRECURSOR_QUANTITY,
)


class DuckDB:
    """Base class for DuckDB operations."""

    def __init__(
        self,
        database_path: Optional[Union[Path, str]] = None,
        max_memory: str = "16GB",
        worker_threads: int = 8,
    ):
        """Initialize DuckDB connection.

        Args:
            database_path: Optional path to DuckDB database file
            max_memory: Maximum memory to use
            worker_threads: Number of worker threads
        """
        self.logger = logging.getLogger(f"{__name__}.{self.__class__.__name__}")
        self.logger.info(f"duckdb uses {worker_threads} threads.")
        self.logger.info(f"duckdb uses {max_memory} of memory.")

        # Store database path for cleanup
        self._database_path = str(database_path) if database_path else None

        start_time = time.time()
        self._duckdb = duckdb.connect(
            database=self._database_path or ":memory:",
            config={
                "memory_limit": max_memory,
                "threads": worker_threads,
            },
        )
        self.logger.info(
            f"Time to initialize duckdb {(time.time() - start_time):.2f} seconds"
        )

    def destroy_database(self):
        """Close DuckDB connection and clean up database file."""
        self.cleanup_duckdb()

    def cleanup_duckdb(self):
        """Check if DuckDB connection is closed, then delete the database file."""
        # Close connection if it is still open
        if self._duckdb:
            try:
                self._duckdb.close()
                self.logger.info("[Check] DuckDB connection closed.")
            except Exception as e:
                self.logger.info(f"Failed to close DuckDB connection: {e}")
            finally:
                self._duckdb = None

        db_file = Path(self._database_path)
        # Delete the database file using pathlib
        if db_file.exists():
            try:
                db_file.unlink()
                self.logger.info(f"[CleanUp] Database file deleted: {db_file}")
            except Exception as e:
                self.logger.info(f"Failed to delete database file: {e}")

    def query_to_df(self, query: str) -> pd.DataFrame:
        """Execute query and return result as DataFrame."""
        return self._duckdb.execute(query).df()

    def secure_query_to_df(
        self, query: str, params: Optional[Dict[str, Any]] = None
    ) -> pd.DataFrame:
        """Execute parameterized query safely and return result as DataFrame.

        Args:
            query: SQL query with placeholders (e.g., 'SELECT * FROM table WHERE id = ?')
            params: Dictionary of parameters to safely substitute

        Returns:
            DataFrame with query results

        Example:
            >>> secure_query_to_df("SELECT * FROM proteins WHERE accession = ?", {"accession": "P12345"})
        """
        if params:
            # DuckDB uses positional parameters with ?
            param_values = list(params.values())
            return self._duckdb.execute(query, param_values).df()
        else:
            return self._duckdb.execute(query).df()

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

    def check_tables_in_db(self, tables: List[str]) -> bool:
        """Check if the database connection is valid."""
        tables = [t[0] for t in duckdb.execute("SHOW TABLES").fetchall()]
        for table in tables:
            if table not in tables:
                logging.info(f"Table {table} not in database")
                return False
        return True

    @staticmethod
    def check_tables_in_parquet_db(
        database_path: Union[Path, str], tables: List[str]
    ) -> bool:
        """Check if parquet files exist in the database directory.

        Args:
            database_path: Path to the parquet database directory
            tables: List of table names to check (without .parquet extension)

        Returns:
            True if all required parquet files exist, False otherwise
        """
        database_path = Path(database_path)
        if not database_path.exists() or not database_path.is_dir():
            return False

        for table in tables:
            parquet_file = database_path / f"{table}.parquet"
            if not parquet_file.exists():
                logging.info(f"Parquet file {parquet_file} not found")
                return False
        return True

    def get_query_performance_stats(
        self, query: str, iterations: int = 5
    ) -> Dict[str, Any]:
        """
        Measure query performance statistics.

        This method executes a query multiple times and measures execution time
        to provide performance statistics. Useful for benchmarking and optimizing
        database queries.

        Args:
            query: SQL query to test
            iterations: Number of iterations to run (default: 5)

        Returns:
            Dictionary with performance statistics including:
            - query: The tested query
            - iterations: Number of iterations run
            - min_time: Minimum execution time in seconds
            - max_time: Maximum execution time in seconds
            - avg_time: Average execution time in seconds
            - std_time: Standard deviation of execution times
            - error: Error message if query fails

        Example:
            >>> stats = db.get_query_performance_stats("SELECT COUNT(*) FROM proteins")
            >>> print(f"Average time: {stats['avg_time']:.4f}s")
        """
        import time

        times = []
        for i in range(iterations):
            start_time = time.time()
            try:
                result = self.query_to_df(query)
                end_time = time.time()
                execution_time = end_time - start_time
                times.append(execution_time)
                self.logger.debug(f"Iteration {i+1}: {execution_time:.4f}s")
            except Exception as e:
                self.logger.error(f"Query failed on iteration {i+1}: {e}")
                return {"error": str(e)}

        if times:
            return {
                "query": query,
                "iterations": iterations,
                "min_time": min(times),
                "max_time": max(times),
                "avg_time": sum(times) / len(times),
                "std_time": (
                    sum((t - sum(times) / len(times)) ** 2 for t in times) / len(times)
                )
                ** 0.5,
            }
        else:
            return {"error": "No successful iterations"}

    def analyze_tables(self):
        """
        Run ANALYZE on all tables to update statistics for query optimization.

        This method updates table statistics that DuckDB uses for query planning
        and optimization. Running ANALYZE can significantly improve query performance
        by providing the query planner with accurate table statistics.
        """
        if not self._duckdb:
            self.logger.warning("Database not initialized")
            return

        tables = self._duckdb.execute("SHOW TABLES").fetchall()
        for table in tables:
            table_name = table[0]
            try:
                self._duckdb.execute(f"ANALYZE {table_name}")
                self.logger.info(f"Analyzed table {table_name}")
            except Exception as e:
                self.logger.warning(f"Failed to analyze table {table_name}: {e}")

    def vacuum_database(self):
        """
        Run VACUUM to reclaim storage and optimize the database.

        VACUUM reclaims storage by removing deleted rows and optimizing
        the database file structure. This can reduce file size and improve
        overall database performance.
        """
        if not self._duckdb:
            self.logger.warning("Database not initialized")
            return

        try:
            self._duckdb.execute("VACUUM")
            self.logger.info("Database vacuum completed")
        except Exception as e:
            self.logger.warning(f"Failed to vacuum database: {e}")


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
        database_name = create_uuid_filename("diann-report", ".duckdb")
        super().__init__(database_name, max_memory, worker_threads)
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
