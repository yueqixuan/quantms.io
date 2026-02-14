"""DuckDB connection management, configuration, and resource settings."""

import duckdb
from pathlib import Path


class DuckDBEngine:
    """
    Manages a DuckDB connection with configurable resources.

    Used by Dataset (shared connection) and standalone structures
    (own connection). All DuckDB tuning lives here.
    """

    def __init__(
        self,
        memory_limit: str = "16GB",
        threads: int = 4,
        temp_directory: str | None = None,
    ):
        self._conn = duckdb.connect(":memory:")
        self._conn.execute(f"SET memory_limit='{memory_limit}'")
        self._conn.execute(f"SET threads={threads}")
        if temp_directory:
            self._conn.execute(f"SET temp_directory='{temp_directory}'")
        self._conn.execute("SET enable_progress_bar=true")

    @property
    def connection(self) -> duckdb.DuckDBPyConnection:
        return self._conn

    def register_parquet(self, name: str, file_path: str | Path) -> None:
        """Register a Parquet file as a lazy DuckDB view."""
        self._conn.execute(
            f"CREATE OR REPLACE VIEW {name} AS "
            f"SELECT * FROM read_parquet('{file_path}')"
        )

    def register_partitioned_parquet(self, name: str, directory: str | Path) -> None:
        """Register a Hive-partitioned Parquet directory as a DuckDB view."""
        glob_pattern = str(Path(directory) / "**" / "*.parquet")
        self._conn.execute(
            f"CREATE OR REPLACE VIEW {name} AS "
            f"SELECT * FROM read_parquet('{glob_pattern}', hive_partitioning=true)"
        )

    def execute(self, sql: str, params=None):
        if params:
            return self._conn.execute(sql, params)
        return self._conn.execute(sql)

    def close(self):
        self._conn.close()

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()


def create_engine(**kwargs) -> DuckDBEngine:
    """Factory for creating a DuckDB engine with default settings."""
    return DuckDBEngine(**kwargs)
