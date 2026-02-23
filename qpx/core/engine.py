"""DuckDB connection management, configuration, and resource settings."""

from __future__ import annotations

import logging
import re

import duckdb
from pathlib import Path

logger = logging.getLogger(__name__)

_SAFE_SET_VALUE = re.compile(r"^[\w./:@\-]+$")


def _validate_set_value(value: str, name: str) -> str:
    """Validate a value used in a DuckDB SET statement."""
    if not _SAFE_SET_VALUE.match(value):
        raise ValueError(f"Invalid character in DuckDB config '{name}': {value!r}")
    return value


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
        s3_config: dict | None = None,
    ):
        self._conn = duckdb.connect(":memory:")
        self._conn.execute(
            f"SET memory_limit='{_validate_set_value(memory_limit, 'memory_limit')}'"
        )
        self._conn.execute(f"SET threads={int(threads)}")
        if temp_directory:
            self._conn.execute(
                f"SET temp_directory='{_validate_set_value(temp_directory, 'temp_directory')}'"
            )
        self._conn.execute("SET enable_progress_bar=true")
        if s3_config is not None:
            self._setup_s3(s3_config)

    def _setup_s3(self, config: dict) -> None:
        """Configure DuckDB httpfs for S3 access."""
        self._conn.execute("INSTALL httpfs")
        self._conn.execute("LOAD httpfs")
        if "region" in config:
            self._conn.execute(
                f"SET s3_region='{_validate_set_value(config['region'], 's3_region')}'"
            )
        if "access_key_id" in config and "secret_access_key" in config:
            self._conn.execute(
                f"SET s3_access_key_id='{_validate_set_value(config['access_key_id'], 's3_access_key_id')}'"
            )
            self._conn.execute(
                f"SET s3_secret_access_key='{_validate_set_value(config['secret_access_key'], 's3_secret_access_key')}'"
            )
        if "endpoint" in config:
            self._conn.execute(
                f"SET s3_endpoint='{_validate_set_value(config['endpoint'], 's3_endpoint')}'"
            )
        if config.get("anonymous", False):
            self._conn.execute("SET s3_url_style='path'")

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

    def register_s3_parquet(self, name: str, s3_path: str) -> None:
        """Register an S3 Parquet file as a DuckDB view."""
        self._conn.execute(
            f"CREATE OR REPLACE VIEW {name} AS "
            f"SELECT * FROM read_parquet('{s3_path}')"
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


def create_converter_connection(
    memory_limit: str = "16GB",
    threads: int = 4,
) -> duckdb.DuckDBPyConnection:
    """Create a validated DuckDB in-memory connection for converter use.

    Unifies DuckDB configuration across converters. Uses the same validation
    as DuckDBEngine to prevent injection via SET values.
    """
    conn = duckdb.connect(":memory:")
    conn.execute(
        f"SET memory_limit='{_validate_set_value(memory_limit, 'memory_limit')}'"
    )
    conn.execute(f"SET threads={int(threads)}")
    conn.execute("SET enable_progress_bar=true")
    return conn
