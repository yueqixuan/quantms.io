"""BaseConverter -- abstract base for all tool-specific converters.

Each converter creates its own DuckDB in-memory connection, loads
tool-specific output (TSV, mzTab, etc.), transforms via SQL into the
QPX schema, and pipes rows into the appropriate Writer.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
import logging
from pathlib import Path

import duckdb

from qpx.core.scores import score_ontology_entries


def resolve_columns(
    field_mappings: dict[str, list[str]],
    available_columns: set[str],
) -> dict[str, str]:
    """Resolve QPX fields to actual tool column names.

    For each QPX field, tries candidates in order and returns the first match.
    Skips fields where no candidate matches the available columns.

    Args:
        field_mappings: QPX field name -> ordered list of candidate column names.
        available_columns: Set of column names present in the input data.

    Returns:
        Dict mapping QPX field name -> resolved tool column name.
    """
    resolved = {}
    for qpx_field, candidates in field_mappings.items():
        for candidate in candidates:
            if candidate in available_columns:
                resolved[qpx_field] = candidate
                break
    return resolved


class BaseConverter(ABC):
    """Abstract base class for converter adapters.

    Subclasses implement ``convert()`` which performs:
        1. Load tool output into DuckDB tables
        2. SQL transform into QPX schema columns
        3. Stream results into a Writer (FeatureWriter, PsmWriter, etc.)

    After conversion, call ``write_score_ontology()`` to emit ontology
    entries for all discovered score names.
    """

    def __init__(
        self,
        duckdb_memory: str = "16GB",
        duckdb_threads: int = 4,
        conn: duckdb.DuckDBPyConnection | None = None,
    ):
        self.logger = logging.getLogger(
            f"{self.__class__.__module__}.{self.__class__.__name__}"
        )
        if conn is not None:
            self._conn = conn
            self._owns_conn = False
        else:
            self._conn = duckdb.connect(":memory:")
            self._conn.execute(f"SET memory_limit='{duckdb_memory}'")
            self._conn.execute(f"SET threads={duckdb_threads}")
            self._owns_conn = True
        # Accumulate discovered score names during conversion
        self._discovered_scores: set[str] = set()

    # ------------------------------------------------------------------
    # Abstract interface
    # ------------------------------------------------------------------

    @abstractmethod
    def convert(self, **kwargs) -> None:
        """Run the full conversion pipeline and write output files."""
        ...

    # ------------------------------------------------------------------
    # Helpers shared by all converters
    # ------------------------------------------------------------------

    def _table_exists(self, name: str) -> bool:
        """Return True if *name* is already registered as a DuckDB table."""
        tables = {t[0] for t in self._conn.execute("SHOW TABLES").fetchall()}
        return name in tables

    def _load_tsv(self, table_name: str, path: str, **kwargs) -> None:
        """Load a TSV file into a DuckDB table.

        Args:
            table_name: Target table name in DuckDB.
            path: Path to the TSV file.
            **kwargs: Extra ``read_csv_auto`` options (e.g. ``columns``).
        """
        opts = ", delim='\\t', header=true, auto_detect=true"
        for k, v in kwargs.items():
            if isinstance(v, str):
                opts += f", {k}='{v}'"
            elif isinstance(v, list):
                cols = ", ".join(f"'{c}'" for c in v)
                opts += f", {k}=[{cols}]"
            else:
                opts += f", {k}={v}"
        sql = f"CREATE OR REPLACE TABLE {table_name} AS SELECT * FROM read_csv_auto('{path}'{opts})"
        self._conn.execute(sql)
        count = self._conn.execute(f"SELECT COUNT(*) FROM {table_name}").fetchone()[0]
        self.logger.info(f"Loaded {count:,} rows into '{table_name}' from {path}")

    def _load_parquet(self, table_name: str, path: str) -> None:
        """Load a Parquet file into a DuckDB table."""
        sql = f"CREATE OR REPLACE TABLE {table_name} AS SELECT * FROM read_parquet('{path}')"
        self._conn.execute(sql)
        count = self._conn.execute(f"SELECT COUNT(*) FROM {table_name}").fetchone()[0]
        self.logger.info(f"Loaded {count:,} rows into '{table_name}' from {path}")

    def _query_batched(self, sql: str, batch_size: int = 500_000):
        """Execute *sql* and yield result batches as Arrow record batches.

        This avoids materialising the full result set in Python.
        """
        result = self._conn.execute(sql)
        while True:
            batch = result.fetch_record_batch(batch_size)
            if batch.num_rows == 0:
                break
            yield batch

    def _query_to_arrow(self, sql: str):
        """Execute *sql* and return the full result as a PyArrow Table."""
        return self._conn.execute(sql).fetch_arrow_table()

    # ------------------------------------------------------------------
    # Score tracking and ontology writing
    # ------------------------------------------------------------------

    def _track_scores(self, records: list[dict]) -> None:
        """Collect score names from a batch of records into ``_discovered_scores``."""
        for rec in records:
            scores = rec.get("additional_scores")
            if scores:
                for s in scores:
                    name = s.get("score_name")
                    if name:
                        self._discovered_scores.add(name)

    def write_score_ontology(
        self,
        output_path: str | Path,
        view: str = "psm",
        extra_entries: list[dict] | None = None,
    ) -> Path | None:
        """Write ontology.parquet with CV term entries for discovered scores.

        Args:
            output_path: Path for the ontology Parquet file.
            view: QPX view name (``"psm"``, ``"feature"``, ``"pg"``).
            extra_entries: Additional ontology entries to include.

        Returns:
            Path to the written file, or ``None`` if no entries to write.
        """
        entries = score_ontology_entries(self._discovered_scores, view=view)
        if extra_entries:
            entries.extend(extra_entries)
        if not entries:
            return None

        from qpx.writers.ontology import OntologyWriter

        output_path = Path(output_path)
        with OntologyWriter(output_path, creator="qpx") as writer:
            writer.write_batch(entries)
        self.logger.info(
            "Wrote %d ontology entries (%d scores) to %s",
            len(entries),
            len(self._discovered_scores),
            output_path,
        )
        return output_path

    # ------------------------------------------------------------------
    # Lifecycle
    # ------------------------------------------------------------------

    def close(self):
        """Close the DuckDB connection (only if we created it)."""
        if self._conn and self._owns_conn:
            try:
                self._conn.close()
            except Exception:
                pass
            self._conn = None

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()
