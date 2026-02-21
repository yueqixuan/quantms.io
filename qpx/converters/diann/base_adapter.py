"""Shared base adapter for DIA-NN converters.

Provides common data-loading helpers used by both
:class:`DiannFeatureAdapter` and :class:`DiannPgAdapter`.
"""

from __future__ import annotations

from qpx.converters.base import BaseConverter


class DiaNNBaseAdapter(BaseConverter):
    """Base adapter with DIA-NN-specific loading utilities.

    Subclasses inherit ``_load_diann_report()`` so the logic is defined
    once rather than duplicated across the feature and PG adapters.
    """

    def _load_diann_report(self, path: str) -> None:
        """Load a DIA-NN ``report.tsv`` into the DuckDB ``report`` table.

        Uses DuckDB ``read_csv_auto`` with tab-delimited auto-detection.
        If the ``report`` table already exists it is replaced.

        Args:
            path: Filesystem path to the DIA-NN ``report.tsv``.
        """
        if self._table_exists("report"):
            self.logger.debug("report table already loaded -- skipping reload")
            return
        self._conn.execute(
            f"""
            CREATE TABLE report AS
            SELECT * FROM read_csv_auto('{path}',
                delim='\\t', header=true, auto_detect=true)
            """
        )
        count = self._conn.execute("SELECT COUNT(*) FROM report").fetchone()[0]
        self.logger.info(f"Loaded {count:,} DIA-NN report rows")
