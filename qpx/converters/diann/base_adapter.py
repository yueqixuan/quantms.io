"""Shared base adapter for DIA-NN converters.

Provides common data-loading helpers used by both
:class:`DiannFeatureAdapter` and :class:`DiannPgAdapter`.
"""

from __future__ import annotations

import pandas as pd

from qpx.converters.base import BaseConverter
from qpx.converters.channel_labels import normalize_label, read_sdrf_labels
from qpx.core.sql import escape_path, sql_build, validate_identifier


class DiaNNBaseAdapter(BaseConverter):
    """Base adapter with DIA-NN-specific loading utilities.

    Subclasses inherit ``_load_diann_report()`` so the logic is defined
    once rather than duplicated across the feature and PG adapters.
    """

    def _load_diann_report(self, path: str) -> None:
        """Create a DuckDB view over a DIA-NN report file.

        Uses ``CREATE VIEW`` so DuckDB reads from the file lazily with
        column pruning and predicate pushdown — avoiding loading the
        entire dataset into memory (critical for multi-GB reports).

        Supports both TSV (``report.tsv``) and Parquet (``report.parquet``)
        formats.  If the ``report`` view already exists it is skipped.

        Args:
            path: Filesystem path to the DIA-NN report file.
        """
        if self._table_exists("report"):
            self.logger.debug("report view already loaded -- skipping reload")
            return
        safe_path = escape_path(path)
        if path.endswith(".parquet"):
            self._conn.execute(
                sql_build(
                    "CREATE VIEW report AS SELECT * FROM read_parquet('$path')",
                    path=safe_path,
                )
            )
        else:
            self._conn.execute(
                sql_build(
                    """CREATE VIEW report AS
                SELECT * FROM read_csv_auto('$path',
                    delim='\t', header=true, auto_detect=true,
                    null_padding=true)""",
                    path=safe_path,
                )
            )
        count = self._conn.execute("SELECT COUNT(*) FROM report").fetchone()[0]
        self.logger.info(f"DIA-NN report view created ({count:,} rows)")

    def _register_channel_labels(
        self,
        report_columns: set[str],
        sdrf_path: str | None,
    ) -> str | None:
        """Register canonical DIA-NN channel labels and return the source column."""
        channel_col = next((col for col in ("Channel", "Label") if col in report_columns), None)
        if channel_col is None:
            return None

        quoted_channel = validate_identifier(channel_col)
        raw_channels = [
            row[0]
            for row in self._conn.execute(
                sql_build(
                    """
                    SELECT DISTINCT CAST($channel AS VARCHAR)
                    FROM report
                    ORDER BY CAST($channel AS VARCHAR)
                    """,
                    channel=quoted_channel,
                )
            ).fetchall()
        ]
        if any(channel is None or not str(channel).strip() for channel in raw_channels):
            raise ValueError(f"DIA-NN {channel_col} contains an empty channel identifier")

        declared_labels = read_sdrf_labels(sdrf_path)
        declared_by_key = {normalize_label(label).casefold(): normalize_label(label) for label in declared_labels or ()}

        rows: list[tuple[str, str]] = []
        canonical_to_raw: dict[str, str] = {}
        for raw_channel in raw_channels:
            raw_text = str(raw_channel).strip()
            normalized = normalize_label(raw_text)
            if declared_by_key:
                key = normalized.casefold()
                if key not in declared_by_key:
                    raise ValueError(
                        f"DIA-NN channel {raw_text!r} is not declared in SDRF "
                        "comment[label]; channel identifiers must match so "
                        "intensities can resolve to run samples"
                    )
                normalized = declared_by_key[key]

            canonical_key = normalized.casefold()
            previous = canonical_to_raw.get(canonical_key)
            if previous is not None and previous != raw_text:
                raise ValueError(f"DIA-NN channels {previous!r} and {raw_text!r} both canonicalize to {normalized!r}")
            canonical_to_raw[canonical_key] = raw_text
            rows.append((raw_text, normalized))

        self._conn.execute("DROP TABLE IF EXISTS diann_channel_labels")
        self._conn.from_df(
            pd.DataFrame(
                rows,
                columns=["channel_value", "canonical_label"],
            )
        ).create("diann_channel_labels")
        return channel_col
