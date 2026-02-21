"""QuantMS PG adapter -- aggregate feature.parquet to pg.parquet.

Reads a feature.parquet file (produced by QuantmsFeatureAdapter) and the
mzTab proteins table, aggregates feature-level intensities per
(anchor_protein, run_file_name) group, and writes the results through
``PgWriter``.
"""

from __future__ import annotations

import logging
import re
from typing import Optional

import pandas as pd

from qpx.converters.base import BaseConverter
from qpx.converters.utils import safe_float
from qpx.converters.mztab import load_mztab_sections
from qpx.writers.pg import PgWriter

logger = logging.getLogger(__name__)


class QuantmsPgAdapter(BaseConverter):
    """Convert QuantMS features + mzTab proteins to ``pg.parquet``.

    Reads the already-produced feature.parquet and the mzTab proteins
    table to create protein-group level quantification.

    Usage::

        with QuantmsPgAdapter() as adapter:
            adapter.convert(
                mztab_path="results.mzTab",
                feature_path="output.feature.parquet",
                output_path="output.pg.parquet",
            )
    """

    def convert(
        self,
        mztab_path: str,
        feature_path: str,
        output_path: str,
        creator: str = "quantms",
        query_batch_size: int = 50_000,
        write_batch_size: int = 10_000,
    ) -> None:
        """Run the feature.parquet + mzTab proteins -> pg.parquet conversion.

        Args:
            mztab_path: Path to the mzTab file (for protein group metadata).
            feature_path: Path to the feature.parquet file.
            output_path: Destination Parquet path.
            creator: Creator tag in Parquet metadata.
        """
        # Step 1: Load mzTab protein group metadata (two-pass: single + group)
        if not self._table_exists("proteins"):
            load_mztab_sections(self._conn, mztab_path)
        single_meta, group_meta = self._build_protein_meta()

        # Step 2: Stream sorted feature rows via DuckDB to avoid full materialization.
        feature_path_sql = feature_path.replace("'", "''")
        self.logger.info(f"Reading features from {feature_path}")
        total_features = self._conn.execute(
            f"SELECT COUNT(*) FROM read_parquet('{feature_path_sql}')"
        ).fetchone()[0]
        self.logger.info(f"Loaded {total_features} features")

        # Step 3: Aggregate features per (anchor_protein, run_file_name) while streaming.
        self.logger.info("Aggregating protein groups from features ...")
        grouped_sql = (
            "SELECT * FROM read_parquet('{path}') "
            "ORDER BY anchor_protein, run_file_name"
        ).format(path=feature_path_sql)

        records_buffer: list[dict] = []
        total_records = 0
        current_anchor: Optional[str] = None
        current_run: Optional[str] = None
        current_rows: list[dict] = []
        total_groups = 0
        processed_groups = 0
        skipped_groups = 0
        failed_groups = 0
        max_bad_group_ratio = 0.20
        min_groups_for_ratio_failure = 25

        def _flush_current_group() -> None:
            nonlocal current_rows, total_records
            nonlocal total_groups, processed_groups, skipped_groups, failed_groups
            if not current_rows or current_anchor is None or current_run is None:
                return
            total_groups += 1
            try:
                rec = self._build_single_pg(
                    anchor_protein=current_anchor,
                    run_file_name=current_run,
                    features=pd.DataFrame.from_records(current_rows),
                    single_meta=single_meta,
                    group_meta=group_meta,
                )
                if rec:
                    records_buffer.append(rec)
                    total_records += 1
                    processed_groups += 1
                else:
                    skipped_groups += 1
                    self.logger.debug(
                        "Skipping PG group (%s, %s): no record built",
                        current_anchor,
                        current_run,
                    )
            except Exception as e:
                failed_groups += 1
                self.logger.debug(
                    "Failed PG group (%s, %s): %s",
                    current_anchor,
                    current_run,
                    e,
                )
            finally:
                current_rows = []

        with PgWriter(output_path, creator=creator) as writer:
            for batch in self._query_batched(grouped_sql, batch_size=query_batch_size):
                batch_df = batch.to_pandas()
                for row in batch_df.to_dict("records"):
                    anchor = str(row["anchor_protein"])
                    run_file = str(row["run_file_name"])
                    if current_anchor is None:
                        current_anchor, current_run = anchor, run_file

                    if anchor != current_anchor or run_file != current_run:
                        _flush_current_group()
                        current_anchor, current_run = anchor, run_file

                    current_rows.append(row)

                    if len(records_buffer) >= write_batch_size:
                        self._track_scores(records_buffer)
                        writer.write_batch(records_buffer)
                        records_buffer.clear()

            _flush_current_group()
            if records_buffer:
                self._track_scores(records_buffer)
                writer.write_batch(records_buffer)

        bad_groups = skipped_groups + failed_groups
        self.logger.info(
            "PG group aggregation summary: total=%d processed=%d skipped=%d failed=%d",
            total_groups,
            processed_groups,
            skipped_groups,
            failed_groups,
        )
        if bad_groups > 0:
            self.logger.warning(
                "PG conversion encountered problematic groups: %d/%d "
                "(skipped=%d, failed=%d)",
                bad_groups,
                total_groups,
                skipped_groups,
                failed_groups,
            )
        if total_groups > 0 and processed_groups == 0:
            raise ValueError(
                "PG conversion failed: all protein groups were skipped or failed "
                f"(total={total_groups}, skipped={skipped_groups}, failed={failed_groups})."
            )
        if (
            total_groups >= min_groups_for_ratio_failure
            and total_groups > 0
            and (bad_groups / total_groups) > max_bad_group_ratio
        ):
            raise ValueError(
                "PG conversion failed: too many problematic protein groups "
                f"({bad_groups}/{total_groups}, skipped={skipped_groups}, failed={failed_groups}, "
                f"threshold={max_bad_group_ratio:.0%})."
            )

        self.logger.info(
            "PG conversion complete -> %s (%d records, streamed)",
            output_path,
            total_records,
        )

    # ------------------------------------------------------------------
    # Protein group metadata
    # ------------------------------------------------------------------

    def _build_protein_meta(
        self,
    ) -> tuple[dict[str, dict], dict[str, dict]]:
        """Build protein metadata from mzTab proteins table (two-pass).

        Pass 1 — ``single_protein`` rows: keyed by individual accession.
        Pass 2 — ``indistinguishable_protein_group`` rows: keyed by
        sorted, semicolon-joined accessions for order-independent lookup.

        Returns:
            (single_meta, group_meta) tuple of dicts.
        """
        single_meta: dict[str, dict] = {}
        group_meta: dict[str, dict] = {}
        try:
            df = self._conn.execute("SELECT * FROM proteins").df()
            for row in df.to_dict("records"):
                result_type = str(row.get("opt_global_result_type", "single_protein"))
                accession = str(row.get("accession", ""))
                if not accession or accession == "null":
                    continue
                if result_type not in (
                    "single_protein",
                    "indistinguishable_protein_group",
                ):
                    continue

                entry = self._parse_protein_row(row, accession)

                if result_type == "single_protein":
                    single_meta[accession] = entry
                else:
                    # Normalize key: sort accessions for order-independent match
                    sorted_key = ";".join(
                        sorted(self._normalize_pg_accessions(accession))
                    )
                    if not sorted_key:
                        continue
                    group_meta[sorted_key] = entry
        except Exception as e:
            self.logger.error(f"Error loading protein metadata: {e}")

        return single_meta, group_meta

    def _parse_protein_row(self, row: dict, accession: str) -> dict:
        """Extract metadata fields from a single mzTab protein row."""
        description = row.get("description", "")
        description = (
            description if pd.notna(description) and description != "null" else None
        )

        # Gene names from description
        gg_accessions = []
        if description:
            gn = re.search(r"GN=([^\s]+)", str(description))
            if gn:
                gg_accessions = [gn.group(1)]

        # Global q-value
        global_qvalue = safe_float(
            row.get(
                "best_search_engine_score[1]",
                row.get("best_search_engine_score_1"),
            )
        )

        # Decoy
        is_decoy_raw = row.get(
            "opt_global_cv_pride:0000303_decoy_hit",
            row.get("opt_global_cv_PRIDE:0000303_decoy_hit", 0),
        )
        is_decoy = str(is_decoy_raw).strip() in ("1", "true", "True")

        # Sequence coverage
        seq_cov = safe_float(row.get("protein_coverage")) or 0.0

        # Protein names
        pg_names = self._extract_protein_names(accession)

        return {
            "pg_names": pg_names,
            "gg_accessions": gg_accessions or None,
            "global_qvalue": global_qvalue,
            "is_decoy": is_decoy,
            "sequence_coverage": seq_cov,
        }

    # ------------------------------------------------------------------
    # Record building
    # ------------------------------------------------------------------

    def _build_single_pg(
        self,
        anchor_protein: str,
        run_file_name: str,
        features: pd.DataFrame,
        single_meta: dict[str, dict],
        group_meta: dict[str, dict],
    ) -> Optional[dict]:
        """Build a single PG record from aggregated features.

        Merges metadata from both single_protein and group rows:
        - q-value from single_protein (per-protein score)
        - decoy/coverage from whichever has them (prefer single_protein)
        """

        # Protein group accessions from first feature's pg_accessions.
        # Support list<dict>, list[str], semicolon-delimited strings, and null values.
        first_pg = features["pg_accessions"].iloc[0]
        pg_accessions = self._normalize_pg_accessions(first_pg)
        if not pg_accessions:
            self.logger.warning(
                "Empty/invalid pg_accessions for anchor=%s run=%s; "
                "falling back to anchor_protein",
                anchor_protein,
                run_file_name,
            )
            pg_accessions = [anchor_protein]

        # Lookup: sorted key for group, anchor for single protein
        sorted_key = ";".join(sorted(pg_accessions))
        gmeta = group_meta.get(sorted_key, {})
        smeta = single_meta.get(anchor_protein, {})

        # Merge: prefer single_protein, fall back to group
        pg_names_str = smeta.get("pg_names") or gmeta.get("pg_names")
        pg_names = pg_names_str.split(";") if pg_names_str else None

        gg_accessions = smeta.get("gg_accessions") or gmeta.get("gg_accessions")
        global_qvalue = (
            smeta.get("global_qvalue")
            if smeta.get("global_qvalue") is not None
            else gmeta.get("global_qvalue")
        )
        is_decoy = smeta.get("is_decoy", gmeta.get("is_decoy", False))
        seq_coverage = (
            smeta.get("sequence_coverage")
            if smeta.get("sequence_coverage") is not None
            else gmeta.get("sequence_coverage")
        )

        # Fallback: if neither dict had data, use feature-level info
        if not smeta and not gmeta:
            is_decoy_raw = features["is_decoy"].iloc[0]
            is_decoy = bool(is_decoy_raw) if pd.notna(is_decoy_raw) else False
            if "gg_accessions" in features.columns:
                first_gg = features["gg_accessions"].iloc[0]
                if isinstance(first_gg, list) and len(first_gg) > 0:
                    gg_accessions = first_gg

        # Aggregate intensities: sum per label across all features
        intensities = self._aggregate_intensities(features)

        # Peptide and feature counts
        unique_sequences = features["sequence"].nunique()
        total_sequences = unique_sequences

        feature_keys = set(zip(features["peptidoform"], features["charge"]))
        unique_features = len(feature_keys)
        total_features = len(features)

        # Peptides per protein
        peptides = [
            {"protein_name": acc, "peptide_count": unique_sequences}
            for acc in pg_accessions
        ]

        # Additional scores
        additional_scores = []
        if seq_coverage is not None:
            additional_scores.append(
                {
                    "score_name": "sequence_coverage_percent",
                    "score_value": seq_coverage,
                    "higher_better": True,
                }
            )
        additional_scores.append(
            {
                "score_name": "peptide_count",
                "score_value": float(unique_sequences),
                "higher_better": True,
            }
        )

        return {
            "pg_accessions": pg_accessions,
            "pg_names": pg_names,
            "gg_accessions": (
                gg_accessions if isinstance(gg_accessions, list) else None
            ),
            "gg_names": None,
            "anchor_protein": anchor_protein,
            "run_file_name": run_file_name,
            "global_qvalue": global_qvalue,
            "pg_qvalue": None,
            "intensities": intensities or None,
            "additional_intensities": None,
            "is_decoy": is_decoy,
            "contaminant": False,
            "peptides": peptides,
            "peptide_counts": {
                "unique_sequences": unique_sequences,
                "total_sequences": total_sequences,
            },
            "feature_counts": {
                "unique_features": unique_features,
                "total_features": total_features,
            },
            "sequence_coverage": seq_coverage,
            "molecular_weight": None,
            "additional_scores": additional_scores or None,
            "cv_params": None,
        }

    def _aggregate_intensities(self, features: pd.DataFrame) -> list[dict]:
        """Sum feature intensities per label.

        Feature intensities are stored as list<intensity> where
        intensity = {label: string, intensity: float32}.
        """
        label_sums: dict[str, float] = {}
        for intensities in features["intensities"]:
            if intensities is None:
                continue
            for entry in intensities:
                label = entry.get("label", "LFQ")
                value = float(entry.get("intensity", 0.0) or 0.0)
                label_sums[label] = label_sums.get(label, 0.0) + value

        return [
            {"label": label, "intensity": total}
            for label, total in sorted(label_sums.items())
        ]

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _normalize_pg_accessions(raw_accessions) -> list[str]:
        """Normalize mixed pg_accessions shapes to a clean list of strings."""
        normalized: list[str] = []

        def _add(value) -> None:
            if value is None:
                return
            if pd.isna(value):
                return
            text = str(value).strip()
            if not text or text.lower() == "null":
                return
            normalized.append(text)

        items = raw_accessions if isinstance(raw_accessions, list) else [raw_accessions]
        for item in items:
            value = item.get("accession") if isinstance(item, dict) else item
            if isinstance(value, str):
                for token in value.split(";"):
                    _add(token)
            else:
                _add(value)
        return normalized

    @staticmethod
    def _extract_protein_names(accession: str) -> str:
        """Extract protein names from UniProt-style accession."""
        parts_list = []
        for acc in accession.split(";"):
            if "|" in acc:
                pieces = acc.split("|")
                if len(pieces) >= 3:
                    parts_list.append(pieces[2])
                elif len(pieces) >= 2:
                    parts_list.append(pieces[1])
                else:
                    parts_list.append(acc)
            else:
                parts_list.append(acc)
        return ";".join(parts_list)
