"""MaxQuant PG adapter -- proteinGroups.txt to pg.parquet.

Loads ``proteinGroups.txt`` into DuckDB, transforms protein-group rows
into ``PgSchema``, and writes through ``PgWriter``.

Key schema changes:
    - ``reference_file_name`` -> ``run_file_name``
    - ``is_decoy`` (int32, 0/1) -> ``is_decoy`` (bool)
    - Intensities struct: ``{sample_accession, channel, intensity}`` ->
      ``{label, intensity}``
"""

from __future__ import annotations

import logging
from typing import Optional

import pandas as pd

from qpx.converters.base import BaseConverter
from qpx.converters.utils import mq_flag_to_bool, safe_float
from qpx.writers.pg import PgWriter

logger = logging.getLogger(__name__)


class MaxQuantPgAdapter(BaseConverter):
    """Convert MaxQuant ``proteinGroups.txt`` to ``pg.parquet``.

    Usage::

        with MaxQuantPgAdapter() as adapter:
            adapter.convert(
                protein_groups_path="proteinGroups.txt",
                output_path="pg.parquet",
                sdrf_path="sdrf.tsv",
            )
    """

    def convert(
        self,
        protein_groups_path: str,
        output_path: str,
        sdrf_path: Optional[str] = None,
        chunksize: int = 100_000,
        creator: str = "maxquant",
    ) -> None:
        """Run the proteinGroups.txt -> pg.parquet conversion.

        Args:
            protein_groups_path: Path to MaxQuant ``proteinGroups.txt``.
            output_path: Destination Parquet path.
            sdrf_path: Optional SDRF for sample mapping.
            chunksize: Rows per batch.
            creator: Creator tag in Parquet metadata.
        """
        # Step 1: Load proteinGroups.txt into DuckDB
        self._load_protein_groups(protein_groups_path)

        # Step 2: Load SDRF mapping
        sample_map, experiment_type, tmt_channels = self._load_sdrf(sdrf_path)

        # Step 3: Detect intensity columns in the data
        intensity_cols = self._detect_intensity_columns(experiment_type)

        # Step 4: Stream and transform
        self.logger.info("Transforming MaxQuant protein groups ...")

        total = self._conn.execute("SELECT COUNT(*) FROM protein_groups").fetchone()[0]

        with PgWriter(output_path, creator=creator) as writer:
            offset = 0
            while offset < total:
                df = self._conn.execute(
                    f"SELECT * FROM protein_groups LIMIT {chunksize} OFFSET {offset}"
                ).df()
                if df.empty:
                    break
                records = self._transform_batch(
                    df, sample_map, experiment_type, intensity_cols
                )
                if records:
                    self._track_scores(records)
                    writer.write_batch(records)
                offset += chunksize

        self.logger.info(f"MaxQuant PG conversion complete -> {output_path}")

    # ------------------------------------------------------------------
    # Data loading
    # ------------------------------------------------------------------

    def _load_protein_groups(self, path: str) -> None:
        """Load proteinGroups.txt into DuckDB."""
        self._conn.execute(
            f"""
            CREATE TABLE protein_groups AS
            SELECT * FROM read_csv_auto('{path}',
                delim='\\t', header=true, auto_detect=true)
            """
        )
        count = self._conn.execute("SELECT COUNT(*) FROM protein_groups").fetchone()[0]
        self.logger.info(f"Loaded {count:,} MaxQuant protein groups")

    def _load_sdrf(
        self, sdrf_path: Optional[str]
    ) -> tuple[dict, str, list]:
        if not sdrf_path:
            return {}, "LFQ", []

        from qpx.core.sdrf import SDRFHandler
        handler = SDRFHandler(sdrf_path)
        sample_map = handler.get_sample_map_run()
        experiment_type = handler.get_experiment_type_from_sdrf()

        tmt_channels: list[str] = []
        if experiment_type and "TMT" in experiment_type.upper():
            labels = handler.sdrf_table.get("comment[label]")
            if labels is not None:
                tmt_labels = [l for l in labels.unique() if l and "TMT" in str(l).upper()]
                tmt_channels = sorted(tmt_labels)

        return sample_map, experiment_type or "LFQ", tmt_channels

    def _detect_intensity_columns(self, experiment_type: str) -> dict[str, list[str]]:
        """Detect sample-specific intensity columns in proteinGroups.txt."""
        cols = self._conn.execute(
            "SELECT column_name FROM information_schema.columns WHERE table_name='protein_groups'"
        ).fetchall()
        col_names = [c[0] for c in cols]

        # Find Intensity <sample> columns
        intensity_cols = [c for c in col_names if c.startswith("Intensity ") and c != "Intensity"]
        lfq_cols = [c for c in col_names if c.startswith("LFQ intensity ")]
        ibaq_cols = [c for c in col_names if c.startswith("iBAQ ") and c != "iBAQ"]

        return {
            "intensity": intensity_cols,
            "lfq": lfq_cols,
            "ibaq": ibaq_cols,
        }

    # ------------------------------------------------------------------
    # Transform
    # ------------------------------------------------------------------

    def _transform_batch(
        self,
        df: pd.DataFrame,
        sample_map: dict,
        experiment_type: str,
        intensity_cols: dict[str, list[str]],
    ) -> list[dict]:
        records: list[dict] = []
        for row in df.to_dict("records"):
            try:
                recs = self._transform_row(
                    row, sample_map, experiment_type, intensity_cols
                )
                records.extend(recs)
            except Exception as e:
                self.logger.debug(f"Skipping MaxQuant PG row: {e}")
        return records

    def _transform_row(
        self,
        row,
        sample_map: dict,
        experiment_type: str,
        intensity_cols: dict[str, list[str]],
    ) -> list[dict]:
        """Transform a single proteinGroups.txt row.

        MaxQuant PG rows contain intensities for ALL runs in one row, so
        we expand them into one QPX PG record per run.
        """
        records: list[dict] = []

        # Protein group identity
        pg_acc_raw = str(row.get("Protein IDs", ""))
        pg_accessions = pg_acc_raw.split(";") if pg_acc_raw else []

        pg_names_raw = row.get("Protein names")
        pg_names = (
            str(pg_names_raw).split(";")
            if pd.notna(pg_names_raw) and pg_names_raw
            else None
        )

        gg_raw = row.get("Gene names")
        gg_accessions = (
            str(gg_raw).split(";")
            if pd.notna(gg_raw) and gg_raw
            else None
        )

        # Anchor protein (first in Majority protein IDs, or first PG accession)
        majority_raw = row.get("Majority protein IDs")
        if pd.notna(majority_raw) and majority_raw:
            anchor_protein = str(majority_raw).split(";")[0]
        elif pg_accessions:
            anchor_protein = pg_accessions[0]
        else:
            anchor_protein = ""

        # Quality metrics
        global_qvalue = safe_float(row.get("Q-value"))
        is_decoy = mq_flag_to_bool(row.get("Reverse", ""))
        contaminant_val = 1 if mq_flag_to_bool(row.get("Potential contaminant", "")) else 0

        # Sequence coverage and molecular weight
        seq_coverage = safe_float(row.get("Sequence coverage [%]"))
        mol_weight = safe_float(row.get("Mol. weight [kDa]"))

        # Peptide counts
        peptide_count_total = int(row.get("Peptides", 0) or 0)
        peptide_count_unique = int(row.get("Unique peptides", 0) or 0)
        peptide_count_razor = int(row.get("Razor + unique peptides", 0) or 0)

        # Additional scores
        andromeda = safe_float(row.get("Score"))
        additional_scores = []
        if andromeda is not None:
            additional_scores.append(
                {"score_name": "andromeda_score", "score_value": andromeda, "higher_better": True}
            )

        # Peptides per protein
        peptides = [
            {"protein_name": acc, "peptide_count": peptide_count_total}
            for acc in pg_accessions
        ]

        # Expand per run using intensity columns
        for intensity_col in intensity_cols.get("intensity", []):
            # Extract run name from column "Intensity <run_name>"
            run_name = intensity_col.replace("Intensity ", "")
            intensity_val = safe_float(row.get(intensity_col)) or 0.0

            if intensity_val <= 0:
                continue

            # Intensities (new schema: {label, intensity})
            label = "LFQ"
            intensities = [{"label": label, "intensity": float(intensity_val)}]

            # Additional intensities (LFQ, iBAQ)
            additional_intensities = []
            lfq_col = f"LFQ intensity {run_name}"
            lfq_val = safe_float(row.get(lfq_col))
            ibaq_col = f"iBAQ {run_name}"
            ibaq_val = safe_float(row.get(ibaq_col))

            extra_vals = []
            if lfq_val is not None:
                extra_vals.append({"intensity_name": "lfq", "intensity_value": float(lfq_val)})
            if ibaq_val is not None:
                extra_vals.append({"intensity_name": "ibaq", "intensity_value": float(ibaq_val)})

            if extra_vals:
                additional_intensities.append(
                    {"label": label, "intensities": extra_vals}
                )

            rec = {
                "pg_accessions": pg_accessions,
                "pg_names": pg_names,
                "gg_accessions": gg_accessions,
                "gg_names": None,
                "anchor_protein": anchor_protein,
                "run_file_name": run_name,
                "global_qvalue": global_qvalue,
                "pg_qvalue": None,
                "intensities": intensities,
                "additional_intensities": additional_intensities or None,
                "is_decoy": is_decoy,
                "contaminant": contaminant_val,
                "peptides": peptides,
                "peptide_counts": {
                    "unique_sequences": peptide_count_unique,
                    "total_sequences": peptide_count_total,
                },
                "feature_counts": {
                    "unique_features": peptide_count_unique,
                    "total_features": peptide_count_total,
                },
                "sequence_coverage": seq_coverage,
                "molecular_weight": mol_weight,
                "additional_scores": additional_scores or None,
                "cv_params": None,
            }
            records.append(rec)

        # If no per-sample intensity columns, emit one record with total Intensity
        if not records:
            total_intensity = safe_float(row.get("Intensity")) or 0.0
            if total_intensity > 0:
                records.append(
                    {
                        "pg_accessions": pg_accessions,
                        "pg_names": pg_names,
                        "gg_accessions": gg_accessions,
                        "gg_names": None,
                        "anchor_protein": anchor_protein,
                        "run_file_name": "unknown",
                        "global_qvalue": global_qvalue,
                        "pg_qvalue": None,
                        "intensities": [{"label": "LFQ", "intensity": float(total_intensity)}],
                        "additional_intensities": None,
                        "is_decoy": is_decoy,
                        "contaminant": contaminant_val,
                        "peptides": peptides,
                        "peptide_counts": {
                            "unique_sequences": peptide_count_unique,
                            "total_sequences": peptide_count_total,
                        },
                        "feature_counts": {
                            "unique_features": peptide_count_unique,
                            "total_features": peptide_count_total,
                        },
                        "sequence_coverage": seq_coverage,
                        "molecular_weight": mol_weight,
                        "additional_scores": additional_scores or None,
                        "cv_params": None,
                    }
                )

        return records
