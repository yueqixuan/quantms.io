"""DIA-NN PG adapter -- report.tsv + pg_matrix.tsv to pg.parquet.

Loads a DIA-NN ``report.tsv`` and optionally the ``pg_matrix.tsv`` into
DuckDB, aggregates protein-group intensities via SQL, and writes
through ``PgWriter``.

Key schema changes:
    - ``reference_file_name`` -> ``run_file_name``
    - ``is_decoy`` (int, 0/1) -> ``is_decoy`` (bool)
    - Intensities struct: ``{sample_accession, channel, intensity}`` ->
      ``{label, intensity}``
"""

from __future__ import annotations

import logging
from typing import Optional

import pandas as pd

from qpx.converters.base import BaseConverter
from qpx.converters.diann.constants import FIELD_MAPPINGS
from qpx.converters.utils import safe_float
from qpx.writers.pg import PgWriter

logger = logging.getLogger(__name__)

# Extra columns needed for PG aggregation but not in FIELD_MAPPINGS
_PG_EXTRA_COLS = [
    ('"Proteotypic"', "proteotypic"),
    ('"Stripped.Sequence"', "stripped_sequence"),
    ('"Precursor.Id"', "precursor_id"),
]


class DiannPgAdapter(BaseConverter):
    """Convert DIA-NN report + PG matrix to ``pg.parquet``.

    Usage::

        with DiannPgAdapter() as adapter:
            adapter.convert(
                diann_report="report.tsv",
                pg_matrix_path="pg_matrix.tsv",
                output_path="pg.parquet",
                sdrf_path="sdrf.tsv",
            )
    """

    def convert(
        self,
        diann_report: str,
        pg_matrix_path: str,
        output_path: str,
        sdrf_path: Optional[str] = None,
        file_num: int = 20,
        creator: str = "diann",
    ) -> None:
        """Run the DIA-NN report+matrix -> pg.parquet conversion.

        Args:
            diann_report: Path to the DIA-NN ``report.tsv``.
            pg_matrix_path: Path to the DIA-NN ``pg_matrix.tsv``.
            output_path: Destination Parquet path.
            sdrf_path: Optional SDRF file for sample mapping.
            file_num: Number of runs to process per batch.
            creator: Creator tag in Parquet metadata.
        """
        # Step 1: Load report into DuckDB
        self._load_diann_report(diann_report)

        # Step 2: Load PG matrix
        pg_matrix = self._load_pg_matrix(pg_matrix_path)

        # Step 3: Load SDRF mapping
        sample_map: dict[str, str] = {}
        if sdrf_path:
            from qpx.core.sdrf import SDRFHandler
            handler = SDRFHandler(sdrf_path)
            sample_map = handler.get_sample_map_run()

        # Step 4: Get unique runs and process in batches
        runs = self._get_unique_runs()

        with PgWriter(output_path, creator=creator) as writer:
            for i in range(0, len(runs), file_num):
                batch_runs = runs[i : i + file_num]
                self.logger.info(
                    f"Processing PG runs {i+1}-{min(i+file_num, len(runs))} of {len(runs)}"
                )
                records = self._process_batch(
                    batch_runs, pg_matrix, sample_map
                )
                if records:
                    self._track_scores(records)
                    writer.write_batch(records)

        self.logger.info(f"DIA-NN PG conversion complete -> {output_path}")

    # ------------------------------------------------------------------
    # Data loading
    # ------------------------------------------------------------------

    def _load_diann_report(self, path: str) -> None:
        """Load DIA-NN report into DuckDB."""
        self._conn.execute(
            f"""
            CREATE TABLE report AS
            SELECT * FROM read_csv_auto('{path}',
                delim='\\t', header=true, auto_detect=true)
            """
        )
        count = self._conn.execute("SELECT COUNT(*) FROM report").fetchone()[0]
        self.logger.info(f"Loaded {count:,} DIA-NN report rows for PG")

    def _load_pg_matrix(self, path: str) -> pd.DataFrame:
        """Load the DIA-NN PG matrix TSV."""
        pg_map = FIELD_MAPPINGS["pg"]
        pg_col = pg_map["pg_accessions"][0]   # "Protein.Group"
        names_col = pg_map["pg_names"][0]      # "Protein.Names"
        genes_col = pg_map["gg_accessions"][0] # "Genes"

        header = pd.read_csv(path, sep="\t", nrows=0).columns.tolist()
        mzml_cols = [c for c in header if c.endswith(".mzML")]
        usecols = [pg_col, names_col, genes_col] + mzml_cols
        df = pd.read_csv(path, sep="\t", usecols=usecols)
        df.rename(
            columns={
                pg_col: "pg_accessions",
                names_col: "pg_names",
                genes_col: "gg_accessions",
            },
            inplace=True,
        )
        # Strip .mzML from column names
        df.columns = [c.replace(".mzML", "") for c in df.columns]
        return df

    def _get_unique_runs(self) -> list[str]:
        """Get sorted list of unique Run values from the report."""
        run_col = FIELD_MAPPINGS["pg"]["run_file_name"][0]
        rows = self._conn.execute(
            f'SELECT DISTINCT "{run_col}" FROM report ORDER BY "{run_col}"'
        ).fetchall()
        return [r[0] for r in rows]

    # ------------------------------------------------------------------
    # Batch processing
    # ------------------------------------------------------------------

    def _process_batch(
        self,
        runs: list[str],
        pg_matrix: pd.DataFrame,
        sample_map: dict[str, str],
    ) -> list[dict]:
        """Process a batch of runs for PG quantification."""
        records: list[dict] = []

        # Build SQL SELECT clause from FIELD_MAPPINGS
        pg_map = FIELD_MAPPINGS["pg"]
        select_parts = []
        for qpx_field, candidates in pg_map.items():
            col = candidates[0]
            select_parts.append(f'"{col}" AS {qpx_field}')
        # Add extra columns needed for aggregation
        for src, alias in _PG_EXTRA_COLS:
            select_parts.append(f'{src} AS {alias}')

        select_clause = ",\n                ".join(select_parts)

        # Use constants-derived column names for filtering
        run_col = pg_map["run_file_name"][0]
        pg_col = pg_map["pg_accessions"][0]

        placeholders = ", ".join(["?" for _ in runs])
        report_df = self._conn.execute(
            f"""
            SELECT
                {select_clause}
            FROM report
            WHERE "{run_col}" IN ({placeholders})
              AND "{pg_col}" IS NOT NULL
            """,
            runs,
        ).df()

        if report_df.empty:
            return []

        # Strip extension from run file names
        report_df["run_file_name"] = (
            report_df["run_file_name"].astype(str).str.replace(r"\.mzML$", "", regex=True)
        )

        for run_name in runs:
            run_name_clean = run_name.replace(".mzML", "")
            run_report = report_df[report_df["run_file_name"] == run_name_clean].copy()
            if run_report.empty:
                continue

            # Aggregate per protein group
            pg_groups = run_report.groupby(
                ["pg_accessions", "pg_names", "gg_accessions", "run_file_name"],
                dropna=False,
            )

            for (pg_acc, pg_nm, gg_acc, ref), group in pg_groups:
                rec = self._build_pg_record(
                    pg_acc=str(pg_acc),
                    pg_names_raw=pg_nm,
                    gg_acc_raw=gg_acc,
                    run_file_name=str(ref),
                    group=group,
                    pg_matrix=pg_matrix,
                    sample_map=sample_map,
                )
                if rec:
                    records.append(rec)

        return records

    def _build_pg_record(
        self,
        pg_acc: str,
        pg_names_raw,
        gg_acc_raw,
        run_file_name: str,
        group: pd.DataFrame,
        pg_matrix: pd.DataFrame,
        sample_map: dict[str, str],
    ) -> Optional[dict]:
        """Build a single PG record."""

        pg_accessions = pg_acc.split(";")
        pg_names = (
            str(pg_names_raw).split(";")
            if pd.notna(pg_names_raw) and pg_names_raw
            else None
        )
        gg_accessions = (
            str(gg_acc_raw).split(";")
            if pd.notna(gg_acc_raw) and gg_acc_raw
            else None
        )

        anchor_protein = pg_accessions[0] if pg_accessions else ""
        global_qvalue = safe_float(group["global_qvalue"].iloc[0])

        # Peptide/feature counts
        total_sequences = group["stripped_sequence"].nunique()
        proteotypic_mask = group["proteotypic"].astype(str).isin(["1", "1.0", "True"])
        unique_sequences = group.loc[proteotypic_mask, "stripped_sequence"].nunique()
        total_features = group["precursor_id"].nunique() if "precursor_id" in group.columns else total_sequences
        unique_features = (
            group.loc[proteotypic_mask, "precursor_id"].nunique()
            if "precursor_id" in group.columns
            else unique_sequences
        )

        # PG quantity from matrix
        pg_quantity = 0.0
        if run_file_name in pg_matrix.columns:
            match = pg_matrix[pg_matrix["pg_accessions"] == pg_acc]
            if not match.empty:
                pg_quantity = safe_float(match[run_file_name].iloc[0]) or 0.0

        # Intensities (new schema: {label, intensity})
        label = "LFQ"
        intensities = [{"label": label, "intensity": float(pg_quantity)}]

        # Additional intensities (LFQ)
        lfq_val = safe_float(group["lfq"].iloc[0])
        additional_intensities = None
        if lfq_val is not None:
            additional_intensities = [
                {
                    "label": label,
                    "intensities": [
                        {"intensity_name": "lfq", "intensity_value": float(lfq_val)},
                    ],
                }
            ]

        # Peptides per protein
        peptide_count = max(total_sequences, 1)
        peptides = [
            {"protein_name": acc, "peptide_count": peptide_count}
            for acc in pg_accessions
        ]

        # Additional scores
        qvalue_val = safe_float(group["qvalue"].iloc[0])
        additional_scores = []
        if qvalue_val is not None:
            additional_scores.append(
                {"score_name": "qvalue", "score_value": qvalue_val, "higher_better": False}
            )

        return {
            "pg_accessions": pg_accessions,
            "pg_names": pg_names,
            "gg_accessions": gg_accessions,
            "gg_names": None,
            "anchor_protein": anchor_protein,
            "run_file_name": run_file_name,
            "global_qvalue": global_qvalue,
            "pg_qvalue": None,
            "intensities": intensities,
            "additional_intensities": additional_intensities,
            "is_decoy": False,
            "contaminant": None,
            "peptides": peptides,
            "peptide_counts": {
                "unique_sequences": unique_sequences,
                "total_sequences": total_sequences,
            },
            "feature_counts": {
                "unique_features": unique_features,
                "total_features": total_features,
            },
            "sequence_coverage": None,
            "molecular_weight": None,
            "additional_scores": additional_scores or None,
            "cv_params": None,
        }
