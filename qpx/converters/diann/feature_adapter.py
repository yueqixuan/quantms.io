"""DIA-NN Feature adapter -- report.tsv to feature.parquet.

Loads a DIA-NN ``report.tsv`` into DuckDB, joins with mzML scan info,
transforms via SQL into ``FeatureSchema``, and writes through ``FeatureWriter``.

Key schema changes handled:
    - ``reference_file_name`` -> ``run_file_name``
    - ``precursor_charge`` (int32) -> ``charge`` (int16)
    - ``is_decoy`` (0/1 int) -> ``is_decoy`` (bool)
    - ``scan`` (string) -> ``scan`` (list<int32>)
    - ``start_ion_mobility``/``stop_ion_mobility`` -> ``ion_mobility_start``/``ion_mobility_stop``
    - Intensities struct: ``{sample_accession, channel, intensity}`` ->
      ``{label, intensity}``
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional

import pandas as pd

from qpx.converters.base import BaseConverter
from qpx.converters.utils import safe_float
from qpx.writers.feature import FeatureWriter

logger = logging.getLogger(__name__)

# DIA-NN report columns to load (matching DIANN_MAP from core/common.py)
_DIANN_REPORT_COLS = [
    "Precursor.Quantity",
    "RT.Start",
    "RT.Stop",
    "RT",
    "Predicted.RT",
    "Protein.Group",
    "Protein.Ids",
    "PEP",
    "Global.Q.Value",
    "Global.PG.Q.Value",
    "Q.Value",
    "PG.Q.Value",
    "Precursor.Normalised",
    "PG.MaxLFQ",
    "Quantity.Quality",
    "Precursor.Charge",
    "Stripped.Sequence",
    "Modified.Sequence",
    "Genes",
    "Run",
]


class DiannFeatureAdapter(BaseConverter):
    """Convert DIA-NN report to ``feature.parquet``.

    Usage::

        with DiannFeatureAdapter() as adapter:
            adapter.convert(
                diann_report="report.tsv",
                mzml_info_folder="ms_info/",
                output_path="feature.parquet",
                sdrf_path="sdrf.tsv",
            )
    """

    def convert(
        self,
        diann_report: str,
        mzml_info_folder: str,
        output_path: str,
        sdrf_path: Optional[str] = None,
        qvalue_threshold: float = 0.01,
        file_num: int = 50,
        creator: str = "diann",
    ) -> None:
        """Run the DIA-NN report -> feature.parquet conversion.

        Args:
            diann_report: Path to the DIA-NN ``report.tsv``.
            mzml_info_folder: Folder containing ``*_ms_info.parquet`` files.
            output_path: Destination Parquet path.
            sdrf_path: Optional SDRF file for sample mapping.
            qvalue_threshold: Filter threshold for Q.Value.
            file_num: Number of runs to process per batch.
            creator: Creator tag in Parquet metadata.
        """
        # Step 1: Load report into DuckDB
        self._load_diann_report(diann_report)

        # Step 2: Load SDRF for sample mapping
        sample_map: dict[str, str] = {}
        if sdrf_path:
            sample_map = self._load_sdrf_sample_map(sdrf_path)

        # Step 3: Get available MS info files
        info_files = list(Path(mzml_info_folder).glob("*_ms_info.parquet"))
        run_names = [f.stem.replace("_ms_info", "") for f in info_files]

        self.logger.info(f"Found {len(run_names)} MS info files for {len(run_names)} runs")

        # Step 4: Process in batches
        with FeatureWriter(output_path, creator=creator) as writer:
            for i in range(0, len(run_names), file_num):
                batch_runs = run_names[i : i + file_num]
                self.logger.info(
                    f"Processing runs {i+1}-{min(i+file_num, len(run_names))} of {len(run_names)}"
                )
                records = self._process_batch(
                    batch_runs, mzml_info_folder, qvalue_threshold, sample_map
                )
                if records:
                    self._track_scores(records)
                    writer.write_batch(records)

        self.logger.info(f"DIA-NN feature conversion complete -> {output_path}")

    # ------------------------------------------------------------------
    # Data loading
    # ------------------------------------------------------------------

    def _load_diann_report(self, path: str) -> None:
        """Load DIA-NN report TSV into DuckDB ``report`` table."""
        self._conn.execute(
            f"""
            CREATE TABLE report AS
            SELECT * FROM read_csv_auto('{path}',
                delim='\\t', header=true, auto_detect=true)
            """
        )
        count = self._conn.execute("SELECT COUNT(*) FROM report").fetchone()[0]
        self.logger.info(f"Loaded {count:,} DIA-NN report rows")

    def _load_sdrf_sample_map(self, sdrf_path: str) -> dict[str, str]:
        """Load SDRF and create run_file -> sample_accession mapping."""
        from qpx.core.sdrf import SDRFHandler

        handler = SDRFHandler(sdrf_path)
        return handler.get_sample_map_run()

    # ------------------------------------------------------------------
    # Batch processing
    # ------------------------------------------------------------------

    def _process_batch(
        self,
        runs: list[str],
        mzml_info_folder: str,
        qvalue_threshold: float,
        sample_map: dict[str, str],
    ) -> list[dict]:
        """Process a batch of runs and return feature records."""
        # Fetch and filter report data for these runs
        placeholders = ", ".join(["?" for _ in runs])
        report_df = self._conn.execute(
            f"""
            SELECT
                "Precursor.Quantity" AS intensity,
                "RT.Start" AS rt_start,
                "RT.Stop" AS rt_stop,
                "RT" AS rt,
                "Predicted.RT" AS predicted_rt,
                "Protein.Group" AS pg_accessions,
                "Protein.Ids" AS mp_accessions,
                "PEP" AS posterior_error_probability,
                "Global.Q.Value" AS global_qvalue,
                "Global.PG.Q.Value" AS pg_global_qvalue,
                "Q.Value" AS qvalue,
                "PG.Q.Value" AS pg_qvalue,
                "Precursor.Normalised" AS normalize_intensity,
                "PG.MaxLFQ" AS lfq,
                "Quantity.Quality" AS precursor_quantification_score,
                "Precursor.Charge" AS charge,
                "Stripped.Sequence" AS sequence,
                "Modified.Sequence" AS peptidoform,
                "Genes" AS gg_names,
                "Run" AS run_file_name
            FROM report
            WHERE "Run" IN ({placeholders})
              AND "Q.Value" < {qvalue_threshold}
              AND "Protein.Group" IS NOT NULL
            """,
            runs,
        ).df()

        if report_df.empty:
            return []

        # Strip .mzML extension from run names
        report_df["run_file_name"] = (
            report_df["run_file_name"].astype(str).str.replace(r"\.mzML$", "", regex=True)
        )

        # Integrate scan info from ms_info files
        records: list[dict] = []
        for run_name in runs:
            run_data = report_df[report_df["run_file_name"] == run_name].copy()
            if run_data.empty:
                continue

            # Load scan info
            ms_info_path = next(
                Path(mzml_info_folder).glob(f"*{run_name}_ms_info.parquet"), None
            )
            if ms_info_path:
                run_data = self._merge_scan_info(run_data, ms_info_path)

            # Build feature records
            for row in run_data.to_dict("records"):
                rec = self._build_feature_record(row, sample_map)
                if rec:
                    records.append(rec)

        return records

    def _merge_scan_info(self, run_data: pd.DataFrame, ms_info_path: Path) -> pd.DataFrame:
        """Merge MS scan info (scan, observed_mz) with report data."""
        target = pd.read_parquet(
            ms_info_path, columns=["rt", "scan", "precursor_mz"]
        )
        target = target.rename(columns={"precursor_mz": "observed_mz"})
        target["rt"] = target["rt"] / 60  # Convert to minutes

        run_data = run_data.sort_values("rt")
        result = pd.merge_asof(
            run_data, target, on="rt", direction="nearest", suffixes=("", "_scan")
        )
        return result

    # ------------------------------------------------------------------
    # Record building
    # ------------------------------------------------------------------

    def _build_feature_record(self, row, sample_map: dict) -> Optional[dict]:
        """Build a single feature record from a DIA-NN report row."""
        run_file_name = str(row.get("run_file_name", ""))
        sequence = str(row.get("sequence", ""))
        peptidoform = str(row.get("peptidoform", ""))
        charge = int(row.get("charge", 0))

        # Scan (list<int32>)
        scan_raw = row.get("scan")
        if pd.notna(scan_raw):
            try:
                scan = [int(scan_raw)]
            except (ValueError, TypeError):
                scan = []
        else:
            scan = []

        # Protein accessions
        pg_acc_raw = str(row.get("pg_accessions", ""))
        pg_accessions = pg_acc_raw.split(";") if pg_acc_raw else []
        anchor_protein = pg_accessions[0] if pg_accessions else ""

        # Is decoy (bool) -- DIA-NN typically does not have decoys in output
        is_decoy = False

        # Unique peptide
        unique = 0 if ";" in pg_acc_raw else 1

        # Intensities (new schema: {label, intensity})
        intensity_val = safe_float(row.get("intensity")) or 0.0
        label = "LFQ"
        intensities = [{"label": label, "intensity": float(intensity_val)}]

        # Additional intensities
        lfq_val = safe_float(row.get("lfq"))
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

        # Additional scores
        additional_scores = []
        for score_name, col_name, lower_better in [
            ("qvalue", "qvalue", True),
            ("pg_qvalue", "pg_qvalue", True),
            ("global_qvalue", "global_qvalue", True),
        ]:
            val = safe_float(row.get(col_name))
            if val is not None:
                additional_scores.append(
                    {
                        "score_name": score_name,
                        "score_value": val,
                        "higher_better": not lower_better,
                    }
                )

        # CV params
        pqs = row.get("precursor_quantification_score")
        cv_params = None
        if pd.notna(pqs):
            cv_params = [
                {
                    "cv_name": "precursor_quantification_score",
                    "cv_value": str(pqs),
                }
            ]

        # Gene names
        gg_raw = row.get("gg_names")
        gg_names = (
            str(gg_raw).split(",") if pd.notna(gg_raw) and gg_raw else None
        )

        # Observed m/z (from scan merge or None)
        observed_mz = safe_float(row.get("observed_mz")) or 0.0

        return {
            "sequence": sequence,
            "peptidoform": peptidoform,
            "modifications": None,
            "charge": charge,
            "posterior_error_probability": safe_float(row.get("posterior_error_probability")),
            "is_decoy": is_decoy,
            "calculated_mz": 0.0,  # Needs PyOpenMS to compute; left as 0
            "observed_mz": observed_mz,
            "additional_scores": additional_scores or None,
            "predicted_rt": safe_float(row.get("predicted_rt")),
            "run_file_name": run_file_name,
            "cv_params": cv_params,
            "scan": scan,
            "rt": safe_float(row.get("rt")),
            "ion_mobility": None,
            "intensities": intensities,
            "additional_intensities": additional_intensities,
            "pg_accessions": pg_accessions or None,
            "anchor_protein": anchor_protein,
            "unique": unique,
            "pg_global_qvalue": safe_float(row.get("pg_global_qvalue")),
            "pg_positions": None,
            "ion_mobility_start": None,
            "ion_mobility_stop": None,
            "gg_accessions": None,
            "gg_names": gg_names,
            "id_run_file_name": None,
            "rt_start": safe_float(row.get("rt_start")),
            "rt_stop": safe_float(row.get("rt_stop")),
        }
