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
from qpx.converters.diann.constants import FIELD_MAPPINGS
from qpx.converters.diann.constants import to_modifications, to_proforma
from qpx.converters.ptm import compute_precursor_mz
from qpx.converters.utils import safe_float
from qpx.writers.feature import FeatureWriter

logger = logging.getLogger(__name__)

# Derive report columns from FIELD_MAPPINGS
_DIANN_REPORT_COLS = list({
    candidates[0]
    for candidates in FIELD_MAPPINGS["feature"].values()
})

# Fields to skip in the SQL SELECT (sourced elsewhere or duplicate mapping)
_FEATURE_SQL_SKIP = {"lfq", "observed_mz"}

# Alias overrides: constants key -> SQL alias expected by _build_feature_record
_FEATURE_ALIAS_OVERRIDES = {"lfq_maxlfq": "lfq"}


class DiannFeatureAdapter(BaseConverter):
    """Convert DIA-NN report to ``feature.parquet``.

    Usage::

        with DiannFeatureAdapter() as adapter:
            adapter.convert(
                diann_report="report.tsv",
                output_path="feature.parquet",
                mzml_info_folder="ms_info/",
                sdrf_path="sdrf.tsv",
            )
    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        # Cache for computed m/z values: (modified_seq, charge) -> float
        self._mz_cache: dict[tuple[str, int], float | None] = {}

    def convert(
        self,
        diann_report: str,
        output_path: str,
        mzml_info_folder: Optional[str] = None,
        sdrf_path: Optional[str] = None,
        qvalue_threshold: float = 0.01,
        file_num: int = 50,
        creator: str = "diann",
    ) -> None:
        """Run the DIA-NN report -> feature.parquet conversion.

        Args:
            diann_report: Path to the DIA-NN ``report.tsv``.
            output_path: Destination Parquet path.
            mzml_info_folder: Optional folder with ``*_ms_info.parquet`` files.
                When provided, scan numbers and observed m/z are merged from
                these files.  When absent, runs are discovered directly from
                the report and scan/observed_mz fields are left empty.
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

        # Step 3: Discover runs — prefer ms_info files, fall back to report
        if mzml_info_folder:
            info_files = list(Path(mzml_info_folder).glob("*_ms_info.parquet"))
        else:
            info_files = []

        if info_files:
            run_names = [f.stem.replace("_ms_info", "") for f in info_files]
            self.logger.info(
                f"Found {len(run_names)} MS info files — using them for scan/mz merge"
            )
        else:
            # Get unique run names directly from the DIA-NN report
            run_col = FIELD_MAPPINGS["feature"]["run_file_name"][0]
            rows = self._conn.execute(
                f'SELECT DISTINCT "{run_col}" FROM report ORDER BY "{run_col}"'
            ).fetchall()
            run_names = [
                r[0].replace(".mzML", "").replace(".raw", "")
                for r in rows
            ]
            self.logger.info(
                f"No MS info files — discovered {len(run_names)} runs from report"
            )

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
        mzml_info_folder: Optional[str],
        qvalue_threshold: float,
        sample_map: dict[str, str],
    ) -> list[dict]:
        """Process a batch of runs and return feature records."""
        # Build SQL SELECT clause from FIELD_MAPPINGS
        feature_map = FIELD_MAPPINGS["feature"]
        select_parts = []
        for qpx_field, candidates in feature_map.items():
            if qpx_field in _FEATURE_SQL_SKIP:
                continue
            col = candidates[0]
            alias = _FEATURE_ALIAS_OVERRIDES.get(qpx_field, qpx_field)
            select_parts.append(f'"{col}" AS {alias}')

        select_clause = ",\n                ".join(select_parts)

        # Use constants-derived column names for filtering
        run_col = feature_map["run_file_name"][0]
        qvalue_col = feature_map["qvalue"][0]
        pg_col = feature_map["pg_accessions"][0]

        placeholders = ", ".join(["?" for _ in runs])
        report_df = self._conn.execute(
            f"""
            SELECT
                {select_clause}
            FROM report
            WHERE "{run_col}" IN ({placeholders})
              AND "{qvalue_col}" < {qvalue_threshold}
              AND "{pg_col}" IS NOT NULL
            """,
            runs,
        ).df()

        if report_df.empty:
            return []

        # Strip file extensions from run names
        report_df["run_file_name"] = (
            report_df["run_file_name"]
            .astype(str)
            .str.replace(r"\.(mzML|raw|d)$", "", regex=True)
        )

        # Integrate scan info from ms_info files (when available)
        records: list[dict] = []
        for run_name in runs:
            run_data = report_df[report_df["run_file_name"] == run_name].copy()
            if run_data.empty:
                continue

            # Load scan info only if mzml_info_folder is available
            if mzml_info_folder:
                ms_info_path = next(
                    Path(mzml_info_folder).glob(f"*{run_name}_ms_info.parquet"),
                    None,
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
        modified_seq = str(row.get("modified_sequence", ""))
        charge = int(row.get("charge", 0))

        # Peptidoform — convert DIA-NN notation to ProForma
        peptidoform = to_proforma(modified_seq)

        # Modifications — parse structured mods from DIA-NN notation
        modifications = to_modifications(modified_seq, sequence)

        # Calculated m/z — compute from peptidoform + charge via pyOpenMS
        # Use cache on (modified_seq, charge) to avoid recomputation
        cache_key = (modified_seq, charge)
        if cache_key not in self._mz_cache:
            self._mz_cache[cache_key] = compute_precursor_mz(
                modified_seq, charge
            )
        calculated_mz = self._mz_cache[cache_key] or 0.0

        # Observed m/z — only available when mzML info was merged
        observed_mz = safe_float(row.get("observed_mz")) or 0.0

        # Scan — from MS2.Scan column (DIA-NN report) or from mzML merge
        scan = []
        ms2_scan = row.get("ms2_scan")
        if pd.notna(ms2_scan):
            try:
                scan = [int(ms2_scan)]
            except (ValueError, TypeError):
                pass
        # Fallback to scan from mzML merge
        if not scan:
            scan_raw = row.get("scan")
            if pd.notna(scan_raw):
                try:
                    scan = [int(scan_raw)]
                except (ValueError, TypeError):
                    pass

        # Protein accessions (list of pg_protein structs)
        pg_acc_raw = str(row.get("pg_accessions", ""))
        pg_acc_list = pg_acc_raw.split(";") if pg_acc_raw else []
        pg_accessions = [
            {"accession": acc, "start": None, "end": None} for acc in pg_acc_list
        ] or None
        anchor_protein = pg_acc_list[0] if pg_acc_list else ""

        # Multi-protein accessions (all mapped proteins)
        mp_acc_raw = row.get("mp_accessions")
        mp_accessions = (
            str(mp_acc_raw).split(";")
            if pd.notna(mp_acc_raw) and mp_acc_raw
            else None
        )

        # Is decoy (bool) — DIA-NN typically filters decoys from output
        is_decoy = False

        # Unique peptide — proteotypic if only one protein group
        unique = 0 if len(pg_acc_list) > 1 else 1

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

        # Gene names — DIA-NN uses semicolons for multi-gene entries
        gg_raw = row.get("gg_names")
        gg_names = (
            str(gg_raw).split(";") if pd.notna(gg_raw) and gg_raw else None
        )

        return {
            "sequence": sequence,
            "peptidoform": peptidoform,
            "modifications": modifications,
            "charge": charge,
            "posterior_error_probability": safe_float(row.get("posterior_error_probability")),
            "is_decoy": is_decoy,
            "calculated_mz": calculated_mz,
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
            "pg_accessions": pg_accessions,
            "anchor_protein": anchor_protein,
            "unique": unique,
            "pg_global_qvalue": safe_float(row.get("pg_global_qvalue")),
            "ion_mobility_start": None,
            "ion_mobility_stop": None,
            "gg_accessions": None,
            "gg_names": gg_names,
            "id_run_file_name": None,
            "rt_start": safe_float(row.get("rt_start")),
            "rt_stop": safe_float(row.get("rt_stop")),
        }
