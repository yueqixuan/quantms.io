"""MaxQuant Feature adapter -- evidence.txt to feature.parquet.

Loads ``evidence.txt`` into DuckDB, optionally merges SDRF for sample
mapping and TMT channels, transforms into ``FeatureSchema``, and writes
through ``FeatureWriter``.

Key schema changes:
    - ``reference_file_name`` -> ``run_file_name``
    - ``precursor_charge`` (int32) -> ``charge`` (int16)
    - ``is_decoy`` (int32, 0/1) -> ``is_decoy`` (bool)
    - ``scan`` (string) -> ``scan`` (list<int32>)
    - ``start_ion_mobility``/``stop_ion_mobility`` -> ``ion_mobility_start``/``ion_mobility_stop``
    - Intensities struct: ``{sample_accession, channel, intensity}`` ->
      ``{label, intensity}``
"""

from __future__ import annotations

import logging
import re
from typing import Optional

import pandas as pd

from qpx.converters.base import resolve_columns
from qpx.converters.maxquant.base_adapter import MaxQuantBaseAdapter
from qpx.converters.maxquant.constants import FIELD_MAPPINGS
from qpx.converters.maxquant.constants import to_proforma
from qpx.converters.ptm import from_proforma
from qpx.converters.utils import mq_flag_to_bool, safe_float
from qpx.writers.feature import FeatureWriter

logger = logging.getLogger(__name__)

# Derive field map from constants
_FEATURE_MAP = FIELD_MAPPINGS["feature"]


class MaxQuantFeatureAdapter(MaxQuantBaseAdapter):
    """Convert MaxQuant ``evidence.txt`` to ``feature.parquet``.

    Usage::

        with MaxQuantFeatureAdapter() as adapter:
            adapter.convert(
                evidence_path="evidence.txt",
                output_path="feature.parquet",
                sdrf_path="sdrf.tsv",
            )
    """

    def convert(
        self,
        evidence_path: str,
        output_path: str,
        sdrf_path: Optional[str] = None,
        chunksize: int = 500_000,
        creator: str = "maxquant",
    ) -> None:
        """Run the evidence.txt -> feature.parquet conversion.

        Args:
            evidence_path: Path to MaxQuant ``evidence.txt``.
            output_path: Destination Parquet path.
            sdrf_path: Optional SDRF file for sample/channel mapping.
            chunksize: Rows per batch.
            creator: Creator tag in Parquet metadata.
        """
        # Step 1: Load evidence into DuckDB
        self._load_evidence(evidence_path)

        # Step 2: Resolve column mappings against actual input columns
        actual_cols = {
            c[0]
            for c in self._conn.execute(
                "SELECT column_name FROM information_schema.columns "
                "WHERE table_name='evidence'"
            ).fetchall()
        }
        self._resolved = resolve_columns(_FEATURE_MAP, actual_cols)

        # Step 3: Load SDRF for sample mapping
        sample_map, experiment_type, tmt_channels = self._load_sdrf(sdrf_path)

        # Step 4: Stream and transform
        self.logger.info("Transforming MaxQuant features ...")

        with FeatureWriter(
            output_path, creator=creator, compression=self._compression
        ) as writer:
            for batch in self._query_batched("SELECT * FROM evidence", chunksize):
                df = batch.to_pandas()
                records = self._transform_batch(
                    df, sample_map, experiment_type, tmt_channels
                )
                if records:
                    self._track_scores(records)
                    writer.write_batch(records)

        self.logger.info(f"MaxQuant feature conversion complete -> {output_path}")

    # ------------------------------------------------------------------
    # Data loading
    # ------------------------------------------------------------------

    def _load_evidence(self, path: str) -> None:
        """Load evidence.txt into DuckDB."""
        self._conn.execute(f"""
            CREATE TABLE evidence AS
            SELECT * FROM read_csv_auto('{path}',
                delim='\\t', header=true, auto_detect=true)
            """)
        count = self._conn.execute("SELECT COUNT(*) FROM evidence").fetchone()[0]
        self.logger.info(f"Loaded {count:,} MaxQuant evidence rows")

    # ------------------------------------------------------------------
    # Transform
    # ------------------------------------------------------------------

    def _transform_batch(
        self,
        df: pd.DataFrame,
        sample_map: dict,
        experiment_type: str,
        tmt_channels: list[str],
    ) -> list[dict]:
        """Transform a batch of evidence.txt rows into QPX feature records."""
        records: list[dict] = []
        skipped = 0
        # Pre-extract column arrays for faster per-row access than to_dict("records")
        col_arrays = {col: df[col].values for col in df.columns}
        n_rows = len(df)
        for i in range(n_rows):
            try:
                row = {col: vals[i] for col, vals in col_arrays.items()}
                rec = self._transform_row(
                    row, sample_map, experiment_type, tmt_channels
                )
                if rec:
                    records.append(rec)
            except Exception as e:
                skipped += 1
                self.logger.debug(f"Skipping MaxQuant feature row: {e}")
        if skipped:
            total = skipped + len(records)
            self.logger.warning(
                "Skipped %d / %d rows (%.1f%%) in batch",
                skipped,
                total,
                100 * skipped / total if total else 0,
            )
        return records

    def _transform_row(
        self,
        row,
        sample_map: dict,
        experiment_type: str,
        tmt_channels: list[str],
    ) -> Optional[dict]:
        """Transform a single evidence.txt row."""
        r = self._resolved  # shorthand for resolved column mappings

        sequence = str(row.get(r.get("sequence", "Sequence"), ""))
        peptidoform = to_proforma(
            str(row.get(r.get("modified_sequence", "Modified sequence"), "")),
        )
        modifications = from_proforma(peptidoform, sequence) if peptidoform else None
        charge = int(row.get(r.get("charge", "Charge"), 0))
        run_file_name = str(row.get(r.get("run_file_name", "Raw file"), ""))

        # is_decoy (bool)
        is_decoy = mq_flag_to_bool(row.get(r.get("is_decoy", "Reverse"), ""))

        # Scan (list<int32>)
        scan_raw = row.get(r.get("scan", "MS/MS scan number"))
        scan = []
        if pd.notna(scan_raw):
            try:
                scan = [int(scan_raw)]
            except (ValueError, TypeError):
                pass

        # Calculated m/z from neutral mass and charge
        # calculated_mz = (mass + charge * proton_mass) / charge
        _PROTON_MASS = 1.00727646677
        neutral_mass = safe_float(row.get(r.get("mass", "Mass")))
        if neutral_mass and charge:
            calculated_mz = (neutral_mass + charge * _PROTON_MASS) / charge
        else:
            calculated_mz = 0.0

        # m/z
        observed_mz = safe_float(row.get(r.get("observed_mz", "m/z"))) or 0.0

        # RT
        rt = safe_float(row.get(r.get("rt", "Calibrated retention time")))
        rt_start = safe_float(
            row.get(r.get("rt_start", "Calibrated retention time start"))
        )
        rt_stop = safe_float(
            row.get(r.get("rt_stop", "Calibrated retention time finish"))
        )

        # PEP
        pep = safe_float(row.get(r.get("posterior_error_probability", "PEP")))

        # Ion mobility
        ion_mobility = safe_float(row.get(r.get("ion_mobility", "1/K0")))

        # Protein accessions (list of pg_protein structs)
        pg_acc_raw = str(row.get(r.get("pg_accessions", "Leading proteins"), ""))
        pg_acc_list = pg_acc_raw.split(";") if pg_acc_raw else []
        anchor_protein = str(
            row.get(r.get("anchor_protein", "Leading razor protein"), "")
        )
        if not anchor_protein and pg_acc_list:
            anchor_protein = pg_acc_list[0]
        pg_accessions = [
            {"accession": acc, "start": None, "end": None, "pre": None, "post": None}
            for acc in pg_acc_list
        ] or None

        # Unique peptide indicator
        unique = len(pg_acc_list) <= 1

        # Gene names
        gg_raw = row.get(r.get("gg_names", "Gene names"))
        gg_names = str(gg_raw).split(";") if pd.notna(gg_raw) and gg_raw else None

        # Mass error (ppm) — direct column from evidence.txt
        mass_error_ppm = safe_float(
            row.get(r.get("mass_error_ppm", "Mass error [ppm]"))
        )

        # Missed cleavages (dedicated field from evidence.txt)
        mc_raw = row.get(r.get("missed_cleavages", "Missed cleavages"))
        missed_cleavages = int(mc_raw) if pd.notna(mc_raw) else None

        # MBR detection: Type == "MULTI-MATCH" means transferred (MBR) feature
        evidence_type = str(row.get("Type", "")).strip().upper()
        id_run_file_name = None if evidence_type == "MULTI-MATCH" else run_file_name

        # Build intensities (new schema: {label, intensity})
        intensities, additional_intensities = self._build_intensities(
            row, run_file_name, sample_map, experiment_type, tmt_channels
        )

        # Additional scores
        additional_scores = []
        andromeda = safe_float(row.get(r.get("andromeda_score", "Score")))
        if andromeda is not None:
            additional_scores.append(
                {
                    "score_name": "andromeda_score",
                    "score_value": andromeda,
                    "higher_better": True,
                }
            )
        delta = safe_float(row.get(r.get("andromeda_delta_score", "Delta score")))
        if delta is not None:
            additional_scores.append(
                {
                    "score_name": "andromeda_delta_score",
                    "score_value": delta,
                    "higher_better": True,
                }
            )

        return {
            "sequence": sequence,
            "peptidoform": peptidoform,
            "modifications": modifications,
            "charge": charge,
            "posterior_error_probability": pep,
            "is_decoy": is_decoy,
            "calculated_mz": calculated_mz,
            "observed_mz": observed_mz,
            "mass_error_ppm": mass_error_ppm,
            "additional_scores": additional_scores or None,
            "predicted_rt": None,
            "run_file_name": run_file_name,
            "cv_params": None,
            "scan": scan,
            "rt": rt,
            "ion_mobility": ion_mobility,
            "missed_cleavages": missed_cleavages,
            "intensities": intensities or None,
            "additional_intensities": additional_intensities,
            "pg_accessions": pg_accessions,
            "anchor_protein": anchor_protein,
            "unique": unique,
            "pg_global_qvalue": None,
            "ion_mobility_start": None,
            "ion_mobility_stop": None,
            "gg_accessions": gg_names,
            "gg_names": gg_names,
            "id_run_file_name": id_run_file_name,
            "rt_start": rt_start,
            "rt_stop": rt_stop,
        }

    def _build_intensities(
        self,
        row,
        run_file_name: str,
        sample_map: dict,
        experiment_type: str,
        tmt_channels: list[str],
    ) -> tuple[list[dict], Optional[list[dict]]]:
        """Build intensities and additional_intensities from evidence row."""
        intensities: list[dict] = []
        additional_intensities: list[dict] = []

        exp_type = re.sub(r"\d+", "", experiment_type).upper()

        if exp_type in ("TMT", "ITRAQ") and tmt_channels:
            # TMT/iTRAQ: one intensity per channel
            for i, channel_name in enumerate(tmt_channels):
                # Try reporter intensity columns
                for col_name in [
                    f"Reporter intensity {i}",
                    f"Reporter intensity {i+1}",
                ]:
                    val = row.get(col_name)
                    if pd.notna(val) and float(val) > 0:
                        intensities.append(
                            {"label": channel_name, "intensity": float(val)}
                        )

                        # Corrected intensity
                        corr_col = col_name.replace(
                            "Reporter intensity", "Reporter intensity corrected"
                        )
                        corr_val = row.get(corr_col)
                        if pd.notna(corr_val):
                            additional_intensities.append(
                                {
                                    "label": channel_name,
                                    "intensities": [
                                        {
                                            "intensity_name": "corrected_reporter_intensity",
                                            "intensity_value": float(corr_val),
                                        }
                                    ],
                                }
                            )
                        break
        else:
            # LFQ: single intensity
            int_col = _FEATURE_MAP["intensity"][0]
            intensity_val = safe_float(row.get(int_col)) or 0.0
            label = "LFQ"
            intensities.append({"label": label, "intensity": float(intensity_val)})

        return intensities, additional_intensities or None
