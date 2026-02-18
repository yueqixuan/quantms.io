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

from qpx.converters.base import BaseConverter
from qpx.converters.utils import clean_peptidoform, mq_flag_to_bool, safe_float
from qpx.writers.feature import FeatureWriter

logger = logging.getLogger(__name__)


class MaxQuantFeatureAdapter(BaseConverter):
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

        # Step 2: Load SDRF for sample mapping
        sample_map, experiment_type, tmt_channels = self._load_sdrf(sdrf_path)

        # Step 3: Stream and transform
        self.logger.info("Transforming MaxQuant features ...")

        total = self._conn.execute("SELECT COUNT(*) FROM evidence").fetchone()[0]

        with FeatureWriter(output_path, creator=creator) as writer:
            offset = 0
            while offset < total:
                df = self._conn.execute(
                    f"SELECT * FROM evidence LIMIT {chunksize} OFFSET {offset}"
                ).df()
                if df.empty:
                    break
                records = self._transform_batch(
                    df, sample_map, experiment_type, tmt_channels
                )
                if records:
                    self._track_scores(records)
                    writer.write_batch(records)
                offset += chunksize

        self.logger.info(f"MaxQuant feature conversion complete -> {output_path}")

    # ------------------------------------------------------------------
    # Data loading
    # ------------------------------------------------------------------

    def _load_evidence(self, path: str) -> None:
        """Load evidence.txt into DuckDB."""
        self._conn.execute(
            f"""
            CREATE TABLE evidence AS
            SELECT * FROM read_csv_auto('{path}',
                delim='\\t', header=true, auto_detect=true)
            """
        )
        count = self._conn.execute("SELECT COUNT(*) FROM evidence").fetchone()[0]
        self.logger.info(f"Loaded {count:,} MaxQuant evidence rows")

    def _load_sdrf(
        self, sdrf_path: Optional[str]
    ) -> tuple[dict, str, list]:
        """Load SDRF and return (sample_map, experiment_type, tmt_channels)."""
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
        for row in df.to_dict("records"):
            try:
                rec = self._transform_row(
                    row, sample_map, experiment_type, tmt_channels
                )
                if rec:
                    records.append(rec)
            except Exception as e:
                self.logger.debug(f"Skipping MaxQuant feature row: {e}")
        return records

    def _transform_row(
        self,
        row,
        sample_map: dict,
        experiment_type: str,
        tmt_channels: list[str],
    ) -> Optional[dict]:
        """Transform a single evidence.txt row."""

        sequence = str(row.get("Sequence", ""))
        peptidoform = clean_peptidoform(str(row.get("Modified sequence", "")))
        charge = int(row.get("Charge", 0))
        run_file_name = str(row.get("Raw file", ""))

        # is_decoy (bool)
        is_decoy = mq_flag_to_bool(row.get("Reverse", ""))

        # Scan (list<int32>)
        scan_raw = row.get("MS/MS scan number")
        scan = []
        if pd.notna(scan_raw):
            try:
                scan = [int(scan_raw)]
            except (ValueError, TypeError):
                pass

        # m/z
        observed_mz = safe_float(row.get("m/z")) or 0.0

        # RT
        rt = safe_float(row.get("Calibrated retention time"))
        rt_start = safe_float(row.get("Calibrated retention time start"))
        rt_stop = safe_float(row.get("Calibrated retention time finish"))

        # PEP
        pep = safe_float(row.get("PEP"))

        # Ion mobility
        ion_mobility = safe_float(row.get("1/K0"))

        # Protein accessions
        pg_acc_raw = str(row.get("Leading proteins", ""))
        pg_accessions = pg_acc_raw.split(";") if pg_acc_raw else []
        anchor_protein = str(row.get("Leading razor protein", ""))
        if not anchor_protein and pg_accessions:
            anchor_protein = pg_accessions[0]

        # Unique peptide indicator
        unique = 0 if len(pg_accessions) > 1 else 1

        # Gene names
        gg_raw = row.get("Gene names")
        gg_names = str(gg_raw).split(";") if pd.notna(gg_raw) and gg_raw else None

        # Build intensities (new schema: {label, intensity})
        intensities, additional_intensities = self._build_intensities(
            row, run_file_name, sample_map, experiment_type, tmt_channels
        )

        # Additional scores
        additional_scores = []
        andromeda = safe_float(row.get("Score"))
        if andromeda is not None:
            additional_scores.append(
                {"score_name": "andromeda_score", "score_value": andromeda, "higher_better": True}
            )
        delta = safe_float(row.get("Delta score"))
        if delta is not None:
            additional_scores.append(
                {"score_name": "andromeda_delta_score", "score_value": delta, "higher_better": True}
            )

        return {
            "sequence": sequence,
            "peptidoform": peptidoform,
            "modifications": None,
            "charge": charge,
            "posterior_error_probability": pep,
            "is_decoy": is_decoy,
            "calculated_mz": 0.0,
            "observed_mz": observed_mz,
            "additional_scores": additional_scores or None,
            "predicted_rt": None,
            "run_file_name": run_file_name,
            "cv_params": None,
            "scan": scan,
            "rt": rt,
            "ion_mobility": ion_mobility,
            "intensities": intensities or None,
            "additional_intensities": additional_intensities,
            "pg_accessions": pg_accessions or None,
            "anchor_protein": anchor_protein,
            "unique": unique,
            "pg_global_qvalue": None,
            "pg_positions": None,
            "ion_mobility_start": None,
            "ion_mobility_stop": None,
            "gg_accessions": None,
            "gg_names": gg_names,
            "id_run_file_name": None,
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
                for col_name in [f"Reporter intensity {i}", f"Reporter intensity {i+1}"]:
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
            intensity_val = safe_float(row.get("Intensity")) or 0.0
            label = "LFQ"
            intensities.append({"label": label, "intensity": float(intensity_val)})

        return intensities, additional_intensities or None
