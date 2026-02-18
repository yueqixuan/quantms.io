"""MaxQuant PSM adapter -- msms.txt to psm.parquet.

Loads ``msms.txt`` into DuckDB, transforms via SQL into ``PsmSchema``,
and streams results through ``PsmWriter``.

Key schema changes:
    - ``reference_file_name`` -> ``run_file_name``  (Raw file column)
    - ``precursor_charge`` (int32) -> ``charge`` (int16)
    - ``is_decoy`` (int32, 0/1) -> ``is_decoy`` (bool)
    - ``scan`` (string) -> ``scan`` (list<int32>)
"""

from __future__ import annotations

import logging
from typing import Optional

import pandas as pd

from qpx.converters.base import BaseConverter
from qpx.converters.utils import clean_peptidoform, mq_flag_to_bool, safe_float
from qpx.writers.psm import PsmWriter

logger = logging.getLogger(__name__)

# MaxQuant columns to read from msms.txt (matching MAXQUANT_PSM_MAP)
_MQ_PSM_USECOLS = [
    "Sequence",
    "Proteins",
    "PEP",
    "Modified sequence",
    "Reverse",
    "m/z",
    "Scan number",
    "Retention time",
    "Charge",
    "Raw file",
    "Score",
    "Delta score",
    "PIF",
    "Masses",
    "Intensities",
    "Number of matches",
    "1/K0",
]

# Optional spectral data columns
_MQ_PSM_SPECTRAL_COLS = [
    "Matches",
    "Fragmentation",
    "Mass analyzer",
    "Type",
]


class MaxQuantPsmAdapter(BaseConverter):
    """Convert MaxQuant ``msms.txt`` to ``psm.parquet``.

    Usage::

        with MaxQuantPsmAdapter() as adapter:
            adapter.convert(
                msms_path="msms.txt",
                output_path="psm.parquet",
            )
    """

    def convert(
        self,
        msms_path: str,
        output_path: str,
        chunksize: int = 500_000,
        spectral_data: bool = False,
        creator: str = "maxquant",
    ) -> None:
        """Run the msms.txt -> psm.parquet conversion.

        Args:
            msms_path: Path to MaxQuant ``msms.txt``.
            output_path: Destination Parquet path.
            chunksize: Rows per batch.
            spectral_data: Whether to include spectral arrays.
            creator: Creator tag in Parquet metadata.
        """
        # Step 1: Load msms.txt into DuckDB
        self._load_msms(msms_path)

        # Step 2: Stream and transform
        self.logger.info("Transforming MaxQuant PSMs ...")

        total = self._conn.execute("SELECT COUNT(*) FROM msms").fetchone()[0]

        with PsmWriter(output_path, creator=creator) as writer:
            offset = 0
            while offset < total:
                df = self._conn.execute(
                    f"SELECT * FROM msms LIMIT {chunksize} OFFSET {offset}"
                ).df()
                if df.empty:
                    break
                records = self._transform_batch(df, spectral_data)
                if records:
                    self._track_scores(records)
                    writer.write_batch(records)
                offset += chunksize

        self.logger.info(f"MaxQuant PSM conversion complete -> {output_path}")

    # ------------------------------------------------------------------
    # Data loading
    # ------------------------------------------------------------------

    def _load_msms(self, path: str) -> None:
        """Load msms.txt into DuckDB."""
        self._conn.execute(
            f"""
            CREATE TABLE msms AS
            SELECT * FROM read_csv_auto('{path}',
                delim='\\t', header=true, auto_detect=true)
            """
        )
        count = self._conn.execute("SELECT COUNT(*) FROM msms").fetchone()[0]
        self.logger.info(f"Loaded {count:,} MaxQuant PSM rows from msms.txt")

    # ------------------------------------------------------------------
    # Transform
    # ------------------------------------------------------------------

    def _transform_batch(
        self, df: pd.DataFrame, spectral_data: bool
    ) -> list[dict]:
        """Transform a batch of msms.txt rows into QPX PSM records."""
        records: list[dict] = []
        for row in df.to_dict("records"):
            try:
                rec = self._transform_row(row, spectral_data)
                if rec:
                    records.append(rec)
            except Exception as e:
                self.logger.debug(f"Skipping MaxQuant PSM row: {e}")
        return records

    def _transform_row(self, row, spectral_data: bool) -> Optional[dict]:
        """Transform a single msms.txt row into a QPX PSM record."""

        sequence = str(row.get("Sequence", ""))
        peptidoform = clean_peptidoform(str(row.get("Modified sequence", "")))
        charge = int(row.get("Charge", 0))

        # is_decoy (bool) -- MaxQuant uses '+' for Reverse
        is_decoy = mq_flag_to_bool(row.get("Reverse", ""))

        # Scan (list<int32>)
        scan_raw = row.get("Scan number")
        scan = []
        if pd.notna(scan_raw):
            try:
                scan = [int(scan_raw)]
            except (ValueError, TypeError):
                pass

        # Run file name (strip extension)
        run_file_name = str(row.get("Raw file", ""))

        # m/z values
        observed_mz = safe_float(row.get("m/z")) or 0.0

        # Calculated m/z -- needs to be computed from peptidoform + charge
        # We set 0.0 here; downstream can recompute with PyOpenMS if needed
        calculated_mz = 0.0

        # RT
        rt = safe_float(row.get("Retention time"))

        # PEP
        pep = safe_float(row.get("PEP"))

        # Protein accessions
        proteins_raw = str(row.get("Proteins", ""))
        protein_accessions = (
            proteins_raw.split(";") if proteins_raw else []
        )

        # Ion mobility
        ion_mobility = safe_float(row.get("1/K0"))

        # Additional scores
        additional_scores = []
        andromeda_score = safe_float(row.get("Score"))
        if andromeda_score is not None:
            additional_scores.append(
                {"score_name": "andromeda_score", "score_value": andromeda_score, "higher_better": True}
            )
        delta_score = safe_float(row.get("Delta score"))
        if delta_score is not None:
            additional_scores.append(
                {"score_name": "andromeda_delta_score", "score_value": delta_score, "higher_better": True}
            )
        pif = safe_float(row.get("PIF"))
        if pif is not None:
            additional_scores.append(
                {"score_name": "parent_ion_fraction", "score_value": pif, "higher_better": True}
            )

        # CV params (from Fragmentation, Mass analyzer, Type)
        cv_params = []
        for cv_col in ["Fragmentation", "Mass analyzer", "Type"]:
            val = row.get(cv_col)
            if pd.notna(val) and val:
                cv_params.append({"cv_name": cv_col.lower(), "cv_value": str(val)})

        # Spectral arrays
        mz_array = None
        intensity_array = None
        charge_array = None
        ion_type_array = None
        num_peaks = None

        if spectral_data:
            num_peaks_raw = row.get("Number of matches")
            if pd.notna(num_peaks_raw):
                try:
                    num_peaks = int(num_peaks_raw)
                except (ValueError, TypeError):
                    pass

            masses_raw = row.get("Masses")
            if pd.notna(masses_raw) and masses_raw:
                try:
                    mz_array = [float(x) for x in str(masses_raw).split(";") if x.strip()]
                except ValueError:
                    pass

            intensities_raw = row.get("Intensities")
            if pd.notna(intensities_raw) and intensities_raw:
                try:
                    intensity_array = [float(x) for x in str(intensities_raw).split(";") if x.strip()]
                except ValueError:
                    pass

        return {
            "sequence": sequence,
            "peptidoform": peptidoform,
            "modifications": None,
            "charge": charge,
            "posterior_error_probability": pep,
            "is_decoy": is_decoy,
            "calculated_mz": calculated_mz,
            "observed_mz": observed_mz,
            "additional_scores": additional_scores or None,
            "predicted_rt": None,
            "run_file_name": run_file_name,
            "cv_params": cv_params or None,
            "scan": scan,
            "rt": rt,
            "ion_mobility": ion_mobility,
            "protein_accessions": protein_accessions or None,
            "mz_array": mz_array,
            "intensity_array": intensity_array,
            "charge_array": charge_array,
            "ion_type_array": ion_type_array,
            "ion_mobility_array": None,
        }
