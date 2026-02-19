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

from qpx.converters.base import BaseConverter, resolve_columns
from qpx.converters.maxquant.constants import FIELD_MAPPINGS
from qpx.converters.utils import clean_peptidoform, mq_flag_to_bool, safe_float
from qpx.writers.psm import PsmWriter

logger = logging.getLogger(__name__)

# Derive expected columns from FIELD_MAPPINGS
_PSM_MAP = FIELD_MAPPINGS["psm"]

# MaxQuant columns to read from msms.txt (derived from constants)
_MQ_PSM_USECOLS = list({
    candidates[0]
    for candidates in _PSM_MAP.values()
}) + [
    "Proteins",       # protein accessions (not in FIELD_MAPPINGS)
    "PIF",            # parent ion fraction score
    "Masses",         # spectral masses
    "Intensities",    # spectral intensities
    "Number of matches",
    "1/K0",           # ion mobility
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

        # Step 2: Resolve column mappings against actual input columns
        actual_cols = {
            c[0] for c in self._conn.execute(
                "SELECT column_name FROM information_schema.columns "
                "WHERE table_name='msms'"
            ).fetchall()
        }
        self._resolved = resolve_columns(_PSM_MAP, actual_cols)

        # Step 3: Stream and transform
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
        r = self._resolved  # shorthand for resolved column mappings

        sequence = str(row.get(r.get("sequence", "Sequence"), ""))
        peptidoform = clean_peptidoform(
            str(row.get(r.get("modified_sequence", "Modified sequence"), ""))
        )
        charge = int(row.get(r.get("charge", "Charge"), 0))

        # is_decoy (bool) -- MaxQuant uses '+' for Reverse
        is_decoy = mq_flag_to_bool(
            row.get(r.get("is_decoy", "Reverse"), "")
        )

        # Scan (list<int32>)
        scan_raw = row.get(r.get("scan", "Scan number"))
        scan = []
        if pd.notna(scan_raw):
            try:
                scan = [int(scan_raw)]
            except (ValueError, TypeError):
                pass

        # Run file name (strip extension)
        run_file_name = str(row.get(r.get("run_file_name", "Raw file"), ""))

        # m/z values
        observed_mz = safe_float(
            row.get(r.get("observed_mz", "m/z"))
        ) or 0.0

        # Calculated m/z -- needs to be computed from peptidoform + charge
        # We set 0.0 here; downstream can recompute with PyOpenMS if needed
        calculated_mz = 0.0

        # RT
        rt = safe_float(row.get(r.get("rt", "Retention time")))

        # PEP
        pep = safe_float(
            row.get(r.get("posterior_error_probability", "PEP"))
        )

        # Protein accessions
        proteins_raw = str(row.get("Proteins", ""))
        protein_accessions = (
            proteins_raw.split(";") if proteins_raw else []
        )

        # Ion mobility
        ion_mobility = safe_float(row.get("1/K0"))

        # Additional scores
        additional_scores = []
        andromeda_score = safe_float(
            row.get(r.get("andromeda_score", "Score"))
        )
        if andromeda_score is not None:
            additional_scores.append(
                {"score_name": "andromeda_score", "score_value": andromeda_score, "higher_better": True}
            )
        delta_score = safe_float(
            row.get(r.get("andromeda_delta_score", "Delta score"))
        )
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
