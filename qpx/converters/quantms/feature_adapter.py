"""QuantMS Feature adapter -- mzTab + MSstats to feature.parquet.

Loads an mzTab file and its companion MSstats input into DuckDB, joins
them via SQL, and streams the result through ``FeatureWriter``.

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

from qpx.converters.base import BaseConverter, resolve_columns
from qpx.converters.ptm import from_proforma
from qpx.converters.quantms.constants import FIELD_MAPPINGS
from qpx.converters.utils import safe_float, parse_scan_numbers, resolve_run_file, get_cv_value
from qpx.converters.mztab import (
    load_mztab_sections,
    load_msstats,
    extract_ms_runs,
    extract_modifications,
    extract_score_names,
)
from qpx.core.cv_terms import CV_PEPTIDOFORM_SEQUENCE, CV_DECOY_PEPTIDE
from qpx.writers.feature import FeatureWriter

logger = logging.getLogger(__name__)

# Derive field map from constants
_FEATURE_MAP = FIELD_MAPPINGS["feature"]


class QuantmsFeatureAdapter(BaseConverter):
    """Convert QuantMS mzTab + MSstats data to ``feature.parquet``.

    Usage::

        with QuantmsFeatureAdapter() as adapter:
            adapter.convert(
                mztab_path="results.mzTab",
                msstats_path="msstats_in.csv",
                output_path="feature.parquet",
            )
    """

    def convert(
        self,
        mztab_path: str,
        msstats_path: str,
        output_path: str,
        sdrf_path: Optional[str] = None,
        file_batch_size: int = 10,
        chunksize: int = 500_000,
        creator: str = "quantms",
    ) -> None:
        """Run the mzTab+MSstats -> feature.parquet conversion.

        Args:
            mztab_path: Path to the mzTab file.
            msstats_path: Path to the MSstats input file.
            output_path: Destination Parquet path.
            sdrf_path: Optional SDRF path for sample mapping.
            file_batch_size: Number of run files to batch.
            chunksize: Rows per streaming batch.
            creator: Creator tag in Parquet metadata.
        """
        # Step 1: Load data into DuckDB (skip if already loaded)
        if not self._table_exists("psms"):
            load_mztab_sections(self._conn, mztab_path)
        if not self._table_exists("msstats"):
            load_msstats(self._conn, msstats_path)

        # Step 2: Extract auxiliary lookups
        ms_runs = extract_ms_runs(self._conn)
        self._modifications_meta = extract_modifications(self._conn)
        score_names = extract_score_names(self._conn)

        # Step 3: Build PSM lookup for merging best-PSM info with features
        psm_lookup = self._build_psm_lookup(ms_runs)

        # Step 4: Build protein q-value map
        protein_qvalue_map = self._build_protein_qvalue_map()

        # Step 5: Determine experiment type (LFQ or TMT)
        experiment_type = self._detect_experiment_type()

        # Step 6: Stream aggregated features and write
        self.logger.info("Aggregating MSstats data and writing features ...")

        with FeatureWriter(output_path, creator=creator) as writer:
            for batch_df in self._iter_feature_batches(file_batch_size):
                records = self._transform_batch(
                    batch_df,
                    psm_lookup,
                    protein_qvalue_map,
                    experiment_type,
                    ms_runs,
                )
                if records:
                    self._track_scores(records)
                    writer.write_batch(records)

        self.logger.info(f"Feature conversion complete -> {output_path}")

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    # Standard TMT channel index → label mappings
    _TMT_CHANNEL_MAP: dict[int, str] = {
        1: "TMT126", 2: "TMT127N", 3: "TMT127C",
        4: "TMT128N", 5: "TMT128C", 6: "TMT129N",
        7: "TMT129C", 8: "TMT130N", 9: "TMT130C",
        10: "TMT131N", 11: "TMT131C",
        # TMT16plex extensions
        12: "TMT132N", 13: "TMT132C", 14: "TMT133N",
        15: "TMT133C", 16: "TMT134N",
        # TMT18plex extensions
        17: "TMT134C", 18: "TMT135N",
    }

    def _detect_experiment_type(self) -> str:
        """Detect experiment type from MSstats Channel column."""
        try:
            channels = self._conn.execute(
                "SELECT DISTINCT \"Channel\" FROM msstats LIMIT 20"
            ).fetchall()
            channel_vals = [str(c[0]).upper() for c in channels if c[0]]
            if any("TMT" in c for c in channel_vals):
                return "TMT"
            if any("ITRAQ" in c for c in channel_vals):
                return "iTRAQ"
            # Numeric channels with >1 distinct value indicate isobaric labeling
            if len(channel_vals) > 1:
                try:
                    int_vals = [int(float(v)) for v in channel_vals]
                    if all(1 <= v <= 18 for v in int_vals):
                        return "TMT"
                except (ValueError, TypeError):
                    pass
        except Exception:
            pass
        return "LFQ"

    def _build_psm_lookup(
        self, ms_runs: dict[int, str]
    ) -> dict[tuple, dict]:
        """Build a lookup from (run_file_name, peptidoform, charge) -> PSM info."""
        lookup: dict[tuple, dict] = {}
        try:
            psm_df = self._conn.execute("SELECT * FROM psms").df()

            for row in psm_df.to_dict("records"):
                spectra_ref = str(row.get("spectra_ref", ""))
                run_file = resolve_run_file(spectra_ref, ms_runs) or ""

                # mzTab column embeds CV_PEPTIDOFORM_SEQUENCE (MS:1000889)
                peptidoform = str(get_cv_value(row, CV_PEPTIDOFORM_SEQUENCE, "peptidoform_sequence", ""))
                charge = str(row.get("charge", "0"))

                key = (run_file, peptidoform, charge)
                pep = safe_float(
                    row.get(
                        "opt_global_posterior_error_probability_score",
                        row.get("opt_global_Posterior_Error_Probability_score"),
                    )
                )
                calc_mz = safe_float(row.get("calc_mass_to_charge"))
                obs_mz = safe_float(row.get("exp_mass_to_charge"))
                # mzTab column embeds CV_DECOY_PEPTIDE (MS:1002217)
                is_decoy_raw = str(get_cv_value(row, CV_DECOY_PEPTIDE, "decoy_peptide", "0")).strip()
                is_decoy = is_decoy_raw == "1"
                accession = str(row.get("accession", ""))

                scan = parse_scan_numbers(spectra_ref)

                if key not in lookup:
                    lookup[key] = {
                        "pep": pep,
                        "calculated_mz": calc_mz,
                        "observed_mz": obs_mz,
                        "is_decoy": is_decoy,
                        "accession": accession,
                        "scan": scan,
                        "id_run_file_name": run_file,
                    }
        except Exception as e:
            self.logger.warning(f"Could not build PSM lookup: {e}")

        return lookup

    def _build_protein_qvalue_map(self) -> dict[str, float]:
        """Build protein accession -> global q-value."""
        try:
            rows = self._conn.execute(
                """
                SELECT accession, best_search_engine_score_1
                FROM proteins
                WHERE accession IS NOT NULL
                  AND best_search_engine_score_1 IS NOT NULL
                """
            ).fetchall()
            return {str(acc): float(qval) for acc, qval in rows}
        except Exception:
            return {}

    def _iter_feature_batches(self, file_batch_size: int):
        """Yield DataFrames of MSstats data grouped by reference files."""
        try:
            ref_col = self._detect_ref_column()
            refs = self._conn.execute(
                f'SELECT DISTINCT "{ref_col}" FROM msstats ORDER BY "{ref_col}"'
            ).fetchall()
            ref_list = [r[0] for r in refs]

            for i in range(0, len(ref_list), file_batch_size):
                batch_refs = ref_list[i : i + file_batch_size]
                placeholders = ", ".join(["?" for _ in batch_refs])
                df = self._conn.execute(
                    f'SELECT * FROM msstats WHERE "{ref_col}" IN ({placeholders})',
                    batch_refs,
                ).df()
                if not df.empty:
                    yield df
        except Exception as e:
            self.logger.error(f"Error iterating feature batches: {e}")

    def _detect_ref_column(self) -> str:
        """Detect the reference/run column name in MSstats table.

        Uses FIELD_MAPPINGS candidates if available, falls back to common names.
        """
        cols = self._conn.execute(
            "SELECT column_name FROM information_schema.columns WHERE table_name='msstats'"
        ).fetchall()
        col_names = {c[0] for c in cols}
        # Try FIELD_MAPPINGS candidates first
        for candidate in _FEATURE_MAP.get("run_file_name", []):
            if candidate in col_names:
                return candidate
        # Fallback
        for candidate in ["Reference", "reference", "Run", "run"]:
            if candidate in col_names:
                return candidate
        return "Reference"

    def _transform_batch(
        self,
        df: pd.DataFrame,
        psm_lookup: dict,
        protein_qvalue_map: dict,
        experiment_type: str,
        ms_runs: dict,
    ) -> list[dict]:
        """Transform an MSstats batch DataFrame into QPX feature records."""
        records: list[dict] = []

        # Resolve column names against actual DataFrame columns
        col_map = resolve_columns(_FEATURE_MAP, set(df.columns))

        # Group by (peptidoform, pg_accessions, run_file, charge)
        peptidoform_col = col_map.get("peptidoform", "PeptideSequence")
        protein_col = col_map.get("pg_accessions", "ProteinName")
        ref_col = col_map.get("run_file_name", "Reference")
        charge_col = col_map.get("charge", "Charge")
        intensity_col = col_map.get("intensity", "Intensity")
        channel_col = col_map.get("channel", "Channel")

        # Ensure charge column exists
        if charge_col not in df.columns:
            df[charge_col] = 0

        grouping = [peptidoform_col, protein_col, ref_col]
        if charge_col in df.columns:
            grouping.append(charge_col)

        for group_key, group_data in df.groupby(grouping, dropna=False):
            try:
                rec = self._build_feature_record(
                    group_key,
                    group_data,
                    grouping,
                    col_map,
                    psm_lookup,
                    protein_qvalue_map,
                    experiment_type,
                )
                if rec:
                    records.append(rec)
            except Exception as e:
                self.logger.debug(f"Skipping feature group: {e}")

        return records

    def _build_feature_record(
        self,
        group_key,
        group_data: pd.DataFrame,
        grouping: list[str],
        col_map: dict,
        psm_lookup: dict,
        protein_qvalue_map: dict,
        experiment_type: str,
    ) -> Optional[dict]:
        """Build a single feature record from a grouped MSstats batch."""
        if len(grouping) == 4:
            peptidoform, protein_name, run_file_name, charge = group_key
        else:
            peptidoform, protein_name, run_file_name = group_key
            charge = 0

        peptidoform = str(peptidoform) if peptidoform else ""
        protein_name = str(protein_name) if protein_name else ""
        run_file_name = str(run_file_name).split(".")[0] if run_file_name else ""
        charge = int(float(charge)) if charge not in (None, "", "null") else 0

        # Extract plain sequence from peptidoform
        sequence = re.sub(r"[^A-Z]", "", peptidoform.upper()) if peptidoform else ""

        # Parse modifications from peptidoform (may contain ProForma-like tags)
        modifications = from_proforma(
            peptidoform, sequence, meta=self._modifications_meta
        ) if peptidoform and peptidoform != sequence else None

        # Build intensities (new schema: {label, intensity})
        intensity_col = col_map.get("intensity", "Intensity")
        channel_col = col_map.get("channel", "Channel")

        intensities = []
        for row in group_data.to_dict("records"):
            channel_raw = row.get(channel_col, "")
            if experiment_type == "TMT" and pd.notna(channel_raw) and channel_raw:
                # Map numeric channel indices to TMT label names
                try:
                    ch_int = int(float(channel_raw))
                    label = self._TMT_CHANNEL_MAP.get(ch_int, str(channel_raw))
                except (ValueError, TypeError):
                    label = str(channel_raw)
            else:
                label = str(channel_raw) if pd.notna(channel_raw) and channel_raw else "LFQ"
            intensity_val = safe_float(row.get(intensity_col, 0.0)) or 0.0

            intensities.append({"label": label, "intensity": float(intensity_val)})

        # PSM lookup for additional info
        psm_key = (run_file_name, peptidoform, str(charge))
        psm_info = psm_lookup.get(psm_key, {})

        # is_decoy (bool)
        is_decoy = psm_info.get("is_decoy", False)

        # Protein accessions (list<pg_protein>; MSstats does not provide start/end)
        acc_list = protein_name.split(";") if protein_name else []
        pg_accessions = (
            [{"accession": acc, "start": None, "end": None} for acc in acc_list]
            if acc_list
            else None
        )
        anchor_protein = acc_list[0] if acc_list else ""

        # Unique peptide indicator
        unique = len(acc_list) <= 1

        # Protein global q-value
        pg_global_qvalue = protein_qvalue_map.get(anchor_protein)

        # Scan from PSM lookup
        scan = psm_info.get("scan", [])

        # RT from first row
        rt_col = col_map.get("rt", "RetentionTime")
        first_row = group_data.iloc[0]
        rt = safe_float(first_row.get(rt_col))

        return {
            "sequence": sequence,
            "peptidoform": peptidoform,
            "modifications": modifications,
            "charge": charge,
            "posterior_error_probability": psm_info.get("pep"),
            "is_decoy": is_decoy,
            "calculated_mz": psm_info.get("calculated_mz") or 0.0,
            "observed_mz": psm_info.get("observed_mz") or 0.0,
            "additional_scores": None,
            "predicted_rt": None,
            "run_file_name": run_file_name,
            "cv_params": None,
            "scan": scan,
            "rt": rt,
            "ion_mobility": None,
            "intensities": intensities or None,
            "additional_intensities": None,
            "pg_accessions": pg_accessions,
            "anchor_protein": anchor_protein,
            "unique": unique,
            "pg_global_qvalue": pg_global_qvalue,
            "pg_positions": None,
            "ion_mobility_start": None,
            "ion_mobility_stop": None,
            "gg_accessions": None,
            "gg_names": None,
            "id_run_file_name": psm_info.get("id_run_file_name"),
            "rt_start": None,
            "rt_stop": None,
        }

    @staticmethod
    def _detect_msstats_columns(df: pd.DataFrame) -> dict[str, str]:
        """Detect the actual MSstats column names in a DataFrame.

        .. deprecated:: Use ``resolve_columns(FIELD_MAPPINGS["feature"], ...)`` instead.
           Kept for backward compatibility with external callers.
        """
        return resolve_columns(_FEATURE_MAP, set(df.columns))
