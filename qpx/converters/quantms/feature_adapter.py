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
from qpx.converters.utils import (
    safe_float,
    parse_scan_numbers,
    resolve_run_file,
)
from qpx.core.cleavage import count_missed_cleavages
from qpx.converters.mztab import (
    load_mztab_sections,
    load_msstats,
    extract_ms_runs,
    extract_modifications,
    extract_score_names,
)
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
        enzyme_name: str | None = None,
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
            enzyme_name: Enzyme name for computing missed cleavages.
        """
        self._enzyme_name = enzyme_name
        # Step 1: Load data into DuckDB (skip if already loaded)
        if not self._table_exists("psms"):
            load_mztab_sections(self._conn, mztab_path)
        if not self._table_exists("msstats"):
            load_msstats(self._conn, msstats_path)

        # Step 2: Extract auxiliary lookups
        ms_runs = extract_ms_runs(self._conn)
        self._modifications_meta = extract_modifications(self._conn)
        _score_names = extract_score_names(self._conn)

        # Step 3: Build PSM lookup for merging best-PSM info with features
        psm_lookup = self._build_psm_lookup(ms_runs)

        # Step 4: Build protein q-value map
        protein_qvalue_map = self._build_protein_qvalue_map()

        # Step 4b: Build protein -> gene map
        protein_gene_map = self._build_protein_gene_map()

        # Step 5: Determine experiment type (LFQ or TMT)
        experiment_type = self._detect_experiment_type()

        # Step 6: Stream aggregated features and write
        self.logger.info("Aggregating MSstats data and writing features ...")

        with FeatureWriter(output_path, creator=creator, compression=self._compression) as writer:
            for batch_df in self._iter_feature_batches(file_batch_size):
                records = self._transform_batch(
                    batch_df,
                    psm_lookup,
                    protein_qvalue_map,
                    protein_gene_map,
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
        1: "TMT126",
        2: "TMT127N",
        3: "TMT127C",
        4: "TMT128N",
        5: "TMT128C",
        6: "TMT129N",
        7: "TMT129C",
        8: "TMT130N",
        9: "TMT130C",
        10: "TMT131N",
        11: "TMT131C",
        # TMT16plex extensions
        12: "TMT132N",
        13: "TMT132C",
        14: "TMT133N",
        15: "TMT133C",
        16: "TMT134N",
        # TMT18plex extensions
        17: "TMT134C",
        18: "TMT135N",
    }

    def _detect_experiment_type(self) -> str:
        """Detect experiment type from MSstats Channel column."""
        try:
            channels = self._conn.execute(
                'SELECT DISTINCT "Channel" FROM msstats LIMIT 20'
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

    def _build_psm_lookup(self, ms_runs: dict[int, str]) -> dict[tuple, dict]:
        """Build a lookup from (run_file_name, peptidoform, charge) -> PSM info.

        Uses a DuckDB SQL pre-aggregation to extract only the columns and
        first-occurrence rows needed, avoiding a full Python-side iteration
        over all 500K+ PSM rows.
        """
        lookup: dict[tuple, dict] = {}

        # Detect which CV column names are present for peptidoform / decoy
        actual_cols = {
            c[0]
            for c in self._conn.execute(
                "SELECT column_name FROM information_schema.columns "
                "WHERE table_name='psms'"
            ).fetchall()
        }

        # Peptidoform column (CV_PEPTIDOFORM_SEQUENCE = MS:1000889)
        pf_lo, pf_hi = (
            "opt_global_cv_ms:1000889_peptidoform_sequence",
            "opt_global_cv_MS:1000889_peptidoform_sequence",
        )
        if pf_lo in actual_cols:
            pf_col = pf_lo
        elif pf_hi in actual_cols:
            pf_col = pf_hi
        else:
            pf_col = None

        # Decoy column (CV_DECOY_PEPTIDE = MS:1002217)
        dec_lo, dec_hi = (
            "opt_global_cv_ms:1002217_decoy_peptide",
            "opt_global_cv_MS:1002217_decoy_peptide",
        )
        if dec_lo in actual_cols:
            dec_col = dec_lo
        elif dec_hi in actual_cols:
            dec_col = dec_hi
        else:
            dec_col = None

        # PEP column
        pep_col = None
        for candidate in [
            "opt_global_posterior_error_probability_score",
            "opt_global_Posterior_Error_Probability_score",
        ]:
            if candidate in actual_cols:
                pep_col = candidate
                break

        try:
            # Build SQL to extract exactly the needed columns
            pf_expr = f'CAST("{pf_col}" AS VARCHAR)' if pf_col else "''"
            dec_expr = f'CAST("{dec_col}" AS VARCHAR)' if dec_col else "'0'"
            pep_expr = f'TRY_CAST("{pep_col}" AS DOUBLE)' if pep_col else "NULL"

            sql = f"""
                SELECT
                    spectra_ref,
                    {pf_expr} AS peptidoform,
                    CAST(charge AS VARCHAR) AS charge,
                    {pep_expr} AS pep,
                    TRY_CAST(calc_mass_to_charge AS DOUBLE) AS calc_mz,
                    TRY_CAST(exp_mass_to_charge AS DOUBLE) AS obs_mz,
                    {dec_expr} AS is_decoy_raw,
                    CAST(accession AS VARCHAR) AS accession
                FROM psms
            """
            rows = self._conn.execute(sql).fetchall()

            for (
                spectra_ref,
                peptidoform,
                charge,
                pep,
                calc_mz,
                obs_mz,
                is_decoy_raw,
                accession,
            ) in rows:
                spectra_ref = str(spectra_ref) if spectra_ref else ""
                run_file_raw = resolve_run_file(spectra_ref, ms_runs) or ""
                run_file = run_file_raw.split(".")[0] if run_file_raw else ""
                peptidoform = str(peptidoform) if peptidoform else ""
                charge = str(charge) if charge else "0"

                key = (run_file, peptidoform, charge)
                if key not in lookup:
                    is_decoy = (
                        str(is_decoy_raw).strip() == "1" if is_decoy_raw else False
                    )
                    scan = parse_scan_numbers(spectra_ref)
                    lookup[key] = {
                        "pep": float(pep) if pep is not None else None,
                        "calculated_mz": (
                            float(calc_mz) if calc_mz is not None else None
                        ),
                        "observed_mz": float(obs_mz) if obs_mz is not None else None,
                        "is_decoy": is_decoy,
                        "accession": str(accession) if accession else "",
                        "scan": scan,
                        "id_run_file_name": run_file,
                    }
        except Exception as e:
            self.logger.warning(f"Could not build PSM lookup: {e}")

        return lookup

    def _build_protein_qvalue_map(self) -> dict[str, float]:
        """Build protein accession -> global q-value."""
        try:
            # Detect the actual score column name (may contain brackets)
            cols = {
                c[0]
                for c in self._conn.execute(
                    "SELECT column_name FROM information_schema.columns "
                    "WHERE table_name='proteins'"
                ).fetchall()
            }
            score_col = None
            for candidate in [
                "best_search_engine_score[1]",
                "best_search_engine_score_1",
            ]:
                if candidate in cols:
                    score_col = candidate
                    break
            if not score_col:
                return {}
            rows = self._conn.execute(f"""
                SELECT accession, "{score_col}"
                FROM proteins
                WHERE accession IS NOT NULL
                  AND "{score_col}" IS NOT NULL
                """).fetchall()
            return {str(acc): float(qval) for acc, qval in rows}
        except Exception:
            return {}

    def _build_protein_gene_map(self) -> dict[str, list[str]]:
        """Build protein accession -> gene symbol(s) from mzTab proteins table.

        Gene names are extracted from the UniProt description field using
        the ``GN=`` tag, mirroring the logic in
        ``QuantmsPgAdapter._parse_protein_row``.
        """
        gene_map: dict[str, list[str]] = {}
        try:
            rows = self._conn.execute(
                "SELECT accession, description FROM proteins "
                "WHERE accession IS NOT NULL"
            ).fetchall()
            for acc, desc in rows:
                acc_str = str(acc).strip()
                if not acc_str or acc_str == "null":
                    continue
                if desc and str(desc) != "null":
                    gn = re.search(r"GN=([^\s]+)", str(desc))
                    if gn:
                        gene_map[acc_str] = [gn.group(1)]
        except Exception:
            pass
        return gene_map

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
        protein_gene_map: dict,
        experiment_type: str,
        ms_runs: dict,
    ) -> list[dict]:
        """Transform an MSstats batch DataFrame into QPX feature records.

        Dispatches to a fast LFQ path (direct row iteration, no groupby)
        or an isobaric path (groupby to aggregate channels per feature).
        """
        col_map = resolve_columns(_FEATURE_MAP, set(df.columns))
        charge_col = col_map.get("charge", "Charge")
        if charge_col not in df.columns:
            df[charge_col] = 0

        if experiment_type == "LFQ":
            return self._transform_batch_lfq(
                df, col_map, psm_lookup, protein_qvalue_map, protein_gene_map
            )
        return self._transform_batch_isobaric(
            df, col_map, psm_lookup, protein_qvalue_map, protein_gene_map, experiment_type
        )

    # ------ LFQ fast path (1 row per feature, no groupby) ------

    def _transform_batch_lfq(
        self,
        df: pd.DataFrame,
        col_map: dict,
        psm_lookup: dict,
        protein_qvalue_map: dict,
        protein_gene_map: dict,
    ) -> list[dict]:
        """Build feature records for LFQ data by iterating rows directly.

        For LFQ, each MSstats row is a unique (peptidoform, protein, run, charge)
        combination, so groupby is pure overhead.  We use numpy column arrays
        for fast element access and avoid ``to_dict("records")`` entirely.
        """
        records: list[dict] = []

        pf_col = col_map.get("peptidoform", "PeptideSequence")
        prot_col = col_map.get("pg_accessions", "ProteinName")
        ref_col = col_map.get("run_file_name", "Reference")
        charge_col = col_map.get("charge", "Charge")
        intensity_col = col_map.get("intensity", "Intensity")
        rt_col = col_map.get("rt", "RetentionTime")

        # Pre-extract numpy arrays for fast element access
        pf_arr = df[pf_col].values
        prot_arr = df[prot_col].values
        ref_arr = df[ref_col].values
        charge_arr = df[charge_col].values
        int_arr = df[intensity_col].values
        rt_arr = df[rt_col].values if rt_col in df.columns else None

        _sub = re.sub
        _from_proforma = from_proforma
        _safe_float = safe_float
        mods_meta = self._modifications_meta

        n = len(df)
        for i in range(n):
            try:
                peptidoform = (
                    str(pf_arr[i])
                    if pf_arr[i] is not None and pd.notna(pf_arr[i])
                    else ""
                )
                protein_name = (
                    str(prot_arr[i])
                    if prot_arr[i] is not None and pd.notna(prot_arr[i])
                    else ""
                )
                run_file_name = (
                    str(ref_arr[i]).split(".")[0]
                    if ref_arr[i] is not None and pd.notna(ref_arr[i])
                    else ""
                )

                charge_raw = charge_arr[i]
                charge = (
                    int(float(charge_raw))
                    if charge_raw is not None
                    and charge_raw != ""
                    and pd.notna(charge_raw)
                    else 0
                )

                sequence = (
                    _sub(r"[^A-Z]", "", peptidoform.upper()) if peptidoform else ""
                )

                modifications = (
                    _from_proforma(peptidoform, sequence, meta=mods_meta)
                    if peptidoform and peptidoform != sequence
                    else None
                )

                intensity_val = _safe_float(int_arr[i]) or 0.0
                intensities = [{"label": "LFQ", "intensity": float(intensity_val)}]

                rt = _safe_float(rt_arr[i]) if rt_arr is not None else None

                psm_key = (run_file_name, peptidoform, str(charge))
                psm_info = psm_lookup.get(psm_key, {})

                _calc = psm_info.get("calculated_mz")
                _obs = psm_info.get("observed_mz")
                mass_error_ppm = (
                    1e6 * (_obs - _calc) / _calc
                    if _calc and _obs
                    else None
                )

                acc_list = protein_name.split(";") if protein_name else []
                anchor_protein = acc_list[0] if acc_list else ""

                records.append(
                    {
                        "sequence": sequence,
                        "peptidoform": peptidoform,
                        "modifications": modifications,
                        "charge": charge,
                        "posterior_error_probability": psm_info.get("pep"),
                        "is_decoy": psm_info.get("is_decoy", False),
                        "calculated_mz": _calc or 0.0,
                        "observed_mz": _obs or 0.0,
                        "mass_error_ppm": mass_error_ppm,
                        "missed_cleavages": count_missed_cleavages(sequence, self._enzyme_name) if self._enzyme_name else None,
                        "additional_scores": None,
                        "predicted_rt": None,
                        "run_file_name": run_file_name,
                        "cv_params": None,
                        "scan": psm_info.get("scan", []),
                        "rt": rt,
                        "ion_mobility": None,
                        "intensities": intensities,
                        "additional_intensities": None,
                        "pg_accessions": (
                            [
                                {"accession": a, "start": None, "end": None, "pre": None, "post": None}
                                for a in acc_list
                            ]
                            if acc_list
                            else None
                        ),
                        "anchor_protein": anchor_protein,
                        "unique": len(acc_list) <= 1,
                        "pg_global_qvalue": protein_qvalue_map.get(anchor_protein),
                        "pg_positions": None,
                        "ion_mobility_start": None,
                        "ion_mobility_stop": None,
                        "gg_accessions": protein_gene_map.get(anchor_protein),
                        "gg_names": protein_gene_map.get(anchor_protein),
                        "id_run_file_name": psm_info.get("id_run_file_name"),
                        "rt_start": None,
                        "rt_stop": None,
                    }
                )
            except Exception as e:
                self.logger.debug(f"Skipping feature row {i}: {e}")

        return records

    # ------ Isobaric path (TMT / iTRAQ -- groupby to aggregate channels) ------

    def _transform_batch_isobaric(
        self,
        df: pd.DataFrame,
        col_map: dict,
        psm_lookup: dict,
        protein_qvalue_map: dict,
        protein_gene_map: dict,
        experiment_type: str,
    ) -> list[dict]:
        """Build feature records for isobaric (TMT/iTRAQ) data.

        Uses groupby to aggregate channel intensities per feature, but
        accesses column values via ``.values`` arrays instead of
        ``to_dict("records")``.
        """
        records: list[dict] = []

        pf_col = col_map.get("peptidoform", "PeptideSequence")
        prot_col = col_map.get("pg_accessions", "ProteinName")
        ref_col = col_map.get("run_file_name", "Reference")
        charge_col = col_map.get("charge", "Charge")
        intensity_col = col_map.get("intensity", "Intensity")
        channel_col = col_map.get("channel", "Channel")
        rt_col = col_map.get("rt", "RetentionTime")

        grouping = [pf_col, prot_col, ref_col, charge_col]

        _sub = re.sub
        _from_proforma = from_proforma
        _safe_float = safe_float
        _tmt_map = self._TMT_CHANNEL_MAP
        mods_meta = self._modifications_meta
        is_tmt = experiment_type == "TMT"

        for group_key, group_data in df.groupby(grouping, dropna=False):
            try:
                peptidoform, protein_name, run_file_name, charge = group_key

                peptidoform = str(peptidoform) if peptidoform else ""
                protein_name = str(protein_name) if protein_name else ""
                run_file_name = (
                    str(run_file_name).split(".")[0] if run_file_name else ""
                )
                charge = int(float(charge)) if charge not in (None, "", "null") else 0

                sequence = (
                    _sub(r"[^A-Z]", "", peptidoform.upper()) if peptidoform else ""
                )

                modifications = (
                    _from_proforma(peptidoform, sequence, meta=mods_meta)
                    if peptidoform and peptidoform != sequence
                    else None
                )

                # Build intensities from channel values (no to_dict)
                ch_vals = group_data[channel_col].values
                int_vals = group_data[intensity_col].values
                intensities = []
                for j in range(len(ch_vals)):
                    ch_raw = ch_vals[j]
                    if (
                        is_tmt
                        and ch_raw is not None
                        and pd.notna(ch_raw)
                        and ch_raw != ""
                    ):
                        try:
                            label = _tmt_map.get(int(float(ch_raw)), str(ch_raw))
                        except (ValueError, TypeError):
                            label = str(ch_raw)
                    else:
                        label = (
                            str(ch_raw)
                            if ch_raw is not None and pd.notna(ch_raw) and ch_raw != ""
                            else "LFQ"
                        )
                    iv = _safe_float(int_vals[j]) or 0.0
                    intensities.append({"label": label, "intensity": float(iv)})

                rt = (
                    _safe_float(group_data[rt_col].values[0])
                    if rt_col in group_data.columns
                    else None
                )

                psm_key = (run_file_name, peptidoform, str(charge))
                psm_info = psm_lookup.get(psm_key, {})

                _calc = psm_info.get("calculated_mz")
                _obs = psm_info.get("observed_mz")
                mass_error_ppm = (
                    1e6 * (_obs - _calc) / _calc
                    if _calc and _obs
                    else None
                )

                acc_list = protein_name.split(";") if protein_name else []
                anchor_protein = acc_list[0] if acc_list else ""

                records.append(
                    {
                        "sequence": sequence,
                        "peptidoform": peptidoform,
                        "modifications": modifications,
                        "charge": charge,
                        "posterior_error_probability": psm_info.get("pep"),
                        "is_decoy": psm_info.get("is_decoy", False),
                        "calculated_mz": _calc or 0.0,
                        "observed_mz": _obs or 0.0,
                        "mass_error_ppm": mass_error_ppm,
                        "missed_cleavages": count_missed_cleavages(sequence, self._enzyme_name) if self._enzyme_name else None,
                        "additional_scores": None,
                        "predicted_rt": None,
                        "run_file_name": run_file_name,
                        "cv_params": None,
                        "scan": psm_info.get("scan", []),
                        "rt": rt,
                        "ion_mobility": None,
                        "intensities": intensities or None,
                        "additional_intensities": None,
                        "pg_accessions": (
                            [
                                {"accession": a, "start": None, "end": None, "pre": None, "post": None}
                                for a in acc_list
                            ]
                            if acc_list
                            else None
                        ),
                        "anchor_protein": anchor_protein,
                        "unique": len(acc_list) <= 1,
                        "pg_global_qvalue": protein_qvalue_map.get(anchor_protein),
                        "pg_positions": None,
                        "ion_mobility_start": None,
                        "ion_mobility_stop": None,
                        "gg_accessions": protein_gene_map.get(anchor_protein),
                        "gg_names": protein_gene_map.get(anchor_protein),
                        "id_run_file_name": psm_info.get("id_run_file_name"),
                        "rt_start": None,
                        "rt_stop": None,
                    }
                )
            except Exception as e:
                self.logger.debug(f"Skipping feature group: {e}")

        return records

    @staticmethod
    def _detect_msstats_columns(df: pd.DataFrame) -> dict[str, str]:
        """Detect the actual MSstats column names in a DataFrame.

        .. deprecated:: Use ``resolve_columns(FIELD_MAPPINGS["feature"], ...)`` instead.
           Kept for backward compatibility with external callers.
        """
        return resolve_columns(_FEATURE_MAP, set(df.columns))
