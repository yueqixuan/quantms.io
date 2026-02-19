"""FragPipe PSM adapter -- psm.tsv to psm.parquet.

Loads a FragPipe ``psm.tsv`` into DuckDB, transforms via SQL into
``PsmSchema``, and streams results through ``PsmWriter``.

Key schema changes:
    - ``precursor_charge`` (int32) -> ``charge`` (int16)
    - ``is_decoy`` (int32, 0/1) -> ``is_decoy`` (bool)
    - ``scan`` (string) -> ``scan`` (list<int32>)
    - ``source_file`` -> ``run_file_name``
"""

from __future__ import annotations

import logging
from typing import Optional, Tuple

import pandas as pd

from qpx.converters.base import BaseConverter, resolve_columns
from qpx.converters.fragpipe.constants import FIELD_MAPPINGS
from qpx.converters.utils import safe_float
from qpx.core.scores import normalize_score_name
from qpx.writers.psm import PsmWriter

logger = logging.getLogger(__name__)

# Derive field map from constants
_PSM_MAP = FIELD_MAPPINGS["psm"]


def _parse_spectrum_id(identifier: str) -> Tuple[str, int]:
    """Parse FragPipe Spectrum identifier into (source_file, scan_number).

    Format: ``<source_file>.<scan>.<scan>.<charge>``
    """
    tokens = str(identifier).split(".")
    source_file = tokens[0] if tokens else ""
    scan_number = 0
    if len(tokens) >= 2:
        try:
            scan_number = int(tokens[1].lstrip("0") or "0")
        except ValueError:
            pass
    return source_file, scan_number


class FragPipePsmAdapter(BaseConverter):
    """Convert FragPipe ``psm.tsv`` to ``psm.parquet``.

    Usage::

        with FragPipePsmAdapter() as adapter:
            adapter.convert(
                psm_path="psm.tsv",
                output_path="psm.parquet",
            )
    """

    def convert(
        self,
        psm_path: str,
        output_path: str,
        chunksize: int = 500_000,
        creator: str = "fragpipe",
    ) -> None:
        """Run the psm.tsv -> psm.parquet conversion.

        Args:
            psm_path: Path to FragPipe ``psm.tsv``.
            output_path: Destination Parquet path.
            chunksize: Rows per batch.
            creator: Creator tag in Parquet metadata.
        """
        # Step 1: Load psm.tsv into DuckDB
        self._load_psm_tsv(psm_path)

        # Step 2: Resolve column mappings against actual input columns
        actual_cols = {
            c[0] for c in self._conn.execute(
                "SELECT column_name FROM information_schema.columns "
                "WHERE table_name='fragpipe_psms'"
            ).fetchall()
        }
        self._resolved = resolve_columns(_PSM_MAP, actual_cols)

        # Step 3: Stream and transform
        self.logger.info("Transforming FragPipe PSMs ...")

        total = self._conn.execute("SELECT COUNT(*) FROM fragpipe_psms").fetchone()[0]

        with PsmWriter(output_path, creator=creator) as writer:
            offset = 0
            while offset < total:
                df = self._conn.execute(
                    f"SELECT * FROM fragpipe_psms LIMIT {chunksize} OFFSET {offset}"
                ).df()
                if df.empty:
                    break
                records = self._transform_batch(df)
                if records:
                    self._track_scores(records)
                    writer.write_batch(records)
                offset += chunksize

        self.logger.info(f"FragPipe PSM conversion complete -> {output_path}")

    # ------------------------------------------------------------------
    # Data loading
    # ------------------------------------------------------------------

    def _load_psm_tsv(self, path: str) -> None:
        """Load FragPipe psm.tsv into DuckDB."""
        self._conn.execute(
            f"""
            CREATE TABLE fragpipe_psms AS
            SELECT * FROM read_csv_auto('{path}',
                delim='\\t', header=true, auto_detect=true)
            """
        )
        count = self._conn.execute("SELECT COUNT(*) FROM fragpipe_psms").fetchone()[0]
        self.logger.info(f"Loaded {count:,} FragPipe PSM rows")

    # ------------------------------------------------------------------
    # Transform
    # ------------------------------------------------------------------

    def _transform_batch(self, df: pd.DataFrame) -> list[dict]:
        records: list[dict] = []
        for row in df.to_dict("records"):
            try:
                rec = self._transform_row(row)
                if rec:
                    records.append(rec)
            except Exception as e:
                self.logger.debug(f"Skipping FragPipe PSM row: {e}")
        return records

    def _transform_row(self, row) -> Optional[dict]:
        """Transform a single FragPipe psm.tsv row."""
        r = self._resolved  # shorthand for resolved column mappings

        sequence = str(row.get(r.get("sequence", "Peptide"), ""))
        charge = int(row.get(r.get("charge", "Charge"), 0))

        # Parse Spectrum column for source file and scan
        spectrum_raw = str(row.get("Spectrum", ""))
        source_file, scan_number = _parse_spectrum_id(spectrum_raw)

        # run_file_name is the source file stem
        run_file_name = source_file

        # scan (list<int32>)
        scan = [scan_number] if scan_number > 0 else []

        # m/z
        observed_mz = safe_float(row.get("Calibrated Observed M/Z",
                                            row.get(r.get("observed_mz", "Observed M/Z")))) or 0.0
        calculated_mz = safe_float(
            row.get(r.get("calculated_mz", "Calculated M/Z"))
        ) or 0.0

        # RT
        rt = safe_float(row.get(r.get("rt", "Retention")))

        # PEP (FragPipe uses PeptideProphet Probability; convert to 1-prob)
        pp_prob = safe_float(row.get("PeptideProphet Probability"))
        pep = None
        if pp_prob is not None:
            pep = 1.0 - pp_prob if pp_prob <= 1.0 else pp_prob

        # Is decoy (bool) -- FragPipe typically marks with "rev_" prefix or
        # "Is Unique" column; we default to False for now
        is_decoy = False

        # Protein accessions
        protein_id = str(row.get(r.get("pg_accessions", "Protein"), ""))
        protein_accessions = [protein_id] if protein_id else []

        # Peptidoform -- FragPipe's "Modified Peptide" or the peptide itself
        peptidoform = str(
            row.get(r.get("modified_sequence", "Modified Peptide"), sequence)
        )

        # Ion mobility
        ion_mobility = safe_float(row.get("Ion Mobility"))

        # Additional scores
        additional_scores = []
        hyperscore = safe_float(row.get("Hyperscore"))
        if hyperscore is not None:
            additional_scores.append(
                {"score_name": normalize_score_name("MSFragger:Hyperscore"), "score_value": hyperscore, "higher_better": True}
            )
        expectation = safe_float(row.get("Expectation"))
        if expectation is not None:
            additional_scores.append(
                {"score_name": normalize_score_name("MSFragger:Expectation"), "score_value": expectation, "higher_better": False}
            )

        # CV params (missed cleavages, enzymatic termini)
        cv_params = []
        missed = row.get("Number of Missed Cleavages")
        if pd.notna(missed):
            cv_params.append(
                {"cv_name": "number_of_missed_cleavages", "cv_value": str(int(missed))}
            )
        termini = row.get("Number of Enzymatic Termini")
        if pd.notna(termini):
            cv_params.append(
                {"cv_name": "number_of_enzymatic_termini", "cv_value": str(int(termini))}
            )

        # Modifications -- parse Assigned Modifications if present
        modifications = None
        assigned_mods = row.get("Assigned Modifications")
        if pd.notna(assigned_mods) and assigned_mods:
            modifications = self._parse_assigned_modifications(str(assigned_mods))

        return {
            "sequence": sequence,
            "peptidoform": peptidoform,
            "modifications": modifications,
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
            "mz_array": None,
            "intensity_array": None,
            "charge_array": None,
            "ion_type_array": None,
            "ion_mobility_array": None,
        }

    @staticmethod
    def _parse_assigned_modifications(text: str) -> Optional[list[dict]]:
        """Parse FragPipe ``Assigned Modifications`` field.

        Format: ``<pos><AA>(<mass>), ...``
        Example: ``7M(15.9949), 2C(57.0215)``
        """
        if not text:
            return None

        mods: dict[str, dict] = {}

        for token in text.split(","):
            token = token.strip()
            if not token:
                continue
            try:
                # Split "7M(15.9949)" into position+aa and mass
                paren_idx = token.index("(")
                pos_aa = token[:paren_idx].strip()
                mass = float(token[paren_idx + 1 : -1])

                # Determine position
                if pos_aa.lower() == "n-term":
                    position = 0
                    aa = None
                else:
                    aa = pos_aa[-1] if pos_aa else None
                    position = int(pos_aa[:-1]) if len(pos_aa) > 1 else 0

                mod_name = f"CHEMMOD:{mass:+.4f}"
                if mod_name not in mods:
                    mods[mod_name] = {
                        "name": mod_name,
                        "accession": None,
                        "positions": [],
                    }
                mods[mod_name]["positions"].append(
                    {"position": position, "amino_acid": aa, "scores": None}
                )
            except (ValueError, IndexError):
                continue

        return list(mods.values()) if mods else None
