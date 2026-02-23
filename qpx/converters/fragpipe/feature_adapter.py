"""FragPipe Feature adapter -- combined_ion.tsv or combined_peptide.tsv to feature.parquet.

Loads a FragPipe ``combined_ion.tsv`` or ``combined_peptide.tsv`` into DuckDB,
transforms feature rows into ``FeatureSchema``, and writes through ``FeatureWriter``.

Auto-detects the input format:
    - ``combined_ion.tsv``: has a ``Charge`` (singular) column — one row per precursor
    - ``combined_peptide.tsv``: has a ``Charges`` (plural, comma-separated) column — one
      row per peptide sequence
"""

from __future__ import annotations

import logging
from typing import Optional

import pandas as pd

from qpx.converters.base import BaseConverter, resolve_columns
from qpx.converters.fragpipe.constants import FIELD_MAPPINGS
from qpx.converters.fragpipe.constants import to_modifications, to_proforma
from qpx.converters.utils import safe_float
from qpx.writers.feature import FeatureWriter

logger = logging.getLogger(__name__)

# Derive field map from constants
_FEATURE_MAP = FIELD_MAPPINGS["feature"]


def _extract_anchor_protein(protein_str: str) -> str:
    """Extract the first UniProt accession from a FragPipe Protein field.

    Handles formats like ``sp|P12345|PROT_HUMAN`` or plain accession ``P12345``.
    """
    pg_proteins = _extract_pg_proteins(protein_str, start=None, end=None)
    return pg_proteins[0]["accession"] if pg_proteins else ""


def _extract_pg_proteins(
    protein_str: str,
    start: int | None = None,
    end: int | None = None,
) -> list[dict]:
    """Extract pg_protein structs from a FragPipe Protein field.

    Handles formats like ``sp|P12345|PROT_HUMAN`` or plain accession ``P12345``.
    Returns list of {accession, start, end}. When the tool does not provide
    start/end, pass None and they will be encapsulated as null in the struct.
    """
    if not protein_str:
        return []
    result = []
    for part in protein_str.split(","):
        part = part.strip()
        if not part:
            continue
        acc = part.split("|")[1] if "|" in part and len(part.split("|")) >= 2 else part
        if acc:
            result.append({"accession": acc, "start": start, "end": end})
    return result


class FragPipeFeatureAdapter(BaseConverter):
    """Convert FragPipe ``combined_ion.tsv`` or ``combined_peptide.tsv`` to ``feature.parquet``.

    Usage::

        with FragPipeFeatureAdapter() as adapter:
            adapter.convert(
                feature_path="combined_ion.tsv",
                output_path="feature.parquet",
            )
    """

    def convert(
        self,
        feature_path: str,
        output_path: str,
        sdrf_path: Optional[str] = None,
        chunksize: int = 500_000,
        creator: str = "fragpipe",
    ) -> None:
        """Run the combined_ion/peptide.tsv -> feature.parquet conversion.

        Args:
            feature_path: Path to FragPipe ``combined_ion.tsv`` or ``combined_peptide.tsv``.
            output_path: Destination Parquet path.
            sdrf_path: Optional SDRF file (not used in current impl).
            chunksize: Rows per batch.
            creator: Creator tag in Parquet metadata.
        """
        # Step 1: Load feature file into DuckDB
        self._load_feature_file(feature_path)

        # Step 2: Resolve column mappings against actual input columns
        actual_cols = {
            c[0]
            for c in self._conn.execute(
                "SELECT column_name FROM information_schema.columns "
                "WHERE table_name='fragpipe_features'"
            ).fetchall()
        }
        self._resolved = resolve_columns(_FEATURE_MAP, actual_cols)

        # Step 3: Detect format and experiment columns
        format_type = self._detect_format()
        experiments = self._detect_experiment_columns()
        self.logger.info(f"Detected format: {format_type}, experiments: {experiments}")

        # Step 3: Stream and transform
        with FeatureWriter(output_path, creator=creator, compression=self._compression) as writer:
            for batch in self._query_batched("SELECT * FROM fragpipe_features", chunksize):
                df = batch.to_pandas()
                records = self._transform_batch(df, experiments, format_type)
                if records:
                    self._track_scores(records)
                    writer.write_batch(records)

        self.logger.info(f"FragPipe feature conversion complete -> {output_path}")

    # ------------------------------------------------------------------
    # Data loading
    # ------------------------------------------------------------------

    def _load_feature_file(self, path: str) -> None:
        """Load combined_ion.tsv or combined_peptide.tsv into DuckDB."""
        self._conn.execute(f"""
            CREATE TABLE fragpipe_features AS
            SELECT * FROM read_csv_auto('{path}',
                delim='\\t', header=true, auto_detect=true)
            """)
        count = self._conn.execute("SELECT COUNT(*) FROM fragpipe_features").fetchone()[
            0
        ]
        self.logger.info(f"Loaded {count:,} FragPipe feature rows")

    def _detect_format(self) -> str:
        """Return 'ion' if Charge column exists, 'peptide' if Charges."""
        cols = self._conn.execute(
            "SELECT column_name FROM information_schema.columns "
            "WHERE table_name='fragpipe_features'"
        ).fetchall()
        col_names = [c[0] for c in cols]
        return "ion" if "Charge" in col_names else "peptide"

    def _detect_experiment_columns(self) -> list[str]:
        """Detect per-experiment intensity columns.

        FragPipe feature files have columns like ``<experiment> Intensity``.
        """
        cols = self._conn.execute(
            "SELECT column_name FROM information_schema.columns "
            "WHERE table_name='fragpipe_features'"
        ).fetchall()
        col_names = [c[0] for c in cols]

        experiments = set()
        for col in col_names:
            if col.endswith(" Intensity"):
                experiments.add(col[: -len(" Intensity")])

        return sorted(experiments)

    # ------------------------------------------------------------------
    # Transform
    # ------------------------------------------------------------------

    def _transform_batch(
        self,
        df: pd.DataFrame,
        experiments: list[str],
        format_type: str,
    ) -> list[dict]:
        records: list[dict] = []
        # Pre-extract column arrays for faster per-row access than to_dict("records")
        col_arrays = {col: df[col].values for col in df.columns}
        n_rows = len(df)
        for i in range(n_rows):
            try:
                row = {col: vals[i] for col, vals in col_arrays.items()}
                recs = self._transform_row(row, experiments, format_type)
                records.extend(recs)
            except Exception as e:
                self.logger.debug(f"Skipping FragPipe feature row: {e}")
        return records

    def _transform_row(
        self,
        row,
        experiments: list[str],
        format_type: str,
    ) -> list[dict]:
        """Transform a single row into one or more feature records.

        Expands into one record per experiment (run) with non-zero intensity.
        For peptide format, also expands per charge state.
        """
        records: list[dict] = []
        r = self._resolved  # shorthand for resolved column mappings

        # Sequence
        sequence = str(row.get(r.get("sequence", "Peptide Sequence"), ""))

        # Peptidoform -- build ProForma from sequence + Assigned Modifications
        mods_col = r.get("modifications", "Assigned Modifications")
        assigned_mods_raw = row.get(mods_col, "")
        assigned_mods_str = (
            str(assigned_mods_raw)
            if pd.notna(assigned_mods_raw) and assigned_mods_raw
            else ""
        )
        peptidoform = to_proforma(assigned_mods_str, sequence)

        # Protein mapping (schema: list<pg_protein>; FragPipe does not provide start/end)
        protein_raw = str(row.get(r.get("pg_accessions", "Protein"), ""))
        pg_accessions = _extract_pg_proteins(protein_raw, start=None, end=None) or None
        anchor_protein = pg_accessions[0]["accession"] if pg_accessions else ""

        # Gene
        gene_raw = row.get(r.get("gg_names", "Gene"), "")
        gg_names = [str(gene_raw)] if pd.notna(gene_raw) and gene_raw else None

        # Modifications (reuse assigned_mods_str already extracted for peptidoform)
        modifications = None
        if assigned_mods_str:
            modifications = to_modifications(assigned_mods_str, sequence)

        # M/Z
        mz = safe_float(row.get(r.get("observed_mz", "M/Z"))) or 0.0

        # Determine charge states
        if format_type == "ion":
            charges = [int(row.get(r.get("charge", "Charge"), 0))]
        else:
            charges_col = r.get("charges", "Charges")
            charges_raw = str(row.get(charges_col, "0"))
            charges = [int(c.strip()) for c in charges_raw.split(",") if c.strip()]
            if not charges:
                charges = [0]

        # Expand per experiment and per charge
        for experiment in experiments:
            int_col = f"{experiment} Intensity"
            intensity_val = safe_float(row.get(int_col)) or 0.0
            if intensity_val <= 0:
                continue

            intensities = [{"label": "LFQ", "intensity": float(intensity_val)}]

            for charge in charges:
                rec = {
                    "sequence": sequence,
                    "peptidoform": peptidoform,
                    "modifications": modifications,
                    "charge": charge,
                    "posterior_error_probability": None,
                    "is_decoy": False,
                    "calculated_mz": float(mz),
                    "observed_mz": float(mz),
                    "additional_scores": None,
                    "predicted_rt": None,
                    "run_file_name": experiment,
                    "cv_params": None,
                    "scan": [],
                    "rt": None,
                    "ion_mobility": None,
                    "intensities": intensities,
                    "additional_intensities": None,
                    "pg_accessions": pg_accessions,
                    "anchor_protein": anchor_protein,
                    "unique": len(pg_accessions) <= 1 if pg_accessions else True,
                    "pg_global_qvalue": None,
                    "ion_mobility_start": None,
                    "ion_mobility_stop": None,
                    "gg_accessions": None,
                    "gg_names": gg_names,
                    "id_run_file_name": experiment,
                    "rt_start": None,
                    "rt_stop": None,
                }
                records.append(rec)

        return records
