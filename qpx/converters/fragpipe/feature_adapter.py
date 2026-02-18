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

from qpx.converters.base import BaseConverter
from qpx.converters.fragpipe.psm_adapter import FragPipePsmAdapter
from qpx.converters.utils import safe_float
from qpx.writers.feature import FeatureWriter

logger = logging.getLogger(__name__)


def _extract_anchor_protein(protein_str: str) -> str:
    """Extract the first UniProt accession from a FragPipe Protein field.

    Handles formats like ``sp|P12345|PROT_HUMAN`` or plain accession ``P12345``.
    """
    if not protein_str:
        return ""
    first = protein_str.split(",")[0].strip()
    parts = first.split("|")
    if len(parts) >= 2:
        return parts[1]
    return first


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

        # Step 2: Detect format and experiment columns
        format_type = self._detect_format()
        experiments = self._detect_experiment_columns()
        self.logger.info(
            f"Detected format: {format_type}, experiments: {experiments}"
        )

        # Step 3: Stream and transform
        total = self._conn.execute(
            "SELECT COUNT(*) FROM fragpipe_features"
        ).fetchone()[0]

        with FeatureWriter(output_path, creator=creator) as writer:
            offset = 0
            while offset < total:
                df = self._conn.execute(
                    f"SELECT * FROM fragpipe_features LIMIT {chunksize} OFFSET {offset}"
                ).df()
                if df.empty:
                    break
                records = self._transform_batch(df, experiments, format_type)
                if records:
                    self._track_scores(records)
                    writer.write_batch(records)
                offset += chunksize

        self.logger.info(f"FragPipe feature conversion complete -> {output_path}")

    # ------------------------------------------------------------------
    # Data loading
    # ------------------------------------------------------------------

    def _load_feature_file(self, path: str) -> None:
        """Load combined_ion.tsv or combined_peptide.tsv into DuckDB."""
        self._conn.execute(
            f"""
            CREATE TABLE fragpipe_features AS
            SELECT * FROM read_csv_auto('{path}',
                delim='\\t', header=true, auto_detect=true)
            """
        )
        count = self._conn.execute(
            "SELECT COUNT(*) FROM fragpipe_features"
        ).fetchone()[0]
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
        for row in df.to_dict("records"):
            try:
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

        # Sequence
        sequence = str(row.get("Peptide Sequence", ""))

        # Peptidoform (modified sequence)
        peptidoform = str(
            row.get("Modified Sequence", row.get("Modified Peptide", sequence))
        )

        # Protein mapping
        protein_raw = str(row.get("Protein", row.get("Protein ID", "")))
        anchor_protein = _extract_anchor_protein(protein_raw)
        protein_id = str(row.get("Protein ID", anchor_protein))
        pg_accessions = [protein_id] if protein_id else None

        # Gene
        gene_raw = row.get("Gene", "")
        gg_names = (
            [str(gene_raw)] if pd.notna(gene_raw) and gene_raw else None
        )

        # Modifications
        mods_raw = row.get("Assigned Modifications", "")
        modifications = None
        if pd.notna(mods_raw) and mods_raw:
            modifications = FragPipePsmAdapter._parse_assigned_modifications(
                str(mods_raw)
            )

        # M/Z
        mz = safe_float(row.get("M/Z")) or 0.0

        # Determine charge states
        if format_type == "ion":
            charges = [int(row.get("Charge", 0))]
        else:
            charges_raw = str(row.get("Charges", "0"))
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
                    "scan": [0],
                    "rt": None,
                    "ion_mobility": None,
                    "intensities": intensities,
                    "additional_intensities": None,
                    "pg_accessions": pg_accessions,
                    "anchor_protein": anchor_protein,
                    "unique": 1,
                    "pg_global_qvalue": None,
                    "pg_positions": None,
                    "ion_mobility_start": None,
                    "ion_mobility_stop": None,
                    "gg_accessions": None,
                    "gg_names": gg_names,
                    "id_run_file_name": None,
                    "rt_start": None,
                    "rt_stop": None,
                }
                records.append(rec)

        return records
