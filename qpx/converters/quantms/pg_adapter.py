"""QuantMS PG adapter -- mzTab + MSstats to pg.parquet.

Loads an mzTab file and MSstats input into DuckDB, performs SQL-based
protein-group quantification, and writes the results through ``PgWriter``.

Key schema changes:
    - ``reference_file_name`` -> ``run_file_name``
    - ``is_decoy`` (int, 0/1) -> ``is_decoy`` (bool)
    - Intensities struct: ``{sample_accession, channel, intensity}`` ->
      ``{label, intensity}``
"""

from __future__ import annotations

import logging
import re
from typing import Optional

import pandas as pd

from qpx.converters.base import BaseConverter
from qpx.converters.utils import safe_float
from qpx.converters.mztab import (
    load_mztab_sections,
    load_msstats,
    extract_ms_runs,
)
from qpx.writers.pg import PgWriter

logger = logging.getLogger(__name__)


class QuantmsPgAdapter(BaseConverter):
    """Convert QuantMS mzTab + MSstats data to ``pg.parquet``.

    Usage::

        with QuantmsPgAdapter() as adapter:
            adapter.convert(
                mztab_path="results.mzTab",
                msstats_path="msstats_in.csv",
                output_path="pg.parquet",
            )
    """

    def convert(
        self,
        mztab_path: str,
        msstats_path: str,
        output_path: str,
        sdrf_path: Optional[str] = None,
        compute_ibaq: bool = True,
        file_batch_size: int = 2,
        creator: str = "quantms",
    ) -> None:
        """Run the mzTab+MSstats -> pg.parquet conversion.

        Args:
            mztab_path: Path to the mzTab file.
            msstats_path: Path to the MSstats input file.
            output_path: Destination Parquet path.
            sdrf_path: Optional SDRF for extra sample mapping.
            compute_ibaq: Whether to compute iBAQ intensities.
            file_batch_size: Number of reference files per batch.
            creator: Creator tag in Parquet metadata.
        """
        # Step 1: Load into DuckDB (skip if already loaded)
        if not self._table_exists("psms"):
            load_mztab_sections(self._conn, mztab_path)
        if not self._table_exists("msstats"):
            load_msstats(self._conn, msstats_path)

        # Step 2: Load protein-group info and create SQL join
        protein_groups = self._load_protein_groups()
        self._create_protein_group_table(protein_groups)
        self._create_joined_view()

        # Step 3: Process in batches per reference file
        self.logger.info("Quantifying protein groups ...")

        with PgWriter(output_path, creator=creator) as writer:
            for file_batch in self._iter_file_batches(file_batch_size):
                batch_data = self._fetch_batch(file_batch)
                if batch_data.empty:
                    continue
                records = self._build_pg_records(batch_data, compute_ibaq)
                if records:
                    self._track_scores(records)
                    writer.write_batch(records)

        self.logger.info(f"PG conversion complete -> {output_path}")

    # ------------------------------------------------------------------
    # Protein group loading
    # ------------------------------------------------------------------

    def _load_protein_groups(self) -> list[dict]:
        """Load protein-group metadata from the mzTab proteins table."""
        groups: list[dict] = []
        try:
            df = self._conn.execute("SELECT * FROM proteins").df()
            for row in df.to_dict("records"):
                result_type = str(
                    row.get("opt_global_result_type", "single_protein")
                )
                accession = str(row.get("accession", ""))
                if not accession or accession == "null":
                    continue
                if result_type not in ("single_protein", "indistinguishable_protein_group"):
                    continue

                description = row.get("description", "")
                description = description if pd.notna(description) and description != "null" else None

                # Gene names from description
                gg_accessions = []
                if description:
                    gn = re.search(r"GN=([^\s]+)", str(description))
                    if gn:
                        gg_accessions = [gn.group(1)]

                # Global q-value
                global_qvalue = safe_float(
                    row.get("best_search_engine_score[1]",
                            row.get("best_search_engine_score_1", 0))
                ) or 0.0

                # Peptide count
                peptide_count = 0
                pep_raw = row.get("opt_global_nr_found_peptides")
                if pep_raw not in (None, "", "null"):
                    try:
                        peptide_count = int(float(pep_raw))
                    except (ValueError, TypeError):
                        pass

                # Decoy
                is_decoy_raw = row.get("opt_global_cv_pride:0000303_decoy_hit",
                                       row.get("opt_global_cv_PRIDE:0000303_decoy_hit", 0))
                is_decoy = str(is_decoy_raw).strip() in ("1", "true", "True")

                # Sequence coverage
                seq_cov = safe_float(row.get("protein_coverage")) or 0.0

                # Protein name extraction
                pg_names = self._extract_protein_names(accession)

                groups.append(
                    {
                        "pg_accessions": accession,
                        "anchor_protein": accession.split(";")[0] if accession else "",
                        "pg_names": pg_names,
                        "gg_accessions": ";".join(gg_accessions) if gg_accessions else None,
                        "description": description,
                        "global_qvalue": global_qvalue,
                        "peptide_count": peptide_count,
                        "is_decoy": is_decoy,
                        "sequence_coverage": seq_cov,
                        "sequence_length": safe_float(row.get("sequence_length")) or 0,
                    }
                )
        except Exception as e:
            self.logger.error(f"Error loading protein groups: {e}")

        return groups

    def _create_protein_group_table(self, groups: list[dict]) -> None:
        """Load protein groups into a DuckDB table."""
        if not groups:
            return
        df = pd.DataFrame(groups)
        self._conn.execute("DROP TABLE IF EXISTS protein_groups")
        self._conn.from_df(df).create("protein_groups")
        self._conn.execute(
            "CREATE INDEX IF NOT EXISTS idx_pg_acc ON protein_groups(pg_accessions)"
        )

    def _create_joined_view(self) -> None:
        """Create a joined view between MSstats and protein groups."""
        try:
            # Detect column names
            ref_col = self._detect_column("Reference", "reference", "Run")
            protein_col = self._detect_column("ProteinName", "pg_accessions", "Protein")
            pform_col = self._detect_column("PeptideSequence", "peptidoform", "Peptide")

            self._conn.execute(f"""
                DROP VIEW IF EXISTS msstats_pg;
                CREATE VIEW msstats_pg AS
                SELECT
                    m.*,
                    COALESCE(pg.anchor_protein, m."{protein_col}") AS anchor_protein,
                    pg.pg_names,
                    pg.gg_accessions,
                    pg.global_qvalue,
                    pg.peptide_count,
                    pg.is_decoy,
                    pg.sequence_coverage,
                    pg.sequence_length
                FROM msstats m
                LEFT JOIN protein_groups pg ON m."{protein_col}" = pg.pg_accessions
            """)
        except Exception as e:
            self.logger.error(f"Error creating joined view: {e}")
            raise

    def _detect_column(self, *candidates: str) -> str:
        """Return the first column name that exists in the msstats table."""
        cols = self._conn.execute(
            "SELECT column_name FROM information_schema.columns WHERE table_name='msstats'"
        ).fetchall()
        existing = {c[0] for c in cols}
        for c in candidates:
            if c in existing:
                return c
        return candidates[0]

    # ------------------------------------------------------------------
    # Batched processing
    # ------------------------------------------------------------------

    def _iter_file_batches(self, batch_size: int):
        """Yield lists of reference file names in batches."""
        ref_col = self._detect_column("Reference", "reference", "Run")
        refs = self._conn.execute(
            f'SELECT DISTINCT "{ref_col}" FROM msstats ORDER BY "{ref_col}"'
        ).fetchall()
        ref_list = [r[0] for r in refs]
        for i in range(0, len(ref_list), batch_size):
            yield ref_list[i : i + batch_size]

    def _fetch_batch(self, file_batch: list[str]) -> pd.DataFrame:
        """Fetch joined data for a batch of reference files."""
        ref_col = self._detect_column("Reference", "reference", "Run")
        placeholders = ", ".join(["?" for _ in file_batch])
        return self._conn.execute(
            f"""
            SELECT * FROM msstats_pg
            WHERE "{ref_col}" IN ({placeholders})
              AND anchor_protein IS NOT NULL
            """,
            file_batch,
        ).df()

    # ------------------------------------------------------------------
    # Record building
    # ------------------------------------------------------------------

    def _build_pg_records(
        self,
        batch_data: pd.DataFrame,
        compute_ibaq: bool,
    ) -> list[dict]:
        """Build PG records from a batch of joined MSstats+protein data."""
        records: list[dict] = []

        ref_col = self._detect_batch_ref_col(batch_data)
        channel_col = self._detect_batch_col(batch_data, "Channel", "channel")
        intensity_col = self._detect_batch_col(batch_data, "Intensity", "intensity")
        pform_col = self._detect_batch_col(batch_data, "PeptideSequence", "peptidoform")
        charge_col = self._detect_batch_col(batch_data, "Charge", "charge")

        for (anchor, ref_file), group in batch_data.groupby(
            ["anchor_protein", ref_col]
        ):
            try:
                rec = self._build_single_pg(
                    anchor_protein=str(anchor),
                    run_file_name=str(ref_file).split(".")[0],
                    group=group,
                    channel_col=channel_col,
                    intensity_col=intensity_col,
                    pform_col=pform_col,
                    charge_col=charge_col,
                    compute_ibaq=compute_ibaq,
                )
                if rec:
                    records.append(rec)
            except Exception as e:
                self.logger.debug(f"Skipping PG group: {e}")

        return records

    def _build_single_pg(
        self,
        anchor_protein: str,
        run_file_name: str,
        group: pd.DataFrame,
        channel_col: str,
        intensity_col: str,
        pform_col: str,
        charge_col: str,
        compute_ibaq: bool,
    ) -> Optional[dict]:
        """Build a single PG record from a group."""

        # Protein group info
        pg_accessions = str(group["pg_accessions"].iloc[0] if "pg_accessions" in group.columns
                            else anchor_protein).split(";")
        pg_names_raw = group.get("pg_names")
        pg_names = (
            str(pg_names_raw.iloc[0]).split(";")
            if pg_names_raw is not None and pd.notna(pg_names_raw.iloc[0])
            else None
        )

        gg_acc_raw = group.get("gg_accessions")
        gg_accessions = None
        if gg_acc_raw is not None and pd.notna(gg_acc_raw.iloc[0]):
            val = gg_acc_raw.iloc[0]
            if isinstance(val, str) and val:
                gg_accessions = val.split(";")

        global_qvalue = safe_float(group["global_qvalue"].iloc[0]) if "global_qvalue" in group.columns else None
        is_decoy = bool(group["is_decoy"].iloc[0]) if "is_decoy" in group.columns else False
        seq_coverage = safe_float(group["sequence_coverage"].iloc[0]) if "sequence_coverage" in group.columns else None

        # Peptide and feature counts
        if pform_col in group.columns:
            total_sequences = group[pform_col].nunique()
        else:
            total_sequences = 0
        unique_sequences = total_sequences  # simplified

        if charge_col in group.columns:
            total_features = len(set(zip(group[pform_col], group[charge_col])))
        else:
            total_features = total_sequences
        unique_features = total_features

        # Peptide count from protein group metadata
        pep_count_raw = group.get("peptide_count")
        peptide_count_val = int(pep_count_raw.iloc[0]) if pep_count_raw is not None and pd.notna(pep_count_raw.iloc[0]) else total_sequences

        # Build intensities (new schema: {label, intensity})
        # TMT channel index → label mapping
        _TMT_MAP = {
            1: "TMT126", 2: "TMT127N", 3: "TMT127C",
            4: "TMT128N", 5: "TMT128C", 6: "TMT129N",
            7: "TMT129C", 8: "TMT130N", 9: "TMT130C",
            10: "TMT131N", 11: "TMT131C",
            12: "TMT132N", 13: "TMT132C", 14: "TMT133N",
            15: "TMT133C", 16: "TMT134N", 17: "TMT134C", 18: "TMT135N",
        }
        intensities = []
        additional_intensities = []
        for channel_val, ch_group in group.groupby(channel_col):
            # Map numeric TMT channels to label names
            try:
                ch_int = int(float(channel_val))
                label = _TMT_MAP.get(ch_int, str(channel_val))
            except (ValueError, TypeError):
                label = str(channel_val) if pd.notna(channel_val) else "LFQ"
            sum_intensity = float(ch_group[intensity_col].sum())

            intensities.append({"label": label, "intensity": sum_intensity})

            # Additional intensities (iBAQ if requested)
            extra = []
            if compute_ibaq:
                seq_len = group.get("sequence_length")
                if seq_len is not None and pd.notna(seq_len.iloc[0]) and float(seq_len.iloc[0]) > 0:
                    ibaq = sum_intensity / max(float(seq_len.iloc[0]), 1)
                    extra.append({"intensity_name": "ibaq", "intensity_value": ibaq})

            if extra:
                additional_intensities.append({"label": label, "intensities": extra})

        # Peptides per protein
        peptides = [
            {"protein_name": acc, "peptide_count": peptide_count_val}
            for acc in pg_accessions
        ]

        # Additional scores
        additional_scores = []
        if seq_coverage:
            additional_scores.append(
                {"score_name": "sequence_coverage_percent", "score_value": seq_coverage, "higher_better": True}
            )
        additional_scores.append(
            {"score_name": "peptide_count", "score_value": float(peptide_count_val), "higher_better": True}
        )

        return {
            "pg_accessions": pg_accessions,
            "pg_names": pg_names,
            "gg_accessions": gg_accessions,
            "gg_names": None,
            "anchor_protein": anchor_protein.split(";")[0],
            "run_file_name": run_file_name,
            "global_qvalue": global_qvalue,
            "pg_qvalue": None,
            "intensities": intensities or None,
            "additional_intensities": additional_intensities or None,
            "is_decoy": is_decoy,
            "contaminant": 0,
            "peptides": peptides,
            "peptide_counts": {
                "unique_sequences": unique_sequences,
                "total_sequences": total_sequences,
            },
            "feature_counts": {
                "unique_features": unique_features,
                "total_features": total_features,
            },
            "sequence_coverage": seq_coverage,
            "molecular_weight": None,
            "additional_scores": additional_scores or None,
            "cv_params": None,
        }

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _extract_protein_names(accession: str) -> str:
        """Extract protein names from UniProt-style accession."""
        parts_list = []
        for acc in accession.split(";"):
            if "|" in acc:
                pieces = acc.split("|")
                if len(pieces) >= 3:
                    parts_list.append(pieces[2])
                elif len(pieces) >= 2:
                    parts_list.append(pieces[1])
                else:
                    parts_list.append(acc)
            else:
                parts_list.append(acc)
        return ";".join(parts_list)

    def _detect_batch_ref_col(self, df: pd.DataFrame) -> str:
        for c in ["Reference", "reference", "Run", "run", "reference_file_name"]:
            if c in df.columns:
                return c
        return "Reference"

    def _detect_batch_col(self, df: pd.DataFrame, *candidates: str) -> str:
        for c in candidates:
            if c in df.columns:
                return c
        return candidates[0]
