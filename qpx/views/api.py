"""Computed views — cross-structure projections using shared DuckDB connection."""

from __future__ import annotations
from typing import TYPE_CHECKING

from qpx.views.base import BaseView
from qpx.core.convert import QueryResult

if TYPE_CHECKING:
    from qpx.dataset import Dataset


class ProteinView(BaseView):
    """Protein abundance per sample — joins PG + Run."""

    def abundance(self, q_value_threshold: float = 0.01) -> QueryResult:
        q_value_threshold = float(q_value_threshold)
        sql = """
        SELECT pg.anchor_protein AS protein_accession,
               rs.sample_accession,
               i.label,
               i.intensity AS abundance,
               pg.global_qvalue,
               pg.gg_names AS gene_names
        FROM pg,
             run r,
             UNNEST(r.samples) AS _t1(rs),
             UNNEST(pg.intensities) AS _t2(i)
        WHERE pg.run_file_name = r.run_file_name
          AND i.label = rs.label
          AND (pg.global_qvalue IS NULL OR pg.global_qvalue <= $1)
        """
        return QueryResult(self._engine.execute(sql, [q_value_threshold]))


class PeptideView(BaseView):
    """Peptide abundance per sample — joins Feature + Run."""

    def abundance(self) -> QueryResult:
        sql = """
        SELECT f.sequence, f.peptidoform, f.charge, f.anchor_protein,
               rs.sample_accession, i.label, i.intensity AS abundance
        FROM feature f,
             run r,
             UNNEST(r.samples) AS _t1(rs),
             UNNEST(f.intensities) AS _t2(i)
        WHERE f.run_file_name = r.run_file_name
          AND i.label = rs.label
        """
        return QueryResult(self._engine.execute(sql))


class IdentificationSummaryView(BaseView):
    """Per-run summary: proteins, peptides, features."""

    def summary(self) -> QueryResult:
        sql = """
        SELECT f.run_file_name,
               COUNT(DISTINCT f.anchor_protein) AS n_proteins,
               COUNT(DISTINCT f.sequence) AS n_peptides,
               COUNT(*) AS n_features
        FROM feature f
        GROUP BY f.run_file_name
        """
        return QueryResult(self._engine.execute(sql))

    def plot(self, figsize=(12, 6)):
        """Bar chart of identifications per run."""
        from qpx.views.plotting import grouped_bar_chart

        df = self.summary().to_df()
        return grouped_bar_chart(
            df,
            x="run_file_name",
            y_cols=["n_proteins", "n_peptides"],
            title="Identifications per Run",
            xlabel="Run",
            ylabel="Count",
            figsize=figsize,
        )


class RunSummaryView(BaseView):
    """Per-run statistics with intensity distribution."""

    def summary(self):
        sql = """
        SELECT f.run_file_name,
               COUNT(DISTINCT f.sequence) AS n_peptides,
               COUNT(DISTINCT f.anchor_protein) AS n_proteins,
               COUNT(*) AS n_features
        FROM feature f
        GROUP BY f.run_file_name
        """
        return self._execute_cached("run_summary", sql)

    def plot(self, figsize=(12, 6)):
        """Bar chart of run-level statistics."""
        from qpx.views.plotting import grouped_bar_chart

        df = self.summary().to_df()
        return grouped_bar_chart(
            df,
            x="run_file_name",
            y_cols=["n_peptides", "n_proteins", "n_features"],
            title="Run Summary",
            xlabel="Run",
            ylabel="Count",
            figsize=figsize,
        )


class ModificationView(BaseView):
    """Modification frequency across PSMs."""

    def frequency(self):
        sql = """
        SELECT m.name AS modification_name,
               m.accession AS modification_accession,
               COUNT(*) AS n_psms,
               COUNT(DISTINCT p.sequence) AS n_peptides
        FROM psm p, UNNEST(p.modifications) AS _t(m)
        WHERE p.modifications IS NOT NULL
        GROUP BY m.name, m.accession
        ORDER BY n_psms DESC
        """
        return self._execute_cached("mod_frequency", sql)

    def plot(self, top_n=20, figsize=(10, 6)):
        """Bar chart of top modifications by PSM count."""
        from qpx.views.plotting import bar_chart

        df = self.frequency().to_df().head(top_n)
        return bar_chart(
            df,
            x="modification_name",
            y="n_psms",
            title=f"Top {top_n} Modifications by PSM Count",
            xlabel="Modification",
            ylabel="PSM Count",
            figsize=figsize,
        )


class QualityControlView(BaseView):
    """Dataset-wide QC metrics."""

    def metrics(self):
        sql = """
        SELECT COUNT(DISTINCT f.sequence) AS unique_peptides,
               COUNT(DISTINCT f.anchor_protein) AS unique_proteins,
               COUNT(DISTINCT f.run_file_name) AS n_runs,
               COUNT(*) AS total_features
        FROM feature f
        WHERE f.is_decoy = false
        """
        return self._execute_cached("qc_metrics", sql)

    def plot(self, figsize=(8, 5)):
        """Bar chart of dataset QC summary metrics."""
        import pandas as pd
        from qpx.views.plotting import bar_chart

        df = self.metrics().to_df()
        summary = pd.DataFrame(
            {
                "metric": [
                    "Unique Peptides",
                    "Unique Proteins",
                    "Runs",
                    "Total Features",
                ],
                "value": [
                    df["unique_peptides"].iloc[0],
                    df["unique_proteins"].iloc[0],
                    df["n_runs"].iloc[0],
                    df["total_features"].iloc[0],
                ],
            }
        )
        return bar_chart(
            summary,
            x="metric",
            y="value",
            title="Dataset QC Summary",
            figsize=figsize,
            rotation=0,
        )
