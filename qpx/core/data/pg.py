"""PG data structure — Protein Groups with quantification."""

from qpx.core.convert import QueryResult
from qpx.core.data.base import BaseStructure
from qpx.core.data.loader import load_schema
from qpx.core.query import _escape_sql_string

PgSchema = load_schema("pg")


class PG(BaseStructure):
    """Protein groups with per-run quantification."""

    _schema_class = PgSchema

    def by_protein(self, protein: str) -> "PG":
        """Filter by anchor protein."""
        return self.filter(f"anchor_protein = '{_escape_sql_string(protein)}'")

    def by_run(self, run_file_name: str) -> "PG":
        """Filter by run file."""
        return self.filter(f"run_file_name = '{_escape_sql_string(run_file_name)}'")

    def targets_only(self) -> "PG":
        """Filter to target proteins only (exclude decoys)."""
        return self.filter("is_decoy = false")

    def _intensity_label_field(self) -> str:
        """Detect whether the intensities struct uses 'label' (new) or 'channel' (old)."""
        try:
            row = self._engine.execute(f"SELECT typeof(intensities) FROM {self._query.source} LIMIT 1").fetchone()
            if row:
                type_str = row[0].lower()
                if "channel" in type_str and "label" not in type_str:
                    return "channel"
        except Exception:
            pass
        return "label"

    def protein_intensities(self) -> QueryResult:
        """
        Flatten intensities: one row per (protein, label, intensity).
        """
        ilf = self._intensity_label_field()
        sql = f"""
        SELECT anchor_protein, pg_accessions, gg_names, run_file_name,
               global_qvalue, i.{ilf} AS label, i.intensity
        FROM {self._query.source},
             UNNEST(intensities) AS _t(i)
        """
        return QueryResult(self._engine.execute(sql))
