"""PeptideProteinMapWriter — writes peptide_protein_map.parquet files."""

from qpx.writers.base import BaseWriter
from qpx.core.models.peptide_protein_map import PeptideProteinMapSchema


class PeptideProteinMapWriter(BaseWriter):
    _schema_class = PeptideProteinMapSchema
