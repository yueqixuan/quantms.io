"""Peptide-to-protein mapping structure."""

from qpx.core.data.base import BaseStructure
from qpx.core.models.peptide_protein_map import PeptideProteinMapSchema


class PeptideProteinMap(BaseStructure):
    """Deduplicated peptide-to-protein mapping."""

    _schema_class = PeptideProteinMapSchema

    def by_protein(self, protein: str) -> "PeptideProteinMap":
        """Filter mappings by protein accession."""
        return self.filter(f"protein_accession = '{protein}'")

    def by_peptide(self, sequence: str) -> "PeptideProteinMap":
        """Filter mappings by peptide sequence."""
        return self.filter(f"sequence = '{sequence}'")

    def unique_peptides(self) -> "PeptideProteinMap":
        """Filter to peptides that map to exactly one protein."""
        return self.filter("is_unique = true")

    def proteotypic(self) -> "PeptideProteinMap":
        """Filter to proteotypic peptides (unique to one gene)."""
        return self.filter("is_proteotypic = true")
