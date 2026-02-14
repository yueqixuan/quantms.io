"""Peptide-protein mapping model — loaded from schemas/peptide_protein_map.yaml."""
from qpx.core.models.loader import load_schema

PeptideProteinMapSchema = load_schema("peptide_protein_map")
