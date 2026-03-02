"""Unit tests for converter utility methods added in recent commits.

Tests cover:
- BaseConverter._escape_path (SQL path injection prevention)
- MaxQuantPgAdapter._extract_protein_names / _extract_gene_names
- MaxQuantFeatureAdapter._detect_pg_columns / _extract_qvalue_map / _extract_gene_map
"""

import pandas as pd
import pytest

from qpx.converters.maxquant.pg_adapter import MaxQuantPgAdapter
from qpx.converters.maxquant.feature_adapter import MaxQuantFeatureAdapter
from qpx.converters.base import BaseConverter

# ---------------------------------------------------------------------------
# BaseConverter._escape_path
# ---------------------------------------------------------------------------


class TestEscapePath:
    """Validate SQL path escaping."""

    def test_normal_path(self):
        assert BaseConverter._escape_path("/data/results.txt") == "/data/results.txt"

    def test_path_with_single_quote(self):
        assert BaseConverter._escape_path("O'Brien/data.txt") == "O''Brien/data.txt"

    def test_path_with_multiple_quotes(self):
        assert BaseConverter._escape_path("a'b'c.txt") == "a''b''c.txt"

    def test_windows_path(self):
        assert (
            BaseConverter._escape_path("C:\\Users\\test\\file.txt")
            == "C:\\Users\\test\\file.txt"
        )

    def test_empty_path(self):
        assert BaseConverter._escape_path("") == ""


# ---------------------------------------------------------------------------
# MaxQuantPgAdapter._extract_protein_names
# ---------------------------------------------------------------------------


class TestExtractProteinNames:
    """Validate UniProt-style accession parsing."""

    def test_uniprot_sp_format(self):
        result = MaxQuantPgAdapter._extract_protein_names(
            ["sp|P12345|PROT_HUMAN", "sp|Q99999|ANOT_MOUSE"]
        )
        assert result == ["PROT_HUMAN", "ANOT_MOUSE"]

    def test_bare_accession(self):
        result = MaxQuantPgAdapter._extract_protein_names(["P12345"])
        assert result == ["P12345"]

    def test_two_pipe_parts(self):
        result = MaxQuantPgAdapter._extract_protein_names(["tr|A0A0A0"])
        assert result == ["A0A0A0"]

    def test_empty_list(self):
        result = MaxQuantPgAdapter._extract_protein_names([])
        assert result is None

    def test_mixed_formats(self):
        result = MaxQuantPgAdapter._extract_protein_names(
            ["sp|P12345|PROT_HUMAN", "CONTAM_Q99"]
        )
        assert result == ["PROT_HUMAN", "CONTAM_Q99"]


# ---------------------------------------------------------------------------
# MaxQuantPgAdapter._extract_gene_names
# ---------------------------------------------------------------------------


class TestExtractGeneNames:
    """Validate GN= tag extraction from FASTA headers."""

    def test_single_gene(self):
        result = MaxQuantPgAdapter._extract_gene_names(
            "sp|P12345|PROT_HUMAN Protein X OS=Homo sapiens GN=TP53 PE=1 SV=1"
        )
        assert result == ["TP53"]

    def test_multiple_headers(self):
        result = MaxQuantPgAdapter._extract_gene_names(
            "sp|P12345|X GN=TP53 PE=1;sp|Q99999|Y GN=BRCA1 PE=1"
        )
        assert result == ["TP53", "BRCA1"]

    def test_no_gene_name(self):
        result = MaxQuantPgAdapter._extract_gene_names(
            "sp|P12345|PROT_HUMAN Protein X OS=Homo sapiens PE=1 SV=1"
        )
        assert result is None

    def test_duplicate_gene_names(self):
        result = MaxQuantPgAdapter._extract_gene_names(
            "sp|P12345|X GN=TP53 PE=1;sp|Q99999|Y GN=TP53 PE=1"
        )
        assert result == ["TP53"]

    def test_empty_string(self):
        result = MaxQuantPgAdapter._extract_gene_names("")
        assert result is None


# ---------------------------------------------------------------------------
# MaxQuantFeatureAdapter._detect_pg_columns
# ---------------------------------------------------------------------------


class TestDetectPgColumns:
    """Validate column detection in proteinGroups DataFrame."""

    def test_standard_columns(self):
        df = pd.DataFrame(
            {"Protein IDs": ["P1"], "Q-value": [0.01], "Gene names": ["TP53"]}
        )
        acc, qval, gene = MaxQuantFeatureAdapter._detect_pg_columns(df)
        assert acc == "Protein IDs"
        assert qval == "Q-value"
        assert gene == "Gene names"

    def test_alternative_columns(self):
        df = pd.DataFrame(
            {"Majority protein IDs": ["P1"], "q-value": [0.01], "Gene Names": ["TP53"]}
        )
        acc, qval, gene = MaxQuantFeatureAdapter._detect_pg_columns(df)
        assert acc == "Majority protein IDs"
        assert qval == "q-value"
        assert gene == "Gene Names"

    def test_missing_columns(self):
        df = pd.DataFrame({"some_col": [1]})
        acc, qval, gene = MaxQuantFeatureAdapter._detect_pg_columns(df)
        assert acc is None
        assert qval is None
        assert gene is None


# ---------------------------------------------------------------------------
# MaxQuantFeatureAdapter._extract_qvalue_map
# ---------------------------------------------------------------------------


class TestExtractQvalueMap:
    """Validate q-value map extraction."""

    def test_basic_qvalue_map(self):
        df = pd.DataFrame({"Protein IDs": ["P1;P2", "P3"], "Q-value": [0.01, 0.05]})
        result = MaxQuantFeatureAdapter._extract_qvalue_map(
            df, "Protein IDs", "Q-value"
        )
        assert result["P1"] == 0.01
        assert result["P2"] == 0.01
        assert result["P3"] == 0.05

    def test_first_occurrence_wins(self):
        df = pd.DataFrame({"Protein IDs": ["P1", "P1"], "Q-value": [0.01, 0.99]})
        result = MaxQuantFeatureAdapter._extract_qvalue_map(
            df, "Protein IDs", "Q-value"
        )
        assert result["P1"] == 0.01

    def test_no_qval_col(self):
        df = pd.DataFrame({"Protein IDs": ["P1"]})
        result = MaxQuantFeatureAdapter._extract_qvalue_map(df, "Protein IDs", None)
        assert result == {}


# ---------------------------------------------------------------------------
# MaxQuantFeatureAdapter._extract_gene_map
# ---------------------------------------------------------------------------


class TestExtractGeneMap:
    """Validate gene map extraction with Fasta header fallback."""

    def test_from_gene_column(self):
        df = pd.DataFrame({"Protein IDs": ["P1;P2"], "Gene names": ["TP53;BRCA1"]})
        result = MaxQuantFeatureAdapter._extract_gene_map(
            df, "Protein IDs", "Gene names"
        )
        assert result["P1"] == ["TP53", "BRCA1"]
        assert result["P2"] == ["TP53", "BRCA1"]

    def test_fasta_fallback(self):
        df = pd.DataFrame(
            {
                "Protein IDs": ["P1"],
                "Gene names": [None],
                "Fasta headers": ["sp|P1|X GN=TP53 PE=1"],
            }
        )
        result = MaxQuantFeatureAdapter._extract_gene_map(
            df, "Protein IDs", "Gene names"
        )
        assert result["P1"] == ["TP53"]

    def test_no_gene_info(self):
        df = pd.DataFrame({"Protein IDs": ["P1"], "Gene names": [None]})
        result = MaxQuantFeatureAdapter._extract_gene_map(
            df, "Protein IDs", "Gene names"
        )
        assert result == {}
