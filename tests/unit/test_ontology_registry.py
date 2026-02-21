"""Tests for PublicOntology registry and OBO parser."""

import tempfile
from pathlib import Path

import pytest

from qpx.core.obo import (
    CVTerm,
    _normalize_name,
    extract_obo_version,
    parse_obo,
    terms_to_arrow,
    write_terms_parquet,
    read_parquet_ontology_metadata,
)
from qpx.core.ontology import PublicOntology

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

# Minimal OBO snippet for offline testing
_MINI_OBO = """\
format-version: 1.2
data-version: 4.1.235
ontology: ms

[Term]
id: MS:1001143
name: PSM-level search engine specific statistic
def: "A score or other statistic." [PSI:MS]

[Term]
id: MS:1001492
name: percolator:score
def: "The Percolator score." [PSI:MS]
is_a: MS:1001143 ! PSM-level search engine specific statistic
relationship: has_order MS:1002108 ! higher score is better
synonym: "percolator score" EXACT []

[Term]
id: MS:1001493
name: percolator:PEP
def: "The Percolator posterior error probability." [PSI:MS]
is_a: MS:1001143 ! PSM-level search engine specific statistic
relationship: has_order MS:1002109 ! lower score is better

[Term]
id: MS:1002523
name: Q Exactive HF
def: "Thermo Scientific Q Exactive HF." [PSI:MS]

[Term]
id: MS:9999999
name: obsolete term
def: "This term is obsolete." [PSI:MS]
is_obsolete: true
"""


@pytest.fixture
def mini_obo_terms():
    """Parse the mini OBO snippet."""
    return parse_obo(_MINI_OBO, source="MS")


@pytest.fixture
def mini_parquet(mini_obo_terms, tmp_path):
    """Write mini OBO terms to a temp Parquet file."""
    path = tmp_path / "test_ms.parquet"
    write_terms_parquet(mini_obo_terms, path, source="test_ms", version="4.1.235")
    return path


@pytest.fixture
def mini_ontology(mini_parquet):
    """Create a PublicOntology from the mini Parquet file."""
    return PublicOntology.from_parquet(mini_parquet)


# ---------------------------------------------------------------------------
# OBO parser tests
# ---------------------------------------------------------------------------


class TestNormalizeName:
    def test_basic(self):
        assert _normalize_name("percolator:score") == "percolator_score"

    def test_spaces_and_dots(self):
        assert _normalize_name("Comet.xcorr") == "comet_xcorr"

    def test_mixed_case(self):
        assert _normalize_name("MSFragger:Hyperscore") == "msfragger_hyperscore"

    def test_already_normalized(self):
        assert _normalize_name("andromeda_score") == "andromeda_score"

    def test_special_chars(self):
        assert _normalize_name("x!tandem:expect") == "x_tandem_expect"


class TestExtractOboVersion:
    def test_data_version(self):
        assert extract_obo_version(_MINI_OBO) == "4.1.235"

    def test_remark_version(self):
        text = "remark: version: 1.2.3\n[Term]"
        assert extract_obo_version(text) == "1.2.3"

    def test_no_version(self):
        assert extract_obo_version("[Term]\nid: X:1") == "unknown"


class TestParseObo:
    def test_term_count(self, mini_obo_terms):
        assert len(mini_obo_terms) == 5

    def test_accession(self, mini_obo_terms):
        terms_by_acc = {t.accession: t for t in mini_obo_terms}
        assert "MS:1001492" in terms_by_acc

    def test_name(self, mini_obo_terms):
        terms_by_acc = {t.accession: t for t in mini_obo_terms}
        assert terms_by_acc["MS:1001492"].name == "percolator:score"

    def test_normalized_name(self, mini_obo_terms):
        terms_by_acc = {t.accession: t for t in mini_obo_terms}
        assert terms_by_acc["MS:1001492"].normalized_name == "percolator_score"

    def test_definition(self, mini_obo_terms):
        terms_by_acc = {t.accession: t for t in mini_obo_terms}
        assert "Percolator score" in terms_by_acc["MS:1001492"].definition

    def test_is_score(self, mini_obo_terms):
        terms_by_acc = {t.accession: t for t in mini_obo_terms}
        assert terms_by_acc["MS:1001492"].is_score is True
        assert terms_by_acc["MS:1002523"].is_score is False

    def test_higher_better(self, mini_obo_terms):
        terms_by_acc = {t.accession: t for t in mini_obo_terms}
        assert terms_by_acc["MS:1001492"].higher_better is True
        assert terms_by_acc["MS:1001493"].higher_better is False

    def test_is_a_parents(self, mini_obo_terms):
        terms_by_acc = {t.accession: t for t in mini_obo_terms}
        assert "MS:1001143" in terms_by_acc["MS:1001492"].is_a

    def test_synonyms(self, mini_obo_terms):
        terms_by_acc = {t.accession: t for t in mini_obo_terms}
        assert "percolator score" in terms_by_acc["MS:1001492"].synonyms

    def test_obsolete(self, mini_obo_terms):
        terms_by_acc = {t.accession: t for t in mini_obo_terms}
        assert terms_by_acc["MS:9999999"].is_obsolete is True
        assert terms_by_acc["MS:1001492"].is_obsolete is False


# ---------------------------------------------------------------------------
# Parquet round-trip tests
# ---------------------------------------------------------------------------


class TestParquetRoundTrip:
    def test_write_and_read_metadata(self, mini_parquet):
        meta = read_parquet_ontology_metadata(mini_parquet)
        assert meta["qpx_ontology_source"] == "test_ms"
        assert meta["qpx_ontology_version"] == "4.1.235"

    def test_arrow_table_schema(self, mini_obo_terms):
        table = terms_to_arrow(mini_obo_terms)
        assert "accession" in table.schema.names
        assert "normalized_name" in table.schema.names
        assert "is_score" in table.schema.names
        assert "is_obsolete" in table.schema.names
        assert len(table) == 5


# ---------------------------------------------------------------------------
# PublicOntology API tests
# ---------------------------------------------------------------------------


class TestPublicOntology:
    def test_version(self, mini_ontology):
        assert mini_ontology.version == "4.1.235"

    def test_len(self, mini_ontology):
        # 5 terms total, 1 obsolete → 4 non-obsolete
        assert len(mini_ontology) == 4

    def test_search(self, mini_ontology):
        df = mini_ontology.search("percolator")
        assert len(df) >= 2
        assert "percolator:score" in df["name"].values

    def test_lookup(self, mini_ontology):
        term = mini_ontology.lookup("MS:1001492")
        assert term is not None
        assert term["name"] == "percolator:score"

    def test_lookup_missing(self, mini_ontology):
        assert mini_ontology.lookup("MS:0000000") is None

    def test_lookup_name(self, mini_ontology):
        term = mini_ontology.lookup_name("percolator:score")
        assert term is not None
        assert term["accession"] == "MS:1001492"

    def test_lookup_name_case_insensitive(self, mini_ontology):
        term = mini_ontology.lookup_name("Percolator:Score")
        assert term is not None

    def test_lookup_normalized(self, mini_ontology):
        term = mini_ontology.lookup_normalized("percolator_score")
        assert term is not None
        assert term["accession"] == "MS:1001492"

    def test_validate(self, mini_ontology):
        result = mini_ontology.validate(
            ["percolator_score", "unknown_score", "percolator_pep"]
        )
        assert result["percolator_score"] == True
        assert result["unknown_score"] == False
        assert result["percolator_pep"] == True

    def test_validate_empty(self, mini_ontology):
        result = mini_ontology.validate([])
        assert len(result) == 0

    def test_standardize(self, mini_ontology):
        acc = mini_ontology.standardize("percolator:score")
        assert acc == "MS:1001492"

    def test_standardize_missing(self, mini_ontology):
        assert mini_ontology.standardize("nonexistent_term") is None

    def test_scores(self, mini_ontology):
        scores_df = mini_ontology.scores()
        assert len(scores_df) >= 2
        accessions = scores_df["accession"].tolist()
        assert "MS:1001492" in accessions
        assert "MS:1001493" in accessions

    def test_df(self, mini_ontology):
        df = mini_ontology.df()
        assert len(df) == 4  # excludes obsolete

    def test_df_include_obsolete(self, mini_ontology):
        df = mini_ontology.df(include_obsolete=True)
        assert len(df) == 5

    def test_from_obo(self):
        with tempfile.NamedTemporaryFile(suffix=".obo", mode="w", delete=False) as f:
            f.write(_MINI_OBO)
            f.flush()
            onto = PublicOntology.from_obo(f.name, source="MS")
            assert onto.version == "4.1.235"
            assert len(onto) == 4
            onto.close()

    def test_context_manager(self, mini_parquet):
        with PublicOntology.from_parquet(mini_parquet) as onto:
            assert len(onto) == 4

    def test_repr(self, mini_ontology):
        r = repr(mini_ontology)
        assert "PublicOntology" in r
        assert "4.1.235" in r


# ---------------------------------------------------------------------------
# Three-layer resolution tests
# ---------------------------------------------------------------------------


class TestThreeLayerResolution:
    def test_bundled_fallback(self):
        """Bundled Parquet files should be available for psi_ms and pride_cv."""
        bundled = (
            Path(__file__).parent.parent.parent / "qpx" / "core" / "ontology" / "data"
        )
        assert (bundled / "psi_ms.parquet").is_file()
        assert (bundled / "pride_cv.parquet").is_file()

    def test_repo_checkout_found(self):
        """Repo checkout data/ should be found when running from the repo."""
        repo_data = Path(__file__).parent.parent.parent / "data" / "ontology"
        assert (repo_data / "psi_ms.parquet").is_file()

    def test_auto_update_false_skips_network(self):
        """auto_update=False should not make any network calls."""
        onto = PublicOntology("psi_ms", auto_update=False)
        assert onto.version is not None
        assert len(onto) > 0
        onto.close()


# ---------------------------------------------------------------------------
# Integration with scores.py
# ---------------------------------------------------------------------------


class TestScoresIntegration:
    def test_lookup_score(self):
        from qpx.core.scores import lookup_score

        info = lookup_score("percolator_score")
        assert info is not None
        assert info["ontology_accession"] == "MS:1001492"

    def test_is_higher_better(self):
        from qpx.core.scores import is_higher_better

        assert is_higher_better("percolator_score") is True

    def test_builtin_score(self):
        from qpx.core.scores import lookup_score

        info = lookup_score("diann_qvalue")
        assert info is not None
        assert info["higher_better"] is False

    def test_score_ontology_entries(self):
        from qpx.core.scores import score_ontology_entries

        entries = score_ontology_entries(
            {"percolator_score", "diann_qvalue"}, view="psm"
        )
        assert len(entries) == 2
        names = {e["field_name"] for e in entries}
        assert "percolator_score" in names
        assert "diann_qvalue" in names
        # Check ontology_version is present
        for e in entries:
            assert "ontology_version" in e
