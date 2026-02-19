"""Tests for field_ontology_entries() with source provenance."""
import pytest


def test_field_ontology_entries_with_source_provenance():
    from qpx.core.scores import field_ontology_entries
    resolved = {"intensity": "Precursor.Quantity", "rt": "RT"}
    entries = field_ontology_entries(
        view="feature",
        resolved_mappings=resolved,
        tool_name="DIA-NN",
    )
    # Should have entries for resolved fields that have CV mappings
    rt_entries = [e for e in entries if e["field_name"] == "rt"]
    assert len(rt_entries) == 1
    assert rt_entries[0]["source_column_name"] == "RT"
    assert rt_entries[0]["source_tool"] == "DIA-NN"
    assert rt_entries[0]["ontology_accession"] == "MS:1000016"
    assert rt_entries[0]["view"] == "feature"


def test_field_ontology_entries_missing_cv_still_written():
    from qpx.core.scores import field_ontology_entries
    resolved = {"lfq": "Precursor.Normalised"}
    entries = field_ontology_entries(
        view="feature",
        resolved_mappings=resolved,
        tool_name="DIA-NN",
    )
    # lfq has no CV term but should still produce an entry
    assert len(entries) >= 1
    lfq_entries = [e for e in entries if e["field_name"] == "lfq"]
    assert len(lfq_entries) == 1
    assert lfq_entries[0]["source_column_name"] == "Precursor.Normalised"
    assert lfq_entries[0]["source_tool"] == "DIA-NN"
    # ontology_accession should be None (no CV term)
    assert lfq_entries[0]["ontology_accession"] is None


def test_field_ontology_entries_backward_compat():
    """Old call signature (no resolved_mappings) still works."""
    from qpx.core.scores import field_ontology_entries
    entries = field_ontology_entries(view="psm")
    # Should still produce entries for _FIELD_CV_MAP fields
    field_names = {e["field_name"] for e in entries}
    assert "posterior_error_probability" in field_names
    assert "rt" in field_names
    # Old entries should have source_column_name = None
    for e in entries:
        assert "source_column_name" in e
        assert "source_tool" in e
