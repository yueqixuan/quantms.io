"""Tests for DIA-NN converter constants."""


def test_diann_constants():
    """DIA-NN constants: structure, field mappings are lists, key fields present."""
    from qpx.converters.diann.constants import TOOL_NAME, TOOL_VERSIONS, FIELD_MAPPINGS

    assert TOOL_NAME == "DIA-NN"
    assert isinstance(TOOL_VERSIONS, str)
    assert "feature" in FIELD_MAPPINGS
    assert "pg" in FIELD_MAPPINGS

    # Field mappings are lists with at least one candidate
    for view, fields in FIELD_MAPPINGS.items():
        for qpx_field, candidates in fields.items():
            assert isinstance(candidates, list), (
                f"FIELD_MAPPINGS['{view}']['{qpx_field}'] must be a list"
            )
            assert len(candidates) > 0

    # Key fields present
    feature = FIELD_MAPPINGS["feature"]
    assert "intensity" in feature
    assert "posterior_error_probability" in feature
    assert "rt" in feature
    assert "pg_accessions" in feature
