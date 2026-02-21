"""Tests for DIA-NN converter constants."""


def test_diann_constants_structure():
    from qpx.converters.diann.constants import TOOL_NAME, TOOL_VERSIONS, FIELD_MAPPINGS

    assert TOOL_NAME == "DIA-NN"
    assert isinstance(TOOL_VERSIONS, str)
    assert "feature" in FIELD_MAPPINGS
    assert "pg" in FIELD_MAPPINGS


def test_diann_field_mappings_are_lists():
    from qpx.converters.diann.constants import FIELD_MAPPINGS

    for view, fields in FIELD_MAPPINGS.items():
        for qpx_field, candidates in fields.items():
            assert isinstance(
                candidates, list
            ), f"FIELD_MAPPINGS['{view}']['{qpx_field}'] must be a list"
            assert len(candidates) > 0


def test_diann_key_fields_present():
    from qpx.converters.diann.constants import FIELD_MAPPINGS

    feature = FIELD_MAPPINGS["feature"]
    assert "intensity" in feature
    assert "posterior_error_probability" in feature
    assert "rt" in feature
    assert "pg_accessions" in feature
