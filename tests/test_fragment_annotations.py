"""
Test for the new fragment annotation arrays in PSM schema.
"""

import pytest

try:
    import pyarrow as pa
    from qpx.core.format import PSM_UNIQUE_FIELDS
    from qpx.core.common import PSM_SCHEMA

    PYARROW_AVAILABLE = True
except ImportError:
    PYARROW_AVAILABLE = False


@pytest.mark.skipif(not PYARROW_AVAILABLE, reason="PyArrow not available")
def test_psm_schema_contains_fragment_annotation_fields():
    """Test that PSM schema contains the new fragment annotation fields."""

    # Check that the new fields are present
    field_names = [field.name for field in PSM_UNIQUE_FIELDS]
    required_new_fields = ["charge_array", "ion_type_array", "ion_mobility_array"]

    for field_name in required_new_fields:
        assert field_name in field_names, f"{field_name} not found in PSM_UNIQUE_FIELDS"

    # Check that the PSM schema can be created successfully
    schema = PSM_SCHEMA
    assert schema is not None
    assert len(schema) > 0

    # Verify the new fields are in the full schema
    full_field_names = [field.name for field in schema]
    for field_name in required_new_fields:
        assert field_name in full_field_names, f"{field_name} not found in PSM_SCHEMA"


@pytest.mark.skipif(not PYARROW_AVAILABLE, reason="PyArrow not available")
def test_fragment_annotation_fields_types():
    """Test that the new fragment annotation fields have correct types."""

    # Get the new fields from PSM_UNIQUE_FIELDS
    field_dict = {field.name: field for field in PSM_UNIQUE_FIELDS}

    # Test charge_array field
    charge_field = field_dict["charge_array"]
    assert charge_field.type == pa.list_(
        pa.int32()
    ), "charge_array should be list of int32"
    assert charge_field.nullable, "charge_array should be nullable"

    # Test ion_type_array field
    ion_type_field = field_dict["ion_type_array"]
    assert ion_type_field.type == pa.list_(
        pa.string()
    ), "ion_type_array should be list of string"
    assert ion_type_field.nullable, "ion_type_array should be nullable"

    # Test ion_mobility_array field
    ion_mobility_field = field_dict["ion_mobility_array"]
    assert ion_mobility_field.type == pa.list_(
        pa.float32()
    ), "ion_mobility_array should be list of float32"
    assert ion_mobility_field.nullable, "ion_mobility_array should be nullable"


@pytest.mark.skipif(not PYARROW_AVAILABLE, reason="PyArrow not available")
def test_fragment_annotation_fields_have_descriptions():
    """Test that the new fragment annotation fields have proper descriptions."""

    field_dict = {field.name: field for field in PSM_UNIQUE_FIELDS}

    # Check that all new fields have descriptions
    for field_name in ["charge_array", "ion_type_array", "ion_mobility_array"]:
        field = field_dict[field_name]
        assert field.metadata is not None, f"{field_name} should have metadata"
        assert (
            b"description" in field.metadata
        ), f"{field_name} should have description in metadata"
        description = field.metadata[b"description"].decode()
        assert len(description) > 0, f"{field_name} should have non-empty description"
        assert (
            "fragment" in description.lower()
        ), f"{field_name} description should mention 'fragment'"


@pytest.mark.skipif(not PYARROW_AVAILABLE, reason="PyArrow not available")
def test_create_table_with_fragment_annotations():
    """Test creating a PyArrow table with the new fragment annotation fields."""

    # Create sample data with all required fields including the new ones
    sample_data = {
        # Basic required fields
        "sequence": ["PEPTIDE"],
        "peptidoform": ["PEPTIDE"],
        "precursor_charge": [2],
        "posterior_error_probability": [0.01],
        "is_decoy": [0],
        "calculated_mz": [450.2],
        "observed_mz": [450.21],
        "mp_accessions": [["P12345"]],
        "predicted_rt": [None],
        "reference_file_name": ["sample1"],
        "scan": ["1000"],
        "rt": [120.5],
        "ion_mobility": [None],
        # Array fields
        "number_peaks": [3],
        "mz_array": [[100.1, 200.2, 300.3]],
        "intensity_array": [[1000.0, 2000.0, 3000.0]],
        # New fragment annotation fields
        "charge_array": [[1, 1, 2]],
        "ion_type_array": [["b", "y", "b"]],
        "ion_mobility_array": [[0.8, 0.9, 1.0]],
    }

    # Create a minimal schema for testing
    test_fields = [
        pa.field("sequence", pa.string()),
        pa.field("peptidoform", pa.string()),
        pa.field("precursor_charge", pa.int32()),
        pa.field("posterior_error_probability", pa.float32()),
        pa.field("is_decoy", pa.int32()),
        pa.field("calculated_mz", pa.float32()),
        pa.field("observed_mz", pa.float32()),
        pa.field("mp_accessions", pa.list_(pa.string())),
        pa.field("predicted_rt", pa.float32()),
        pa.field("reference_file_name", pa.string()),
        pa.field("scan", pa.string()),
        pa.field("rt", pa.float32()),
        pa.field("ion_mobility", pa.float32()),
        pa.field("number_peaks", pa.int32()),
        pa.field("mz_array", pa.list_(pa.float32())),
        pa.field("intensity_array", pa.list_(pa.float32())),
        pa.field("charge_array", pa.list_(pa.int32())),
        pa.field("ion_type_array", pa.list_(pa.string())),
        pa.field("ion_mobility_array", pa.list_(pa.float32())),
    ]

    test_schema = pa.schema(test_fields)

    # Create table - this should not raise an exception
    table = pa.table(sample_data, schema=test_schema)

    # Verify table structure
    assert table.num_rows == 1
    assert table.num_columns == len(test_fields)

    # Verify the new fields contain expected data
    assert table.column("charge_array").to_pylist() == [[1, 1, 2]]
    assert table.column("ion_type_array").to_pylist() == [["b", "y", "b"]]

    # For float32 arrays, use approximate comparison due to precision
    ion_mobility_values = table.column("ion_mobility_array").to_pylist()[0]
    expected_values = [0.8, 0.9, 1.0]
    assert len(ion_mobility_values) == len(expected_values)
    for actual, expected in zip(ion_mobility_values, expected_values):
        assert abs(actual - expected) < 1e-6, f"Expected {expected}, got {actual}"


@pytest.mark.skipif(not PYARROW_AVAILABLE, reason="PyArrow not available")
def test_nullable_fragment_annotation_fields():
    """Test that fragment annotation fields can contain null values."""

    # Create sample data where fragment annotation fields are null
    sample_data = {
        "sequence": ["PEPTIDE"],
        "peptidoform": ["PEPTIDE"],
        "precursor_charge": [2],
        "mz_array": [[100.1, 200.2, 300.3]],
        "intensity_array": [[1000.0, 2000.0, 3000.0]],
        "number_peaks": [3],
        # New fragment annotation fields as null
        "charge_array": [None],
        "ion_type_array": [None],
        "ion_mobility_array": [None],
    }

    test_fields = [
        pa.field("sequence", pa.string()),
        pa.field("peptidoform", pa.string()),
        pa.field("precursor_charge", pa.int32()),
        pa.field("mz_array", pa.list_(pa.float32())),
        pa.field("intensity_array", pa.list_(pa.float32())),
        pa.field("number_peaks", pa.int32()),
        pa.field("charge_array", pa.list_(pa.int32())),
        pa.field("ion_type_array", pa.list_(pa.string())),
        pa.field("ion_mobility_array", pa.list_(pa.float32())),
    ]

    test_schema = pa.schema(test_fields)

    # Create table with null values - this should not raise an exception
    table = pa.table(sample_data, schema=test_schema)

    # Verify null values are handled correctly
    assert table.column("charge_array").to_pylist() == [None]
    assert table.column("ion_type_array").to_pylist() == [None]
    assert table.column("ion_mobility_array").to_pylist() == [None]
