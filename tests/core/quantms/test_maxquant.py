"""
Optimized MaxQuant test module focused on core functionality testing.
"""

import os
import tempfile
from pathlib import Path

import pandas as pd
import pytest
import pyarrow.parquet as pq

from quantmsio.core.maxquant.maxquant import (
    MaxQuant,
    parse_modifications_from_peptidoform,
    read_evidence,
    read_msms,
    read_protein_groups,
)

# Test data paths
TEST_DATA_ROOT = Path(__file__).parent.parent.parent / "examples" / "maxquant"
TEST_FILES = {
    "msms": TEST_DATA_ROOT / "maxquant_simple/msms.txt",
    "evidence": TEST_DATA_ROOT / "maxquant_simple/evidence.txt",
    "sdrf": TEST_DATA_ROOT / "maxquant_simple/sdrf.tsv",
    "protein_groups": TEST_DATA_ROOT / "maxquant_full/proteinGroups.txt",
}


# Test data validation
@pytest.fixture(scope="module")
def validate_test_data():
    """Validate that all required test files exist."""
    missing_files = []
    for name, path in TEST_FILES.items():
        if not path.exists():
            missing_files.append(f"{name}: {path}")

    if missing_files:
        pytest.skip(f"Missing test files: {', '.join(missing_files)}")

    return TEST_FILES


def test_core_functionality(validate_test_data):
    """Test core MaxQuant functionality."""
    test_files = validate_test_data
    processor = MaxQuant()

    # Test basic file reading
    msms_df = processor.read_msms(str(test_files["msms"]))
    evidence_df = processor.read_evidence(str(test_files["evidence"]))

    assert msms_df is not None and len(msms_df) > 0, "MSMS file should contain data"
    assert (
        evidence_df is not None and len(evidence_df) > 0
    ), "Evidence file should contain data"

    # Test processing functionality
    psm_table = processor.process_msms_to_psm_table(msms_df)
    feature_table = processor.process_evidence_to_feature_table(evidence_df)

    assert psm_table is not None, "PSM table processing should succeed"
    assert feature_table is not None, "Feature table processing should succeed"


def test_backward_compatibility(validate_test_data):
    """Test backward compatibility functions."""
    test_files = validate_test_data

    # Test backward compatibility functions
    msms_df_compat = read_msms(str(test_files["msms"]))
    evidence_df_compat = read_evidence(str(test_files["evidence"]))
    pg_df_compat = read_protein_groups(str(test_files["protein_groups"]))

    assert (
        msms_df_compat is not None and len(msms_df_compat) > 0
    ), "Backward compatible MSMS reading should work"
    assert (
        evidence_df_compat is not None and len(evidence_df_compat) > 0
    ), "Backward compatible evidence reading should work"
    assert (
        pg_df_compat is not None and len(pg_df_compat) > 0
    ), "Backward compatible protein groups reading should work"


def test_batch_processing(validate_test_data):
    """Test batch processing functionality."""
    test_files = validate_test_data
    processor = MaxQuant()

    # Test batch processing functionality
    chunk_count = 0
    for df_chunk in processor.iter_batch(str(test_files["msms"]), chunksize=100):
        assert (
            df_chunk is not None and len(df_chunk) > 0
        ), f"Batch chunk {chunk_count} should contain data"
        chunk_count += 1
        if chunk_count >= 3:  # Limit test to first few batches for performance
            break

    assert chunk_count > 0, "Should process at least one batch"


def test_modification_parsing():
    """Test modification parsing functionality."""
    # Test modification parsing
    test_cases = [
        ("TESTPEPTIDE", None, "Unmodified peptide should return None"),
        ("TEST(Oxidation)PEPTIDE", list, "Modified peptide should return list"),
        ("(Acetyl)PEPTIDE", list, "N-terminal modification should be parsed"),
        ("PEPTIDE(Amidated)", list, "C-terminal modification should be parsed"),
    ]

    for peptidoform, expected_type, description in test_cases:
        modifications = parse_modifications_from_peptidoform(peptidoform)

        if expected_type is None:
            assert modifications is None, description
        else:
            assert isinstance(modifications, expected_type), description
            for mod in modifications:
                assert "name" in mod, "Modification should have 'name' field"
                assert "positions" in mod, "Modification should have 'positions' field"
                assert isinstance(mod["positions"], list), "Positions should be a list"

    # Test edge cases for modification parsing
    edge_cases = [("", None), (None, None)]
    for peptidoform, expected in edge_cases:
        result = parse_modifications_from_peptidoform(peptidoform)
        assert result == expected, f"Edge case failed for {peptidoform}"


def test_protein_groups_processing():
    """Test protein groups processing functionality."""
    # Test protein groups processing
    sample_data = {
        "Protein IDs": ["P12345;Q67890", "P11111", "P22222"],
        "Fasta headers": ["Protein A;Protein B", "Protein C", "Protein D"],
        "Gene names": ["GENEA;GENEB", "GENEC", "GENED"],
        "Q-value": [0.001, 0.01, 0.05],
        "Intensity": [1000.0, 2000.0, 3000.0],
        "LFQ intensity": [800.0, 1800.0, 2800.0],
        "Razor + unique peptides": [5, 3, 7],
        "Unique peptides": [3, 2, 5],
        "Sequence coverage [%]": [15.5, 25.0, 30.0],
        "Mol. weight [kDa]": [50.5, 75.2, 100.1],
        "Score": [85.2, 92.1, 78.3],
        "Reverse": ["", "+", ""],
        "Potential contaminant": ["", "", "+"],
        "MS/MS count": [10, 15, 12],
    }

    df = pd.DataFrame(sample_data)
    processor = MaxQuant()
    table = processor.process_protein_groups_to_pg_table(df)
    result_df = table.to_pandas()

    # Test basic conversion
    required_columns = [
        "pg_accessions",
        "pg_names",
        "gg_accessions",
        "intensities",
        "additional_intensities",
        "anchor_protein",
        "additional_scores",
    ]
    for col in required_columns:
        assert col in result_df.columns, f"Required column '{col}' missing"

    # Test protein accession splitting and flags
    pg_accessions = result_df.iloc[0]["pg_accessions"]
    if hasattr(pg_accessions, "to_pylist"):
        pg_accessions = pg_accessions.to_pylist()
    elif not isinstance(pg_accessions, list):
        pg_accessions = list(pg_accessions) if pg_accessions is not None else []

    assert (
        isinstance(pg_accessions, list) and len(pg_accessions) == 2
    ), "Protein accessions should be correctly split"
    assert (
        result_df.iloc[0]["anchor_protein"] == "P12345"
    ), "Anchor protein should be first"
    assert result_df.iloc[1]["is_decoy"] == 1, "Reverse protein should be decoy"
    assert result_df.iloc[2]["contaminant"] == 1, "Contaminant should be flagged"


def test_error_handling(validate_test_data):
    """Test error handling and edge cases."""
    processor = MaxQuant()

    # Test invalid file paths
    invalid_files = ["nonexistent_file.txt", "/invalid/path/file.txt", ""]

    for invalid_file in invalid_files:
        with pytest.raises(Exception, match=".*"):
            processor.read_msms(invalid_file)

    # Test empty DataFrame handling
    empty_df = pd.DataFrame()
    empty_data_tests = [
        (processor.process_msms_to_psm_table, "PSM processing"),
        (processor.process_evidence_to_feature_table, "Feature processing"),
        (processor.process_protein_groups_to_pg_table, "Protein groups processing"),
    ]

    for process_func, description in empty_data_tests:
        try:
            result = process_func(empty_df)
            assert (
                result is not None
            ), f"{description} should return valid result or raise exception"
        except Exception as e:
            assert isinstance(
                e, (ValueError, KeyError, IndexError)
            ), f"{description} should raise appropriate exception type"


def test_data_validation(validate_test_data):
    """Test data validation and integrity checks."""
    processor = MaxQuant()
    test_files = validate_test_data

    # Test file size validation
    for name, filepath in test_files.items():
        size = os.path.getsize(filepath)
        assert size > 100, f"Test file too small: {name} at {filepath} ({size} bytes)"

    # Test PSM data validation
    msms_df = processor.read_msms(str(test_files["msms"]))
    psm_table = processor.process_msms_to_psm_table(msms_df.head(50).copy())
    psm_df = psm_table.to_pandas()

    # Check required PSM fields
    required_psm_fields = [
        "sequence",
        "peptidoform",
        "precursor_charge",
        "scan",
        "rt",
        "calculated_mz",
        "observed_mz",
        "is_decoy",
        "protein_accessions",
    ]
    for field in required_psm_fields:
        assert field in psm_df.columns, f"Missing PSM field: {field}"

    # Validate required fields have no null values
    required_non_null_psm = ["sequence", "peptidoform", "precursor_charge", "scan"]
    for field in required_non_null_psm:
        assert (
            not psm_df[field].isna().any()
        ), f"PSM field '{field}' should not contain null values"

    # Test Feature data validation
    evidence_df = processor.read_evidence(str(test_files["evidence"]))
    feature_table = processor.process_evidence_to_feature_table(
        evidence_df.head(50).copy()
    )
    feature_df = feature_table.to_pandas()

    # Validate charge values and m/z values are positive
    assert (
        feature_df["precursor_charge"] > 0
    ).all(), "Precursor charges should be positive"
    assert (feature_df["calculated_mz"] > 0).all(), "Calculated m/z should be positive"
    assert (feature_df["observed_mz"] > 0).all(), "Observed m/z should be positive"


def test_file_writing(validate_test_data):
    """Test file writing functionality and comprehensive workflow."""
    test_files = validate_test_data
    processor = MaxQuant()

    with tempfile.TemporaryDirectory() as temp_dir:
        try:
            # Define output file paths
            output_files = {
                "msms": os.path.join(temp_dir, "msms_output.parquet"),
                "evidence": os.path.join(temp_dir, "evidence_output.parquet"),
                "pg": os.path.join(temp_dir, "pg_output.parquet"),
            }

            # Test direct processing and writing
            msms_df = processor.read_msms(str(test_files["msms"]))
            evidence_df = processor.read_evidence(str(test_files["evidence"]))
            pg_df = processor.read_protein_groups(str(test_files["protein_groups"]))

            psm_table = processor.process_msms_to_psm_table(msms_df)
            feature_table = processor.process_evidence_to_feature_table(evidence_df)
            pg_table = processor.process_protein_groups_to_pg_table(pg_df)

            # Write files
            pq.write_table(psm_table, output_files["msms"])
            pq.write_table(feature_table, output_files["evidence"])
            pq.write_table(pg_table, output_files["pg"])

            # Validate direct writing outputs
            for file_type, filepath in output_files.items():
                assert os.path.exists(
                    filepath
                ), f"Direct output file {file_type} not created"
                table = pq.read_table(filepath)
                assert len(table) > 0, f"Direct output file {file_type} is empty"
                size = os.path.getsize(filepath)
                assert (
                    size > 100
                ), f"Direct output file {file_type} too small ({size} bytes)"

            # Test comprehensive workflow using high-level methods
            workflow_files = {
                "psm": os.path.join(temp_dir, "workflow_msms.parquet"),
                "feature": os.path.join(temp_dir, "workflow_evidence.parquet"),
                "pg": os.path.join(temp_dir, "workflow_pg.parquet"),
            }

            processor.write_psm_to_file(str(test_files["msms"]), workflow_files["psm"])
            processor.write_feature_to_file(
                str(test_files["evidence"]),
                str(test_files["sdrf"]),
                workflow_files["feature"],
            )
            processor.write_protein_groups_to_file(
                str(test_files["protein_groups"]),
                str(test_files["sdrf"]),
                workflow_files["pg"],
            )

            # Validate workflow outputs
            for file_type, filepath in workflow_files.items():
                assert os.path.exists(
                    filepath
                ), f"Workflow output file {file_type} not created"
                table = pq.read_table(filepath)
                assert len(table) > 0, f"Workflow output file {file_type} is empty"
                df = table.to_pandas()
                assert (
                    len(df.columns) > 5
                ), f"Workflow output file {file_type} should have multiple columns"

        except Exception as e:
            pytest.fail(f"File writing and workflow test failed: {str(e)}")


if __name__ == "__main__":
    # Direct execution for development/debugging
    import sys

    # Mock test data validation for direct execution
    class MockTestData:
        def __init__(self):
            self.data = TEST_FILES

        def __getitem__(self, key):
            return self.data[key]

        def items(self):
            return self.data.items()

    mock_data = MockTestData()

    try:
        print("Running MaxQuant tests (8 test nodes)...")

        # Run the 8 test functions
        test_core_functionality(mock_data)
        print("Core functionality tests passed!")

        test_backward_compatibility(mock_data)
        print("Backward compatibility tests passed!")

        test_batch_processing(mock_data)
        print("Batch processing tests passed!")

        test_modification_parsing()
        print("Modification parsing tests passed!")

        test_protein_groups_processing()
        print("Protein groups processing tests passed!")

        test_error_handling(mock_data)
        print("Error handling tests passed!")

        test_data_validation(mock_data)
        print("Data validation tests passed!")

        test_file_writing(mock_data)
        print("File writing tests passed!")

        print("All 8 MaxQuant test nodes passed!")

    except Exception as e:
        print(f"Test failed: {str(e)}")
        sys.exit(1)
