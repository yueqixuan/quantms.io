import logging
import tempfile
import os
from pathlib import Path
import pytest

from qpx.core.idxml_utils.idxml import IdXmlPsm

# Test data path
TEST_DATA_ROOT = Path(__file__).parents[2] / "examples"
IDXML_TEST_PATH = (
    TEST_DATA_ROOT
    / "idxml/SF_200217_pPeptideLibrary_pool1_HCDnlETcaD_OT_rep2_consensus_fdr_pep_luciphor.idXML"
)

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


def test_idxml_psm_initialization():
    """Test IdXmlPsm initialization."""
    # Skip if test data doesn't exist
    if not IDXML_TEST_PATH.exists():
        pytest.skip(f"Test data not found: {IDXML_TEST_PATH}")

    # Test initialization
    idxml_psm = IdXmlPsm(str(IDXML_TEST_PATH))

    # Check that PSM data can be loaded
    assert idxml_psm.get_psm_count() > 0
    assert idxml_psm.search_metadata is not None
    logger.info(f"Successfully loaded {idxml_psm.get_psm_count()} PSMs")


def test_iter_psm_table():
    """Test PSM table iteration functionality."""
    # Skip if test data doesn't exist
    if not IDXML_TEST_PATH.exists():
        pytest.skip(f"Test data not found: {IDXML_TEST_PATH}")

    idxml_psm = IdXmlPsm(str(IDXML_TEST_PATH))

    # Test iterating over PSM table
    psm_count = 0
    for table in idxml_psm.iter_psm_table(chunksize=1000):
        assert table is not None
        psm_count += len(table)

    assert psm_count > 0
    assert psm_count == idxml_psm.get_psm_count()
    logger.info(f"Processed {psm_count} PSMs in chunks")


def test_convert_to_parquet():
    """Test converting idXML PSMs to parquet format."""
    # Skip if test data doesn't exist
    if not IDXML_TEST_PATH.exists():
        pytest.skip(f"Test data not found: {IDXML_TEST_PATH}")

    # Create temporary output file
    temp_output_file = tempfile.NamedTemporaryFile(suffix=".parquet", delete=False)
    temp_output_path = temp_output_file.name
    temp_output_file.close()

    try:
        logger.info("=== DEBUG: Starting test_convert_to_parquet ===")
        logger.info(f"Input file: {IDXML_TEST_PATH}")
        logger.info(f"Output file: {temp_output_path}")

        # Initialize IdXmlPsm with verbose logging
        logger.info("=== DEBUG: Initializing IdXmlPsm ===")
        idxml_psm = IdXmlPsm(str(IDXML_TEST_PATH))
        logger.info(f"=== DEBUG: Loaded {idxml_psm.get_psm_count()} PSMs ===")

        # Debug: Check first few PSM records
        if idxml_psm.psm_records:
            logger.info("=== DEBUG: First PSM record sample ===")
            first_psm = idxml_psm.psm_records[0]
            logger.info(f"Sequence: {first_psm.get('sequence', 'N/A')}")
            logger.info(f"Charge: {first_psm.get('charge', 'N/A')}")
            logger.info(f"RT: {first_psm.get('rt', 'N/A')}")
            logger.info(f"MZ: {first_psm.get('mz', 'N/A')}")
            logger.info(
                f"Protein accessions: {first_psm.get('protein_accessions', [])}"
            )
            logger.info(f"Modifications: {first_psm.get('modifications', [])}")
            logger.info(f"Peptidoform: {first_psm.get('peptidoform', 'N/A')}")
            logger.info(f"Additional scores: {first_psm.get('additional_scores', [])}")

        # Test converting to parquet with verbose logging
        logger.info("=== DEBUG: Starting parquet conversion ===")
        idxml_psm.convert_to_parquet(temp_output_path, chunksize=1000)
        logger.info("=== DEBUG: Parquet conversion completed ===")

        # Verify the file was created and has content
        logger.info("=== DEBUG: Verifying output file ===")
        assert os.path.exists(
            temp_output_path
        ), f"Output file not created: {temp_output_path}"
        file_size = os.path.getsize(temp_output_path)
        assert file_size > 0, f"Output file is empty: {temp_output_path}"
        logger.info(f"=== DEBUG: Output file size: {file_size} bytes ===")
        logger.info(f"Successfully wrote PSM data to {temp_output_path}")

        # Try to read the parquet file to verify it's valid
        logger.info("=== DEBUG: Reading parquet file for verification ===")
        import pyarrow.parquet as pq

        table = pq.read_table(temp_output_path)
        assert len(table) > 0, f"Parquet file contains no data: {len(table)} rows"
        logger.info(f"=== DEBUG: Parquet file contains {len(table)} PSMs ===")
        logger.info(f"=== DEBUG: Parquet schema: {table.schema} ===")
        logger.info(f"Verified parquet file contains {len(table)} PSMs")

        # Debug: Check first few rows of parquet data
        if len(table) > 0:
            logger.info("=== DEBUG: First parquet row sample ===")
            first_row = table.to_pylist()[0]
            for key, value in first_row.items():
                logger.info(f"  {key}: {value}")

        logger.info("=== DEBUG: test_convert_to_parquet completed successfully ===")

    except Exception as e:
        logger.error(f"=== DEBUG: Error in test_convert_to_parquet: {e} ===")
        logger.error(f"=== DEBUG: Exception type: {type(e)} ===")
        import traceback

        logger.error(f"=== DEBUG: Traceback: {traceback.format_exc()} ===")
        raise
    finally:
        # Clean up temporary file
        if os.path.exists(temp_output_path):
            os.unlink(temp_output_path)
            logger.info(f"=== DEBUG: Cleaned up temporary file: {temp_output_path} ===")


def test_psm_record_transformation():
    """Test PSM record transformation functionality."""
    # Skip if test data doesn't exist
    if not IDXML_TEST_PATH.exists():
        pytest.skip(f"Test data not found: {IDXML_TEST_PATH}")

    idxml_psm = IdXmlPsm(str(IDXML_TEST_PATH))

    # Get first batch of PSMs
    for table in idxml_psm.iter_psm_table(chunksize=1):
        if len(table) > 0:
            transformed = table.to_pylist()[0]

            assert transformed is not None

            # Check required fields are present
            required_fields = [
                "sequence",
                "peptidoform",
                "precursor_charge",
                "protein_accessions",
                "scan",
                "reference_file_name",
                "is_decoy",
            ]

            for field in required_fields:
                assert field in transformed, f"Required field '{field}' missing"

            # Check data types
            assert isinstance(transformed["sequence"], str)
            assert isinstance(transformed["peptidoform"], str)
            assert isinstance(transformed["precursor_charge"], int)
            assert isinstance(transformed["protein_accessions"], list)
            assert isinstance(transformed["is_decoy"], int)

            logger.info("PSM record transformation test passed")
            break


def test_psm_schema_compliance():
    """Test that generated PSM data complies with PSM schema."""
    # Skip if test data doesn't exist
    if not IDXML_TEST_PATH.exists():
        pytest.skip(f"Test data not found: {IDXML_TEST_PATH}")

    idxml_psm = IdXmlPsm(str(IDXML_TEST_PATH))

    # Get a sample table
    for table in idxml_psm.iter_psm_table(chunksize=10):
        # Check that the table has the expected schema
        from qpx.core.common import PSM_SCHEMA

        # Verify schema compatibility
        assert table.schema.equals(PSM_SCHEMA, check_metadata=False)
        logger.info("PSM schema compliance test passed")
        break


def test_error_handling():
    """Test error handling for invalid inputs."""
    # Test with non-existent file
    with pytest.raises(Exception):
        IdXmlPsm("non_existent_file.idXML")

    logger.info("Error handling test passed")


def test_metadata_extraction():
    """Test metadata and modification extraction."""
    # Skip if test data doesn't exist
    if not IDXML_TEST_PATH.exists():
        pytest.skip(f"Test data not found: {IDXML_TEST_PATH}")

    idxml_psm = IdXmlPsm(str(IDXML_TEST_PATH))

    # Check that metadata is properly extracted
    assert idxml_psm.search_metadata is not None
    assert idxml_psm.psm_records is not None

    # Check first few PSMs
    for table in idxml_psm.iter_psm_table(chunksize=5):
        for psm in table.to_pylist():  # Check first 5 PSMs
            # Check that basic fields are present
            assert "sequence" in psm
            assert "precursor_charge" in psm
            assert "protein_accessions" in psm

            # Check modifications field
            assert "modifications" in psm
            assert isinstance(psm["modifications"], list)

            # Check peptidoform field
            assert "peptidoform" in psm
            assert isinstance(psm["peptidoform"], str)
        break  # Only check first batch

    logger.info("Metadata extraction test passed")


def test_openms_integration():
    """Test integration with OpenMS functionality."""
    from qpx.core.openms import OpenMSHandler

    # Skip if test data doesn't exist
    if not IDXML_TEST_PATH.exists():
        pytest.skip(f"Test data not found: {IDXML_TEST_PATH}")

    # Test that IdXmlPsm uses OpenMS functionality correctly
    idxml_psm = IdXmlPsm(str(IDXML_TEST_PATH))

    # Check that OpenMS handler is initialized
    assert idxml_psm.openms_handler is not None
    assert isinstance(idxml_psm.openms_handler, OpenMSHandler)

    # Check that PSM data was processed using OpenMS
    assert len(idxml_psm.psm_records) > 0

    # Check that each PSM has expected fields processed by OpenMS
    for psm in idxml_psm.psm_records[:3]:  # Check first 3 PSMs
        assert "sequence" in psm
        assert "charge" in psm
        assert "protein_accessions" in psm
        assert "modifications" in psm
        assert "peptidoform" in psm

    # Check that search metadata was extracted using OpenMS
    assert idxml_psm.search_metadata is not None
    assert isinstance(idxml_psm.search_metadata, dict)

    logger.info("OpenMS integration test passed")


if __name__ == "__main__":
    # Run tests manually if needed
    test_idxml_psm_initialization()
    test_iter_psm_table()
    test_convert_to_parquet()
    test_psm_record_transformation()
    test_psm_schema_compliance()
    test_error_handling()
    test_metadata_extraction()
    test_openms_integration()
    print("All tests passed!")
