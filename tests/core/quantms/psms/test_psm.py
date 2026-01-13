import logging
from pathlib import Path
import tempfile
import os
import time
import pandas as pd
import pyarrow as pa

from qpx.core.quantms.psm import Psm
from qpx.core.quantms.mztab import MzTabIndexer
from qpx.utils.mztab_utils import (
    fetch_ms_runs_from_mztab_line,
    extract_ms_runs_from_metadata,
    extract_ms_runs_from_lines,
    _clean_file_path,
    _is_ms_run_location_line,
    _extract_ms_run_id_from_key,
)
from qpx.utils.logger import get_logger

logger = get_logger(__name__)

TEST_DATA_ROOT = Path(__file__).parents[3] / "examples"

# Global test data path
MZTAB_TEST_PATH = (
    TEST_DATA_ROOT
    / "quantms/dda-lfq-full/PXD007683-LFQ.sdrf_openms_design_openms.mzTab.gz"
)


def cleanup_test_resources(indexer=None, temp_db_path=None, temp_output_path=None):
    """
    Clean up test resources including indexer connection and temporary files.

    Args:
        indexer: MzTabIndexer instance to close
        temp_db_path: Path to temporary database file to delete
        temp_output_path: Path to temporary output file to delete
    """
    # Close indexer connection first
    if indexer:
        try:
            indexer.close()
        except Exception as e:
            logger.warning(f"Failed to close indexer: {e}")

    # Clean up temporary files with retry on Windows
    for temp_path in [temp_db_path, temp_output_path]:
        if temp_path and os.path.exists(temp_path):
            try:
                os.unlink(temp_path)
            except PermissionError:
                # On Windows, sometimes we need a brief delay
                time.sleep(0.1)
                try:
                    os.unlink(temp_path)
                except Exception as e:
                    logger.warning(f"Failed to delete temporary file {temp_path}: {e}")


MZTAB_TEST_TMT_PATH = (
    TEST_DATA_ROOT
    / "quantms/dda-plex-full/PXD007683TMT.sdrf_openms_design_openms.mzTab.gz"
)

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


def test_generate_report():
    """Test PSM report generation functionality."""
    # Skip if test data doesn't exist
    if not MZTAB_TEST_PATH.exists():
        print(f"Test data not found: {MZTAB_TEST_PATH}")
        return

    # Create a temporary database path (don't create the file, let DuckDB do it)
    temp_db_file = tempfile.NamedTemporaryFile(suffix=".duckdb", delete=False)
    temp_db_path = temp_db_file.name
    temp_db_file.close()
    if os.path.exists(temp_db_path):
        os.unlink(temp_db_path)

    indexer = None
    try:
        # Create the MzTabIndexer first
        indexer = MzTabIndexer.create(
            mztab_path=MZTAB_TEST_PATH, database_path=temp_db_path
        )

        # Pass the indexer to PSM
        psm = Psm(indexer)
        logger.info(f"Number of psms: {psm._indexer._get_num_psms()}")

    finally:
        cleanup_test_resources(indexer=indexer, temp_db_path=temp_db_path)


def test_iter_psm_table():
    """Test PSM table iteration functionality."""
    # Skip if test data doesn't exist
    if not MZTAB_TEST_PATH.exists():
        print(f"Test data not found: {MZTAB_TEST_PATH}")
        return

    # Create a temporary database path (don't create the file, let DuckDB do it)
    temp_db_file = tempfile.NamedTemporaryFile(suffix=".duckdb", delete=False)
    temp_db_path = temp_db_file.name
    temp_db_file.close()
    if os.path.exists(temp_db_path):
        os.unlink(temp_db_path)

    indexer = None
    try:
        # Create the MzTabIndexer first
        indexer = MzTabIndexer.create(
            mztab_path=MZTAB_TEST_PATH, database_path=temp_db_path
        )

        # Pass the indexer to PSM
        psm = Psm(indexer)

        # Test iterating over PSM table
        psm_count = 0
        for df in psm.iter_psm_table(chunksize=1000):
            assert isinstance(df, pa.Table)
            psm_count += len(df)
        assert psm_count > 0
        print(f"Processed {psm_count} PSMs")

    finally:
        cleanup_test_resources(indexer=indexer, temp_db_path=temp_db_path)


def test_write_to_file():
    """Test PSM file writing functionality."""
    # Skip if test data doesn't exist
    if not MZTAB_TEST_PATH.exists():
        print(f"Test data not found: {MZTAB_TEST_PATH}")
        return

    # Create temporary file paths (don't create the files, let the libraries do it)
    temp_db_file = tempfile.NamedTemporaryFile(suffix=".duckdb", delete=False)
    temp_db_path = temp_db_file.name
    temp_db_file.close()
    if os.path.exists(temp_db_path):
        os.unlink(temp_db_path)
    temp_output_file = tempfile.NamedTemporaryFile(suffix=".parquet", delete=False)
    temp_output_path = temp_output_file.name
    temp_output_file.close()

    indexer = None
    try:
        # Create the MzTabIndexer first
        indexer = MzTabIndexer.create(
            mztab_path=MZTAB_TEST_PATH, database_path=temp_db_path
        )

        # Pass the indexer to PSM
        psm = Psm(indexer)

        # Test writing PSM data to file
        psm.convert_to_parquet(temp_output_path, chunksize=1000)

        # Verify the file was created and has content
        assert os.path.exists(temp_output_path)
        assert os.path.getsize(temp_output_path) > 0
        print(f"Successfully wrote PSM data to {temp_output_path}")

    finally:
        cleanup_test_resources(
            indexer=indexer,
            temp_db_path=temp_db_path,
            temp_output_path=temp_output_path,
        )


def test_extract_ms_runs():
    """Test MS run extraction with both DuckDB and parquet backends."""

    # Skip if test data doesn't exist
    if not MZTAB_TEST_PATH.exists():
        logger.error(f"Test data not found: {MZTAB_TEST_PATH}")
        return

    logger.info(f"Test data found: {MZTAB_TEST_PATH}")

    # Create a temporary database path for DuckDB
    temp_db_file = tempfile.NamedTemporaryFile(suffix=".duckdb", delete=False)
    temp_db_path = temp_db_file.name
    temp_db_file.close()
    if os.path.exists(temp_db_path):
        os.unlink(temp_db_path)
    # else:
    #     # Create a temporary directory for parquet backend
    #     temp_db_path = tempfile.mktemp(suffix="_parquet")
    #     os.makedirs(temp_db_path, exist_ok=True)

    indexer = None
    try:
        # Create the MzTabIndexer with real data
        backend = "duckdb"
        logger.info(
            f"Creating MzTabIndexer with {backend} backend, path: {temp_db_path}"
        )
        indexer = MzTabIndexer.create(
            mztab_path=MZTAB_TEST_PATH,
            database_path=temp_db_path,
        )
        logger.info(f"MzTabIndexer created with {backend} backend")

        # Create PSM instance
        psm = Psm(indexer)
        logger.info(f"PSM created with {backend} backend")

        # Test the _extract_ms_runs function
        ms_runs = psm._extract_ms_runs()

        # Verify basic structure
        assert isinstance(
            ms_runs, dict
        ), f"MS runs should be a dictionary for {backend} backend"
        assert (
            len(ms_runs) > 0
        ), f"Should have at least one MS run for {backend} backend"

        # Print the ms_runs to console for inspection
        print(f"\nExtracted {len(ms_runs)} MS runs with {backend} backend:")
        for run_id, file_name in sorted(ms_runs.items()):
            logger.debug(f"  ms_run[{run_id}]: {file_name}")

        # Verify all keys and values are properly formatted
        for run_id, file_name in ms_runs.items():
            assert isinstance(
                run_id, int
            ), f"Run ID should be integer for {backend} backend, got {type(run_id)}"
            assert isinstance(
                file_name, str
            ), f"File name should be string for {backend} backend, got {type(file_name)}"
            assert run_id > 0, f"Run ID should be positive for {backend} backend"
            assert (
                len(file_name) > 0
            ), f"File name should not be empty for {backend} backend"

            # Verify file_name format (should be a clean filename without extension)
            assert (
                "." not in file_name
            ), f"File name should not have extension for {backend} backend, got '{file_name}'"

        # Print results for inspection
        logger.debug(f"\nExtracted {len(ms_runs)} MS runs with {backend} backend:")
        for run_id, file_name in sorted(ms_runs.items()):
            logger.debug(f"  ms_run[{run_id}]: {file_name}")

        # Additional assertions based on expected data structure
        # The test data should have multiple MS runs
        assert (
            len(ms_runs) >= 1
        ), f"Expected at least 1 MS run for {backend} backend, got {len(ms_runs)}"

        # All file names should be unique
        file_names = list(ms_runs.values())
        unique_file_names = set(file_names)
        assert len(file_names) == len(
            unique_file_names
        ), f"All file names should be unique for {backend} backend"

        # Run IDs should be sequential starting from 1
        run_ids = list(ms_runs.keys())
        run_ids.sort()
        assert (
            run_ids[0] == 1
        ), f"First run ID should be 1 for {backend} backend, got {run_ids[0]}"

        # Check for gaps in run IDs (optional - depends on data)
        for i in range(len(run_ids) - 1):
            assert (
                run_ids[i + 1] == run_ids[i] + 1
            ), f"Run IDs should be sequential for {backend} backend, found gap between {run_ids[i]} and {run_ids[i+1]}"

        logger.info(f"{backend} backend test passed successfully")

    finally:
        cleanup_test_resources(indexer=indexer, temp_db_path=temp_db_path)


def test_ms_run_utils():
    """Test the MS run utility functions."""

    # Test _clean_file_path
    assert _clean_file_path("file:///path/to/file.mzML") == "file"
    assert _clean_file_path("/data/sample_001.mzML") == "sample_001"
    assert _clean_file_path("relative/path/file.txt") == "file"
    assert _clean_file_path("simple.mzML") == "simple"

    # Test _is_ms_run_location_line
    assert (
        _is_ms_run_location_line(
            ["MTD", "ms_run[1]-location", "file:///path/file.mzML"]
        )
        == True
    )
    assert _is_ms_run_location_line(["MTD", "ms_run[1]-format", "mzML"]) == False
    assert _is_ms_run_location_line(["PRT", "accession", "P12345"]) == False

    # Test _extract_ms_run_id_from_key
    assert _extract_ms_run_id_from_key("ms_run[1]-location") == 1
    assert _extract_ms_run_id_from_key("ms_run[10]-location") == 10

    # Test fetch_ms_runs_from_mztab_line
    ms_runs = {}
    line = "MTD\tms_run[1]-location\tfile:///path/to/file1.mzML"
    ms_runs = fetch_ms_runs_from_mztab_line(line, ms_runs)
    assert ms_runs == {1: "file1"}

    line2 = "MTD\tms_run[2]-location\tfile:///data/sample_001.mzML"
    ms_runs = fetch_ms_runs_from_mztab_line(line2, ms_runs)
    assert ms_runs == {1: "file1", 2: "sample_001"}

    # Test extract_ms_runs_from_metadata
    metadata_df = pd.DataFrame(
        {
            "key": ["ms_run[1]-location", "ms_run[2]-location", "other-key"],
            "value": [
                "file:///path/to/file1.mzML",
                "file:///path/to/file2.mzML",
                "other-value",
            ],
        }
    )
    result = extract_ms_runs_from_metadata(metadata_df)
    assert result == {1: "file1", 2: "file2"}

    # Test extract_ms_runs_from_lines
    lines = [
        "MTD\tms_run[1]-location\tfile:///path/to/file1.mzML",
        "MTD\tms_run[2]-location\tfile:///path/to/file2.mzML",
        "MTD\tother-key\tother-value",
    ]
    result = extract_ms_runs_from_lines(lines)
    assert result == {1: "file1", 2: "file2"}


def test_ms_run_extraction():
    """Test MS run extraction utility functions."""

    # Skip if test data doesn't exist
    if not MZTAB_TEST_PATH.exists():
        logger.error(f"Test data not found: {MZTAB_TEST_PATH}")
        return

    logger.info(f"Test data found: {MZTAB_TEST_PATH}")

    # Create a temporary database path for DuckDB
    temp_db_file = tempfile.NamedTemporaryFile(suffix=".duckdb", delete=False)
    temp_db_path = temp_db_file.name
    temp_db_file.close()
    if os.path.exists(temp_db_path):
        os.unlink(temp_db_path)

    indexer = None
    try:
        # Create the MzTabIndexer with real data
        backend = "duckdb"
        logger.info(
            f"Creating MzTabIndexer with {backend} backend, path: {temp_db_path}"
        )
        indexer = MzTabIndexer.create(
            mztab_path=MZTAB_TEST_PATH,
            database_path=temp_db_path,
        )
        logger.info(f"MzTabIndexer created with {backend} backend")

        # Test the utility function directly
        metadata_df = indexer.get_metadata()
        ms_runs = extract_ms_runs_from_metadata(metadata_df)

        # Verify basic structure
        assert isinstance(
            ms_runs, dict
        ), f"MS runs should be a dictionary for {backend} backend"
        assert (
            len(ms_runs) > 0
        ), f"Should have at least one MS run for {backend} backend"

        # Print the ms_runs to console for inspection
        print(
            f"\nExtracted {len(ms_runs)} MS runs with {backend} backend using utility function:"
        )
        for run_id, file_name in sorted(ms_runs.items()):
            logger.debug(f"  ms_run[{run_id}]: {file_name}")

        # Verify all keys and values are properly formatted
        for run_id, file_name in ms_runs.items():
            assert isinstance(
                run_id, int
            ), f"Run ID should be integer for {backend} backend, got {type(run_id)}"
            assert isinstance(
                file_name, str
            ), f"File name should be string for {backend} backend, got {type(file_name)}"
            assert run_id > 0, f"Run ID should be positive for {backend} backend"
            assert (
                len(file_name) > 0
            ), f"File name should not be empty for {backend} backend"

            # Verify file_name format (should be a clean filename without extension)
            assert (
                "." not in file_name
            ), f"File name should not have extension for {backend} backend, got '{file_name}'"

        # Additional assertions based on expected data structure
        # The test data should have multiple MS runs
        assert (
            len(ms_runs) >= 1
        ), f"Expected at least 1 MS run for {backend} backend, got {len(ms_runs)}"

        # All file names should be unique
        file_names = list(ms_runs.values())
        unique_file_names = set(file_names)
        assert len(file_names) == len(
            unique_file_names
        ), f"All file names should be unique for {backend} backend"

        # Run IDs should be sequential starting from 1
        run_ids = list(ms_runs.keys())
        run_ids.sort()
        assert (
            run_ids[0] == 1
        ), f"First run ID should be 1 for {backend} backend, got {run_ids[0]}"

        logger.info(f"{backend} backend utility function test passed successfully")

    finally:
        cleanup_test_resources(indexer=indexer, temp_db_path=temp_db_path)


def test_get_protein_qvalue():
    """Test protein q-value extraction functionality."""

    # Skip if test data doesn't exist
    if not MZTAB_TEST_PATH.exists():
        logger.error(f"Test data not found: {MZTAB_TEST_PATH}")
        return

    logger.info(f"Test data found: {MZTAB_TEST_PATH}")

    # Create a temporary database path for DuckDB
    temp_db_file = tempfile.NamedTemporaryFile(suffix=".duckdb", delete=False)
    temp_db_path = temp_db_file.name
    temp_db_file.close()
    if os.path.exists(temp_db_path):
        os.unlink(temp_db_path)

    indexer = None
    try:
        # Create the MzTabIndexer with real data
        backend = "duckdb"
        logger.info(
            f"Creating MzTabIndexer with {backend} backend, path: {temp_db_path}"
        )
        indexer = MzTabIndexer.create(
            mztab_path=MZTAB_TEST_PATH,
            database_path=temp_db_path,
        )
        logger.info(f"MzTabIndexer created with {backend} backend")

        # Create PSM instance
        psm = Psm(indexer)
        logger.info(f"PSM created with {backend} backend")

        # Test the _get_protein_qvalue function
        protein_qvalue_map = psm._get_protein_qvalue()

        # Verify basic structure
        assert isinstance(
            protein_qvalue_map, dict
        ), f"Protein q-value map should be a dictionary for {backend} backend"

        # Count the number of proteins
        protein_count = len(protein_qvalue_map)
        logger.info(
            f"Found {protein_count} proteins with q-values in {backend} backend"
        )

        # Print the first 10 proteins with lowest q-values
        if protein_count > 0:
            # Sort proteins by q-value (lowest first)
            sorted_proteins = sorted(
                protein_qvalue_map.items(),
                key=lambda x: float(x[1]) if x[1] is not None else float("inf"),
            )

            logger.info(f"First 10 proteins with lowest q-values ({backend} backend):")
            logger.info(f"{'Rank':<4} {'Protein Accession':<20} {'Q-Value':<10}")
            logger.info("-" * 40)

            for i, (protein_acc, qvalue) in enumerate(sorted_proteins[:10], 1):
                qvalue_str = str(qvalue) if qvalue is not None else "None"
                logger.info(f"{i:<4} {protein_acc:<20} {qvalue_str:<10}")

            # Additional analysis
            qvalues = [float(q) for q in protein_qvalue_map.values() if q is not None]
            if qvalues:
                min_qvalue = min(qvalues)
                max_qvalue = max(qvalues)
                avg_qvalue = sum(qvalues) / len(qvalues)

                logger.info(f"Q-value statistics ({backend} backend):")
                logger.info(f"  Minimum q-value: {min_qvalue:.6f}")
                logger.info(f"  Maximum q-value: {max_qvalue:.6f}")
                logger.info(f"  Average q-value: {avg_qvalue:.6f}")
                logger.info(f"  Proteins with q-values: {len(qvalues)}/{protein_count}")

            # Verify that q-values are reasonable (between 0 and 1)
            for protein_acc, qvalue in protein_qvalue_map.items():
                if qvalue is not None:
                    qvalue_float = float(qvalue)
                    assert (
                        0 <= qvalue_float <= 1
                    ), f"Q-value should be between 0 and 1 for {backend} backend, got {qvalue_float} for protein {protein_acc}"

            # Check that we have some proteins with low q-values (good quality)
            low_qvalue_count = sum(
                1
                for q in protein_qvalue_map.values()
                if q is not None and float(q) < 0.01
            )
            logger.info(
                f"Proteins with q-value < 0.01: {low_qvalue_count}/{protein_count}"
            )

            # Verify that protein accessions are not empty
            for protein_acc in protein_qvalue_map.keys():
                assert isinstance(
                    protein_acc, str
                ), f"Protein accession should be string for {backend} backend"
                assert (
                    len(protein_acc) > 0
                ), f"Protein accession should not be empty for {backend} backend"

            # Analyze protein groups with multiple accessions
            logger.info(f"\n{'='*60}")
            logger.info(
                f"ANALYZING PROTEIN GROUPS WITH MULTIPLE ACCESSIONS ({backend} backend)"
            )
            logger.info(f"{'='*60}")

            # Filter proteins that have semicolons (indicating protein groups with multiple accessions)
            protein_groups = {
                protein: qvalue
                for protein, qvalue in protein_qvalue_map.items()
                if ";" in protein
            }
            logger.info(
                f"Found {len(protein_groups)} protein groups with multiple accessions"
            )

            if protein_groups:
                # Sort by q-value (lowest first) and show top 10
                sorted_groups = sorted(
                    protein_groups.items(),
                    key=lambda x: float(x[1]) if x[1] is not None else float("inf"),
                )

                logger.info(f"\nTop 10 protein groups with lowest q-values:")
                logger.info(f"{'Accession':<50} {'Protein Count':<15} {'Q-value':<10}")
                logger.info(f"{'-'*50} {'-'*15} {'-'*10}")

                for i, (accession, qvalue) in enumerate(sorted_groups[:10], 1):
                    protein_count = len(accession.split(";"))
                    qvalue_str = (
                        f"{float(qvalue):.6f}" if qvalue is not None else "None"
                    )
                    logger.info(f"{accession:<50} {protein_count:<15} {qvalue_str:<10}")
            else:
                logger.info("No protein groups with multiple accessions found")

                # Debug: Show some examples of what accessions look like
                logger.info(
                    f"\nDebug: Sample of protein accessions from {backend} backend:"
                )
                sample_accessions = list(protein_qvalue_map.keys())[:5]
                for acc in sample_accessions:
                    logger.info(f"  '{acc}' (length: {len(acc)})")
        else:
            logger.warning(f"No proteins with q-values found for {backend} backend")

        logger.info(f"{backend} backend protein q-value test passed successfully")

    finally:
        cleanup_test_resources(indexer=indexer, temp_db_path=temp_db_path)


def test_get_metadata_modifications():
    """Test the _get_metadata_modifications method using a real mzTab file."""

    # Skip if test data doesn't exist
    if not MZTAB_TEST_TMT_PATH.exists():
        logger.error(f"Test data not found: {MZTAB_TEST_TMT_PATH}")
        return

    logger.info(f"Testing modification parsing from real file: {MZTAB_TEST_TMT_PATH}")

    # Use a temporary directory for the backend database
    temp_db_file = tempfile.NamedTemporaryFile(suffix=".duckdb", delete=False)
    temp_db_path = temp_db_file.name
    temp_db_file.close()
    if os.path.exists(temp_db_path):
        os.unlink(temp_db_path)

    indexer = None
    try:
        # 1. Create the MzTabIndexer with the real mzTab file
        indexer = MzTabIndexer.create(
            mztab_path=MZTAB_TEST_TMT_PATH,
            database_path=temp_db_path,
        )

        # 2. Instantiate the Psm processor
        psm_processor = Psm(indexer)

        # 3. Access the modifications parsed during __init__
        modifications = psm_processor._modifications

        # Print all modifications
        logger.info(f"Number of modifications: {len(modifications)}")
        for key, (name, sites, positions) in modifications.items():
            logger.info(f"Acession: {key}")
            logger.info(f"Name: {name}")
            logger.info(f"Sites {sites}")
            logger.info(f"Postions: {positions}")

        # 4. Assert the results based on expected content of the test file
        assert isinstance(modifications, dict), "The result should be a dictionary."
        assert len(modifications) > 0, "Should parse at least one modification."

        logger.info(
            f"Successfully parsed {len(modifications)} modifications from the file."
        )
        logger.info("test_get_metadata_modifications_from_file passed successfully.")

    finally:
        cleanup_test_resources(indexer=indexer, temp_db_path=temp_db_path)
