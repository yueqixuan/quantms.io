from pathlib import Path
import tempfile

from quantmsio.core.quantms.feature import Feature
from quantmsio.core.quantms.mztab import MzTabIndexer

TEST_DATA_ROOT = Path(__file__).parent.parent.parent / "examples"

test_data = (
    TEST_DATA_ROOT
    / "quantms/dda-lfq-full/PXD007683-LFQ.sdrf_openms_design_openms.mzTab.gz",
    TEST_DATA_ROOT
    / "quantms/dda-lfq-full/PXD007683-LFQ.sdrf_openms_design_msstats_in.csv.gz",
    TEST_DATA_ROOT / "quantms/dda-lfq-full/PXD007683-LFQ.sdrf.tsv",
)


def test_transform_msstats():
    """Test transforming msstats data."""
    # Resolve file paths
    mztab_file = test_data[0]
    msstats_file = test_data[1]
    sdrf_file = test_data[2]

    # Use a temporary directory for the backend database
    temp_db_path = tempfile.mktemp(suffix=".duckdb")

    # Initialize Feature
    indexer = MzTabIndexer.create(
        mztab_path=mztab_file,
        msstats_path=msstats_file,
        sdrf_path=sdrf_file,
        database_path=temp_db_path,
    )

    # Use new composition pattern
    feature = Feature(
        mztab_indexer=indexer,
    )

    try:
        # Process msstats data
        count = 0
        for msstats in feature.transform_msstats_in():
            # Add assertions to verify the result
            assert msstats is not None
            assert len(msstats) > 0
            assert "peptidoform" in msstats.columns
            count += 1

        # Ensure we got at least one result
        assert count > 0
    finally:
        # Clean up Feature resources
        feature.cleanup()


def test_extract_psm_msg():
    """Test extracting PSM messages."""
    # Resolve file paths
    mztab_file = test_data[0]
    msstats_file = test_data[1]
    sdrf_file = test_data[2]

    # Use a temporary directory for the backend database
    temp_db_path = tempfile.mktemp(suffix=".duckdb")

    # Initialize Feature
    indexer = MzTabIndexer.create(
        mztab_path=mztab_file,
        msstats_path=msstats_file,
        sdrf_path=sdrf_file,
        database_path=temp_db_path,
    )

    # Use new composition pattern
    feature = Feature(
        mztab_indexer=indexer,
    )

    try:
        # Extract PSM messages
        map_dict, pep_dict = feature.extract_psm_msg()

        # Add assertions to verify the result
        assert map_dict is not None
        assert isinstance(map_dict, dict)
        assert len(map_dict) > 0

        assert pep_dict is not None
        assert isinstance(pep_dict, dict)
        assert len(pep_dict) > 0
    finally:
        # Clean up Feature resources
        feature.cleanup()


def test_generate_feature():
    """Test generating features."""
    # Resolve file paths
    mztab_file = test_data[0]
    msstats_file = test_data[1]
    sdrf_file = test_data[2]

    # Use a temporary directory for the backend database
    temp_db_path = tempfile.mktemp(suffix=".duckdb")

    # Initialize Feature
    indexer = MzTabIndexer.create(
        mztab_path=mztab_file,
        msstats_path=msstats_file,
        sdrf_path=sdrf_file,
        database_path=temp_db_path,
    )

    # Use new composition pattern
    feature = Feature(
        mztab_indexer=indexer,
    )

    try:
        # Generate features
        count = 0
        for feature_table in feature.generate_feature():
            # Add assertions to verify the result
            assert feature_table is not None
            assert len(feature_table) > 0
            assert "peptidoform" in feature_table.column_names
            count += 1

        # Ensure we got at least one result
        assert count > 0
    finally:
        # Clean up Feature resources
        feature.cleanup()
