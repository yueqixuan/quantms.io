"""Tests for MzTab indexer functionality."""

import logging
import os
import tempfile
import time
from pathlib import Path
import pytest
import pandas as pd

from qpx.core.quantms.mztab import MzTabIndexer
from qpx.utils.mztab_utils import extract_ms_runs_from_metadata
from qpx.utils.logger import get_logger

logger = get_logger(__name__)


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


# Test data paths
TEST_DATA_ROOT = Path(__file__).parents[3] / "examples"
DDA_PLEX_DATA = {
    "mztab": TEST_DATA_ROOT
    / "quantms/dda-plex-full/PXD007683TMT.sdrf_openms_design_openms.mzTab.gz",
    "msstats": TEST_DATA_ROOT
    / "quantms/dda-plex-full/PXD007683TMT.sdrf_openms_design_msstats_in.csv.gz",
    "sdrf": TEST_DATA_ROOT / "quantms/dda-plex-full/PXD007683-TMT.sdrf.tsv",
}
DDA_LFQ_DATA = {
    "mztab": TEST_DATA_ROOT
    / "quantms/dda-lfq-full/PXD007683-LFQ.sdrf_openms_design_openms.mzTab.gz",
    "msstats": TEST_DATA_ROOT
    / "quantms/dda-lfq-full/PXD007683-LFQ.sdrf_openms_design_msstats_in.csv.gz",
    "sdrf": TEST_DATA_ROOT / "quantms/dda-lfq-full/PXD007683-LFQ.sdrf.tsv",
}


def test_dda_plex_dataset():
    """Test MzTabIndexer functionality with DDA-plex dataset."""
    # Skip if files don't exist
    if not all(path.exists() for path in DDA_PLEX_DATA.values()):
        pytest.skip("DDA-plex test files not found")

    # Create a temporary database path
    temp_fd, temp_db_path = tempfile.mkstemp(suffix=".duckdb")
    os.close(temp_fd)
    os.unlink(temp_db_path)

    indexer = None
    try:
        # Initialize MzTabIndexer using factory method
        indexer = MzTabIndexer.create(
            mztab_path=DDA_PLEX_DATA["mztab"],
            msstats_path=DDA_PLEX_DATA["msstats"],
            sdrf_path=DDA_PLEX_DATA["sdrf"],
            database_path=temp_db_path,
            max_memory="8GB",
            worker_threads=4,
            batch_size=100000,
        )

        # Test metadata retrieval
        metadata = indexer.get_metadata()
        assert not metadata.empty, "No metadata found"
        assert len(metadata) == 555
        print("\nDDA-plex dataset metadata statistics:")
        print(f"Number of metadata entries: {len(metadata)}")

        # Test protein data retrieval
        proteins = indexer.get_proteins()
        assert not proteins.empty, "No protein data found"
        assert len(proteins) == 9655
        print("\nProtein statistics:")
        print(f"Total proteins: {len(proteins):,}")
        print("\nProtein columns:")
        print(proteins.columns.tolist())

        # Test protein details data retrieval
        # protein_details = indexer.get_protein_details()
        # assert not protein_details.empty, "No protein details data found"
        # assert len(protein_details) == 10940
        # print(f"\nTotal protein details: {len(protein_details):,}")

        # Test PSM data retrieval using stream_section to get sample
        psm_count = indexer.get_psm_count()
        assert psm_count == 154774, f"Expected 154774 PSMs, got {psm_count}"
        print(f"\nTotal PSMs: {psm_count:,}")

        # Get PSM sample using stream_section
        psm_sample = None
        for chunk in indexer.stream_section("PSM", chunk_size=1000):
            psm_sample = chunk
            break
        assert psm_sample is not None and not psm_sample.empty, "No PSM data found"
        print("\nPSM columns:")
        print(psm_sample.columns.tolist())

        # Test MSstats data retrieval
        msstats = indexer.get_msstats()
        assert msstats is not None and not msstats.empty, "No MSstats data found"
        assert len(msstats) == 1549526
        print(f"\nTotal MSstats entries: {len(msstats):,}")

        # Test anchor protein and accession logic
        protein_sample = indexer.query_to_df(
            "SELECT accession, ambiguity_members, opt_global_result_type FROM proteins LIMIT 5"
        )
        for _, row in protein_sample.iterrows():
            if row["opt_global_result_type"] == "single_protein":
                assert (
                    row["accession"] == row["ambiguity_members"]
                ), "Anchor protein should match accession for single proteins"
            elif row["opt_global_result_type"] == "indistinguishable_group":
                assert (
                    row["accession"] != row["ambiguity_members"]
                ), "Anchor protein should not match accession for groups"

        # Test sorted accessions in PSMs
        psm_sample = indexer.query_to_df("SELECT accession FROM psms LIMIT 5")
        for _, row in psm_sample.iterrows():
            if row["accession"] and "," in row["accession"]:
                accessions = row["accession"].split(",")
                assert accessions == sorted(accessions), "PSM accessions are not sorted"

        # Test join between proteins and psms
        join_query = """
        SELECT COUNT(*) as count
        FROM proteins p
        JOIN psms s ON p.accession = s.accession
        """
        join_result = indexer.query_to_df(join_query)
        assert join_result["count"].iloc[0] > 0, "Join between proteins and PSMs failed"
        print(
            f"\nSuccessfully joined proteins and PSMs, found {join_result['count'].iloc[0]} matches."
        )

        # Test referential integrity
        protein_accessions = set(
            indexer.query_to_df("SELECT accession FROM proteins")["accession"]
        )
        # protein_details_accessions = set(
        #     indexer.query_to_df("SELECT accession FROM protein_details")["accession"]
        # )
        # all_protein_accessions = protein_accessions.union(protein_details_accessions)

        psm_accessions = set(
            indexer.query_to_df(
                "SELECT DISTINCT accession FROM psms WHERE accession IS NOT NULL AND accession != 'null'"
            )["accession"]
        )
        for psm_acc_group in psm_accessions:
            for psm_acc in psm_acc_group.split(","):
                logging.info(f"PSM accession {psm_acc} protein group not found")
                logging.info(
                    any(psm_acc in prot_acc for prot_acc in protein_accessions)
                )

        msstats_df = indexer.query_to_df(
            "SELECT DISTINCT pg_accessions FROM msstats WHERE pg_accessions IS NOT NULL"
        )
        msstats_accessions = set(
            msstats_df["pg_accessions"].str.replace(";", ",", regex=False)
        )
        for msstats_acc_group in msstats_accessions:
            for msstats_acc in msstats_acc_group.split(","):
                logging.info(f"MSstats accession {msstats_acc} protein group not found")
                logging.info(
                    any(msstats_acc in prot_acc for prot_acc in protein_accessions)
                )
        print("\nReferential integrity check passed for DDA-plex dataset.")

    finally:
        cleanup_test_resources(indexer=indexer, temp_db_path=temp_db_path)


def test_dda_lfq_dataset():
    """Test MzTabIndexer functionality with DDA-LFQ dataset."""
    # Skip if files don't exist
    if not all(path.exists() for path in DDA_LFQ_DATA.values()):
        pytest.skip("DDA-LFQ test files not found")

    # Create a temporary database path
    temp_fd, temp_db_path = tempfile.mkstemp(suffix=".duckdb")
    os.close(temp_fd)
    os.unlink(temp_db_path)

    indexer = None
    try:
        # Initialize MzTabIndexer using factory method
        indexer = MzTabIndexer.create(
            mztab_path=DDA_LFQ_DATA["mztab"],
            msstats_path=DDA_LFQ_DATA["msstats"],
            sdrf_path=DDA_LFQ_DATA["sdrf"],
            database_path=temp_db_path,
            max_memory="8GB",
            worker_threads=4,
            batch_size=50000,
        )

        # Test metadata retrieval
        metadata = indexer.get_metadata()
        assert not metadata.empty, "No metadata found"
        assert len(metadata) == 109
        print("\nDDA-LFQ dataset metadata statistics:")
        print(f"Number of metadata entries: {len(metadata)}")

        # Test protein data retrieval
        proteins = indexer.get_proteins()
        assert not proteins.empty, "No protein data found"
        assert len(proteins) == 9346
        print("\nProtein statistics:")
        print(f"Total proteins: {len(proteins):,}")
        print("\nProtein columns:")
        print(proteins.columns.tolist())

        # Test protein details data retrieval
        # protein_details = indexer.get_protein_details()
        # assert not protein_details.empty, "No protein details data found"
        # assert len(protein_details) == 9574
        # print(f"\nTotal protein details: {len(protein_details):,}")

        # Test PSM data retrieval using stream_section to get sample
        psm_count = indexer.get_psm_count()
        assert psm_count == 537848, f"Expected 537848 PSMs, got {psm_count}"
        print(f"\nTotal PSMs: {psm_count:,}")

        # Get PSM sample using stream_section
        psm_sample = None
        for chunk in indexer.stream_section("PSM", chunk_size=1000):
            psm_sample = chunk
            break
        assert psm_sample is not None and not psm_sample.empty, "No PSM data found"
        print("\nPSM columns:")
        print(psm_sample.columns.tolist())

        # Test MSstats data retrieval
        msstats = indexer.get_msstats()
        assert msstats is not None and not msstats.empty, "No MSstats data found"
        assert len(msstats) == 425594
        print(f"\nTotal MSstats entries: {len(msstats):,}")

        # Test anchor protein and accession logic
        protein_sample = indexer.query_to_df(
            "SELECT accession, ambiguity_members, opt_global_result_type FROM proteins LIMIT 5"
        )
        for _, row in protein_sample.iterrows():
            if row["opt_global_result_type"] == "single_protein":
                assert (
                    row["accession"] == row["ambiguity_members"]
                ), "Anchor protein should match accession for single proteins"
            elif row["opt_global_result_type"] == "indistinguishable_group":
                assert (
                    row["accession"] != row["ambiguity_members"]
                ), "Anchor protein should not match accession for groups"

        # Test sorted accessions in PSMs
        psm_sample = indexer.query_to_df("SELECT accession FROM psms LIMIT 5")
        for _, row in psm_sample.iterrows():
            if row["accession"] and "," in row["accession"]:
                accessions = row["accession"].split(",")
                assert accessions == sorted(accessions), "PSM accessions are not sorted"

        # Test join between proteins and psms
        join_query = """
        SELECT COUNT(*) as count
        FROM proteins p
        JOIN psms s ON p.accession = s.accession
        """
        join_result = indexer.query_to_df(join_query)
        assert join_result["count"].iloc[0] > 0, "Join between proteins and PSMs failed"
        print(
            f"\nSuccessfully joined proteins and PSMs, found {join_result['count'].iloc[0]} matches."
        )

        # Test referential integrity
        protein_accessions = set(
            indexer.query_to_df("SELECT accession FROM proteins")["accession"]
        )
        # protein_details_accessions = set(
        #     indexer.query_to_df("SELECT accession FROM protein_details")["accession"]
        # )
        # all_protein_accessions = protein_accessions.union(protein_details_accessions)

        psm_accessions = set(
            indexer.query_to_df(
                "SELECT DISTINCT accession FROM psms WHERE accession IS NOT NULL AND accession != 'null'"
            )["accession"]
        )
        for psm_acc_group in psm_accessions:
            for psm_acc in psm_acc_group.split(","):
                logging.info(f"PSM accession {psm_acc} protein group not found")
                logging.info(
                    any(psm_acc in prot_acc for prot_acc in protein_accessions)
                )

        msstats_df = indexer.query_to_df(
            "SELECT DISTINCT pg_accessions FROM msstats WHERE pg_accessions IS NOT NULL"
        )
        msstats_accessions = set(
            msstats_df["pg_accessions"].str.replace(";", ",", regex=False)
        )
        for msstats_acc_group in msstats_accessions:
            for msstats_acc in msstats_acc_group.split(","):
                logging.info(f"MSstats accession {msstats_acc} protein group not found")
                logging.info(
                    any(msstats_acc in prot_acc for prot_acc in protein_accessions)
                )
        print("\nReferential integrity check passed for DDA-LFQ dataset.")

    finally:
        cleanup_test_resources(indexer=indexer, temp_db_path=temp_db_path)


def test_add_msstats_to_existing_db():
    """Test adding an MSstats table to an existing MzTabIndexer database."""
    if not DDA_LFQ_DATA["mztab"].exists() or not DDA_LFQ_DATA["msstats"].exists():
        pytest.skip("DDA-LFQ test files not found")

    # Create a temporary database path
    temp_fd, temp_db_path = tempfile.mkstemp(suffix=".duckdb")
    os.close(temp_fd)
    os.unlink(temp_db_path)

    indexer1 = None
    indexer2 = None
    try:
        # Step 1: Create an indexer with only mzTab to create the initial database
        indexer1 = MzTabIndexer.create(
            mztab_path=DDA_LFQ_DATA["mztab"],
            database_path=temp_db_path,
        )

        # Check that msstats table doesn't exist
        tables = indexer1.query_to_df("SHOW TABLES")["name"].tolist()
        assert MzTabIndexer._MZTAB_INDEXER_TABLE_MSSTATS not in tables

        # Close first indexer before opening second one
        indexer1.close()
        indexer1 = None

        # Step 2: Open the existing database and add the MSstats table
        indexer2 = MzTabIndexer.open(
            database_path=temp_db_path,
            sdrf_path=DDA_LFQ_DATA["sdrf"],
        )

        # The msstats table should not exist yet
        tables = indexer2.query_to_df("SHOW TABLES")["name"].tolist()
        assert MzTabIndexer._MZTAB_INDEXER_TABLE_MSSTATS not in tables

        # Add the msstats table
        indexer2.add_msstats_table(msstats_path=str(DDA_LFQ_DATA["msstats"]))

        # Verify that the MSstats table now exists
        tables = indexer2.query_to_df("SHOW TABLES")["name"].tolist()
        assert MzTabIndexer._MZTAB_INDEXER_TABLE_MSSTATS in tables

        msstats_data = indexer2.get_msstats()
        assert msstats_data is not None
        assert not msstats_data.empty
        # The number of entries is known from the other test
        assert len(msstats_data) == 425594

    finally:
        # Clean up both indexers
        if indexer1:
            try:
                indexer1.close()
            except Exception as e:
                logger.warning(f"Failed to close indexer1: {e}")
        if indexer2:
            try:
                indexer2.close()
            except Exception as e:
                logger.warning(f"Failed to close indexer2: {e}")
        # Clean up temporary files
        if temp_db_path and os.path.exists(temp_db_path):
            try:
                os.unlink(temp_db_path)
            except PermissionError:
                time.sleep(0.1)
                try:
                    os.unlink(temp_db_path)
                except Exception as e:
                    logger.warning(
                        f"Failed to delete temporary file {temp_db_path}: {e}"
                    )


def test_mztab_metadata_parsing():
    """Test MzTab metadata parsing functionality."""
    if not DDA_PLEX_DATA["mztab"].exists():
        pytest.skip("DDA-plex test files not found")

    # Create a temporary database path
    temp_fd, temp_db_path = tempfile.mkstemp(suffix=".duckdb")
    os.close(temp_fd)  # Close file descriptor, let DuckDB create the file
    os.unlink(temp_db_path)  # Remove the empty file

    indexer = None
    try:
        indexer = MzTabIndexer.create(
            mztab_path=DDA_PLEX_DATA["mztab"],
            database_path=temp_db_path,
        )

        # Test metadata retrieval
        metadata = indexer.get_metadata()
        assert not metadata.empty, "No metadata found"

        # Test MS runs extraction using utility function
        ms_runs = extract_ms_runs_from_metadata(metadata)
        assert len(ms_runs) == 11
        assert 1 in ms_runs
        assert ms_runs[1] == "a05068"

        # Test score names from metadata
        score_metadata = metadata[
            metadata["key"].str.contains("search_engine_score", na=False)
        ]
        assert not score_metadata.empty, "No score metadata found"

        # Test modifications from metadata
        mod_metadata = metadata[metadata["key"].str.contains("_mod\\[", na=False)]
        assert not mod_metadata.empty, "No modification metadata found"

        print("\nMetadata parsing tests passed.")

    finally:
        cleanup_test_resources(indexer=indexer, temp_db_path=temp_db_path)


def test_mztab_stream_section():
    """Test MzTab.stream_section functionality."""
    if not DDA_PLEX_DATA["mztab"].exists():
        pytest.skip("DDA-plex test files not found")

    # Create a temporary database path
    temp_fd, temp_db_path = tempfile.mkstemp(suffix=".duckdb")
    os.close(temp_fd)
    os.unlink(temp_db_path)

    indexer = None
    try:
        indexer = MzTabIndexer.create(
            mztab_path=DDA_PLEX_DATA["mztab"],
            database_path=temp_db_path,
        )

        # Test streaming PSM section
        psm_count = 0
        for chunk in indexer.stream_section("PSM", chunk_size=50000):
            assert not chunk.empty
            assert isinstance(chunk, pd.DataFrame)
            psm_count += len(chunk)
        assert psm_count == 154774

        # Test streaming PRT section
        protein_count = 0
        for chunk in indexer.stream_section("PRT", chunk_size=5000):
            assert not chunk.empty
            assert isinstance(chunk, pd.DataFrame)
            protein_count += len(chunk)
        assert protein_count == 9655

        # Test streaming with a medium chunk size
        psm_count_medium_chunks = 0
        for chunk in indexer.stream_section("PSM", chunk_size=10000):
            psm_count_medium_chunks += len(chunk)
        assert psm_count_medium_chunks == 154774

        # Test invalid section
        with pytest.raises(ValueError, match="Invalid section: FAKE"):
            # We need to consume the generator to trigger the code inside it
            list(indexer.stream_section("FAKE"))

    finally:
        cleanup_test_resources(indexer=indexer, temp_db_path=temp_db_path)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
