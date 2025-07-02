"""Tests for MzTab indexer functionality."""

import os
from pathlib import Path
import pytest
import pandas as pd

from quantmsio.core.duckdb import MzTabIndexer
from quantmsio.core.quantms.mztab import MzTab

# Test data paths
TEST_DATA_ROOT = Path("tests/examples/quantms")
DDA_PLEX_DATA = {
    "mztab": TEST_DATA_ROOT
    / "dda-plex-full"
    / "PXD007683TMT.sdrf_openms_design_openms.mzTab.gz",
    "msstats": TEST_DATA_ROOT
    / "dda-plex-full"
    / "PXD007683TMT.sdrf_openms_design_msstats_in.csv.gz",
    "sdrf": TEST_DATA_ROOT / "dda-plex-full" / "PXD007683-TMT.sdrf.tsv",
}
DDA_LFQ_DATA = {
    "mztab": TEST_DATA_ROOT
    / "dda-lfq-full"
    / "PXD007683-LFQ.sdrf_openms_design_openms.mzTab.gz",
    "msstats": TEST_DATA_ROOT
    / "dda-lfq-full"
    / "PXD007683-LFQ.sdrf_openms_design_msstats_in.csv.gz",
    "sdrf": TEST_DATA_ROOT / "dda-lfq-full" / "PXD007683-LFQ.sdrf.tsv",
}


def test_dda_plex_dataset():
    """Test MzTabIndexer functionality with DDA-plex dataset."""
    # Skip if files don't exist
    if not all(path.exists() for path in DDA_PLEX_DATA.values()):
        pytest.skip("DDA-plex test files not found")

    # Initialize MzTabIndexer with both mzTab and MSstats files
    with MzTabIndexer(
        mztab_path=DDA_PLEX_DATA["mztab"],
        msstats_path=DDA_PLEX_DATA["msstats"],
        max_memory="8GB",
        worker_threads=4,
        batch_size=100000,
    ) as indexer:
        # Create tables from the files
        indexer.create_tables()

        # Test metadata retrieval
        metadata = indexer.get_metadata()
        assert not metadata.empty, "No metadata found"
        assert len(metadata) == 555
        print("\nDDA-plex dataset metadata statistics:")
        print(f"Number of metadata entries: {len(metadata)}")

        # Test protein data retrieval
        proteins = indexer.get_proteins()
        assert not proteins.empty, "No protein data found"
        assert len(proteins) == 9571
        print("\nProtein statistics:")
        print(f"Total proteins: {len(proteins):,}")
        print("\nProtein columns:")
        print(proteins.columns.tolist())

        # Test protein details data retrieval
        protein_details = indexer.get_protein_details()
        assert not protein_details.empty, "No protein details data found"
        assert len(protein_details) == 10940
        print(f"\nTotal protein details: {len(protein_details):,}")

        # Test PSM data retrieval
        psms = indexer.get_psms()
        assert not psms.empty, "No PSM data found"
        assert len(psms) == 154774
        print(f"\nTotal PSMs: {len(psms):,}")
        print("\nPSM columns:")
        print(psms.columns.tolist())

        # Test MSstats data retrieval
        msstats = indexer.get_msstats()
        assert msstats is not None and not msstats.empty, "No MSstats data found"
        assert len(msstats) == 1549526
        print(f"\nTotal MSstats entries: {len(msstats):,}")

        # Test anchor protein and accession logic
        protein_sample = indexer.query_to_df(
            "SELECT accession, anchor_protein, opt_global_result_type FROM proteins LIMIT 5"
        )
        for _, row in protein_sample.iterrows():
            if row["opt_global_result_type"] == "single_protein":
                assert (
                    row["accession"] == row["anchor_protein"]
                ), "Anchor protein should match accession for single proteins"
            elif row["opt_global_result_type"] == "indistinguishable_group":
                assert (
                    row["accession"] != row["anchor_protein"]
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
        protein_details_accessions = set(
            indexer.query_to_df("SELECT accession FROM protein_details")["accession"]
        )
        all_protein_accessions = protein_accessions.union(protein_details_accessions)

        psm_accessions = set(
            indexer.query_to_df(
                "SELECT DISTINCT accession FROM psms WHERE accession IS NOT NULL AND accession != 'null'"
            )["accession"]
        )
        for psm_acc_group in psm_accessions:
            for psm_acc in psm_acc_group.split(","):
                assert any(
                    psm_acc in prot_acc for prot_acc in all_protein_accessions
                ), f"PSM accession {psm_acc} not in any protein group"

        msstats_df = indexer.query_to_df(
            "SELECT DISTINCT ProteinName FROM msstats WHERE ProteinName IS NOT NULL"
        )
        msstats_accessions = set(
            msstats_df["ProteinName"].str.replace(";", ",", regex=False)
        )
        for msstats_acc_group in msstats_accessions:
            for msstats_acc in msstats_acc_group.split(","):
                assert any(
                    msstats_acc in prot_acc for prot_acc in all_protein_accessions
                ), f"MSstats accession {msstats_acc} not in any protein group"
        print("\nReferential integrity check passed for DDA-plex dataset.")


def test_dda_lfq_dataset():
    """Test MzTabIndexer functionality with DDA-LFQ dataset."""
    # Skip if files don't exist
    if not all(path.exists() for path in DDA_LFQ_DATA.values()):
        pytest.skip("DDA-LFQ test files not found")

    # Initialize MzTabIndexer with both mzTab and MSstats files
    with MzTabIndexer(
        mztab_path=DDA_LFQ_DATA["mztab"],
        msstats_path=DDA_LFQ_DATA["msstats"],
        max_memory="8GB",
        worker_threads=4,
        batch_size=50000,
    ) as indexer:
        # Create tables from the files
        indexer.create_tables()

        # Test metadata retrieval
        metadata = indexer.get_metadata()
        assert not metadata.empty, "No metadata found"
        assert len(metadata) == 109
        print("\nDDA-LFQ dataset metadata statistics:")
        print(f"Number of metadata entries: {len(metadata)}")

        # Test protein data retrieval
        proteins = indexer.get_proteins()
        assert not proteins.empty, "No protein data found"
        assert len(proteins) == 9216
        print("\nProtein statistics:")
        print(f"Total proteins: {len(proteins):,}")
        print("\nProtein columns:")
        print(proteins.columns.tolist())

        # Test protein details data retrieval
        protein_details = indexer.get_protein_details()
        assert not protein_details.empty, "No protein details data found"
        assert len(protein_details) == 9574
        print(f"\nTotal protein details: {len(protein_details):,}")

        # Test PSM data retrieval
        psms = indexer.get_psms()
        assert not psms.empty, "No PSM data found"
        assert len(psms) == 537848
        print(f"\nTotal PSMs: {len(psms):,}")
        print("\nPSM columns:")
        print(psms.columns.tolist())

        # Test MSstats data retrieval
        msstats = indexer.get_msstats()
        assert msstats is not None and not msstats.empty, "No MSstats data found"
        assert len(msstats) == 425594
        print(f"\nTotal MSstats entries: {len(msstats):,}")

        # Test anchor protein and accession logic
        protein_sample = indexer.query_to_df(
            "SELECT accession, anchor_protein, opt_global_result_type FROM proteins LIMIT 5"
        )
        for _, row in protein_sample.iterrows():
            if row["opt_global_result_type"] == "single_protein":
                assert (
                    row["accession"] == row["anchor_protein"]
                ), "Anchor protein should match accession for single proteins"
            elif row["opt_global_result_type"] == "indistinguishable_group":
                assert (
                    row["accession"] != row["anchor_protein"]
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
        protein_details_accessions = set(
            indexer.query_to_df("SELECT accession FROM protein_details")["accession"]
        )
        all_protein_accessions = protein_accessions.union(protein_details_accessions)

        psm_accessions = set(
            indexer.query_to_df(
                "SELECT DISTINCT accession FROM psms WHERE accession IS NOT NULL AND accession != 'null'"
            )["accession"]
        )
        for psm_acc_group in psm_accessions:
            for psm_acc in psm_acc_group.split(","):
                assert any(
                    psm_acc in prot_acc for prot_acc in all_protein_accessions
                ), f"PSM accession {psm_acc} not in any protein group"

        msstats_df = indexer.query_to_df(
            "SELECT DISTINCT ProteinName FROM msstats WHERE ProteinName IS NOT NULL"
        )
        msstats_accessions = set(
            msstats_df["ProteinName"].str.replace(";", ",", regex=False)
        )
        for msstats_acc_group in msstats_accessions:
            for msstats_acc in msstats_acc_group.split(","):
                assert any(
                    msstats_acc in prot_acc for prot_acc in all_protein_accessions
                ), f"MSstats accession {msstats_acc} not in any protein group"
        print("\nReferential integrity check passed for DDA-LFQ dataset.")


def test_add_msstats_to_existing_db(tmp_path):
    """Test adding an MSstats table to an existing MzTabIndexer database."""
    if not DDA_LFQ_DATA["mztab"].exists() or not DDA_LFQ_DATA["msstats"].exists():
        pytest.skip("DDA-LFQ test files not found")

    db_path = tmp_path / "existing.db"

    # Step 1: Create an indexer with only mzTab to create the initial database.
    # We don't use a 'with' statement here to prevent the database from being deleted on exit.
    indexer1 = MzTabIndexer(mztab_path=DDA_LFQ_DATA["mztab"], database_path=db_path)
    db_file_path = indexer1._duckdb_name  # Save the path

    # Check that msstats table doesn't exist
    tables = indexer1.query_to_df("SHOW TABLES")["name"].tolist()
    assert MzTabIndexer._MZTAB_INDEXER_TABLE_MSSTATS not in tables

    # Clean up the first indexer instance without deleting the database file.
    # This is a bit of a workaround to test this specific scenario.
    indexer1._duckdb.close()
    indexer1._duckdb = None
    indexer1._cleanup_temp_dir()

    # Step 2: Create a new indexer from the existing database and add the MSstats table.
    with MzTabIndexer(database_path=db_file_path) as indexer2:
        # The msstats table should not exist yet.
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

    # The DB file will be deleted by indexer2's __exit__ via destroy_database()
    assert not os.path.exists(db_file_path)


def test_mztab_parse_metadata_once():
    """Test MzTab._parse_metadata_once functionality."""
    if not DDA_PLEX_DATA["mztab"].exists():
        pytest.skip("DDA-plex test files not found")

    with MzTab(mztab_path=DDA_PLEX_DATA["mztab"]) as mztab_processor:
        mztab_processor.create_tables()

        # Check that metadata is not parsed yet
        assert not mztab_processor._metadata_parsed
        assert not mztab_processor._metadata_cache

        # Trigger parsing by calling a getter
        ms_runs = mztab_processor.extract_ms_runs()

        # Check that parsing happened and cache is filled
        assert mztab_processor._metadata_parsed
        assert mztab_processor._metadata_cache

        # Print the full metadata cache for inspection
        print("\nMetadata Cache Contents:")
        print("------------------------")
        print("\nMS Runs:")
        for run_id, filename in mztab_processor._metadata_cache["ms_runs"].items():
            print(f"{run_id}: {filename}")

        print("\nScore Names:")
        for score_name, value in mztab_processor._metadata_cache["score_names"].items():
            print(f"{score_name}: {value}")

        print("\nModifications:")
        for mod_id, details in mztab_processor._metadata_cache["modifications"].items():
            print(f"{mod_id}: {details}")

        print("\nMods Map:")
        for mod_name, details in mztab_processor._metadata_cache["mods_map"].items():
            print(f"{mod_name}: {details}")

        # Check MS runs - there are 11 runs in total (a05058 through a05068)
        assert len(ms_runs) == 11
        assert "ms_run[1]" in ms_runs
        assert ms_runs["ms_run[1]"] == "a05068"

        # Check score names - note that X!Tandem:hyperscore overwrites X!Tandem:expect
        score_names = mztab_processor.get_score_names()
        assert len(score_names) == 1
        assert "OpenMS" in score_names
        assert score_names["OpenMS"] == "search_engine_score[1]"

        # Check modifications
        modifications = mztab_processor.get_modifications()
        assert len(modifications) == 3
        assert "UNIMOD:4" in modifications
        assert modifications["UNIMOD:4"] == ["Carbamidomethyl", "1", "C", "Anywhere"]
        assert "UNIMOD:35" in modifications
        assert modifications["UNIMOD:35"] == ["Oxidation", "1", "M", "Anywhere"]
        assert "UNIMOD:737" in modifications
        assert modifications["UNIMOD:737"] == ["TMT6plex", "3", "X", "Any N-term"]

        # Check mods map
        mods_map = mztab_processor.get_mods_map()
        assert (
            len(mods_map) == 6
        )  # Each modification has two entries (name->accession and accession->name)
        assert mods_map["Carbamidomethyl"] == ["UNIMOD:4", "C"]
        assert mods_map["UNIMOD:4"] == ["Carbamidomethyl", "C"]
        assert mods_map["TMT6plex"] == ["UNIMOD:737", "X"]
        assert mods_map["UNIMOD:737"] == ["TMT6plex", "X"]

        # Check idempotency
        ms_runs_again = mztab_processor.extract_ms_runs()
        assert ms_runs == ms_runs_again


def test_mztab_stream_section():
    """Test MzTab.stream_section functionality."""
    if not DDA_PLEX_DATA["mztab"].exists():
        pytest.skip("DDA-plex test files not found")

    with MzTab(mztab_path=DDA_PLEX_DATA["mztab"]) as mztab_processor:
        mztab_processor.create_tables()

        # Test streaming PSM section
        psm_count = 0
        for chunk in mztab_processor.stream_section("PSM", chunk_size=50000):
            assert not chunk.empty
            assert isinstance(chunk, pd.DataFrame)
            psm_count += len(chunk)
        assert psm_count == 154774

        # Test streaming PRT section
        protein_count = 0
        for chunk in mztab_processor.stream_section("PRT", chunk_size=5000):
            assert not chunk.empty
            assert isinstance(chunk, pd.DataFrame)
            protein_count += len(chunk)
        assert protein_count == 9571

        # Test streaming with a medium chunk size
        psm_count_medium_chunks = 0
        for chunk in mztab_processor.stream_section("PSM", chunk_size=10000):
            psm_count_medium_chunks += len(chunk)
        assert psm_count_medium_chunks == 154774

        # Test invalid section
        with pytest.raises(ValueError, match="Invalid section: FAKE"):
            # We need to consume the generator to trigger the code inside it
            list(mztab_processor.stream_section("FAKE"))


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
