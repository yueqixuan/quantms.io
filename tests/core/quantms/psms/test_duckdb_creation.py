import tempfile
import time
from pathlib import Path

import duckdb
import pytest
from qpx.core.quantms.mztab import MzTabIndexer

TEST_DATA_ROOT = Path(__file__).parents[3] / "examples"


# Define test cases to be run, now including the backend type
# All tests use full examples (dda-lfq-full and dda-plex-full)
TEST_CASES = [
    (
        "lfq_duckdb",
        TEST_DATA_ROOT
        / "quantms/dda-lfq-full/PXD007683-LFQ.sdrf_openms_design_openms.mzTab.gz",
        TEST_DATA_ROOT
        / "quantms/dda-lfq-full/PXD007683-LFQ.sdrf_openms_design_msstats_in.csv.gz",
        TEST_DATA_ROOT / "quantms/dda-lfq-full/PXD007683-LFQ.sdrf.tsv",
    ),
    (
        "tmt_duckdb",
        TEST_DATA_ROOT
        / "quantms/dda-plex-full/PXD007683TMT.sdrf_openms_design_openms.mzTab.gz",
        TEST_DATA_ROOT
        / "quantms/dda-plex-full/PXD007683TMT.sdrf_openms_design_msstats_in.csv.gz",
        TEST_DATA_ROOT / "quantms/dda-plex-full/PXD007683-TMT.sdrf.tsv",
    ),
    pytest.param(
        "pxd020192_duckdb",
        Path(
            "/Users/yperez/work/QPX/tissues/PXD020192/PXD020192.sdrf_openms_design_openms.mzTab.gz"
        ),
        Path(
            "/Users/yperez/work/QPX/tissues/PXD020192/PXD020192.sdrf_openms_design_msstats_in.csv.gz"
        ),
        Path("/Users/yperez/work/QPX/tissues/PXD020192/PXD020192.sdrf.tsv"),
        marks=pytest.mark.skip(reason="Local test, big data file"),
    ),
]


@pytest.mark.parametrize("test_name, mztab_path, msstats_path, sdrf_path", TEST_CASES)
def test_create_duckdb_and_parquet(test_name, mztab_path, msstats_path, sdrf_path):
    """
    Test creating a DuckDB database or Parquet store from various mzTab formats.
    This single test function is parameterized to run for each case in TEST_CASES.
    """
    if not mztab_path.exists():
        pytest.skip(f"Test file not found: {mztab_path}")

    with tempfile.TemporaryDirectory() as temp_dir:
        # The output path will be a file for duckdb and a directory for parquet
        output_path = Path(temp_dir) / f"test_{test_name}.duckdb"

        # Start timing
        start_time = time.time()

        # Run the command to create the database or parquet store using the new interface
        mztab = MzTabIndexer.create(
            mztab_path=mztab_path,
            msstats_path=msstats_path,
            sdrf_path=sdrf_path,
            database_path=output_path,
        )

        # End timing and print result
        end_time = time.time()
        elapsed_time = end_time - start_time
        print(f"\nTime to process '{test_name}': {elapsed_time:.2f} seconds")

        # Verification logic depends on the used
        con = None
        try:
            # Close the MzTabIndexer connection first
            cleanup_mztab_indexer(mztab)
            assert output_path.exists(), f"DuckDB file not created: {output_path}"
            con = duckdb.connect(str(output_path), read_only=True)
            tables = [t[0] for t in con.execute("SHOW TABLES").fetchall()]
            assert "metadata" in tables and "proteins" in tables and "psms" in tables
            metadata_count = con.execute("SELECT COUNT(*) FROM metadata").fetchone()[0]
            proteins_count = con.execute("SELECT COUNT(*) FROM proteins").fetchone()[0]
            psms_count = con.execute("SELECT COUNT(*) FROM psms").fetchone()[0]
            if msstats_path:
                msstats_count = con.execute("SELECT COUNT(*) FROM msstats").fetchone()[
                    0
                ]

            # Common assertions for both backends
            assert metadata_count > 0, "No metadata rows found"
            assert proteins_count > 0, "No protein rows found"
            assert psms_count > 0, "No PSM rows found"
            if msstats_path:
                assert msstats_count > 0, "No MSstats rows found"

            # Print final statistics for verification
            # print(f"\n{test_name.upper()} ({backend.upper()}) Statistics:")
            print(f"\n{test_name.upper()} Statistics:")
            print(f"Metadata rows: {metadata_count:,}")
            print(f"Protein rows: {proteins_count:,}")
            print(f"PSM rows: {psms_count:,}")
            if msstats_path:
                print(f"MSstats rows: {msstats_count:,}")

        finally:
            # Clean up connections
            cleanup_mztab_indexer(mztab)
            cleanup_duckdb_connection(con)


def test_two_step_creation_process():
    """
    Test the two-step process:
    1. Create database with mzTab data only
    2. Close and reopen the database
    3. Add MSstats data to the existing database using add_msstats_table
    """
    # Use the LFQ dataset for this test
    mztab_path = (
        TEST_DATA_ROOT
        / "quantms/dda-lfq-full/PXD007683-LFQ.sdrf_openms_design_openms.mzTab.gz"
    )
    msstats_path = (
        TEST_DATA_ROOT
        / "quantms/dda-lfq-full/PXD007683-LFQ.sdrf_openms_design_msstats_in.csv.gz"
    )
    sdrf_path = TEST_DATA_ROOT / "quantms/dda-lfq-full/PXD007683-LFQ.sdrf.tsv"

    if not mztab_path.exists() or not msstats_path.exists():
        pytest.skip("Test files not found")

    with tempfile.TemporaryDirectory() as temp_dir:
        output_path = Path(temp_dir) / "test_two_step.duckdb"

        print(f"\n=== Testing Two-Step Creation Process ===")

        # Step 1: Create database with mzTab data only (no MSstats)
        print("Step 1: Creating database with mzTab data only...")
        start_time = time.time()

        mztab_step1 = MzTabIndexer.create(
            mztab_path=mztab_path,
            database_path=output_path,
        )

        # Verify that only mzTab tables exist, no MSstats table
        tables = [t[0] for t in mztab_step1._duckdb.execute("SHOW TABLES").fetchall()]
        assert "metadata" in tables
        assert "proteins" in tables
        assert "psms" in tables
        assert "msstats" not in tables, "MSstats table should not exist yet"

        # Get initial counts
        metadata_count = mztab_step1._duckdb.execute(
            "SELECT COUNT(*) FROM metadata"
        ).fetchone()[0]
        proteins_count = mztab_step1._duckdb.execute(
            "SELECT COUNT(*) FROM proteins"
        ).fetchone()[0]
        psms_count = mztab_step1._duckdb.execute(
            "SELECT COUNT(*) FROM psms"
        ).fetchone()[0]

        print(f"Initial database created with:")
        print(f"  Metadata rows: {metadata_count:,}")
        print(f"  Protein rows: {proteins_count:,}")
        print(f"  PSM rows: {psms_count:,}")

        # Close the database connection
        mztab_step1._duckdb.close()
        mztab_step1._duckdb = None

        step1_time = time.time() - start_time
        print(f"Step 1 completed in {step1_time:.2f} seconds")

        # Step 2: Reopen the existing database and add MSstats data
        print("\nStep 2: Reopening database and adding MSstats data...")
        start_time = time.time()

        mztab_step2 = MzTabIndexer.open(
            database_path=output_path,
            sdrf_path=sdrf_path,
        )

        # Verify the existing tables are still there
        tables = [t[0] for t in mztab_step2._duckdb.execute("SHOW TABLES").fetchall()]
        assert "metadata" in tables
        assert "proteins" in tables
        assert "psms" in tables
        assert "msstats" not in tables, "MSstats table should not exist yet"

        # Add MSstats table to the existing database
        mztab_step2.add_msstats_table(str(msstats_path))

        # Verify that MSstats table now exists
        tables = [t[0] for t in mztab_step2._duckdb.execute("SHOW TABLES").fetchall()]
        assert "msstats" in tables, "MSstats table should now exist"

        # Get final counts including MSstats
        final_metadata_count = mztab_step2._duckdb.execute(
            "SELECT COUNT(*) FROM metadata"
        ).fetchone()[0]
        final_proteins_count = mztab_step2._duckdb.execute(
            "SELECT COUNT(*) FROM proteins"
        ).fetchone()[0]
        final_psms_count = mztab_step2._duckdb.execute(
            "SELECT COUNT(*) FROM psms"
        ).fetchone()[0]
        msstats_count = mztab_step2._duckdb.execute(
            "SELECT COUNT(*) FROM msstats"
        ).fetchone()[0]

        # Verify counts are consistent
        assert (
            final_metadata_count == metadata_count
        ), "Metadata count should not change"
        assert final_proteins_count == proteins_count, "Protein count should not change"
        assert final_psms_count == psms_count, "PSM count should not change"
        assert msstats_count > 0, "MSstats table should have data"

        step2_time = time.time() - start_time
        print(f"Step 2 completed in {step2_time:.2f} seconds")

        print(f"\nFinal database contains:")
        print(f"  Metadata rows: {final_metadata_count:,}")
        print(f"  Protein rows: {final_proteins_count:,}")
        print(f"  PSM rows: {final_psms_count:,}")
        print(f"  MSstats rows: {msstats_count:,}")

        total_time = step1_time + step2_time
        print(f"\nTotal two-step process time: {total_time:.2f} seconds")

        # Clean up
        mztab_step2._duckdb.close()


def test_two_step_creation_process_tmt():
    """
    Test the two-step process with TMT data:
    1. Create database with mzTab data only
    2. Close and reopen the database
    3. Add MSstats data to the existing database using add_msstats_table
    """
    # Use the TMT dataset for this test
    mztab_path = (
        TEST_DATA_ROOT
        / "quantms/dda-plex-full/PXD007683TMT.sdrf_openms_design_openms.mzTab.gz"
    )
    msstats_path = (
        TEST_DATA_ROOT
        / "quantms/dda-plex-full/PXD007683TMT.sdrf_openms_design_msstats_in.csv.gz"
    )
    sdrf_path = TEST_DATA_ROOT / "quantms/dda-plex-full/PXD007683-TMT.sdrf.tsv"

    if not mztab_path.exists() or not msstats_path.exists():
        pytest.skip("Test files not found")

    with tempfile.TemporaryDirectory() as temp_dir:
        output_path = Path(temp_dir) / "test_two_step_tmt.duckdb"

        print(f"\n=== Testing Two-Step Creation Process (TMT) ===")

        # Step 1: Create database with mzTab data only (no MSstats)
        print("Step 1: Creating database with mzTab data only...")
        start_time = time.time()

        mztab_step1 = MzTabIndexer.create(
            mztab_path=mztab_path,
            msstats_path=None,  # No MSstats in first step
            database_path=output_path,
        )

        # Verify that only mzTab tables exist, no MSstats table
        tables = [t[0] for t in mztab_step1._duckdb.execute("SHOW TABLES").fetchall()]
        assert "metadata" in tables
        assert "proteins" in tables
        assert "psms" in tables
        assert "msstats" not in tables, "MSstats table should not exist yet"

        # Get initial counts
        metadata_count = mztab_step1._duckdb.execute(
            "SELECT COUNT(*) FROM metadata"
        ).fetchone()[0]
        proteins_count = mztab_step1._duckdb.execute(
            "SELECT COUNT(*) FROM proteins"
        ).fetchone()[0]
        psms_count = mztab_step1._duckdb.execute(
            "SELECT COUNT(*) FROM psms"
        ).fetchone()[0]

        print(f"Initial database created with:")
        print(f"  Metadata rows: {metadata_count:,}")
        print(f"  Protein rows: {proteins_count:,}")
        print(f"  PSM rows: {psms_count:,}")

        # Close the database connection
        mztab_step1._duckdb.close()
        mztab_step1._duckdb = None

        step1_time = time.time() - start_time
        print(f"Step 1 completed in {step1_time:.2f} seconds")

        # Step 2: Reopen the existing database and add MSstats data
        print("\nStep 2: Reopening database and adding MSstats data...")
        start_time = time.time()

        mztab_step2 = MzTabIndexer.open(database_path=output_path, sdrf_path=sdrf_path)

        # Verify the existing tables are still there
        tables = [t[0] for t in mztab_step2._duckdb.execute("SHOW TABLES").fetchall()]
        assert "metadata" in tables
        assert "proteins" in tables
        assert "psms" in tables
        assert "msstats" not in tables, "MSstats table should not exist yet"

        # Add MSstats table to the existing database
        mztab_step2.add_msstats_table(str(msstats_path))

        # Verify that MSstats table now exists
        tables = [t[0] for t in mztab_step2._duckdb.execute("SHOW TABLES").fetchall()]
        assert "msstats" in tables, "MSstats table should now exist"

        # Get final counts including MSstats
        final_metadata_count = mztab_step2._duckdb.execute(
            "SELECT COUNT(*) FROM metadata"
        ).fetchone()[0]
        final_proteins_count = mztab_step2._duckdb.execute(
            "SELECT COUNT(*) FROM proteins"
        ).fetchone()[0]
        final_psms_count = mztab_step2._duckdb.execute(
            "SELECT COUNT(*) FROM psms"
        ).fetchone()[0]
        msstats_count = mztab_step2._duckdb.execute(
            "SELECT COUNT(*) FROM msstats"
        ).fetchone()[0]

        # Verify counts are consistent
        assert (
            final_metadata_count == metadata_count
        ), "Metadata count should not change"
        assert final_proteins_count == proteins_count, "Protein count should not change"
        assert final_psms_count == psms_count, "PSM count should not change"
        assert msstats_count > 0, "MSstats table should have data"

        step2_time = time.time() - start_time
        print(f"Step 2 completed in {step2_time:.2f} seconds")

        print(f"\nFinal database contains:")
        print(f"  Metadata rows: {final_metadata_count:,}")
        print(f"  Protein rows: {final_proteins_count:,}")
        print(f"  PSM rows: {final_psms_count:,}")
        print(f"  MSstats rows: {msstats_count:,}")

        total_time = step1_time + step2_time
        print(f"\nTotal two-step process time (TMT): {total_time:.2f} seconds")

        # Clean up
        mztab_step2._duckdb.close()


def test_create_vs_open_database():
    """
    Test creating a new database vs opening an existing one using class methods.
    This demonstrates the flexibility of the MzTabIndexer class with clear API.
    """
    # Use the LFQ dataset for this test
    mztab_path = (
        TEST_DATA_ROOT
        / "quantms/dda-lfq-full/PXD007683-LFQ.sdrf_openms_design_openms.mzTab.gz"
    )
    msstats_path = (
        TEST_DATA_ROOT
        / "quantms/dda-lfq-full/PXD007683-LFQ.sdrf_openms_design_msstats_in.csv.gz"
    )
    sdrf_path = TEST_DATA_ROOT / "quantms/dda-lfq-full/PXD007683-LFQ.sdrf.tsv"

    if not mztab_path.exists() or not msstats_path.exists():
        pytest.skip("Test files not found")

    with tempfile.TemporaryDirectory() as temp_dir:
        output_path = Path(temp_dir) / "test_create_vs_open.duckdb"

        print(f"\n=== Testing Create vs Open Database with Class Methods ===")

        # Test 1: Create a new database using create()
        print("\nTest 1: Creating a new database using MzTabIndexer.create()...")
        start_time = time.time()

        mztab_new = MzTabIndexer.create(
            mztab_path=mztab_path,
            msstats_path=msstats_path,
            sdrf_path=sdrf_path,
            database_path=output_path,
        )

        # Verify all tables exist
        tables = [t[0] for t in mztab_new._duckdb.execute("SHOW TABLES").fetchall()]
        assert "metadata" in tables
        assert "proteins" in tables
        assert "psms" in tables
        assert "msstats" in tables

        # Get counts
        metadata_count = mztab_new._duckdb.execute(
            "SELECT COUNT(*) FROM metadata"
        ).fetchone()[0]
        proteins_count = mztab_new._duckdb.execute(
            "SELECT COUNT(*) FROM proteins"
        ).fetchone()[0]
        psms_count = mztab_new._duckdb.execute("SELECT COUNT(*) FROM psms").fetchone()[
            0
        ]
        msstats_count = mztab_new._duckdb.execute(
            "SELECT COUNT(*) FROM msstats"
        ).fetchone()[0]

        print(f"New database created with:")
        print(f"  Metadata rows: {metadata_count:,}")
        print(f"  Protein rows: {proteins_count:,}")
        print(f"  PSM rows: {psms_count:,}")
        print(f"  MSstats rows: {msstats_count:,}")

        # Close and delete the object to ensure clean connection release
        mztab_new._duckdb.close()
        del mztab_new

        create_time = time.time() - start_time
        print(f"Database creation completed in {create_time:.2f} seconds")

        # Test 2: Open the existing database using open()
        print("\nTest 2: Opening existing database using MzTabIndexer.open()...")
        start_time = time.time()

        mztab_existing = MzTabIndexer.open(
            database_path=output_path,
            sdrf_path=sdrf_path,
        )

        # Verify all tables still exist
        tables = [
            t[0] for t in mztab_existing._duckdb.execute("SHOW TABLES").fetchall()
        ]
        assert "metadata" in tables
        assert "proteins" in tables
        assert "psms" in tables
        assert "msstats" in tables

        # Get counts and verify they match
        final_metadata_count = mztab_existing._duckdb.execute(
            "SELECT COUNT(*) FROM metadata"
        ).fetchone()[0]
        final_proteins_count = mztab_existing._duckdb.execute(
            "SELECT COUNT(*) FROM proteins"
        ).fetchone()[0]
        final_psms_count = mztab_existing._duckdb.execute(
            "SELECT COUNT(*) FROM psms"
        ).fetchone()[0]
        final_msstats_count = mztab_existing._duckdb.execute(
            "SELECT COUNT(*) FROM msstats"
        ).fetchone()[0]

        # Verify counts are consistent
        assert (
            final_metadata_count == metadata_count
        ), "Metadata count should not change"
        assert final_proteins_count == proteins_count, "Protein count should not change"
        assert final_psms_count == psms_count, "PSM count should not change"
        assert final_msstats_count == msstats_count, "MSstats count should not change"

        print(f"Existing database opened with:")
        print(f"  Metadata rows: {final_metadata_count:,}")
        print(f"  Protein rows: {final_proteins_count:,}")
        print(f"  PSM rows: {final_psms_count:,}")
        print(f"  MSstats rows: {final_msstats_count:,}")

        open_time = time.time() - start_time
        print(f"Database opening completed in {open_time:.2f} seconds")

        # Test 3: Demonstrate that we can still add more data to the existing database
        print("\nTest 3: Adding additional data to existing database...")
        start_time = time.time()

        # Add the same MSstats data again (this will replace the existing table)
        mztab_existing.add_msstats_table(str(msstats_path))

        # Verify the table still exists and has the same count
        final_msstats_count_after_add = mztab_existing._duckdb.execute(
            "SELECT COUNT(*) FROM msstats"
        ).fetchone()[0]
        assert (
            final_msstats_count_after_add == msstats_count
        ), "MSstats count should remain the same after re-adding"

        add_time = time.time() - start_time
        print(f"Additional data added in {add_time:.2f} seconds")

        print(f"\nSummary:")
        print(f"  Create new database: {create_time:.2f} seconds")
        print(f"  Open existing database: {open_time:.2f} seconds")
        print(f"  Add data to existing: {add_time:.2f} seconds")

        # Clean up
        mztab_existing._duckdb.close()


def test_create_vs_open_database_tmt():
    """
    Test creating a new database vs opening an existing one using class methods with TMT data.
    This demonstrates the flexibility of the MzTabIndexer class with clear API.
    """
    # Use the TMT dataset for this test
    mztab_path = (
        TEST_DATA_ROOT
        / "quantms/dda-plex-full/PXD007683TMT.sdrf_openms_design_openms.mzTab.gz"
    )
    msstats_path = (
        TEST_DATA_ROOT
        / "quantms/dda-plex-full/PXD007683TMT.sdrf_openms_design_msstats_in.csv.gz"
    )
    sdrf_path = TEST_DATA_ROOT / "quantms/dda-plex-full/PXD007683-TMT.sdrf.tsv"

    if not mztab_path.exists() or not msstats_path.exists():
        pytest.skip("Test files not found")

    with tempfile.TemporaryDirectory() as temp_dir:
        output_path = Path(temp_dir) / "test_create_vs_open_tmt.duckdb"

        print(f"\n=== Testing Create vs Open Database with Class Methods (TMT) ===")

        # Test 1: Create a new database using create()
        print("\nTest 1: Creating a new database using MzTabIndexer.create()...")
        start_time = time.time()

        mztab_new = MzTabIndexer.create(
            mztab_path=mztab_path,
            msstats_path=msstats_path,
            sdrf_path=sdrf_path,
            database_path=output_path,
        )

        # Verify all tables exist
        tables = [t[0] for t in mztab_new._duckdb.execute("SHOW TABLES").fetchall()]
        assert "metadata" in tables
        assert "proteins" in tables
        assert "psms" in tables
        assert "msstats" in tables

        # Get counts
        metadata_count = mztab_new._duckdb.execute(
            "SELECT COUNT(*) FROM metadata"
        ).fetchone()[0]
        proteins_count = mztab_new._duckdb.execute(
            "SELECT COUNT(*) FROM proteins"
        ).fetchone()[0]
        psms_count = mztab_new._duckdb.execute("SELECT COUNT(*) FROM psms").fetchone()[
            0
        ]
        msstats_count = mztab_new._duckdb.execute(
            "SELECT COUNT(*) FROM msstats"
        ).fetchone()[0]

        print(f"New database created with:")
        print(f"  Metadata rows: {metadata_count:,}")
        print(f"  Protein rows: {proteins_count:,}")
        print(f"  PSM rows: {psms_count:,}")
        print(f"  MSstats rows: {msstats_count:,}")

        # Close and delete the object to ensure clean connection release
        mztab_new._duckdb.close()
        del mztab_new

        create_time = time.time() - start_time
        print(f"Database creation completed in {create_time:.2f} seconds")

        # Test 2: Open the existing database using open()
        print("\nTest 2: Opening existing database using MzTabIndexer.open()...")
        start_time = time.time()

        mztab_existing = MzTabIndexer.open(
            database_path=output_path,
            sdrf_path=sdrf_path,
        )

        # Verify all tables still exist
        tables = [
            t[0] for t in mztab_existing._duckdb.execute("SHOW TABLES").fetchall()
        ]
        assert "metadata" in tables
        assert "proteins" in tables
        assert "psms" in tables
        assert "msstats" in tables

        # Get counts and verify they match
        final_metadata_count = mztab_existing._duckdb.execute(
            "SELECT COUNT(*) FROM metadata"
        ).fetchone()[0]
        final_proteins_count = mztab_existing._duckdb.execute(
            "SELECT COUNT(*) FROM proteins"
        ).fetchone()[0]
        final_psms_count = mztab_existing._duckdb.execute(
            "SELECT COUNT(*) FROM psms"
        ).fetchone()[0]
        final_msstats_count = mztab_existing._duckdb.execute(
            "SELECT COUNT(*) FROM msstats"
        ).fetchone()[0]

        # Verify counts are consistent
        assert (
            final_metadata_count == metadata_count
        ), "Metadata count should not change"
        assert final_proteins_count == proteins_count, "Protein count should not change"
        assert final_psms_count == psms_count, "PSM count should not change"
        assert final_msstats_count == msstats_count, "MSstats count should not change"

        print(f"Existing database opened with:")
        print(f"  Metadata rows: {final_metadata_count:,}")
        print(f"  Protein rows: {final_proteins_count:,}")
        print(f"  PSM rows: {final_psms_count:,}")
        print(f"  MSstats rows: {final_msstats_count:,}")

        open_time = time.time() - start_time
        print(f"Database opening completed in {open_time:.2f} seconds")

        # Test 3: Demonstrate that we can still add more data to the existing database
        print("\nTest 3: Adding additional data to existing database...")
        start_time = time.time()

        # Add the same MSstats data again (this will replace the existing table)
        mztab_existing.add_msstats_table(str(msstats_path))

        # Verify the table still exists and has the same count
        final_msstats_count_after_add = mztab_existing._duckdb.execute(
            "SELECT COUNT(*) FROM msstats"
        ).fetchone()[0]
        assert (
            final_msstats_count_after_add == msstats_count
        ), "MSstats count should remain the same after re-adding"

        add_time = time.time() - start_time
        print(f"Additional data added in {add_time:.2f} seconds")

        print(f"\nSummary:")
        print(f"  Create new database: {create_time:.2f} seconds")
        print(f"  Open existing database: {open_time:.2f} seconds")
        print(f"  Add data to existing: {add_time:.2f} seconds")

        # Clean up
        mztab_existing._duckdb.close()


def cleanup_mztab_indexer(mztab_indexer):
    """Clean up MzTabIndexer instance and its connections."""
    if mztab_indexer and hasattr(mztab_indexer, "_duckdb") and mztab_indexer._duckdb:
        try:
            mztab_indexer._duckdb.close()
        except (AttributeError, RuntimeError, OSError):
            # Ignore cleanup errors - connection may already be closed or invalid
            pass


def cleanup_duckdb_connection(con):
    """Clean up DuckDB connection."""
    if con:
        try:
            con.close()
        except (AttributeError, RuntimeError, OSError):
            # Ignore cleanup errors - connection may already be closed or invalid
            pass
