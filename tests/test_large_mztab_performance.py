"""Performance tests for MzTabIndexer with large datasets."""

import os
import time
from pathlib import Path

import pytest
from quantmsio.core.duckdb import MzTabIndexer

# Define the path to the large dataset
LARGE_MZTAB_DATA_ROOT = Path("tissues/PXD020192")
LARGE_DATASET = {
    "mztab": LARGE_MZTAB_DATA_ROOT / "PXD020192.sdrf_openms_design_openms.mzTab.gz",
    "msstats": LARGE_MZTAB_DATA_ROOT / "PXD020192.sdrf_openms_design_msstats_in.csv.gz",
}


# This marker skips the test if the 'CI' environment variable is set
# Most CI systems (GitHub Actions, Travis CI, etc.) set this variable by default
@pytest.mark.skipif(
    os.environ.get("CI") == "true", reason="Skipping performance test in CI environment"
)
def test_large_mztab_performance():
    """
    Tests the performance of MzTabIndexer on a large mzTab dataset.
    This test is intended for local execution and will be skipped in CI.
    """
    # Check if the large dataset files exist
    if not all(path.exists() for path in LARGE_DATASET.values()):
        pytest.skip("Large mzTab test files not found. Skipping performance test.")

    print("\n--- Starting Large mzTab Performance Test ---")

    start_time = time.time()

    with MzTabIndexer(
        mztab_path=LARGE_DATASET["mztab"],
        msstats_path=LARGE_DATASET["msstats"],
        max_memory="16GB",  # Adjust as needed for your system
        worker_threads=4,  # Adjust as needed for your system
        batch_size=100000,  # Using a larger batch size for a larger file
    ) as indexer:
        # Time to create tables
        t0 = time.time()
        indexer.create_tables()
        t1 = time.time()
        print(f"Time to create tables: {t1 - t0:.2f} seconds")

        # Get and print counts
        protein_count = indexer.get_protein_count()
        protein_details_count = indexer.get_protein_details_count()
        psm_count = indexer.get_psm_count()

        print(f"\nTotal proteins: {protein_count:,}")
        print(f"Total protein details: {protein_details_count:,}")
        print(f"Total PSMs: {psm_count:,}")

        msstats = indexer.get_msstats()
        if msstats is not None:
            print(f"Total MSstats entries: {len(msstats):,}")

    end_time = time.time()
    print(f"\n--- Total test duration: {end_time - start_time:.2f} seconds ---")
