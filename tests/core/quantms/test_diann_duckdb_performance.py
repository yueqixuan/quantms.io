"""Performance tests for DIA-NN DuckDB functionality."""

import os
import time
from pathlib import Path

import pandas as pd
import pytest

from qpx.core.diann.diann import DiannDuckDB

# Test data paths
TEST_DATA_ROOT = Path(__file__).parents[2] / "examples" / "diann"
SMALL_DATA = {
    "report": TEST_DATA_ROOT / "small" / "diann_report.tsv",
    "sdrf": TEST_DATA_ROOT / "small" / "PXD019909-DIA.sdrf.tsv",
}
FULL_DATA = {
    "report": TEST_DATA_ROOT / "full" / "diann_report.tsv.gz",
    "sdrf": TEST_DATA_ROOT / "full" / "PXD036609.sdrf.tsv",
    "pg_matrix": TEST_DATA_ROOT / "full" / "diann_report.pg_matrix.tsv",
}


def measure_time(func):
    """Decorator to measure execution time."""

    def wrapper(*args, **kwargs):
        start_time = time.time()
        result = func(*args, **kwargs)
        end_time = time.time()
        execution_time = end_time - start_time
        print(f"\nExecution time for {func.__name__}: {execution_time:.2f} seconds")
        return result

    return wrapper


@measure_time
def test_small_dataset_performance():
    """Test DIA-NN DuckDB performance with small dataset."""
    # Skip if files don't exist
    if not all(path.exists() for path in SMALL_DATA.values()):
        pytest.skip("Small test files not found")

    print("\nSmall Dataset Performance Test")
    print("-----------------------------")

    # Initialize DiannDuckDB and measure initialization time
    init_start = time.time()
    diann_db = DiannDuckDB(
        diann_report_path=SMALL_DATA["report"], max_memory="4GB", worker_threads=2
    )
    init_time = time.time() - init_start
    print(f"Database initialization time: {init_time:.2f} seconds")

    try:
        # Test query performance
        queries = [
            ("Basic count", "SELECT COUNT(*) FROM report"),
            (
                "Unique protein groups",
                'SELECT COUNT(DISTINCT "Protein.Group") FROM report',
            ),
            (
                "Q-value filtering",
                """
                SELECT COUNT(*) 
                FROM report 
                WHERE CAST("Q.Value" AS FLOAT) <= 0.01
            """,
            ),
            (
                "Complex aggregation",
                """
                SELECT "Protein.Group",
                    COUNT(*) as psm_count,
                    AVG(CAST("Q.Value" AS FLOAT)) as avg_qvalue,
                    MAX(CAST("Precursor.Quantity" AS FLOAT)) as max_intensity
                FROM report
                GROUP BY "Protein.Group"
                ORDER BY psm_count DESC
                LIMIT 100
            """,
            ),
        ]

        for query_name, query in queries:
            start_time = time.time()
            result = diann_db.query_to_df(query)
            query_time = time.time() - start_time
            print(f"\n{query_name}:")
            print(f"- Execution time: {query_time:.2f} seconds")
            print(f"- Result size: {len(result)} rows")

    finally:
        diann_db.destroy_database()


@pytest.mark.large_data
@measure_time
def test_full_dataset_performance():
    """Test DIA-NN DuckDB performance with full dataset."""
    # Skip if files don't exist
    if not all(path.exists() for path in FULL_DATA.values()):
        pytest.skip("Full test files not found")

    print("\nFull Dataset Performance Test")
    print("---------------------------")

    # Initialize DiannDuckDB and measure initialization time
    init_start = time.time()
    diann_db = DiannDuckDB(
        diann_report_path=FULL_DATA["report"],
        max_memory="16GB",
        worker_threads=4,
        pg_matrix_path=FULL_DATA["pg_matrix"],
    )
    init_time = time.time() - init_start
    print(f"Database initialization time: {init_time:.2f} seconds")

    try:
        # Test query performance
        queries = [
            ("Basic count", "SELECT COUNT(*) FROM report"),
            (
                "Unique protein groups",
                'SELECT COUNT(DISTINCT "Protein.Group") FROM report',
            ),
            (
                "Q-value filtering",
                """
                SELECT COUNT(*) 
                FROM report 
                WHERE CAST("Q.Value" AS FLOAT) <= 0.01
            """,
            ),
            (
                "Complex aggregation",
                """
                SELECT "Protein.Group",
                    COUNT(*) as psm_count,
                    AVG(CAST("Q.Value" AS FLOAT)) as avg_qvalue,
                    MAX(CAST("Precursor.Quantity" AS FLOAT)) as max_intensity
                FROM report
                GROUP BY "Protein.Group"
                ORDER BY psm_count DESC
                LIMIT 100
            """,
            ),
            (
                "Join with protein groups matrix",
                """
                SELECT r."Protein.Group",
                    COUNT(DISTINCT r."Run") as num_runs,
                    COUNT(*) as psm_count,
                    COUNT(DISTINCT COLUMNS(*)) as matrix_columns
                FROM report r
                LEFT JOIN pg_matrix pg ON r."Protein.Group" = pg."Protein.Group"
                GROUP BY r."Protein.Group"
                ORDER BY psm_count DESC
                LIMIT 100
            """,
            ),
        ]

        for query_name, query in queries:
            start_time = time.time()
            result = diann_db.query_to_df(query)
            query_time = time.time() - start_time
            print(f"\n{query_name}:")
            print(f"- Execution time: {query_time:.2f} seconds")
            print(f"- Result size: {len(result)} rows")

        # Test protein groups matrix performance
        print("\nProtein Groups Matrix Performance:")

        # Full matrix retrieval
        start_time = time.time()
        matrix = diann_db.get_protein_group_matrix()
        matrix_time = time.time() - start_time
        print(f"- Full matrix retrieval: {matrix_time:.2f} seconds")
        print(f"- Matrix size: {len(matrix)} rows × {len(matrix.columns)} columns")

        # Filtered matrix retrieval
        protein_groups = diann_db.get_unique_values("report", "Protein.Group")[:10]
        start_time = time.time()
        filtered_matrix = diann_db.get_protein_group_matrix(protein_groups)
        filtered_time = time.time() - start_time
        print(f"- Filtered matrix retrieval (10 proteins): {filtered_time:.2f} seconds")
        print(f"- Filtered matrix size: {len(filtered_matrix)} rows")

    finally:
        diann_db.destroy_database()


if __name__ == "__main__":
    pytest.main([__file__, "-v", "-s"])
