"""Tests for DIA-NN DuckDB functionality."""

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


def test_small_dataset():
    """Test DIA-NN DuckDB functionality with small dataset."""
    # Skip if files don't exist
    if not all(path.exists() for path in SMALL_DATA.values()):
        pytest.skip("Small test files not found")

    # Initialize DiannDuckDB
    diann_db = DiannDuckDB(
        diann_report_path=SMALL_DATA["report"], max_memory="4GB", worker_threads=2
    )

    try:
        # Test basic statistics
        unique_runs = diann_db.get_unique_values("report", "Run")
        assert len(unique_runs) > 0, "No runs found in the database"
        print(f"\nSmall dataset statistics:")
        print(f"Number of runs: {len(unique_runs)}")

        # Test protein groups
        protein_groups = diann_db.get_unique_values("report", "Protein.Group")
        assert len(protein_groups) > 0, "No protein groups found"
        print(f"Number of protein groups: {len(protein_groups)}")

        # Test peptide sequences
        peptides = diann_db.get_unique_values("report", "Modified.Sequence")
        assert len(peptides) > 0, "No peptide sequences found"
        print(f"Number of unique peptides: {len(peptides)}")

        # Test Q-value distribution
        q_value_stats = diann_db.query_to_df(
            """
            SELECT 
                COUNT(*) as total_psms,
                COUNT(CASE WHEN CAST("Q.Value" AS FLOAT) <= 0.01 THEN 1 END) as psms_q01,
                COUNT(CASE WHEN CAST("Q.Value" AS FLOAT) <= 0.05 THEN 1 END) as psms_q05
            FROM report
        """
        )
        assert q_value_stats["total_psms"].iloc[0] > 0, "No PSMs found"
        print("\nQ-value statistics:")
        print(f"Total PSMs: {q_value_stats['total_psms'].iloc[0]:,}")
        print(f"PSMs at 1% FDR: {q_value_stats['psms_q01'].iloc[0]:,}")
        print(f"PSMs at 5% FDR: {q_value_stats['psms_q05'].iloc[0]:,}")

        # Test protein group Q-value distribution
        pg_qvalue_stats = diann_db.query_to_df(
            """
            SELECT 
                COUNT(DISTINCT "Protein.Group") as total_proteins,
                COUNT(DISTINCT CASE WHEN CAST("PG.Q.Value" AS FLOAT) <= 0.01 THEN "Protein.Group" END) as proteins_q01,
                COUNT(DISTINCT CASE WHEN CAST("PG.Q.Value" AS FLOAT) <= 0.05 THEN "Protein.Group" END) as proteins_q05
            FROM report
        """
        )
        assert pg_qvalue_stats["total_proteins"].iloc[0] > 0, "No protein groups found"
        print("\nProtein group Q-value statistics:")
        print(f"Total protein groups: {pg_qvalue_stats['total_proteins'].iloc[0]:,}")
        print(f"Protein groups at 1% FDR: {pg_qvalue_stats['proteins_q01'].iloc[0]:,}")
        print(f"Protein groups at 5% FDR: {pg_qvalue_stats['proteins_q05'].iloc[0]:,}")

        # Test intensity distribution
        intensity_stats = diann_db.query_to_df(
            """
            SELECT 
                MIN(CAST("Precursor.Quantity" AS FLOAT)) as min_intensity,
                MAX(CAST("Precursor.Quantity" AS FLOAT)) as max_intensity,
                AVG(CAST("Precursor.Quantity" AS FLOAT)) as avg_intensity
            FROM report
            WHERE "Precursor.Quantity" IS NOT NULL
        """
        )
        assert not intensity_stats.empty, "No intensity values found"
        print("\nIntensity statistics:")
        print(f"Min intensity: {intensity_stats['min_intensity'].iloc[0]:.2e}")
        print(f"Max intensity: {intensity_stats['max_intensity'].iloc[0]:.2e}")
        print(f"Average intensity: {intensity_stats['avg_intensity'].iloc[0]:.2e}")

    finally:
        diann_db.destroy_database()


@pytest.mark.large_data
def test_full_dataset():
    """Test DIA-NN DuckDB functionality with full dataset."""
    # Skip if files don't exist
    if not all(path.exists() for path in FULL_DATA.values()):
        pytest.skip("Full test files not found")

    # Initialize DiannDuckDB with protein groups matrix
    diann_db = DiannDuckDB(
        diann_report_path=FULL_DATA["report"],
        max_memory="16GB",
        worker_threads=4,
        pg_matrix_path=FULL_DATA["pg_matrix"],
    )

    try:
        # Test basic statistics
        unique_runs = diann_db.get_unique_values("report", "Run")
        assert len(unique_runs) > 0, "No runs found in the database"
        print(f"\nFull dataset statistics:")
        print(f"Number of runs: {len(unique_runs)}")

        # Test protein groups
        protein_groups = diann_db.get_unique_values("report", "Protein.Group")
        assert len(protein_groups) > 0, "No protein groups found"
        print(f"Number of protein groups: {len(protein_groups)}")

        # Test peptide sequences
        peptides = diann_db.get_unique_values("report", "Modified.Sequence")
        assert len(peptides) > 0, "No peptide sequences found"
        print(f"Number of unique peptides: {len(peptides)}")

        # Test Q-value distribution
        q_value_stats = diann_db.query_to_df(
            """
            SELECT 
                COUNT(*) as total_psms,
                COUNT(CASE WHEN CAST("Q.Value" AS FLOAT) <= 0.01 THEN 1 END) as psms_q01,
                COUNT(CASE WHEN CAST("Q.Value" AS FLOAT) <= 0.05 THEN 1 END) as psms_q05
            FROM report
        """
        )
        assert q_value_stats["total_psms"].iloc[0] > 0, "No PSMs found"
        print("\nQ-value statistics:")
        print(f"Total PSMs: {q_value_stats['total_psms'].iloc[0]:,}")
        print(f"PSMs at 1% FDR: {q_value_stats['psms_q01'].iloc[0]:,}")
        print(f"PSMs at 5% FDR: {q_value_stats['psms_q05'].iloc[0]:,}")

        # Test protein group Q-value distribution
        pg_qvalue_stats = diann_db.query_to_df(
            """
            SELECT 
                COUNT(DISTINCT "Protein.Group") as total_proteins,
                COUNT(DISTINCT CASE WHEN CAST("PG.Q.Value" AS FLOAT) <= 0.01 THEN "Protein.Group" END) as proteins_q01,
                COUNT(DISTINCT CASE WHEN CAST("PG.Q.Value" AS FLOAT) <= 0.05 THEN "Protein.Group" END) as proteins_q05
            FROM report
        """
        )
        assert pg_qvalue_stats["total_proteins"].iloc[0] > 0, "No protein groups found"
        print("\nProtein group Q-value statistics:")
        print(f"Total protein groups: {pg_qvalue_stats['total_proteins'].iloc[0]:,}")
        print(f"Protein groups at 1% FDR: {pg_qvalue_stats['proteins_q01'].iloc[0]:,}")
        print(f"Protein groups at 5% FDR: {pg_qvalue_stats['proteins_q05'].iloc[0]:,}")

        # Test intensity distribution
        intensity_stats = diann_db.query_to_df(
            """
            SELECT 
                MIN(CAST("Precursor.Quantity" AS FLOAT)) as min_intensity,
                MAX(CAST("Precursor.Quantity" AS FLOAT)) as max_intensity,
                AVG(CAST("Precursor.Quantity" AS FLOAT)) as avg_intensity
            FROM report
            WHERE "Precursor.Quantity" IS NOT NULL
        """
        )
        assert not intensity_stats.empty, "No intensity values found"
        print("\nIntensity statistics:")
        print(f"Min intensity: {intensity_stats['min_intensity'].iloc[0]:.2e}")
        print(f"Max intensity: {intensity_stats['max_intensity'].iloc[0]:.2e}")
        print(f"Average intensity: {intensity_stats['avg_intensity'].iloc[0]:.2e}")

        # Additional statistics for full dataset
        gene_stats = diann_db.query_to_df(
            """
            SELECT 
                COUNT(DISTINCT "Genes") as total_genes,
                COUNT(DISTINCT CASE WHEN "Proteotypic" = 1 THEN "Genes" END) as proteotypic_genes
            FROM report
            WHERE "Genes" IS NOT NULL
        """
        )
        assert gene_stats["total_genes"].iloc[0] > 0, "No genes found"
        print("\nGene statistics:")
        print(f"Total genes: {gene_stats['total_genes'].iloc[0]:,}")
        print(f"Proteotypic genes: {gene_stats['proteotypic_genes'].iloc[0]:,}")

        # Test MaxLFQ values
        maxlfq_stats = diann_db.query_to_df(
            """
            SELECT 
                COUNT(DISTINCT "Protein.Group") as proteins_with_lfq,
                MIN(CAST("PG.MaxLFQ" AS FLOAT)) as min_lfq,
                MAX(CAST("PG.MaxLFQ" AS FLOAT)) as max_lfq,
                AVG(CAST("PG.MaxLFQ" AS FLOAT)) as avg_lfq
            FROM report
            WHERE "PG.MaxLFQ" IS NOT NULL
        """
        )
        assert not maxlfq_stats.empty, "No MaxLFQ values found"
        print("\nMaxLFQ statistics:")
        print(f"Proteins with MaxLFQ: {maxlfq_stats['proteins_with_lfq'].iloc[0]:,}")
        print(f"Min MaxLFQ: {maxlfq_stats['min_lfq'].iloc[0]:.2e}")
        print(f"Max MaxLFQ: {maxlfq_stats['max_lfq'].iloc[0]:.2e}")
        print(f"Average MaxLFQ: {maxlfq_stats['avg_lfq'].iloc[0]:.2e}")

        # Test protein groups matrix
        pg_matrix = diann_db.get_protein_group_matrix()
        assert not pg_matrix.empty, "No protein groups matrix data found"
        print("\nProtein groups matrix statistics:")
        print(f"Number of rows in matrix: {len(pg_matrix)}")
        print(f"Number of columns in matrix: {len(pg_matrix.columns)}")

        # Test protein groups matrix with specific proteins
        test_proteins = protein_groups[:5]  # Get first 5 protein groups
        filtered_matrix = diann_db.get_protein_group_matrix(test_proteins)
        assert not filtered_matrix.empty, "No filtered protein groups matrix data found"
        print(f"\nFiltered protein groups matrix statistics (5 proteins):")
        print(f"Number of rows in filtered matrix: {len(filtered_matrix)}")

        # Test referential integrity between report and pg_matrix
        report_protein_groups = set(
            diann_db.get_unique_values("report", "Protein.Group")
        )
        pg_matrix_protein_groups = set(
            diann_db.get_unique_values("pg_matrix", "Protein.Group")
        )
        assert pg_matrix_protein_groups.issubset(
            report_protein_groups
        ), "Found protein groups in pg_matrix not present in the report"
        print("\nReferential integrity check passed for DIA-NN full dataset.")

    finally:
        diann_db.destroy_database()


def test_caching():
    """Test caching functionality."""
    # Skip if files don't exist
    if not all(path.exists() for path in FULL_DATA.values()):
        pytest.skip("Full test files not found")

    # Initialize DiannDuckDB with protein groups matrix
    diann_db = DiannDuckDB(
        diann_report_path=FULL_DATA["report"],
        max_memory="16GB",
        worker_threads=4,
        pg_matrix_path=FULL_DATA["pg_matrix"],
        cache_size=128,
    )

    try:
        # Test caching of unique values
        start_time = time.time()
        unique_runs_1 = diann_db.get_unique_values("report", "Run")
        first_query_time = time.time() - start_time

        start_time = time.time()
        unique_runs_2 = diann_db.get_unique_values("report", "Run")
        cached_query_time = time.time() - start_time

        assert unique_runs_1 == unique_runs_2, "Cached results don't match"
        assert cached_query_time < first_query_time, "Caching not effective"
        print(f"\nCaching performance:")
        print(f"First query time: {first_query_time:.3f}s")
        print(f"Cached query time: {cached_query_time:.3f}s")

        # Test caching of statistics
        start_time = time.time()
        stats_1 = diann_db.get_statistics(0.01)
        first_stats_time = time.time() - start_time

        start_time = time.time()
        stats_2 = diann_db.get_statistics(0.01)
        cached_stats_time = time.time() - start_time

        assert stats_1 == stats_2, "Cached statistics don't match"
        assert cached_stats_time < first_stats_time, "Statistics caching not effective"
        print(f"First statistics time: {first_stats_time:.3f}s")
        print(f"Cached statistics time: {cached_stats_time:.3f}s")

        # Test cache clearing
        diann_db.clear_cache()
        start_time = time.time()
        stats_3 = diann_db.get_statistics(0.01)
        cleared_cache_time = time.time() - start_time

        assert stats_1 == stats_3, "Results after cache clear don't match"
        assert (
            cleared_cache_time > cached_stats_time
        ), "Cache was not cleared effectively"
        print(f"Query time after cache clear: {cleared_cache_time:.3f}s")

    finally:
        diann_db.destroy_database()


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
