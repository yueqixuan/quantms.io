"""
Tests for MSstats analysis functionality in MzTabIndexer.

This module tests the advanced MSstats analysis capabilities that were migrated
from MsstatsIN to MzTabIndexer, demonstrating enhanced performance and functionality
for both LFQ and TMT experiments.
"""

import gzip
import tempfile
from pathlib import Path

from qpx.core.quantms.mztab import MzTabIndexer
from qpx.core.project import create_uuid_filename

# Test data paths
TEST_DATA_ROOT = Path(__file__).parents[2] / "examples"


def _cleanup_mztab_temp_files(indexer):
    """Clean up temporary files created during MzTab processing.

    This function calls the existing cleanup method to remove the mztab_indexer_
    temporary folder if it exists. It's designed for test use only.

    Args:
        indexer: MzTabIndexer instance that might have created temp files
    """
    try:
        # Call the existing cleanup method to remove temporary parquet directory
        if hasattr(indexer, "_cleanup_temp_parquet_dir"):
            indexer._cleanup_temp_parquet_dir()
    except Exception as e:
        # Don't fail tests due to cleanup issues, just warn
        print(f"Warning: Failed to clean up temporary files: {e}")


def test_mztab_indexer_msstats_lfq_full_dataset():
    """Test MzTabIndexer MSstats processing with full LFQ dataset."""
    print("\n" + "=" * 60)
    print("Testing MzTabIndexer MSstats analysis with full LFQ dataset")
    print("=" * 60)

    # File paths
    msstats_gz_file = (
        TEST_DATA_ROOT
        / "quantms/dda-lfq-full/PXD007683-LFQ.sdrf_openms_design_msstats_in.csv.gz"
    )
    sdrf_file = TEST_DATA_ROOT / "quantms/dda-lfq-full/PXD007683-LFQ.sdrf.tsv"

    # Create temporary MzTabIndexer database
    temp_db_path = create_uuid_filename("test_msstats_lfq", ".duckdb")

    try:
        # Create MzTabIndexer with MSstats data
        indexer = MzTabIndexer(
            mztab_path=None,
            database_path=temp_db_path,
            sdrf_path=sdrf_file,
            max_memory="4GB",
            worker_threads=2,
        )

        # Add MSstats data to the indexer
        indexer.add_msstats_table(str(msstats_gz_file))

        # Test enhanced MSstats analysis methods
        experiment_type = indexer.get_msstats_experiment_type()
        print(f"Experiment type: {experiment_type}")
        print(f"Processing file: {msstats_gz_file.name}")
        print(f"SDRF file: {sdrf_file.name}")

        # Test file statistics
        file_stats = indexer.get_msstats_file_statistics()
        if file_stats is not None:
            print(f"\nFile statistics: {len(file_stats)} files")
            print(f"Total features: {file_stats['feature_count'].sum():,}")
            print(f"Total proteins: {file_stats['protein_count'].sum():,}")

        # Test protein summary
        protein_summary = indexer.get_msstats_protein_summary()
        if protein_summary is not None:
            print(f"Protein summary: {len(protein_summary)} proteins")

        # Test peptide summary
        peptide_summary = indexer.get_msstats_peptide_summary()
        if peptide_summary is not None:
            print(f"Peptide summary: {len(peptide_summary)} peptides")

        # Test intensity distribution
        intensity_dist = indexer.get_msstats_intensity_distribution()
        if intensity_dist is not None:
            print(f"Intensity distribution analysis completed")

        # Test iterating through files (legacy compatibility)
        feature_counts = {}
        total_features = 0
        batch_count = 0

        print("\nProcessing msstats data in batches using enhanced methods...")

        for msstats_batch in indexer.iter_msstats_files(file_batch_size=5):
            if msstats_batch is not None and not msstats_batch.empty:
                batch_count += 1
                batch_size = len(msstats_batch)
                total_features += batch_size

                print(f"  Batch {batch_count}: {batch_size:,} features")

                # Count features per reference file in this batch
                if "reference_file_name" in msstats_batch.columns:
                    ref_counts = msstats_batch["reference_file_name"].value_counts()
                    for ref_file, count in ref_counts.items():
                        if ref_file not in feature_counts:
                            feature_counts[ref_file] = 0
                        feature_counts[ref_file] += count

        # Display results
        print(f"\nTotal features processed: {total_features:,}")
        print(f"Total batches: {batch_count}")
        print(f"Unique reference files: {len(feature_counts)}")

        print("\nFeatures per reference file:")
        print("-" * 50)
        for ref_file, count in sorted(feature_counts.items()):
            print(f"  {ref_file}: {count:,} features")

        # Cleanup
        indexer.destroy_database()

        # Assertions
        assert len(feature_counts) > 0, "Should have at least one reference file"
        assert total_features > 0, "Should have at least one feature"
        assert experiment_type == "LFQ", "Should detect LFQ experiment type"

        print(f"\nLFQ test completed successfully!")

    finally:
        # Clean up temporary files
        import os

        # Clean up mztab temporary files
        if "indexer" in locals():
            _cleanup_mztab_temp_files(indexer)

        if os.path.exists(temp_db_path):
            os.unlink(temp_db_path)


def _setup_tmt_test_data():
    """Setup test data for TMT dataset."""
    msstats_gz_file = (
        TEST_DATA_ROOT
        / "quantms/dda-plex-full/PXD007683TMT.sdrf_openms_design_msstats_in.csv.gz"
    )
    sdrf_file = TEST_DATA_ROOT / "quantms/dda-plex-full/PXD007683-TMT.sdrf.tsv"

    return msstats_gz_file, sdrf_file


def _initialize_mztab_indexer(msstats_gz_file, sdrf_file):
    """Initialize MzTabIndexer object for testing."""
    temp_db_path = create_uuid_filename("test_msstats_tmt", ".duckdb")

    indexer = MzTabIndexer(
        mztab_path=None,
        database_path=temp_db_path,
        sdrf_path=sdrf_file,
        max_memory="4GB",
        worker_threads=2,
    )

    # Add MSstats data to the indexer
    indexer.add_msstats_table(str(msstats_gz_file))

    experiment_type = indexer.get_msstats_experiment_type()
    print(f"Experiment type: {experiment_type}")
    print(f"Processing file: {msstats_gz_file.name}")
    print(f"SDRF file: {sdrf_file.name}")

    return indexer, temp_db_path


def _debug_channel_data(indexer):
    """Debug channel data in the dataset."""
    print("\nChecking raw channel values in data...")
    try:
        raw_channels = indexer._duckdb.sql(
            "SELECT DISTINCT Channel FROM msstats ORDER BY Channel"
        ).df()
        print(f"Raw channels in data: {sorted(raw_channels['Channel'].tolist())}")

        print("\nSample raw data before transformation...")
        sample_raw = indexer._duckdb.sql(
            """
            SELECT Channel, Reference, ProteinName, PeptideSequence, Intensity 
            FROM msstats 
            WHERE Reference LIKE '%a05058%' 
            AND ProteinName LIKE '%P55011%' 
            LIMIT 15
        """
        ).df()
        print(sample_raw.to_string(index=False))
    except Exception as e:
        print(f"Debug channel data failed: {e}")
        # Try enhanced methods instead
        exp_type = indexer.get_msstats_experiment_type()
        print(f"Using enhanced analysis - detected experiment type: {exp_type}")


def _debug_batch_data(msstats_batch, batch_count):
    """Debug batch data for the first batch."""
    if batch_count == 1:
        # Handle both old and new column names
        channel_column = "Channel" if "Channel" in msstats_batch.columns else "channel"
        protein_column = (
            "ProteinName" if "ProteinName" in msstats_batch.columns else "pg_accessions"
        )
        ref_column = (
            "Reference_Name"
            if "Reference_Name" in msstats_batch.columns
            else "reference_file_name"
        )
        intensity_column = (
            "Intensity" if "Intensity" in msstats_batch.columns else "intensity"
        )

        if channel_column in msstats_batch.columns:
            unique_channels = sorted(msstats_batch[channel_column].unique())
            print(f"  Channels in batch 1: {unique_channels}")

        if protein_column in msstats_batch.columns:
            debug_protein = msstats_batch[
                msstats_batch[protein_column].str.contains("P55011", na=False)
            ]
            if len(debug_protein) > 0:
                if channel_column in debug_protein.columns:
                    print(
                        f"  Debug protein P55011 channels: {sorted(debug_protein[channel_column].unique())}"
                    )
                print(f"  Debug protein P55011 sample: {len(debug_protein)} rows")

                # Select available columns for display
                display_columns = []
                for col in [
                    ref_column,
                    channel_column,
                    intensity_column,
                    "PeptideSequence",
                ]:
                    if col in debug_protein.columns:
                        display_columns.append(col)

                if display_columns:
                    sample_debug = debug_protein[display_columns].head(3)
                    print("  Sample rows:")
                    print(sample_debug.to_string(index=False, max_colwidth=20))

                if "intensities" in debug_protein.columns:
                    first_intensity_list = debug_protein.iloc[0]["intensities"]
                    print(
                        f"  First peptide intensity list length: {len(first_intensity_list)}"
                    )
                    print("  Sample intensities:")
                    for i, intensity_info in enumerate(first_intensity_list[:3]):
                        print(f"    {i}: {intensity_info}")


def _process_batch_statistics(
    msstats_batch, feature_counts, channel_counts, file_channel_matrix
):
    """Process statistics for a single batch."""
    # Count features per reference file in this batch
    # Handle both old and new column names
    ref_column = (
        "Reference_Name"
        if "Reference_Name" in msstats_batch.columns
        else "reference_file_name"
    )
    if ref_column in msstats_batch.columns:
        ref_counts = msstats_batch[ref_column].value_counts()
        for ref_file, count in ref_counts.items():
            if ref_file not in feature_counts:
                feature_counts[ref_file] = 0
            feature_counts[ref_file] += count

    # Count features per channel (for TMT)
    channel_column = "Channel" if "Channel" in msstats_batch.columns else "channel"
    if channel_column in msstats_batch.columns:
        _process_channel_statistics(
            msstats_batch, channel_counts, file_channel_matrix, channel_column
        )


def _process_channel_column_statistics(msstats_batch, channel_counts, channel_column):
    """Process statistics from the channel column."""
    chan_counts = msstats_batch[channel_column].value_counts()
    for channel, count in chan_counts.items():
        if channel not in channel_counts:
            channel_counts[channel] = 0
        channel_counts[channel] += count


def _process_intensities_statistics(
    row, intensities_channel_counts, file_channel_matrix_from_intensities
):
    """Process statistics from a single row's intensities column."""
    ref_file = row["reference_file_name"]
    if "intensities" not in row or row["intensities"] is None:
        return

    for intensity_info in row["intensities"]:
        channel = intensity_info["channel"]

        # Count channels from intensities
        if channel not in intensities_channel_counts:
            intensities_channel_counts[channel] = 0
        intensities_channel_counts[channel] += 1

        # Count file-channel combinations from intensities
        if ref_file not in file_channel_matrix_from_intensities:
            file_channel_matrix_from_intensities[ref_file] = {}
        if channel not in file_channel_matrix_from_intensities[ref_file]:
            file_channel_matrix_from_intensities[ref_file][channel] = 0
        file_channel_matrix_from_intensities[ref_file][channel] += 1


def _update_file_channel_matrix(
    file_channel_matrix, file_channel_matrix_from_intensities
):
    """Update the main file-channel matrix with intensities data."""
    for ref_file, channels in file_channel_matrix_from_intensities.items():
        if ref_file not in file_channel_matrix:
            file_channel_matrix[ref_file] = {}
        for channel, count in channels.items():
            if channel not in file_channel_matrix[ref_file]:
                file_channel_matrix[ref_file][channel] = 0
            file_channel_matrix[ref_file][channel] += count


def _initialize_missing_channels(channel_counts, intensities_channel_counts):
    """Initialize any missing channels in the main channel counts."""
    for channel in intensities_channel_counts.keys():
        if channel not in channel_counts:
            channel_counts[channel] = 0


def _process_channel_statistics(
    msstats_batch, channel_counts, file_channel_matrix, channel_column="Channel"
):
    """Process channel statistics from batch data."""
    # Count channels from the specified channel column
    _process_channel_column_statistics(msstats_batch, channel_counts, channel_column)

    # Count channels from the 'intensities' column
    intensities_channel_counts = {}
    file_channel_matrix_from_intensities = {}

    for _, row in msstats_batch.iterrows():
        _process_intensities_statistics(
            row, intensities_channel_counts, file_channel_matrix_from_intensities
        )

    # Update the main counting with intensities data
    if intensities_channel_counts:
        _initialize_missing_channels(channel_counts, intensities_channel_counts)
        _update_file_channel_matrix(
            file_channel_matrix, file_channel_matrix_from_intensities
        )


def _display_results(
    feature_counts, channel_counts, file_channel_matrix, total_features, batch_count
):
    """Display test results."""
    print(f"\nTotal features processed: {total_features:,}")
    print(f"Total batches: {batch_count}")
    print(f"Unique reference files: {len(feature_counts)}")
    print(f"Unique channels: {len(channel_counts)}")

    print("\nFeatures per reference file:")
    print("-" * 50)
    for ref_file, count in sorted(feature_counts.items()):
        print(f"  {ref_file}: {count:,} features")

    if channel_counts:
        print("\nFeatures per TMT channel:")
        print("-" * 40)
        for channel, count in sorted(channel_counts.items()):
            print(f"  {channel}: {count:,} features")

    _display_file_channel_matrix(file_channel_matrix)


def _display_file_channel_matrix(file_channel_matrix):
    """Display file x channel matrix."""
    if not file_channel_matrix:
        return

    print("\nFeatures by File × Channel Matrix:")
    print("=" * 80)

    # Get all unique channels and sort them
    all_channels = sorted(
        set(
            channel
            for file_channels in file_channel_matrix.values()
            for channel in file_channels.keys()
        )
    )

    # Print header
    header = f"{'File':<15}"
    for channel in all_channels:
        header += f"{channel:>12}"
    header += f"{'Total':>12}"
    print(header)
    print("-" * len(header))

    # Print data for each file
    for ref_file in sorted(file_channel_matrix.keys()):
        row = f"{ref_file:<15}"
        file_total = 0
        for channel in all_channels:
            count = file_channel_matrix[ref_file].get(channel, 0)
            row += f"{count:>12,}"
            file_total += count
        row += f"{file_total:>12,}"
        print(row)

    # Print channel totals
    totals_row = f"{'Total':<15}"
    grand_total = 0
    for channel in all_channels:
        channel_total = sum(
            file_channels.get(channel, 0)
            for file_channels in file_channel_matrix.values()
        )
        totals_row += f"{channel_total:>12,}"
        grand_total += channel_total
    totals_row += f"{grand_total:>12,}"
    print("-" * len(header))
    print(totals_row)


def test_mztab_indexer_msstats_tmt_full_dataset():
    """Test MzTabIndexer MSstats processing with full TMT dataset."""
    print("\n" + "=" * 60)
    print("Testing MzTabIndexer MSstats analysis with full TMT dataset")
    print("=" * 60)

    # Setup test data
    msstats_gz_file, sdrf_file = _setup_tmt_test_data()

    # Initialize MzTabIndexer
    indexer, temp_db_path = _initialize_mztab_indexer(msstats_gz_file, sdrf_file)

    # Debug channel data (using legacy DuckDB interface)
    _debug_channel_data(indexer)

    # Initialize counters
    feature_counts = {}
    channel_counts = {}
    file_channel_matrix = {}
    total_features = 0
    batch_count = 0

    print("\nProcessing msstats data in batches...")

    # Process batches
    for msstats_batch in indexer.iter_msstats_files(file_batch_size=5):
        if msstats_batch is not None and not msstats_batch.empty:
            batch_count += 1
            batch_size = len(msstats_batch)
            total_features += batch_size

            print(f"  Batch {batch_count}: {batch_size:,} features")

            # Debug batch data
            _debug_batch_data(msstats_batch, batch_count)

            # Process batch statistics
            _process_batch_statistics(
                msstats_batch, feature_counts, channel_counts, file_channel_matrix
            )

    # Display results
    _display_results(
        feature_counts, channel_counts, file_channel_matrix, total_features, batch_count
    )

    # Test enhanced analysis methods
    print("\nTesting enhanced MSstats analysis methods...")
    experiment_type = indexer.get_msstats_experiment_type()
    file_stats = indexer.get_msstats_file_statistics()
    protein_summary = indexer.get_msstats_protein_summary()

    if file_stats is not None:
        print(f"Enhanced file statistics: {len(file_stats)} files")
    if protein_summary is not None:
        print(f"Enhanced protein summary: {len(protein_summary)} proteins")

    # Cleanup
    _cleanup_mztab_temp_files(indexer)
    indexer.destroy_database()
    import os

    if os.path.exists(temp_db_path):
        os.unlink(temp_db_path)

    # Assertions
    assert len(feature_counts) > 0, "Should have at least one reference file"
    assert total_features > 0, "Should have at least one feature"
    assert experiment_type.startswith("TMT"), "Should detect TMT experiment type"
    assert len(channel_counts) > 0, "Should have TMT channels"

    print(f"\nTMT test completed successfully!")


def test_mztab_indexer_msstats_comparison():
    """Compare LFQ and TMT datasets side by side."""
    print("\n" + "=" * 60)
    print("Comparison of LFQ vs TMT datasets")
    print("=" * 60)

    datasets = [
        {
            "name": "LFQ (PXD007683)",
            "msstats_gz": TEST_DATA_ROOT
            / "quantms/dda-lfq-full/PXD007683-LFQ.sdrf_openms_design_msstats_in.csv.gz",
            "sdrf": TEST_DATA_ROOT / "quantms/dda-lfq-full/PXD007683-LFQ.sdrf.tsv",
        },
        {
            "name": "TMT (PXD007683)",
            "msstats_gz": TEST_DATA_ROOT
            / "quantms/dda-plex-full/PXD007683TMT.sdrf_openms_design_msstats_in.csv.gz",
            "sdrf": TEST_DATA_ROOT / "quantms/dda-plex-full/PXD007683-TMT.sdrf.tsv",
        },
    ]

    results = {}

    for dataset in datasets:
        print(f"\nAnalyzing {dataset['name']}...")

        # Create temporary MzTabIndexer database
        temp_db_path = create_uuid_filename(f"{dataset['name']}_msstats", ".duckdb")

        try:
            # Create MzTabIndexer with MSstats data
            indexer = MzTabIndexer(
                mztab_path=None,  # No mzTab file needed
                database_path=temp_db_path,
                sdrf_path=dataset["sdrf"],
                max_memory="4GB",
                worker_threads=2,
            )

            # Add MSstats data to the indexer
            indexer.add_msstats_table(str(dataset["msstats_gz"]))

            # Test enhanced MSstats analysis methods
            experiment_type = indexer.get_msstats_experiment_type()
            print(f"Experiment type: {experiment_type}")
            print(f"Processing file: {dataset['msstats_gz'].name}")
            print(f"SDRF file: {dataset['sdrf'].name}")

            # Test file statistics
            file_stats = indexer.get_msstats_file_statistics()
            if file_stats is not None:
                print(f"\nFile statistics: {len(file_stats)} files")
                print(f"Total features: {file_stats['feature_count'].sum():,}")
                print(f"Total proteins: {file_stats['protein_count'].sum():,}")

            # Test protein summary
            protein_summary = indexer.get_msstats_protein_summary()
            if protein_summary is not None:
                print(f"Protein summary: {len(protein_summary)} proteins")

            # Test peptide summary
            peptide_summary = indexer.get_msstats_peptide_summary()
            if peptide_summary is not None:
                print(f"Peptide summary: {len(peptide_summary)} peptides")

            # Test intensity distribution
            intensity_dist = indexer.get_msstats_intensity_distribution()
            if intensity_dist is not None:
                print(f"Intensity distribution analysis completed")

            # Test iterating through files (legacy compatibility)
            feature_counts = {}
            total_features = 0
            batch_count = 0

            print("\nProcessing msstats data in batches using enhanced methods...")

            for msstats_batch in indexer.iter_msstats_files(file_batch_size=5):
                if msstats_batch is not None and not msstats_batch.empty:
                    batch_count += 1
                    batch_size = len(msstats_batch)
                    total_features += batch_size

                    print(f"  Batch {batch_count}: {batch_size:,} features")

                    # Count features per reference file in this batch
                    if "Reference_Name" in msstats_batch.columns:
                        ref_counts = msstats_batch["Reference_Name"].value_counts()
                        for ref_file, count in ref_counts.items():
                            if ref_file not in feature_counts:
                                feature_counts[ref_file] = 0
                            feature_counts[ref_file] += count

            # Display results
            print(f"\nTotal features processed: {total_features:,}")
            print(f"Total batches: {batch_count}")
            print(f"Unique reference files: {len(feature_counts)}")

            print("\nFeatures per reference file:")
            print("-" * 50)
            for ref_file, count in sorted(feature_counts.items()):
                print(f"  {ref_file}: {count:,} features")

            # Cleanup
            indexer.destroy_database()

            results[dataset["name"]] = {
                "experiment_type": experiment_type,
                "total_features": total_features,
                "reference_files": len(feature_counts),
                "file_size": dataset["msstats_gz"].stat().st_size / (1024 * 1024),  # MB
            }

        finally:
            # Clean up temporary files
            import os

            # Clean up mztab temporary files
            if "indexer" in locals():
                _cleanup_mztab_temp_files(indexer)

            if os.path.exists(temp_db_path):
                os.unlink(temp_db_path)

    # Display comparison
    print(f"\nDataset Comparison:")
    print("-" * 70)
    print(
        f"{'Dataset':<20} {'Type':<8} {'Features':<12} {'Files':<6} {'Size (MB)':<10}"
    )
    print("-" * 70)

    for name, stats in results.items():
        print(
            f"{name:<20} {stats['experiment_type']:<8} {stats['total_features']:<12,} "
            f"{stats['reference_files']:<6} {stats['file_size']:<10.1f}"
        )

    print("\nComparison completed!")


if __name__ == "__main__":
    # Run individual tests
    test_mztab_indexer_msstats_lfq_full_dataset()
    test_mztab_indexer_msstats_tmt_full_dataset()
    test_mztab_indexer_msstats_comparison()
