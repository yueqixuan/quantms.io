"""
Tests for MsstatsIN module using full datasets.

This module tests the MsstatsIN functionality with real, full-sized datasets
for both LFQ and TMT experiments, demonstrating feature counting per reference file.
"""

import gzip
import tempfile
from pathlib import Path

from quantmsio.core.quantms.msstats_in import MsstatsIN

TEST_DATA_ROOT = Path(__file__).parent / "examples"


def test_msstats_in_lfq_full_dataset():
    """Test MsstatsIN processing with full LFQ dataset."""
    print("\n" + "=" * 60)
    print("Testing MsstatsIN with full LFQ dataset")
    print("=" * 60)

    # File paths
    msstats_gz_file = (
        TEST_DATA_ROOT
        / "quantms/dda-lfq-full/PXD007683-LFQ.sdrf_openms_design_msstats_in.csv.gz"
    )
    sdrf_file = TEST_DATA_ROOT / "quantms/dda-lfq-full/PXD007683-LFQ.sdrf.tsv"

    # Extract the gzipped msstats file to a temporary file
    with tempfile.NamedTemporaryFile(suffix=".csv", mode="w", delete=False) as tmp_file:
        with gzip.open(msstats_gz_file, "rt") as gz_file:
            tmp_file.write(gz_file.read())
        msstats_file = tmp_file.name

    try:
        # Initialize MsstatsIN
        msstats_in = MsstatsIN(
            report_path=msstats_file,
            sdrf_path=sdrf_file,
            duckdb_max_memory="4GB",
            duckdb_threads=2,
        )

        print(f"📊 Experiment type: {msstats_in.experiment_type}")
        print(f"📁 Processing file: {msstats_gz_file.name}")
        print(f"🗂️  SDRF file: {sdrf_file.name}")

        # Count features per reference file
        feature_counts = {}
        total_features = 0
        batch_count = 0

        print("\n🔄 Processing msstats data in batches...")

        for msstats_batch in msstats_in.generate_msstats_in(file_num=5):
            batch_count += 1
            batch_size = len(msstats_batch)
            total_features += batch_size

            print(f"  Batch {batch_count}: {batch_size:,} features")

            # Count features per reference file in this batch
            ref_counts = msstats_batch["reference_file_name"].value_counts()
            for ref_file, count in ref_counts.items():
                if ref_file not in feature_counts:
                    feature_counts[ref_file] = 0
                feature_counts[ref_file] += count

        # Display results
        print(f"\n📈 Total features processed: {total_features:,}")
        print(f"📦 Total batches: {batch_count}")
        print(f"🗃️  Unique reference files: {len(feature_counts)}")

        print("\n📋 Features per reference file:")
        print("-" * 50)
        for ref_file, count in sorted(feature_counts.items()):
            print(f"  {ref_file}: {count:,} features")

        # Cleanup
        msstats_in.destroy_duckdb_database()

        # Assertions
        assert len(feature_counts) > 0, "Should have at least one reference file"
        assert total_features > 0, "Should have at least one feature"
        assert msstats_in.experiment_type == "LFQ", "Should detect LFQ experiment type"

        print(f"\n✅ LFQ test completed successfully!")

    finally:
        # Clean up temporary file
        import os

        if os.path.exists(msstats_file):
            os.unlink(msstats_file)


def test_msstats_in_tmt_full_dataset():
    """Test MsstatsIN processing with full TMT dataset."""
    print("\n" + "=" * 60)
    print("Testing MsstatsIN with full TMT dataset")
    print("=" * 60)

    # File paths - DuckDB can read .gz files directly
    msstats_gz_file = (
        TEST_DATA_ROOT
        / "quantms/dda-plex-full/PXD007683TMT.sdrf_openms_design_msstats_in.csv.gz"
    )
    sdrf_file = TEST_DATA_ROOT / "quantms/dda-plex-full/PXD007683-TMT.sdrf.tsv"

    # Initialize MsstatsIN with the gzipped file directly
    msstats_in = MsstatsIN(
        report_path=msstats_gz_file,  # DuckDB handles .gz files natively
        sdrf_path=sdrf_file,
        duckdb_max_memory="4GB",
        duckdb_threads=2,
    )

    print(f"📊 Experiment type: {msstats_in.experiment_type}")
    print(f"📁 Processing file: {msstats_gz_file.name}")
    print(f"🗂️  SDRF file: {sdrf_file.name}")

    # Debug: Check raw channel values in the data
    print("\n🔍 Checking raw channel values in data...")
    raw_channels = msstats_in._duckdb.sql(
        "SELECT DISTINCT Channel FROM report ORDER BY Channel"
    ).df()
    print(f"Raw channels in data: {sorted(raw_channels['Channel'].tolist())}")

    # Debug: Check a sample of raw data before transformation
    print("\n🔍 Sample raw data before transformation...")
    sample_raw = msstats_in._duckdb.sql(
        """
        SELECT Channel, Reference, ProteinName, PeptideSequence, Charge, Intensity 
        FROM report 
        WHERE Reference LIKE '%a05058%' 
        AND ProteinName LIKE '%P55011%' 
        LIMIT 15
    """
    ).df()
    print(sample_raw.to_string(index=False))

    # Count features per reference file and channel
    feature_counts = {}
    channel_counts = {}
    file_channel_matrix = {}  # New: file x channel matrix
    total_features = 0
    batch_count = 0

    print("\n🔄 Processing msstats data in batches...")

    for msstats_batch in msstats_in.generate_msstats_in(file_num=5):
        batch_count += 1
        batch_size = len(msstats_batch)
        total_features += batch_size

        print(f"  Batch {batch_count}: {batch_size:,} features")

        # Debug: Show unique channels in this batch
        if batch_count == 1:
            unique_channels = sorted(msstats_batch["channel"].unique())
            print(f"  Channels in batch 1: {unique_channels}")

            # Debug: Check what happens to a specific protein across channels
            debug_protein = msstats_batch[
                msstats_batch["pg_accessions"].str.contains("P55011", na=False)
            ]
            if len(debug_protein) > 0:
                print(
                    f"  Debug protein P55011 channels: {sorted(debug_protein['channel'].unique())}"
                )
                print(f"  Debug protein P55011 sample: {len(debug_protein)} rows")
                # Show a few sample rows
                sample_debug = debug_protein[
                    ["reference_file_name", "channel", "intensity", "peptidoform"]
                ].head(3)
                print("  Sample rows:")
                print(sample_debug.to_string(index=False, max_colwidth=20))

                # Debug: Check the intensities column to see if channel aggregation worked
                if "intensities" in debug_protein.columns:
                    first_intensity_list = debug_protein.iloc[0]["intensities"]
                    print(
                        f"  First peptide intensity list length: {len(first_intensity_list)}"
                    )
                    print("  Sample intensities:")
                    for i, intensity_info in enumerate(first_intensity_list[:3]):
                        print(f"    {i}: {intensity_info}")

        # Count features per reference file in this batch
        ref_counts = msstats_batch["reference_file_name"].value_counts()
        for ref_file, count in ref_counts.items():
            if ref_file not in feature_counts:
                feature_counts[ref_file] = 0
            feature_counts[ref_file] += count

        # Count features per channel (for TMT)
        if "channel" in msstats_batch.columns:
            # Count channels from the 'channel' column (which only shows surviving row after deduplication)
            chan_counts = msstats_batch["channel"].value_counts()
            for channel, count in chan_counts.items():
                if channel not in channel_counts:
                    channel_counts[channel] = 0
                channel_counts[channel] += count

            # Also count channels from the 'intensities' column (which should have all channels)
            intensities_channel_counts = {}
            file_channel_matrix_from_intensities = {}

            for _, row in msstats_batch.iterrows():
                ref_file = row["reference_file_name"]
                if "intensities" in row and row["intensities"] is not None:
                    for intensity_info in row["intensities"]:
                        channel = intensity_info["channel"]

                        # Count channels from intensities
                        if channel not in intensities_channel_counts:
                            intensities_channel_counts[channel] = 0
                        intensities_channel_counts[channel] += 1

                        # Count file-channel combinations from intensities
                        if ref_file not in file_channel_matrix_from_intensities:
                            file_channel_matrix_from_intensities[ref_file] = {}
                        if (
                            channel
                            not in file_channel_matrix_from_intensities[ref_file]
                        ):
                            file_channel_matrix_from_intensities[ref_file][channel] = 0
                        file_channel_matrix_from_intensities[ref_file][channel] += 1

            # Update the main counting with intensities data
            if intensities_channel_counts:
                print(
                    f"  Channels from intensities in batch {batch_count}: {sorted(intensities_channel_counts.keys())}"
                )
                for channel, count in intensities_channel_counts.items():
                    if channel not in channel_counts:
                        channel_counts[channel] = 0
                    # Note: not adding to avoid double counting, just showing both methods

                # Use intensities-based matrix instead
                for ref_file, channels in file_channel_matrix_from_intensities.items():
                    if ref_file not in file_channel_matrix:
                        file_channel_matrix[ref_file] = {}
                    for channel, count in channels.items():
                        if channel not in file_channel_matrix[ref_file]:
                            file_channel_matrix[ref_file][channel] = 0
                        file_channel_matrix[ref_file][channel] += count

            # Original counting from channel column
            file_channel_counts = msstats_batch.groupby(
                ["reference_file_name", "channel"]
            ).size()
            # Note: Not using this since it only shows surviving channels after deduplication

    # Display results
    print(f"\n📈 Total features processed: {total_features:,}")
    print(f"📦 Total batches: {batch_count}")
    print(f"🗃️  Unique reference files: {len(feature_counts)}")
    print(f"🏷️  Unique channels: {len(channel_counts)}")

    print("\n📋 Features per reference file:")
    print("-" * 50)
    for ref_file, count in sorted(feature_counts.items()):
        print(f"  {ref_file}: {count:,} features")

    if channel_counts:
        print("\n🏷️  Features per TMT channel:")
        print("-" * 40)
        for channel, count in sorted(channel_counts.items()):
            print(f"  {channel}: {count:,} features")

    # New: Display file x channel matrix
    if file_channel_matrix:
        print("\n📊 Features by File × Channel Matrix:")
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

    # Cleanup
    msstats_in.destroy_duckdb_database()

    # Assertions
    assert len(feature_counts) > 0, "Should have at least one reference file"
    assert total_features > 0, "Should have at least one feature"
    assert msstats_in.experiment_type.startswith(
        "TMT"
    ), "Should detect TMT experiment type"
    assert len(channel_counts) > 0, "Should have TMT channels"

    print(f"\n✅ TMT test completed successfully!")


def test_msstats_in_comparison():
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
        print(f"\n🔍 Analyzing {dataset['name']}...")

        # Extract the gzipped msstats file to a temporary file
        with tempfile.NamedTemporaryFile(
            suffix=".csv", mode="w", delete=False
        ) as tmp_file:
            with gzip.open(dataset["msstats_gz"], "rt") as gz_file:
                tmp_file.write(gz_file.read())
            msstats_file = tmp_file.name

        try:
            # Initialize MsstatsIN
            msstats_in = MsstatsIN(
                report_path=msstats_file,
                sdrf_path=dataset["sdrf"],
                duckdb_max_memory="4GB",
                duckdb_threads=2,
            )

            # Quick stats
            total_features = 0
            reference_files = set()

            for msstats_batch in msstats_in.generate_msstats_in(file_num=10):
                total_features += len(msstats_batch)
                reference_files.update(msstats_batch["reference_file_name"].unique())

            results[dataset["name"]] = {
                "experiment_type": msstats_in.experiment_type,
                "total_features": total_features,
                "reference_files": len(reference_files),
                "file_size": dataset["msstats_gz"].stat().st_size / (1024 * 1024),  # MB
            }

            # Cleanup
            msstats_in.destroy_duckdb_database()

        finally:
            # Clean up temporary file
            import os

            if os.path.exists(msstats_file):
                os.unlink(msstats_file)

    # Display comparison
    print(f"\n📊 Dataset Comparison:")
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

    print("\n✅ Comparison completed!")


if __name__ == "__main__":
    # Run individual tests
    test_msstats_in_lfq_full_dataset()
    test_msstats_in_tmt_full_dataset()
    test_msstats_in_comparison()
