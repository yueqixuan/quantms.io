from pathlib import Path

from quantmsio.core.diann.diann import DiaNNConvert
from quantmsio.core.quantms.feature import Feature

TEST_DATA_ROOT = Path(__file__).parent / "examples"

TEST_DATA = (
    TEST_DATA_ROOT / "DIANN/diann_report.tsv",
    TEST_DATA_ROOT / "DIANN/PXD019909-DIA.sdrf.tsv",
    TEST_DATA_ROOT / "DIANN/mzml",
)


def test_transform_feature():
    report_file = TEST_DATA[0]
    sdrf_file = TEST_DATA[1]
    mzml = TEST_DATA[2]
    diann_converter = DiaNNConvert(report_file, sdrf_file)
    try:
        for report in diann_converter.main_report_df(0.05, mzml, 2):
            diann_converter.add_additional_msg(report)
            Feature.convert_to_parquet_format(report)
            Feature.transform_feature(report)
    finally:
        # Clean up DuckDB database
        diann_converter.destroy_duckdb_database()


def test_transform_features():
    report_file = TEST_DATA[0]
    sdrf_file = TEST_DATA[1]
    mzml = TEST_DATA[2]
    diann_converter = DiaNNConvert(report_file, sdrf_file)
    try:
        for report in diann_converter.main_report_df(0.05, mzml, 2):
            diann_converter.add_additional_msg(report)
            Feature.convert_to_parquet_format(report)
            for _, df in Feature.slice(
                report, ["reference_file_name", "precursor_charge"]
            ):
                Feature.transform_feature(df)
    finally:
        # Clean up DuckDB database
        diann_converter.destroy_duckdb_database()


def test_transform_protein_groups():
    """Test transforming DIA-NN protein group data."""
    report_file = TEST_DATA[0]
    sdrf_file = TEST_DATA[1]
    diann_converter = DiaNNConvert(report_file, sdrf_file)

    try:
        # Get some test data for protein groups using the proper SQL format
        refs = diann_converter.get_unique_references("Run")[:1]  # Just test with 1 file
        from quantmsio.core.common import DIANN_PG_MAP
        from quantmsio.core.diann.diann import DIANN_PG_SQL

        report = diann_converter.get_report_from_database(refs, DIANN_PG_SQL)

        # Apply the DIA-NN protein group mapping
        report.rename(columns=DIANN_PG_MAP, inplace=True)
        report.dropna(subset=["pg_accessions"], inplace=True)

        if len(report) > 0:  # Only test if we have data
            # Test one file worth of data
            df = report.head(5).copy()  # Test with 5 protein groups

            if len(df) > 0:
                # Transform the protein group data
                pg_df = diann_converter.generate_pg_matrix(df)

                # Verify the structure
                assert (
                    "intensities" in pg_df.columns
                ), "intensities column should be present"
                assert (
                    "additional_intensities" in pg_df.columns
                ), "additional_intensities column should be present"

                # Check intensities structure
                intensity_sample = pg_df["intensities"].iloc[0]
                assert isinstance(
                    intensity_sample, list
                ), "intensities should be a list"
                assert len(intensity_sample) > 0, "intensities should not be empty"

                intensity_entry = intensity_sample[0]
                assert (
                    "sample_accession" in intensity_entry
                ), "intensity entry should have sample_accession"
                assert (
                    "channel" in intensity_entry
                ), "intensity entry should have channel"
                assert (
                    "intensity" in intensity_entry
                ), "intensity entry should have intensity value"
                assert (
                    intensity_entry["channel"] == "LFQ"
                ), "channel should be LFQ for DIA-NN"

                # Check additional_intensities structure
                additional_intensity_sample = pg_df["additional_intensities"].iloc[0]
                assert isinstance(
                    additional_intensity_sample, list
                ), "additional_intensities should be a list"
                assert (
                    len(additional_intensity_sample) > 0
                ), "additional_intensities should not be empty"

                additional_entry = additional_intensity_sample[0]
                assert (
                    "sample_accession" in additional_entry
                ), "additional_intensity entry should have sample_accession"
                assert (
                    "channel" in additional_entry
                ), "additional_intensity entry should have channel"
                assert (
                    "intensities" in additional_entry
                ), "additional_intensity entry should have intensities array"
                assert (
                    additional_entry["channel"] == "LFQ"
                ), "channel should be LFQ for DIA-NN"

                # Check the additional intensity array structure
                additional_types = additional_entry["intensities"]
                assert isinstance(
                    additional_types, list
                ), "intensities should be a list"
                assert (
                    len(additional_types) == 2
                ), "should have 2 additional intensity types (normalize_intensity and lfq)"

                # Check normalize_intensity and lfq are present
                intensity_names = [item["intensity_name"] for item in additional_types]
                assert (
                    "normalize_intensity" in intensity_names
                ), "normalize_intensity should be present"
                assert "lfq" in intensity_names, "lfq should be present"

                print("Protein group intensity structure test passed!")
    finally:
        # Clean up DuckDB database
        diann_converter.destroy_duckdb_database()
