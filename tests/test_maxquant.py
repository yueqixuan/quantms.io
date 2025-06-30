import os
import tempfile
from pathlib import Path

import pandas as pd
import pytest

from quantmsio.core.maxquant.maxquant import MaxQuant
from quantmsio.core.quantms.feature import Feature
from quantmsio.core.quantms.pg import MzTabProteinGroups
from quantmsio.core.quantms.psm import Psm
from quantmsio.core.sdrf import SDRFHandler

TEST_DATA_ROOT = Path(__file__).parent / "examples"
test_data = (
    TEST_DATA_ROOT / "maxquant/maxquant_simple/msms.txt",
    TEST_DATA_ROOT / "maxquant/maxquant_simple/evidence.txt",
    TEST_DATA_ROOT / "maxquant/maxquant_simple/sdrf.tsv",
)


def test_transform_psm():
    """Test transforming MaxQuant PSM data."""
    # Resolve file path
    msms_file = test_data[0]

    # Initialize MaxQuant
    maxquant = MaxQuant()

    # Process PSM data
    count = 0
    for df in maxquant.iter_batch(str(msms_file), "psm", chunksize=500000):
        # Transform PSM data
        maxquant.transform_psm(df)
        Psm.convert_to_parquet_format(df)
        psm_parquet = Psm.transform_parquet(df)

        # Add assertions to verify the result
        assert psm_parquet is not None
        assert len(psm_parquet) > 0
        assert "peptidoform" in psm_parquet.column_names
        count += 1

    # Ensure we got at least one result
    assert count > 0


def test_transform_feature():
    """Test transforming MaxQuant feature data."""
    # Resolve file paths
    evidence_file = test_data[1]
    sdrf_file = test_data[2]

    # Initialize SDRF handler and MaxQuant
    sdrf = SDRFHandler(sdrf_file)
    maxquant = MaxQuant()
    maxquant.experiment_type = sdrf.get_experiment_type_from_sdrf()
    maxquant._sample_map = sdrf.get_sample_map_run()

    # Process feature data
    count = 0
    for df in maxquant.iter_batch(evidence_file, chunksize=500000):
        # Transform feature data
        maxquant.transform_feature(df)
        Feature.convert_to_parquet_format(df)
        feature = Feature.transform_feature(df)

        # Add assertions to verify the result
        assert feature is not None
        assert len(feature) > 0
        assert "peptidoform" in feature.column_names
        count += 1

    # Ensure we got at least one result
    assert count > 0


def test_transform_features():
    """Test transforming MaxQuant features with slicing."""
    # Resolve file paths
    evidence_file = test_data[1]
    sdrf_file = test_data[2]

    # Initialize SDRF handler and MaxQuant
    sdrf = SDRFHandler(sdrf_file)
    maxquant = MaxQuant()
    maxquant.experiment_type = sdrf.get_experiment_type_from_sdrf()
    maxquant._sample_map = sdrf.get_sample_map_run()

    # Process feature data
    for report in maxquant.iter_batch(evidence_file, chunksize=500000):
        # Transform feature data
        maxquant.transform_feature(report)
        Feature.convert_to_parquet_format(report)

        # Test slicing
        slice_count = 0
        for _, df in Feature.slice(report, ["reference_file_name", "precursor_charge"]):
            feature = Feature.transform_feature(df)

            # Add assertions to verify the result
            assert feature is not None
            assert len(feature) > 0
            assert "peptidoform" in feature.column_names
            slice_count += 1

        # Ensure we got at least one slice
        assert slice_count > 0


def test_maxquant_protein_groups_transform():
    """Test transforming MaxQuant protein group data."""
    # Create sample MaxQuant protein groups data
    sample_data = {
        "Protein IDs": ["P12345;Q67890", "P11111", "P22222"],
        "Protein names": ["Protein A;Protein B", "Protein C", "Protein D"],
        "Gene names": ["GENEA;GENEB", "GENEC", "GENED"],
        "Q-value": [0.001, 0.01, 0.05],
        "Intensity": [1000.0, 2000.0, 3000.0],
        "LFQ intensity": [800.0, 1800.0, 2800.0],
        "iBAQ": [500.0, 1500.0, 2500.0],
        "Razor + unique peptides": [5, 3, 7],
        "Unique peptides": [3, 2, 5],
        "Sequence coverage [%]": [15.5, 25.0, 30.0],
        "Mol. weight [kDa]": [50.5, 75.2, 100.1],
        "Score": [85.2, 92.1, 78.3],
        "Reverse": ["", "+", ""],
        "Potential contaminant": ["", "", "+"],
        "MS/MS count": [10, 15, 12],
    }

    df = pd.DataFrame(sample_data)

    # Initialize MaxQuant converter
    mq = MaxQuant()

    # Transform the data
    mq.transform_protein_groups(df)

    # Test that basic transformations work
    assert "pg_accessions" in df.columns
    assert "pg_names" in df.columns
    assert "gg_accessions" in df.columns
    assert "intensities" in df.columns
    assert "additional_intensities" in df.columns
    assert "peptides" in df.columns
    assert "anchor_protein" in df.columns
    assert "additional_scores" in df.columns

    # Test protein accessions are properly split
    assert isinstance(df.iloc[0]["pg_accessions"], list)
    assert len(df.iloc[0]["pg_accessions"]) == 2
    assert df.iloc[0]["pg_accessions"] == ["P12345", "Q67890"]

    # Test anchor protein is set correctly
    assert df.iloc[0]["anchor_protein"] == "P12345"
    assert df.iloc[1]["anchor_protein"] == "P11111"

    # Test decoy and contaminant flags
    assert df.iloc[0]["is_decoy"] == 0
    assert df.iloc[1]["is_decoy"] == 1  # Had "+" in Reverse
    assert df.iloc[0]["contaminant"] == 0
    assert df.iloc[2]["contaminant"] == 1  # Had "+" in Potential contaminant

    # Test intensities structure
    intensities = df.iloc[0]["intensities"]
    assert isinstance(intensities, list)
    if intensities:  # May be empty if no sample columns
        assert "sample_accession" in intensities[0]
        assert "channel" in intensities[0]
        assert "intensity" in intensities[0]

    # Test peptide count structure
    peptides = df.iloc[0]["peptides"]
    assert isinstance(peptides, list)
    assert len(peptides) == 2  # Should match number of proteins in group
    assert peptides[0]["protein_name"] == "P12345"
    assert peptides[0]["peptide_count"] == 5  # From Razor + unique peptides


def test_maxquant_protein_groups_with_sample_columns():
    """Test MaxQuant protein groups with sample-specific intensity columns."""
    # Create sample data with sample-specific columns
    sample_data = {
        "Protein IDs": ["P12345"],
        "Protein names": ["Protein A"],
        "Gene names": ["GENEA"],
        "Q-value": [0.001],
        "Intensity Sample1": [1000.0],
        "Intensity Sample2": [1200.0],
        "LFQ intensity Sample1": [800.0],
        "LFQ intensity Sample2": [900.0],
        "iBAQ Sample1": [500.0],
        "iBAQ Sample2": [600.0],
        "Razor + unique peptides": [5],
        "Reverse": [""],
        "Potential contaminant": [""],
    }

    df = pd.DataFrame(sample_data)

    # Store sample columns using proper DataFrame attribute setting
    df.attrs["sample_intensity_cols"] = ["Intensity Sample1", "Intensity Sample2"]
    df.attrs["sample_lfq_cols"] = ["LFQ intensity Sample1", "LFQ intensity Sample2"]
    df.attrs["sample_ibaq_cols"] = ["iBAQ Sample1", "iBAQ Sample2"]

    mq = MaxQuant()
    mq.transform_protein_groups(df)

    # Test intensities from sample columns
    intensities = df.iloc[0]["intensities"]
    assert len(intensities) == 2
    assert intensities[0]["sample_accession"] == "Sample1"
    assert intensities[0]["intensity"] == 1000.0
    assert intensities[1]["sample_accession"] == "Sample2"
    assert intensities[1]["intensity"] == 1200.0

    # Test additional intensities
    additional_intensities = df.iloc[0]["additional_intensities"]
    assert len(additional_intensities) >= 2  # LFQ and iBAQ for both samples


def test_maxquant_pg_basic_transformation():
    """Test basic MaxQuant protein groups transformation."""
    # Create sample MaxQuant protein groups data
    sample_data = {
        "Protein IDs": ["P12345;Q67890", "P11111", "P22222"],
        "Protein names": ["Protein A;Protein B", "Protein C", "Protein D"],
        "Gene names": ["GENEA;GENEB", "GENEC", "GENED"],
        "Q-value": [0.001, 0.01, 0.05],
        "Intensity": [1000.0, 2000.0, 3000.0],
        "LFQ intensity": [800.0, 1800.0, 2800.0],
        "iBAQ": [500.0, 1500.0, 2500.0],
        "Razor + unique peptides": [5, 3, 7],
        "Unique peptides": [3, 2, 5],
        "Sequence coverage [%]": [15.5, 25.0, 30.0],
        "Mol. weight [kDa]": [50.5, 75.2, 100.1],
        "Score": [85.2, 92.1, 78.3],
        "Reverse": ["", "+", ""],
        "Potential contaminant": ["", "", "+"],
        "MS/MS count": [10, 15, 12],
    }

    df = pd.DataFrame(sample_data)
    mq = MaxQuant()
    mq.transform_protein_groups(df)

    # Test that basic transformations work
    assert "pg_accessions" in df.columns
    assert "pg_names" in df.columns
    assert "gg_accessions" in df.columns
    assert "intensities" in df.columns
    assert "additional_intensities" in df.columns
    assert "peptides" in df.columns
    assert "anchor_protein" in df.columns
    assert "additional_scores" in df.columns

    # Test protein accessions are properly split
    assert isinstance(df.iloc[0]["pg_accessions"], list)
    assert len(df.iloc[0]["pg_accessions"]) == 2
    assert df.iloc[0]["pg_accessions"] == ["P12345", "Q67890"]

    # Test anchor protein is set correctly
    assert df.iloc[0]["anchor_protein"] == "P12345"
    assert df.iloc[1]["anchor_protein"] == "P11111"

    # Test decoy and contaminant flags
    assert df.iloc[0]["is_decoy"] == 0
    assert df.iloc[1]["is_decoy"] == 1  # Had "+" in Reverse
    assert df.iloc[0]["contaminant"] == 0
    assert df.iloc[2]["contaminant"] == 1  # Had "+" in Potential contaminant

    # Test intensities structure
    intensities = df.iloc[0]["intensities"]
    assert isinstance(intensities, list)
    if intensities:  # May be empty if no sample columns
        assert "sample_accession" in intensities[0]
        assert "channel" in intensities[0]
        assert "intensity" in intensities[0]

    # Test peptide count structure
    peptides = df.iloc[0]["peptides"]
    assert isinstance(peptides, list)
    assert len(peptides) == 2  # Should match number of proteins in group
    assert peptides[0]["protein_name"] == "P12345"
    assert peptides[0]["peptide_count"] == 5  # From Razor + unique peptides


def test_maxquant_pg_sample_specific_columns():
    """Test MaxQuant protein groups with sample-specific columns."""
    # Create sample data with sample-specific columns
    sample_data = {
        "Protein IDs": ["P12345"],
        "Protein names": ["Protein A"],
        "Gene names": ["GENEA"],
        "Q-value": [0.001],
        "Intensity Sample1": [1000.0],
        "Intensity Sample2": [1200.0],
        "LFQ intensity Sample1": [800.0],
        "LFQ intensity Sample2": [900.0],
        "iBAQ Sample1": [500.0],
        "iBAQ Sample2": [600.0],
        "Razor + unique peptides": [5],
        "Reverse": [""],
        "Potential contaminant": [""],
    }

    df = pd.DataFrame(sample_data)

    # Store sample columns using proper DataFrame attribute setting
    df.attrs["sample_intensity_cols"] = ["Intensity Sample1", "Intensity Sample2"]
    df.attrs["sample_lfq_cols"] = ["LFQ intensity Sample1", "LFQ intensity Sample2"]
    df.attrs["sample_ibaq_cols"] = ["iBAQ Sample1", "iBAQ Sample2"]

    mq = MaxQuant()
    mq.transform_protein_groups(df)

    # Test intensities from sample columns
    intensities = df.iloc[0]["intensities"]
    assert len(intensities) == 2
    assert intensities[0]["sample_accession"] == "Sample1"
    assert intensities[0]["intensity"] == 1000.0
    assert intensities[1]["sample_accession"] == "Sample2"
    assert intensities[1]["intensity"] == 1200.0

    # Test additional intensities
    additional_intensities = df.iloc[0]["additional_intensities"]
    assert len(additional_intensities) >= 2  # LFQ and iBAQ for both samples


def test_mztab_pg_gene_extraction():
    """Test gene name extraction from protein descriptions."""
    with tempfile.NamedTemporaryFile(suffix=".mzTab", mode="w", delete=False) as tmp:
        tmp.write("MTD\tversion\t1.0.0\n")
        tmp_path = tmp.name

    try:
        mztab_pg = MzTabProteinGroups(tmp_path)

        # Test various description formats
        test_cases = [
            ("Protein A OS=Homo sapiens GN=GENEA PE=1 SV=1", ["GENEA"]),
            ("No gene name here OS=Homo sapiens PE=1 SV=1", []),
            ("", []),
            ("Multiple GN=GENE1 and GN=GENE2", ["GENE1"]),  # Should extract first
        ]

        for description, expected in test_cases:
            result = mztab_pg._extract_gene_names(description)
            assert result == expected, f"Failed for description: {description}"

    finally:
        # Clean up MzTabProteinGroups resources
        if "mztab_pg" in locals():
            mztab_pg.cleanup()
        # Clean up temporary file
        os.unlink(tmp_path)


if __name__ == "__main__":
    test_maxquant_protein_groups_transform()
    test_maxquant_protein_groups_with_sample_columns()
    test_maxquant_pg_basic_transformation()
    test_maxquant_pg_sample_specific_columns()
    test_mztab_pg_gene_extraction()
    print("All MaxQuant protein group tests passed!")
