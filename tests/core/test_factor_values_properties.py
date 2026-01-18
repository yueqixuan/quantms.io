"""
Property-based tests for multi-factor value support.

These tests verify the correctness properties defined in the design document
for the factor_values feature that replaces the condition field.
"""

import tempfile
from pathlib import Path
from typing import List

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
from hypothesis import given, settings, strategies as st
import pytest

from qpx.core.sdrf import SDRFHandler
from qpx.core.format import IBAQ_SCHEMA
from qpx.operate.report import (
    get_factor_value_by_name,
    get_primary_factor_value,
    add_factor_column,
)

# Test data path
TEST_DATA_ROOT = Path(__file__).parents[1] / "examples"


# Strategies for generating test data
factor_name_strategy = st.text(
    alphabet=st.characters(categories=("L", "N", "P", "S")),
    min_size=1,
    max_size=50,
).filter(lambda x: x.strip() and "[" not in x and "]" not in x)

factor_value_strategy = st.text(
    alphabet=st.characters(categories=("L", "N", "P", "S")),
    min_size=0,
    max_size=100,
)

factor_pair_strategy = st.tuples(factor_name_strategy, factor_value_strategy)

factor_values_list_strategy = st.lists(
    factor_pair_strategy, min_size=0, max_size=5, unique_by=lambda x: x[0]
)


def create_mock_sdrf_with_factors(factor_names: List[str], num_samples: int = 3) -> str:
    """Create a mock SDRF file with specified factor columns."""
    columns = [
        "source name",
        "comment[data file]",
        "comment[fraction identifier]",
        "comment[technical replicate]",
        "comment[label]",
        "characteristics[biological replicate]",
    ]
    columns.extend([f"factor value[{name}]" for name in factor_names])

    data = []
    for i in range(num_samples):
        row = [
            f"sample_{i}",
            f"file_{i}.raw",
            "1",
            "1",
            "label free sample",
            str(i + 1),
        ]
        row.extend([f"value_{i}_{j}" for j in range(len(factor_names))])
        data.append(row)

    df = pd.DataFrame(data, columns=columns)

    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".sdrf.tsv", delete=False, encoding="utf-8"
    ) as f:
        df.to_csv(f, sep="\t", index=False)
        return f.name


class TestFactorNameExtraction:
    """Property tests for factor name extraction completeness."""

    @given(st.lists(factor_name_strategy, min_size=0, max_size=5, unique=True))
    @settings(max_examples=100, deadline=None)
    def test_factor_name_extraction_completeness(self, factor_names: List[str]):
        """
        For any SDRF data, get_factor_names() should return a list containing
        all factor names parsed from factor value[<name>] columns.
        Note: SDRF handler lowercases all column names, so factor names are lowercased.
        """
        if not factor_names:
            return  # Skip empty case as SDRF requires at least some columns

        sdrf_path = create_mock_sdrf_with_factors(factor_names)
        try:
            handler = SDRFHandler(sdrf_path)
            extracted_names = handler.get_factor_names()

            # Property: count should match
            if len(extracted_names) != len(factor_names):
                pytest.fail(
                    f"Expected {len(factor_names)} factor names, got {len(extracted_names)}"
                )

            # Property: all names should be present (lowercased due to SDRF normalization)
            for name in factor_names:
                if name.lower() not in extracted_names:
                    pytest.fail(
                        f"Factor name '{name.lower()}' not found in extracted names"
                    )
        finally:
            Path(sdrf_path).unlink(missing_ok=True)

    def test_empty_sdrf_returns_empty_list(self):
        """SDRF with no factor value columns returns empty list."""
        # Create SDRF without factor columns
        columns = [
            "source name",
            "comment[data file]",
            "comment[fraction identifier]",
            "comment[technical replicate]",
            "comment[label]",
        ]
        data = [["sample_1", "file_1.raw", "1", "1", "label free sample"]]
        df = pd.DataFrame(data, columns=columns)

        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".sdrf.tsv", delete=False, encoding="utf-8"
        ) as f:
            df.to_csv(f, sep="\t", index=False)
            sdrf_path = f.name

        try:
            handler = SDRFHandler(sdrf_path)
            if handler.get_factor_names() != []:
                pytest.fail("Expected empty list for SDRF without factor columns")
        finally:
            Path(sdrf_path).unlink(missing_ok=True)


class TestFactorValuesRoundTrip:
    """Property tests for factor values round-trip consistency."""

    @given(factor_values_list_strategy)
    @settings(max_examples=100, deadline=None)
    def test_factor_values_round_trip(self, factor_pairs):
        """
        For any sample's factor values data, storing to Parquet and reading back
        should produce identical factor names and values.
        """
        if not factor_pairs:
            return

        # Build factor_values structure
        factor_values = [
            {"factor_name": name, "factor_value": value} for name, value in factor_pairs
        ]

        # Create minimal IBAQ-like data
        data = {
            "pg_accessions": [["P12345"]],
            "sequence": ["PEPTIDE"],
            "peptidoform": ["PEPTIDE"],
            "precursor_charge": [2],
            "unique": [1],
            "reference_file_name": ["test_file"],
            "sample_accession": ["sample_1"],
            "run": ["1_1_1"],
            "channel": ["LFQ"],
            "factor_values": [factor_values],
            "fraction": ["1"],
            "biological_replicate": [1],
            "intensity": [1000.0],
        }
        df = pd.DataFrame(data)

        with tempfile.NamedTemporaryFile(suffix=".parquet", delete=False) as f:
            parquet_path = f.name

        try:
            # Write to parquet
            table = pa.Table.from_pandas(df, schema=IBAQ_SCHEMA)
            pq.write_table(table, parquet_path)

            # Read back
            read_df = pd.read_parquet(parquet_path)
            read_factor_values = read_df["factor_values"].iloc[0]

            # Property: round-trip should preserve data
            if len(read_factor_values) != len(factor_values):
                pytest.fail(
                    f"Round-trip length mismatch: {len(read_factor_values)} != {len(factor_values)}"
                )
            for orig, read in zip(factor_values, read_factor_values):
                if orig["factor_name"] != read["factor_name"]:
                    pytest.fail(
                        f"Factor name mismatch: {orig['factor_name']} != {read['factor_name']}"
                    )
                if orig["factor_value"] != read["factor_value"]:
                    pytest.fail(
                        f"Factor value mismatch: {orig['factor_value']} != {read['factor_value']}"
                    )
        finally:
            Path(parquet_path).unlink(missing_ok=True)


class TestGroupingConsistency:
    """Property tests for grouping consistency."""

    @given(
        st.lists(
            st.tuples(
                st.text(min_size=1, max_size=10), st.text(min_size=1, max_size=20)
            ),
            min_size=2,
            max_size=10,
        )
    )
    @settings(max_examples=100, deadline=None)
    def test_grouping_by_factor(self, sample_factor_pairs):
        """Samples with the same factor value should be grouped together."""
        data = []
        for i, (_, factor_val) in enumerate(sample_factor_pairs):
            data.append(
                {
                    "sample_accession": f"sample_{i}",
                    "factor_values": [
                        {"factor_name": "test_factor", "factor_value": factor_val}
                    ],
                }
            )

        df = pd.DataFrame(data)
        df = add_factor_column(df, factor_name="test_factor", column_name="Condition")

        # Property: grouping should be consistent
        for factor_val in df["Condition"].unique():
            group = df[df["Condition"] == factor_val]
            # All samples in group should have same factor value
            for _, row in group.iterrows():
                extracted = get_factor_value_by_name(
                    row["factor_values"], "test_factor"
                )
                if extracted != factor_val:
                    pytest.fail(f"Grouping inconsistency: {extracted} != {factor_val}")

    def test_get_factor_value_by_name_returns_correct_value(self):
        """get_factor_value_by_name returns the correct value for a given factor."""
        factor_values = [
            {"factor_name": "organism part", "factor_value": "brain"},
            {"factor_name": "disease", "factor_value": "Alzheimer"},
        ]

        if get_factor_value_by_name(factor_values, "organism part") != "brain":
            pytest.fail("Expected 'brain' for 'organism part'")
        if get_factor_value_by_name(factor_values, "disease") != "Alzheimer":
            pytest.fail("Expected 'Alzheimer' for 'disease'")
        if get_factor_value_by_name(factor_values, "nonexistent") is not None:
            pytest.fail("Expected None for 'nonexistent'")

    def test_get_primary_factor_value_returns_first(self):
        """get_primary_factor_value returns the first factor's value."""
        factor_values = [
            {"factor_name": "organism part", "factor_value": "brain"},
            {"factor_name": "disease", "factor_value": "Alzheimer"},
        ]

        if get_primary_factor_value(factor_values) != "brain":
            pytest.fail("Expected 'brain' as primary factor value")
        if get_primary_factor_value([]) is not None:
            pytest.fail("Expected None for empty list")
        if get_primary_factor_value(None) is not None:
            pytest.fail("Expected None for None input")


class TestSpecialCharacterHandling:
    """Property tests for special character handling."""

    @given(
        st.text(
            alphabet=st.characters(
                categories=("L", "N", "P", "S"),
                include_characters="|,\"'\\;:<>[]{}",
            ),
            min_size=1,
            max_size=50,
        )
    )
    @settings(max_examples=100, deadline=None)
    def test_special_characters_preserved(self, special_value):
        """
        For any factor value containing special characters, the value should be
        stored and retrieved without data loss.
        """
        factor_values = [{"factor_name": "test", "factor_value": special_value}]

        # Create minimal data
        data = {
            "pg_accessions": [["P12345"]],
            "sequence": ["PEPTIDE"],
            "peptidoform": ["PEPTIDE"],
            "precursor_charge": [2],
            "unique": [1],
            "reference_file_name": ["test_file"],
            "sample_accession": ["sample_1"],
            "run": ["1_1_1"],
            "channel": ["LFQ"],
            "factor_values": [factor_values],
            "fraction": ["1"],
            "biological_replicate": [1],
            "intensity": [1000.0],
        }
        df = pd.DataFrame(data)

        with tempfile.NamedTemporaryFile(suffix=".parquet", delete=False) as f:
            parquet_path = f.name

        try:
            table = pa.Table.from_pandas(df, schema=IBAQ_SCHEMA)
            pq.write_table(table, parquet_path)

            read_df = pd.read_parquet(parquet_path)
            read_value = read_df["factor_values"].iloc[0][0]["factor_value"]

            # Property: special characters should be preserved
            if read_value != special_value:
                pytest.fail(
                    f"Special characters not preserved: {read_value} != {special_value}"
                )
        finally:
            Path(parquet_path).unlink(missing_ok=True)


# Integration tests with real SDRF files
class TestRealSDRFFiles:
    """Integration tests using real SDRF files."""

    def test_maxquant_sdrf_factor_extraction(self):
        """Test factor extraction from MaxQuant SDRF file."""
        sdrf_path = TEST_DATA_ROOT / "maxquant/maxquant_full/PXD001819.sdrf.tsv"
        if not sdrf_path.exists():
            return

        handler = SDRFHandler(sdrf_path)
        factor_names = handler.get_factor_names()

        # Should have factor names
        if len(factor_names) <= 0:
            pytest.fail("Expected at least one factor name")

        # transform_sdrf should produce factor_values column
        df = handler.transform_sdrf()
        if "factor_values" not in df.columns:
            pytest.fail("Expected 'factor_values' column in transformed SDRF")

        # Each row should have factor_values matching factor_names count
        for _, row in df.iterrows():
            if len(row["factor_values"]) != len(factor_names):
                pytest.fail(
                    f"Factor values count mismatch: {len(row['factor_values'])} != {len(factor_names)}"
                )

    def test_quantms_lfq_sdrf_factor_extraction(self):
        """Test factor extraction from quantms LFQ SDRF file."""
        sdrf_path = TEST_DATA_ROOT / "quantms/dda-lfq-full/PXD007683-LFQ.sdrf.tsv"
        if not sdrf_path.exists():
            return

        handler = SDRFHandler(sdrf_path)
        handler.get_factor_names()

        df = handler.transform_sdrf()
        if "factor_values" not in df.columns:
            pytest.fail("Expected 'factor_values' column in transformed SDRF")


class TestMetadataFactorNamesCompleteness:
    """Property tests for metadata factor names completeness."""

    def test_ae_handler_get_factor_names(self):
        """Test that AE handler correctly retrieves factor names from SDRF."""
        from qpx.core.ae import AbsoluteExpressionHander

        sdrf_path = TEST_DATA_ROOT / "maxquant/maxquant_full/PXD001819.sdrf.tsv"
        if not sdrf_path.exists():
            return

        ae_handler = AbsoluteExpressionHander()
        ae_handler.load_sdrf_file(str(sdrf_path))

        factor_names = ae_handler.get_factor_names()

        # Should return factor names from SDRF
        if not isinstance(factor_names, list):
            pytest.fail("Expected factor_names to be a list")

        # Cross-check with SDRFHandler
        sdrf_handler = SDRFHandler(sdrf_path)
        expected_names = sdrf_handler.get_factor_names()

        if factor_names != expected_names:
            pytest.fail(f"Factor names mismatch: {factor_names} != {expected_names}")

    def test_de_handler_get_factor_names(self):
        """Test that DE handler correctly retrieves factor names from SDRF."""
        from qpx.core.de import DifferentialExpressionHandler

        sdrf_path = TEST_DATA_ROOT / "maxquant/maxquant_full/PXD001819.sdrf.tsv"
        if not sdrf_path.exists():
            return

        de_handler = DifferentialExpressionHandler()
        de_handler.load_sdrf_file(str(sdrf_path))

        factor_names = de_handler.get_factor_names()

        # Should return factor names from SDRF
        if not isinstance(factor_names, list):
            pytest.fail("Expected factor_names to be a list")

        # Cross-check with SDRFHandler
        sdrf_handler = SDRFHandler(sdrf_path)
        expected_names = sdrf_handler.get_factor_names()

        if factor_names != expected_names:
            pytest.fail(f"Factor names mismatch: {factor_names} != {expected_names}")

    def test_ae_handler_without_sdrf_returns_empty(self):
        """AE handler without SDRF loaded returns empty list."""
        from qpx.core.ae import AbsoluteExpressionHander

        ae_handler = AbsoluteExpressionHander()
        if ae_handler.get_factor_names() != []:
            pytest.fail("Expected empty list for AE handler without SDRF")

    def test_de_handler_without_sdrf_returns_empty(self):
        """DE handler without SDRF loaded returns empty list."""
        from qpx.core.de import DifferentialExpressionHandler

        de_handler = DifferentialExpressionHandler()
        if de_handler.get_factor_names() != []:
            pytest.fail("Expected empty list for DE handler without SDRF")


class TestEdgeCases:
    """Unit tests for edge cases in factor value handling."""

    def test_single_factor_sdrf(self):
        """Test SDRF with single factor value column."""
        sdrf_path = create_mock_sdrf_with_factors(["organism part"])
        try:
            handler = SDRFHandler(sdrf_path)
            factor_names = handler.get_factor_names()

            if len(factor_names) != 1:
                pytest.fail(f"Expected 1 factor name, got {len(factor_names)}")
            if "organism part" not in factor_names:
                pytest.fail("Expected 'organism part' in factor names")

            df = handler.transform_sdrf()
            if "factor_values" not in df.columns:
                pytest.fail("Expected 'factor_values' column")

            # Each row should have exactly one factor
            for _, row in df.iterrows():
                if len(row["factor_values"]) != 1:
                    pytest.fail(
                        f"Expected 1 factor value, got {len(row['factor_values'])}"
                    )
                if row["factor_values"][0]["factor_name"] != "organism part":
                    pytest.fail("Expected factor_name to be 'organism part'")
        finally:
            Path(sdrf_path).unlink(missing_ok=True)

    def test_multi_factor_sdrf(self):
        """Test SDRF with multiple factor value columns."""
        factor_names = ["organism part", "disease", "cell type"]
        sdrf_path = create_mock_sdrf_with_factors(factor_names)
        try:
            handler = SDRFHandler(sdrf_path)
            extracted_names = handler.get_factor_names()

            if len(extracted_names) != 3:
                pytest.fail(f"Expected 3 factor names, got {len(extracted_names)}")

            df = handler.transform_sdrf()
            if "factor_values" not in df.columns:
                pytest.fail("Expected 'factor_values' column")

            # Each row should have all factors
            for _, row in df.iterrows():
                if len(row["factor_values"]) != 3:
                    pytest.fail(
                        f"Expected 3 factor values, got {len(row['factor_values'])}"
                    )
        finally:
            Path(sdrf_path).unlink(missing_ok=True)

    def test_empty_factor_value(self):
        """Test handling of empty factor values."""
        # Create SDRF with empty factor value
        columns = [
            "source name",
            "comment[data file]",
            "comment[fraction identifier]",
            "comment[technical replicate]",
            "comment[label]",
            "characteristics[biological replicate]",
            "factor value[organism part]",
        ]
        data = [
            ["sample_1", "file_1.raw", "1", "1", "label free sample", "1", ""],
            ["sample_2", "file_2.raw", "1", "1", "label free sample", "2", "brain"],
        ]
        df = pd.DataFrame(data, columns=columns)

        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".sdrf.tsv", delete=False, encoding="utf-8"
        ) as f:
            df.to_csv(f, sep="\t", index=False)
            sdrf_path = f.name

        try:
            handler = SDRFHandler(sdrf_path)
            result_df = handler.transform_sdrf()

            # Empty values should be preserved as empty strings
            if "factor_values" not in result_df.columns:
                pytest.fail("Expected 'factor_values' column")
            for _, row in result_df.iterrows():
                if len(row["factor_values"]) != 1:
                    pytest.fail(
                        f"Expected 1 factor value, got {len(row['factor_values'])}"
                    )
        finally:
            Path(sdrf_path).unlink(missing_ok=True)

    def test_factor_value_with_spaces(self):
        """Test factor names and values with spaces."""
        factor_names = ["organism part", "cell line"]
        sdrf_path = create_mock_sdrf_with_factors(factor_names)
        try:
            handler = SDRFHandler(sdrf_path)
            extracted_names = handler.get_factor_names()

            # Names with spaces should be preserved
            if "organism part" not in extracted_names:
                pytest.fail("Expected 'organism part' in extracted names")
            if "cell line" not in extracted_names:
                pytest.fail("Expected 'cell line' in extracted names")
        finally:
            Path(sdrf_path).unlink(missing_ok=True)

    def test_extract_feature_properties_has_factor_values(self):
        """Test that extract_feature_properties includes factor_values."""
        sdrf_path = TEST_DATA_ROOT / "maxquant/maxquant_full/PXD001819.sdrf.tsv"
        if not sdrf_path.exists():
            return

        handler = SDRFHandler(sdrf_path)
        df = handler.extract_feature_properties()

        if "factor_values" not in df.columns:
            pytest.fail("Expected 'factor_values' column")

        # Verify structure
        for _, row in df.iterrows():
            factor_values = row["factor_values"]
            if not isinstance(factor_values, list):
                pytest.fail("Expected factor_values to be a list")
            for fv in factor_values:
                if "factor_name" not in fv:
                    pytest.fail("Expected 'factor_name' key in factor value")
                if "factor_value" not in fv:
                    pytest.fail("Expected 'factor_value' key in factor value")
