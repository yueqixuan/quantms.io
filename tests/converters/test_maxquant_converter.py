"""MaxQuant converter integration tests using small example dataset.

Runs the full MaxQuant conversion pipeline once (module scope) and validates
that psm.parquet, feature.parquet, and ontology.parquet are written with
correct schemas and plausible values.
"""

import pyarrow.parquet as pq
import pytest
from pathlib import Path

EXAMPLES_DIR = Path(__file__).resolve().parent.parent / "examples" / "maxquant" / "maxquant_simple"

_MSMS = EXAMPLES_DIR / "msms.txt"
_EVIDENCE = EXAMPLES_DIR / "evidence.txt"
_SDRF = EXAMPLES_DIR / "sdrf.tsv"

_PREFIX = "maxquant_test"


# ---------------------------------------------------------------------------
# Module-scoped fixture: run conversion once, share outputs across tests
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def converted_output(tmp_path_factory):
    """Run MaxQuant conversion once for all tests in this module.

    Returns the output directory containing the generated Parquet files.
    If the test data is missing, all tests in this module are skipped.
    """
    from qpx.converters.maxquant.converter import MaxQuantConverter

    if not _MSMS.exists() or not _EVIDENCE.exists():
        pytest.skip(f"Test data not found: {EXAMPLES_DIR}")

    output = tmp_path_factory.mktemp("maxquant_output")

    converter = MaxQuantConverter()
    converter.convert(
        output_folder=output,
        msms_file=_MSMS,
        evidence_file=_EVIDENCE,
        sdrf_file=_SDRF if _SDRF.exists() else None,
        output_prefix=_PREFIX,
    )

    return output


@pytest.fixture(scope="module")
def psm_table(converted_output):
    """Read the psm.parquet produced by the converter."""
    path = converted_output / f"{_PREFIX}.psm.parquet"
    if not path.exists():
        pytest.skip("psm.parquet was not produced")
    return pq.read_table(str(path))


@pytest.fixture(scope="module")
def feature_table(converted_output):
    """Read the feature.parquet produced by the converter."""
    path = converted_output / f"{_PREFIX}.feature.parquet"
    if not path.exists():
        pytest.skip("feature.parquet was not produced")
    return pq.read_table(str(path))


# ---------------------------------------------------------------------------
# PSM conversion tests
# ---------------------------------------------------------------------------


@pytest.mark.integration
class TestMaxQuantPsmConversion:
    """Validate psm.parquet output from MaxQuant conversion."""

    def test_file_exists(self, converted_output):
        path = converted_output / f"{_PREFIX}.psm.parquet"
        assert path.exists(), "psm.parquet was not created"

    def test_has_rows(self, psm_table):
        assert psm_table.num_rows > 0, "psm.parquet is empty"

    def test_key_columns_present(self, psm_table):
        column_names = set(psm_table.column_names)
        expected = {"sequence", "charge", "calculated_mz", "run_file_name", "scan"}
        missing = expected - column_names
        assert not missing, f"Missing columns: {missing}"

    def test_sequence_values_are_nonempty(self, psm_table):
        sequences = psm_table.column("sequence").to_pylist()
        for seq in sequences:
            assert isinstance(seq, str) and len(seq) > 0, (
                f"Expected non-empty string, got {seq!r}"
            )

    def test_charge_values_are_valid(self, psm_table):
        charges = psm_table.column("charge").to_pylist()
        for charge in charges:
            assert 1 <= charge <= 10, f"Charge out of range: {charge}"

    def test_calculated_mz_values_are_nonnegative(self, psm_table):
        mz_values = psm_table.column("calculated_mz").to_pylist()
        for mz in mz_values:
            assert mz is not None and mz >= 0, (
                f"Invalid calculated_mz value: {mz}"
            )

    def test_run_file_names_are_nonempty(self, psm_table):
        run_names = psm_table.column("run_file_name").to_pylist()
        for name in run_names:
            assert isinstance(name, str) and len(name) > 0

    def test_schema_validation(self, psm_table):
        """Verify that the output conforms to the PsmSchema."""
        from qpx.core.data import PsmSchema

        errors = PsmSchema.validate(psm_table)
        assert not errors, f"Schema validation errors: {errors}"


# ---------------------------------------------------------------------------
# Feature conversion tests
# ---------------------------------------------------------------------------


@pytest.mark.integration
class TestMaxQuantFeatureConversion:
    """Validate feature.parquet output from MaxQuant conversion."""

    def test_file_exists(self, converted_output):
        path = converted_output / f"{_PREFIX}.feature.parquet"
        assert path.exists(), "feature.parquet was not created"

    def test_has_rows(self, feature_table):
        assert feature_table.num_rows > 0, "feature.parquet is empty"

    def test_key_columns_present(self, feature_table):
        column_names = set(feature_table.column_names)
        expected = {"sequence", "charge", "calculated_mz", "intensities", "run_file_name"}
        missing = expected - column_names
        assert not missing, f"Missing columns: {missing}"

    def test_sequence_values_are_nonempty(self, feature_table):
        sequences = feature_table.column("sequence").to_pylist()
        for seq in sequences:
            assert isinstance(seq, str) and len(seq) > 0, (
                f"Expected non-empty string, got {seq!r}"
            )

    def test_charge_values_are_valid(self, feature_table):
        charges = feature_table.column("charge").to_pylist()
        for charge in charges:
            assert 1 <= charge <= 10, f"Charge out of range: {charge}"

    def test_intensities_are_nonnegative(self, feature_table):
        rows = feature_table.column("intensities").to_pylist()
        for row_intensities in rows:
            if row_intensities is None:
                continue
            for entry in row_intensities:
                assert entry["intensity"] >= 0, (
                    f"Negative intensity: {entry['intensity']}"
                )

    def test_run_file_names_are_nonempty(self, feature_table):
        run_names = feature_table.column("run_file_name").to_pylist()
        for name in run_names:
            assert isinstance(name, str) and len(name) > 0

    def test_schema_validation(self, feature_table):
        """Verify that the output conforms to the FeatureSchema."""
        from qpx.core.data import FeatureSchema

        errors = FeatureSchema.validate(feature_table)
        assert not errors, f"Schema validation errors: {errors}"


# ---------------------------------------------------------------------------
# Ontology conversion tests
# ---------------------------------------------------------------------------


@pytest.mark.integration
class TestMaxQuantOntologyConversion:
    """Validate ontology.parquet output from MaxQuant conversion.

    The ontology file is only written when scores are discovered during
    conversion. It may legitimately not exist if no scores were tracked.
    """

    def test_ontology_file_exists_or_no_scores(self, converted_output):
        path = converted_output / f"{_PREFIX}.ontology.parquet"
        if not path.exists():
            pytest.skip("ontology.parquet not written (no scores discovered)")

        table = pq.read_table(str(path))
        assert table.num_rows > 0, "ontology.parquet exists but is empty"

    def test_ontology_columns(self, converted_output):
        path = converted_output / f"{_PREFIX}.ontology.parquet"
        if not path.exists():
            pytest.skip("ontology.parquet not written (no scores discovered)")

        table = pq.read_table(str(path))
        column_names = set(table.column_names)
        expected = {"field_name", "view"}
        missing = expected - column_names
        assert not missing, f"Missing ontology columns: {missing}"

    def test_schema_validation(self, converted_output):
        from qpx.core.data import OntologySchema

        path = converted_output / f"{_PREFIX}.ontology.parquet"
        if not path.exists():
            pytest.skip("ontology.parquet not written (no scores discovered)")
        table = pq.read_table(str(path))
        errors = OntologySchema.validate(table)
        assert not errors, f"Schema validation errors: {errors}"
