"""Tests for the QPX validation feature.

Covers: ValidationIssue, ValidationResult dataclasses,
ViewSchema.validate_full(), BaseStructure.validate(),
Dataset.validate(), and the CLI validate command.
"""

import pyarrow as pa
import pytest
from click.testing import CliRunner

from qpx.core.data import (
    ValidationIssue,
    ValidationResult,
    ViewSchema,
    FieldDef,
    FeatureSchema,
    PsmSchema,
    PgSchema,
    MzSchema,
    SampleSchema,
    RunSchema,
    DatasetSchema,
    OntologySchema,
    ProvenanceSchema,
)


# ---------------------------------------------------------------------------
# ValidationResult dataclass tests
# ---------------------------------------------------------------------------


class TestValidationResult:

    def test_is_valid_when_no_issues(self):
        r = ValidationResult(structure="feature")
        assert r.is_valid is True

    def test_is_valid_with_only_warnings(self):
        r = ValidationResult(
            structure="feature",
            issues=[
                ValidationIssue("feature", "null_values", "warning", "col", "msg")
            ],
        )
        assert r.is_valid is True

    def test_is_invalid_with_error(self):
        r = ValidationResult(
            structure="feature",
            issues=[
                ValidationIssue("feature", "missing_column", "error", "col", "msg")
            ],
        )
        assert r.is_valid is False

    def test_errors_property_filters_errors(self):
        r = ValidationResult(
            structure="feature",
            issues=[
                ValidationIssue("feature", "missing_column", "error", "a", "err"),
                ValidationIssue("feature", "null_values", "warning", "b", "warn"),
                ValidationIssue("feature", "type_mismatch", "error", "c", "err2"),
            ],
        )
        assert len(r.errors) == 2
        assert len(r.warnings) == 1

    def test_summary_valid_no_issues(self):
        r = ValidationResult(structure="feature")
        assert r.summary == "feature: VALID"

    def test_summary_valid_with_warnings(self):
        r = ValidationResult(
            structure="feature",
            issues=[
                ValidationIssue("feature", "null_values", "warning", "x", "w"),
            ],
        )
        assert "VALID" in r.summary
        assert "1 warning" in r.summary

    def test_summary_invalid(self):
        r = ValidationResult(
            structure="pg",
            issues=[
                ValidationIssue("pg", "missing_column", "error", "x", "e"),
                ValidationIssue("pg", "null_values", "warning", "y", "w"),
            ],
        )
        assert "INVALID" in r.summary
        assert "1 error" in r.summary
        assert "1 warning" in r.summary


# ---------------------------------------------------------------------------
# ViewSchema.validate_full() tests
# ---------------------------------------------------------------------------


class TestValidateFull:

    def test_valid_table_returns_no_issues(self):
        schema = FeatureSchema.get_arrow_schema()
        arrays = {f.name: pa.nulls(1, type=f.type) for f in schema}
        table = pa.table(arrays, schema=schema)
        result = FeatureSchema.validate_full(table)
        # No errors (may have warnings for nulls in required cols, which is expected)
        assert result.is_valid

    def test_catches_missing_required_column(self):
        schema = FeatureSchema.get_arrow_schema()
        fields_to_keep = [f for f in schema if f.name != "sequence"]
        arrays = {f.name: pa.nulls(1, type=f.type) for f in fields_to_keep}
        table = pa.table(arrays, schema=pa.schema(fields_to_keep))

        result = FeatureSchema.validate_full(table)
        assert not result.is_valid
        missing = [i for i in result.issues if i.check == "missing_column"]
        assert any("sequence" in i.message for i in missing)

    def test_catches_type_mismatch(self):
        schema = FeatureSchema.get_arrow_schema()
        arrays = {}
        wrong_fields = []
        for f in schema:
            if f.name == "charge":
                arrays[f.name] = pa.array(["wrong"], type=pa.string())
                wrong_fields.append(pa.field(f.name, pa.string()))
            else:
                arrays[f.name] = pa.nulls(1, type=f.type)
                wrong_fields.append(f)

        table = pa.table(arrays, schema=pa.schema(wrong_fields))
        result = FeatureSchema.validate_full(table)
        assert not result.is_valid
        type_errors = [i for i in result.issues if i.check == "type_mismatch"]
        assert any("charge" in i.message for i in type_errors)

    def test_allows_optional_column_absent(self):
        schema = FeatureSchema.get_arrow_schema()
        fields_to_keep = [f for f in schema if f.name != "pg_global_qvalue"]
        arrays = {f.name: pa.nulls(1, type=f.type) for f in fields_to_keep}
        table = pa.table(arrays, schema=pa.schema(fields_to_keep))

        result = FeatureSchema.validate_full(table)
        missing_issues = [
            i for i in result.issues
            if i.check == "missing_column" and "pg_global_qvalue" in i.message
        ]
        assert len(missing_issues) == 0

    def test_catches_null_in_non_nullable_column(self):
        """Non-nullable column with null values should produce a warning."""
        schema = FeatureSchema.get_arrow_schema()
        arrays = {}
        for f in schema:
            if f.name == "sequence":
                # Explicitly create column with nulls
                arrays[f.name] = pa.array([None], type=pa.string())
            else:
                arrays[f.name] = pa.nulls(1, type=f.type)
        table = pa.table(arrays, schema=schema)

        result = FeatureSchema.validate_full(table)
        null_issues = [
            i for i in result.issues
            if i.check == "null_values" and i.column == "sequence"
        ]
        assert len(null_issues) == 1
        assert null_issues[0].severity == "warning"

    def test_allows_null_in_nullable_column(self):
        """Nullable columns with null values should not produce a warning."""
        schema = FeatureSchema.get_arrow_schema()
        arrays = {f.name: pa.nulls(1, type=f.type) for f in schema}
        table = pa.table(arrays, schema=schema)

        result = FeatureSchema.validate_full(table)
        # "posterior_error_probability" is nullable (not required), so no warning for it
        pep_null_issues = [
            i for i in result.issues
            if i.check == "null_values"
            and i.column == "posterior_error_probability"
        ]
        assert len(pep_null_issues) == 0

    def test_catches_duplicate_primary_key(self):
        """Duplicate primary key rows should produce a warning."""
        # DatasetSchema has PK = [project_accession] — easy to create dupes
        schema = DatasetSchema.get_arrow_schema()
        row = {f.name: pa.nulls(1, type=f.type) for f in schema}
        # Two rows with same project_accession
        arrays = {}
        for f in schema:
            if f.name == "project_accession":
                arrays[f.name] = pa.array(["PXD001", "PXD001"], type=f.type)
            elif f.name == "creation_date":
                arrays[f.name] = pa.array(["2024-01-01", "2024-01-02"], type=f.type)
            elif f.name == "qpx_version":
                arrays[f.name] = pa.array(["1.0", "1.0"], type=f.type)
            else:
                arrays[f.name] = pa.nulls(2, type=f.type)
        table = pa.table(arrays, schema=schema)

        result = DatasetSchema.validate_full(table)
        pk_issues = [i for i in result.issues if i.check == "duplicate_pk"]
        assert len(pk_issues) == 1
        assert pk_issues[0].severity == "warning"
        assert "duplicate" in pk_issues[0].message.lower()

    def test_no_pk_issue_when_unique(self):
        """Unique primary key rows should not produce a warning."""
        schema = DatasetSchema.get_arrow_schema()
        arrays = {}
        for f in schema:
            if f.name == "project_accession":
                arrays[f.name] = pa.array(["PXD001", "PXD002"], type=f.type)
            elif f.name == "creation_date":
                arrays[f.name] = pa.array(["2024-01-01", "2024-01-02"], type=f.type)
            elif f.name == "qpx_version":
                arrays[f.name] = pa.array(["1.0", "1.0"], type=f.type)
            else:
                arrays[f.name] = pa.nulls(2, type=f.type)
        table = pa.table(arrays, schema=schema)

        result = DatasetSchema.validate_full(table)
        pk_issues = [i for i in result.issues if i.check == "duplicate_pk"]
        assert len(pk_issues) == 0


# ---------------------------------------------------------------------------
# Backward compatibility: validate() returns list[str]
# ---------------------------------------------------------------------------


class TestValidateBackwardCompat:

    def test_validate_returns_list_of_strings(self):
        schema = FeatureSchema.get_arrow_schema()
        fields_to_keep = [f for f in schema if f.name != "sequence"]
        arrays = {f.name: pa.nulls(1, type=f.type) for f in fields_to_keep}
        table = pa.table(arrays, schema=pa.schema(fields_to_keep))

        errors = FeatureSchema.validate(table)
        assert isinstance(errors, list)
        assert all(isinstance(e, str) for e in errors)
        assert any("sequence" in e for e in errors)

    def test_validate_returns_empty_list_for_valid_table(self):
        schema = FeatureSchema.get_arrow_schema()
        arrays = {f.name: pa.nulls(1, type=f.type) for f in schema}
        table = pa.table(arrays, schema=schema)

        errors = FeatureSchema.validate(table)
        assert errors == []


# ---------------------------------------------------------------------------
# validate_full across all schemas
# ---------------------------------------------------------------------------


ALL_SCHEMAS = [
    FeatureSchema, PsmSchema, PgSchema, MzSchema,
    SampleSchema, RunSchema, DatasetSchema, OntologySchema, ProvenanceSchema,
]


@pytest.mark.parametrize("schema_inst", ALL_SCHEMAS, ids=lambda s: s.__name__)
class TestValidateFullAllSchemas:

    def test_valid_table_no_errors(self, schema_inst):
        """A table matching the schema should have no errors."""
        schema = schema_inst.get_arrow_schema()
        arrays = {f.name: pa.nulls(1, type=f.type) for f in schema}
        table = pa.table(arrays, schema=schema)
        result = schema_inst.validate_full(table)
        assert result.is_valid


# ---------------------------------------------------------------------------
# BaseStructure.validate() integration tests
# ---------------------------------------------------------------------------


class TestBaseStructureValidate:

    def test_feature_validates_clean_parquet(self, feature_parquet):
        from qpx.core.data.feature import Feature
        feat = Feature.from_file(feature_parquet)
        result = feat.validate()
        assert isinstance(result, ValidationResult)
        assert result.is_valid

    def test_psm_validates_clean_parquet(self, psm_parquet):
        from qpx.core.data.psm import PSM
        psm = PSM.from_file(psm_parquet)
        result = psm.validate()
        assert result.is_valid

    def test_pg_validates_clean_parquet(self, pg_parquet):
        from qpx.core.data.pg import PG
        pg = PG.from_file(pg_parquet)
        result = pg.validate()
        assert result.is_valid


# ---------------------------------------------------------------------------
# Dataset.validate() integration tests
# ---------------------------------------------------------------------------


class TestDatasetValidate:

    def test_validate_all_structures(self, dataset_dir):
        import qpx
        with qpx.open(dataset_dir) as ds:
            results = ds.validate()
        assert len(results) > 0
        for name, result in results.items():
            assert isinstance(result, ValidationResult)
            assert result.is_valid, f"{name}: {result.summary}"

    def test_validate_specific_structure(self, dataset_dir):
        import qpx
        with qpx.open(dataset_dir) as ds:
            results = ds.validate(structures=["feature"])
        assert "feature" in results
        assert len(results) == 1
        assert results["feature"].is_valid

    def test_validate_missing_structure(self, dataset_dir):
        import qpx
        with qpx.open(dataset_dir) as ds:
            results = ds.validate(structures=["mz"])
        assert "mz" in results
        assert not results["mz"].is_valid
        assert any(i.check == "missing_structure" for i in results["mz"].issues)


# ---------------------------------------------------------------------------
# CLI tests
# ---------------------------------------------------------------------------


class TestCLIValidate:

    def test_validate_dataset(self, dataset_dir):
        from qpx.cli.main import qpx_main
        runner = CliRunner()
        result = runner.invoke(qpx_main, ["validate", "--dataset-path", str(dataset_dir)])
        assert result.exit_code == 0
        assert "VALID" in result.output

    def test_validate_specific_structure(self, dataset_dir):
        from qpx.cli.main import qpx_main
        runner = CliRunner()
        result = runner.invoke(
            qpx_main,
            ["validate", "--dataset-path", str(dataset_dir), "--structure", "feature"],
        )
        assert result.exit_code == 0
        assert "feature" in result.output

    def test_validate_single_file(self, feature_parquet):
        from qpx.cli.main import qpx_main
        runner = CliRunner()
        result = runner.invoke(
            qpx_main, ["validate", "--file", str(feature_parquet)]
        )
        assert result.exit_code == 0
        assert "VALID" in result.output

    def test_validate_no_args_errors(self):
        from qpx.cli.main import qpx_main
        runner = CliRunner()
        result = runner.invoke(qpx_main, ["validate"])
        assert result.exit_code != 0
