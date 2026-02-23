"""Tests for QPX schema layer: FieldDef, ViewSchema, and all schema classes."""

import pyarrow as pa
import pytest

from qpx.core.data import (
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
# FieldDef tests
# ---------------------------------------------------------------------------


class TestFieldDef:
    """Test the FieldDef container that maps to PyArrow fields."""

    def test_required_field_not_nullable(self):
        """A required field should NOT be nullable."""
        f = FieldDef("my_col", pa.string(), nullable=False, doc="test")
        arrow = f.to_arrow_field()
        assert arrow.name == "my_col"
        assert arrow.type == pa.string()
        assert arrow.nullable is False

    def test_nullable_field_is_nullable(self):
        """A field with nullable=True should be nullable."""
        f = FieldDef("score", pa.float64(), nullable=True, doc="test nullable")
        arrow = f.to_arrow_field()
        assert arrow.nullable is True

    def test_optional_field_is_nullable(self):
        """An optional field should also be nullable (column may be absent)."""
        f = FieldDef(
            "opt_col", pa.int32(), nullable=True, optional=True, doc="optional col"
        )
        arrow = f.to_arrow_field()
        assert arrow.nullable is True
        assert f.optional is True

    def test_default_field_is_nullable(self):
        """A field with default flags is nullable=True."""
        f = FieldDef("default_col", pa.string())
        arrow = f.to_arrow_field()
        assert arrow.nullable is True

    def test_field_preserves_arrow_type(self):
        """FieldDef should preserve the exact Arrow type."""
        f = FieldDef("scan", pa.list_(pa.int32()), nullable=False)
        arrow = f.to_arrow_field()
        assert arrow.type == pa.list_(pa.int32())

    def test_field_doc_in_metadata(self):
        """If doc is provided, it should appear in the Arrow field metadata."""
        f = FieldDef(
            "documented_col", pa.string(), nullable=False, doc="My documentation"
        )
        arrow = f.to_arrow_field()
        assert arrow.metadata is not None
        assert b"doc" in arrow.metadata
        assert arrow.metadata[b"doc"] == b"My documentation"

    def test_field_no_doc_no_metadata(self):
        """If no doc is provided, the Arrow field should have no metadata."""
        f = FieldDef("no_doc_col", pa.string(), nullable=False)
        arrow = f.to_arrow_field()
        assert arrow.metadata is None

    def test_field_complex_type(self):
        """FieldDef handles complex nested types (list of struct)."""
        score_type = pa.list_(
            pa.struct(
                [
                    pa.field("score_name", pa.string()),
                    pa.field("score_value", pa.float64()),
                    pa.field("higher_better", pa.bool_(), nullable=True),
                ]
            )
        )
        f = FieldDef("scores", score_type, nullable=True)
        arrow = f.to_arrow_field()
        assert arrow.type == score_type
        assert arrow.nullable is True


# ---------------------------------------------------------------------------
# ViewSchema tests
# ---------------------------------------------------------------------------


class TestViewSchema:
    """Test the ViewSchema base class and its schema generation/validation."""

    def test_get_arrow_schema_returns_pa_schema(self):
        """get_arrow_schema should return a PyArrow Schema object."""
        schema = FeatureSchema.get_arrow_schema()
        assert isinstance(schema, pa.Schema)

    def test_iter_fields_yields_field_def_instances(self):
        """_iter_fields should yield (name, FieldDef) tuples."""
        fields = list(FeatureSchema._iter_fields())
        assert len(fields) > 0
        for name, field in fields:
            assert isinstance(name, str)
            assert isinstance(field, FieldDef)

    def test_schema_field_names_match_iter_fields(self):
        """Schema field names should match what _iter_fields yields."""
        schema = FeatureSchema.get_arrow_schema()
        field_names_from_schema = set(schema.names)
        field_names_from_iter = {name for name, _ in FeatureSchema._iter_fields()}
        assert field_names_from_schema == field_names_from_iter

    def test_validate_valid_table_returns_no_errors(self):
        """A table matching the schema should validate without errors."""
        schema = FeatureSchema.get_arrow_schema()
        arrays = []
        for field in schema:
            arrays.append(pa.nulls(1, type=field.type))
        table = pa.table(
            {field.name: arr for field, arr in zip(schema, arrays)},
            schema=schema,
        )
        errors = FeatureSchema.validate(table)
        assert errors == []

    def test_validate_catches_missing_required_column(self):
        """Validation should catch a missing required column."""
        schema = FeatureSchema.get_arrow_schema()
        fields_to_keep = [f for f in schema if f.name != "sequence"]
        arrays = {f.name: pa.nulls(1, type=f.type) for f in fields_to_keep}
        partial_schema = pa.schema(fields_to_keep)
        table = pa.table(arrays, schema=partial_schema)

        errors = FeatureSchema.validate(table)
        assert len(errors) >= 1
        assert any("sequence" in e for e in errors)

    def test_validate_catches_type_mismatch(self):
        """Validation should catch a type mismatch for a column."""
        schema = FeatureSchema.get_arrow_schema()
        arrays = {}
        for field in schema:
            if field.name == "charge":
                arrays[field.name] = pa.array(["wrong_type"], type=pa.string())
            else:
                arrays[field.name] = pa.nulls(1, type=field.type)

        wrong_schema = pa.schema(
            [pa.field(f.name, pa.string()) if f.name == "charge" else f for f in schema]
        )
        table = pa.table(arrays, schema=wrong_schema)

        errors = FeatureSchema.validate(table)
        assert len(errors) >= 1
        assert any("charge" in e and "int16" in e for e in errors)

    def test_validate_allows_optional_column_absent(self):
        """Validation should NOT report an error for a missing optional column."""
        schema = FeatureSchema.get_arrow_schema()
        assert FeatureSchema._is_optional("pg_global_qvalue") is True

        fields_to_keep = [f for f in schema if f.name != "pg_global_qvalue"]
        arrays = {f.name: pa.nulls(1, type=f.type) for f in fields_to_keep}
        partial_schema = pa.schema(fields_to_keep)
        table = pa.table(arrays, schema=partial_schema)

        errors = FeatureSchema.validate(table)
        assert not any("pg_global_qvalue" in e for e in errors)

    def test_is_optional_returns_true_for_optional_fields(self):
        """_is_optional should return True for fields with optional=True."""
        assert FeatureSchema._is_optional("pg_global_qvalue") is True
        assert FeatureSchema._is_optional("gg_accessions") is True

    def test_is_optional_returns_false_for_required_fields(self):
        """_is_optional should return False for required fields."""
        assert FeatureSchema._is_optional("sequence") is False
        assert FeatureSchema._is_optional("charge") is False

    def test_is_optional_returns_false_for_nonexistent_field(self):
        """_is_optional should return False for a field that doesn't exist."""
        assert FeatureSchema._is_optional("nonexistent_field") is False


# ---------------------------------------------------------------------------
# All schema instances generate non-empty schemas
# ---------------------------------------------------------------------------


ALL_SCHEMAS = [
    FeatureSchema,
    PsmSchema,
    PgSchema,
    MzSchema,
    SampleSchema,
    RunSchema,
    DatasetSchema,
    OntologySchema,
    ProvenanceSchema,
]


@pytest.mark.parametrize("schema_inst", ALL_SCHEMAS, ids=lambda s: s.__name__)
class TestAllSchemas:
    """Parameterized tests that apply to every schema instance."""

    def test_has_primary_key(self, schema_inst):
        """Every schema should have a _primary_key tuple."""
        assert isinstance(schema_inst._primary_key, tuple)
        assert len(schema_inst._primary_key) > 0

    def test_primary_key_fields_exist_in_schema(self, schema_inst):
        """Primary key fields should exist in the generated schema."""
        schema = schema_inst.get_arrow_schema()
        for pk_field in schema_inst._primary_key:
            assert (
                pk_field in schema.names
            ), f"Primary key field '{pk_field}' not found in {schema_inst.__name__} schema"


# ---------------------------------------------------------------------------
# Schema-specific field count sanity checks
# ---------------------------------------------------------------------------


class TestSchemaFieldCounts:
    """Sanity checks on the expected number of fields per schema."""

    def test_feature_schema_field_count(self):
        """FeatureSchema should have a significant number of fields (>=25)."""
        schema = FeatureSchema.get_arrow_schema()
        assert len(schema) >= 25

    def test_psm_schema_has_spectral_fields(self):
        """PsmSchema should have spectral data fields like mz_array."""
        schema = PsmSchema.get_arrow_schema()
        names = schema.names
        assert "mz_array" in names
        assert "intensity_array" in names
        assert "ion_type_array" in names

    def test_pg_schema_has_quantification_fields(self):
        """PgSchema should have quantification fields."""
        schema = PgSchema.get_arrow_schema()
        names = schema.names
        assert "intensities" in names
        assert "peptides" in names
        assert "global_qvalue" in names

    def test_run_schema_has_samples_field(self):
        """RunSchema should have a samples field (list of structs)."""
        schema = RunSchema.get_arrow_schema()
        assert "samples" in schema.names

    def test_sample_schema_has_organism(self):
        """SampleSchema should have an organism field."""
        schema = SampleSchema.get_arrow_schema()
        assert "organism" in schema.names

    def test_mz_schema_has_spectral_arrays(self):
        """MzSchema should have mz and intensity array fields."""
        schema = MzSchema.get_arrow_schema()
        names = schema.names
        assert "mz" in names
        assert "intensity" in names
        assert "ms_level" in names

    def test_dataset_schema_has_project_accession(self):
        """DatasetSchema should have project_accession."""
        schema = DatasetSchema.get_arrow_schema()
        assert "project_accession" in schema.names

    def test_ontology_schema_has_field_name_and_view(self):
        """OntologySchema should have field_name and view."""
        schema = OntologySchema.get_arrow_schema()
        names = schema.names
        assert "field_name" in names
        assert "view" in names

    def test_ontology_schema_has_source_columns(self):
        """Verify ontology schema includes source_column_name and source_tool."""
        schema = OntologySchema.get_arrow_schema()
        field_names = set(schema.names)
        assert "source_column_name" in field_names
        assert "source_tool" in field_names

    def test_provenance_schema_has_step_fields(self):
        """ProvenanceSchema should have step_order and tool_name."""
        schema = ProvenanceSchema.get_arrow_schema()
        names = schema.names
        assert "step_order" in names
        assert "tool_name" in names


# ---------------------------------------------------------------------------
# YAML loader tests
# ---------------------------------------------------------------------------


class TestYAMLLoader:
    """Test the YAML schema loader."""

    def test_load_schema_returns_view_schema(self):
        """load_schema should return a ViewSchema instance."""
        from qpx.core.data.loader import load_schema

        schema = load_schema("feature")
        assert isinstance(schema, ViewSchema)

    def test_type_resolution_nested_structs(self):
        """Verify nested struct types resolve correctly."""
        schema = FeatureSchema.get_arrow_schema()

        # modifications should be list<struct<name, accession, positions>>
        mod_type = schema.field("modifications").type
        assert isinstance(mod_type, pa.ListType)
        inner = mod_type.value_type
        assert isinstance(inner, pa.StructType)
        assert inner.get_field_index("name") >= 0
        assert inner.get_field_index("positions") >= 0

        # intensities should be list<struct<label, intensity>>
        int_type = schema.field("intensities").type
        assert isinstance(int_type, pa.ListType)
        inner = int_type.value_type
        assert isinstance(inner, pa.StructType)
        assert inner.get_field_index("label") >= 0
        assert inner.get_field_index("intensity") >= 0

    def test_sample_extra_columns(self):
        """Verify SampleSchema supports extending with extra string columns."""
        extended = SampleSchema.get_extended_arrow_schema(
            ["bmi", "biological_replicate"]
        )
        assert "bmi" in extended.names
        assert "biological_replicate" in extended.names
        assert extended.field("bmi").type == pa.string()
        assert extended.field("biological_replicate").type == pa.string()
        # Core fields still present
        assert extended.field("organism").type == pa.string()
        assert extended.field("organism_part").type == pa.string()


# ---------------------------------------------------------------------------
# P0-1: PK uniqueness check with list-typed columns (scan field)
# ---------------------------------------------------------------------------

class TestPKUniquenessWithListColumns:
    """PK uniqueness check must work when PK includes list<int32> columns like scan."""

    def test_psm_duplicate_pk_with_list_scan_reports_warning(self):
        """Duplicate PSM rows (same scan list) should produce duplicate_pk warning."""
        schema = PsmSchema.get_arrow_schema()
        row = {f.name: None for f in schema}
        row["sequence"] = "PEPTIDEK"
        row["charge"] = 2
        row["run_file_name"] = "run1"
        row["scan"] = [100]
        row["is_decoy"] = False
        rows = [row, dict(row)]  # exact duplicate

        table = pa.Table.from_pylist(rows, schema=schema)
        result = PsmSchema.validate_full(table)
        checks = [i.check for i in result.issues]
        assert "duplicate_pk" in checks, f"Expected duplicate_pk, got: {checks}"

    def test_psm_distinct_scan_lists_no_duplicate_warning(self):
        """PSM rows with different scan lists should NOT produce duplicate_pk warning."""
        schema = PsmSchema.get_arrow_schema()
        row_a = {f.name: None for f in schema}
        row_a["sequence"] = "PEPTIDEK"
        row_a["charge"] = 2
        row_a["run_file_name"] = "run1"
        row_a["scan"] = [100]
        row_a["is_decoy"] = False
        row_b = dict(row_a)
        row_b["scan"] = [200]  # different scan

        table = pa.Table.from_pylist([row_a, row_b], schema=schema)
        result = PsmSchema.validate_full(table)
        checks = [i.check for i in result.issues]
        assert "duplicate_pk" not in checks, f"Unexpected duplicate_pk with distinct scans"
