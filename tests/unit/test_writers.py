"""Tests for QPX writer layer: BaseWriter subclasses and round-trip validation."""

import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from qpx.writers import (
    FeatureWriter,
    PsmWriter,
    PgWriter,
    SampleWriter,
    RunWriter,
    DatasetWriter,
    OntologyWriter,
    ProvenanceWriter,
)
from qpx.core.data import (
    FeatureSchema,
    PsmSchema,
    PgSchema,
    SampleSchema,
    RunSchema,
    DatasetSchema,
    OntologySchema,
    ProvenanceSchema,
)
from qpx.core.parquet_io import read_parquet_metadata, parquet_row_count

from tests.conftest import (
    make_feature_record,
    make_psm_record,
    make_pg_record,
    make_sample_record,
    make_run_record,
    make_dataset_record,
    make_ontology_record,
    make_provenance_record,
)

# ---------------------------------------------------------------------------
# Writer creates valid Parquet file tests
# ---------------------------------------------------------------------------


class TestFeatureWriter:
    """Test FeatureWriter creates valid Parquet files."""

    def test_creates_valid_parquet(self, tmp_path):
        """FeatureWriter produces a readable Parquet file."""
        path = tmp_path / "test.feature.parquet"
        with FeatureWriter(path) as w:
            w.write_batch([make_feature_record()])

        assert path.exists()
        table = pq.read_table(path)
        assert table.num_rows == 1

    def test_schema_matches(self, tmp_path):
        """Written file schema matches FeatureSchema."""
        path = tmp_path / "test.feature.parquet"
        with FeatureWriter(path) as w:
            w.write_batch([make_feature_record()])

        pf = pq.ParquetFile(path)
        file_schema = pf.schema_arrow
        expected = FeatureSchema.get_arrow_schema()
        for field in expected:
            assert field.name in file_schema.names, f"Missing field: {field.name}"

    def test_multiple_rows(self, tmp_path):
        """Writing multiple records works."""
        path = tmp_path / "test.feature.parquet"
        records = [make_feature_record(sequence=f"SEQ{i}") for i in range(5)]
        with FeatureWriter(path) as w:
            w.write_batch(records)

        assert parquet_row_count(path) == 5

    def test_footer_metadata(self, tmp_path):
        """File footer contains QPX metadata."""
        path = tmp_path / "test.feature.parquet"
        with FeatureWriter(path, creator="test_suite") as w:
            w.write_batch([make_feature_record()])

        meta = read_parquet_metadata(path)
        assert meta["file_type"] == "feature_file"
        assert meta["creator"] == "test_suite"
        assert "qpx_version" in meta
        assert "uuid" in meta
        assert "creation_date" in meta


class TestPsmWriter:
    """Test PsmWriter creates valid Parquet files."""

    def test_creates_valid_parquet(self, tmp_path):
        """PsmWriter produces a readable Parquet file."""
        path = tmp_path / "test.psm.parquet"
        with PsmWriter(path) as w:
            w.write_batch([make_psm_record()])

        assert path.exists()
        table = pq.read_table(path)
        assert table.num_rows == 1

    def test_schema_matches(self, tmp_path):
        """Written file schema matches PsmSchema."""
        path = tmp_path / "test.psm.parquet"
        with PsmWriter(path) as w:
            w.write_batch([make_psm_record()])

        pf = pq.ParquetFile(path)
        file_schema = pf.schema_arrow
        expected = PsmSchema.get_arrow_schema()
        for field in expected:
            assert field.name in file_schema.names

    def test_footer_metadata(self, tmp_path):
        """PSM file footer contains correct file_type."""
        path = tmp_path / "test.psm.parquet"
        with PsmWriter(path) as w:
            w.write_batch([make_psm_record()])

        meta = read_parquet_metadata(path)
        assert meta["file_type"] == "psm_file"


class TestPgWriter:
    """Test PgWriter creates valid Parquet files."""

    def test_creates_valid_parquet(self, tmp_path):
        """PgWriter produces a readable Parquet file."""
        path = tmp_path / "test.pg.parquet"
        with PgWriter(path) as w:
            w.write_batch([make_pg_record()])

        assert path.exists()
        table = pq.read_table(path)
        assert table.num_rows == 1

    def test_schema_matches(self, tmp_path):
        """Written file schema matches PgSchema."""
        path = tmp_path / "test.pg.parquet"
        with PgWriter(path) as w:
            w.write_batch([make_pg_record()])

        pf = pq.ParquetFile(path)
        file_schema = pf.schema_arrow
        expected = PgSchema.get_arrow_schema()
        for field in expected:
            assert field.name in file_schema.names

    def test_footer_metadata(self, tmp_path):
        """PG file footer contains correct file_type."""
        path = tmp_path / "test.pg.parquet"
        with PgWriter(path) as w:
            w.write_batch([make_pg_record()])

        meta = read_parquet_metadata(path)
        assert meta["file_type"] == "pg_file"


class TestSampleWriter:
    """Test SampleWriter creates valid Parquet files."""

    def test_creates_valid_parquet(self, tmp_path):
        """SampleWriter produces a readable Parquet file."""
        path = tmp_path / "test.sample.parquet"
        with SampleWriter(path) as w:
            w.write_batch([make_sample_record()])

        assert path.exists()
        table = pq.read_table(path)
        assert table.num_rows == 1

    def test_schema_matches(self, tmp_path):
        """Written file schema matches SampleSchema."""
        path = tmp_path / "test.sample.parquet"
        with SampleWriter(path) as w:
            w.write_batch([make_sample_record()])

        pf = pq.ParquetFile(path)
        file_schema = pf.schema_arrow
        expected = SampleSchema.get_arrow_schema()
        for field in expected:
            assert field.name in file_schema.names


class TestRunWriter:
    """Test RunWriter creates valid Parquet files."""

    def test_creates_valid_parquet(self, tmp_path):
        """RunWriter produces a readable Parquet file."""
        path = tmp_path / "test.run.parquet"
        with RunWriter(path) as w:
            w.write_batch([make_run_record()])

        assert path.exists()
        table = pq.read_table(path)
        assert table.num_rows == 1

    def test_schema_matches(self, tmp_path):
        """Written file schema matches RunSchema."""
        path = tmp_path / "test.run.parquet"
        with RunWriter(path) as w:
            w.write_batch([make_run_record()])

        pf = pq.ParquetFile(path)
        file_schema = pf.schema_arrow
        expected = RunSchema.get_arrow_schema()
        for field in expected:
            assert field.name in file_schema.names


class TestDatasetWriter:
    """Test DatasetWriter creates valid Parquet files."""

    def test_creates_valid_parquet(self, tmp_path):
        """DatasetWriter produces a readable Parquet file."""
        path = tmp_path / "test.dataset.parquet"
        with DatasetWriter(path) as w:
            w.write_batch([make_dataset_record()])

        assert path.exists()
        table = pq.read_table(path)
        assert table.num_rows == 1

    def test_footer_metadata(self, tmp_path):
        """Dataset file footer contains correct file_type."""
        path = tmp_path / "test.dataset.parquet"
        with DatasetWriter(path) as w:
            w.write_batch([make_dataset_record()])

        meta = read_parquet_metadata(path)
        assert meta["file_type"] == "dataset_file"


class TestOntologyWriter:
    """Test OntologyWriter creates valid Parquet files."""

    def test_creates_valid_parquet(self, tmp_path):
        """OntologyWriter produces a readable Parquet file."""
        path = tmp_path / "test.ontology.parquet"
        with OntologyWriter(path) as w:
            w.write_batch([make_ontology_record()])

        assert path.exists()
        table = pq.read_table(path)
        assert table.num_rows == 1


class TestProvenanceWriter:
    """Test ProvenanceWriter creates valid Parquet files."""

    def test_creates_valid_parquet(self, tmp_path):
        """ProvenanceWriter produces a readable Parquet file."""
        path = tmp_path / "test.provenance.parquet"
        with ProvenanceWriter(path) as w:
            w.write_batch([make_provenance_record()])

        assert path.exists()
        table = pq.read_table(path)
        assert table.num_rows == 1


# ---------------------------------------------------------------------------
# Writer metadata tests (footer, not row data)
# ---------------------------------------------------------------------------


class TestWriterMetadata:
    """Test that writers store metadata in Parquet footer (not in row data)."""

    def test_metadata_in_footer_not_columns(self, feature_parquet):
        """QPX metadata should be in footer, not as columns in the data."""
        table = pq.read_table(feature_parquet)
        data_columns = table.schema.names
        # These should NOT be data columns
        assert "qpx_version" not in data_columns
        assert "creator" not in data_columns
        assert "file_type" not in data_columns
        assert "uuid" not in data_columns

    def test_footer_metadata_is_accessible(self, feature_parquet):
        """Footer metadata is accessible via schema metadata."""
        meta = read_parquet_metadata(feature_parquet)
        assert "qpx_version" in meta
        assert "file_type" in meta
        assert "compression_format" in meta

    def test_software_provider_in_metadata(self, tmp_path):
        """software_provider is stored in footer metadata."""
        path = tmp_path / "test.feature.parquet"
        with FeatureWriter(path, software_provider="my_tool") as w:
            w.write_batch([make_feature_record()])

        meta = read_parquet_metadata(path)
        assert meta["software_provider"] == "my_tool"

    def test_scan_format_in_metadata(self, tmp_path):
        """scan_format is stored in footer metadata when provided."""
        path = tmp_path / "test.feature.parquet"
        with FeatureWriter(path, scan_format="native") as w:
            w.write_batch([make_feature_record()])

        meta = read_parquet_metadata(path)
        assert meta["scan_format"] == "native"


# ---------------------------------------------------------------------------
# Writer schema validation tests
# ---------------------------------------------------------------------------


class TestWriterValidation:
    """Test that writers validate schema on write_table."""

    def test_write_table_rejects_bad_types(self, tmp_path):
        """write_table should raise ValueError for type mismatches."""
        path = tmp_path / "test.feature.parquet"
        writer = FeatureWriter(path)
        try:
            # Create a table with 'charge' as string instead of int16
            schema = FeatureSchema.get_arrow_schema()
            wrong_fields = []
            for f in schema:
                if f.name == "charge":
                    wrong_fields.append(pa.field("charge", pa.string()))
                else:
                    wrong_fields.append(f)
            wrong_schema = pa.schema(wrong_fields)
            arrays = {f.name: pa.nulls(1, type=f.type) for f in wrong_schema}
            # Override charge with a string
            arrays["charge"] = pa.array(["wrong"], type=pa.string())
            table = pa.table(arrays, schema=wrong_schema)

            with pytest.raises(ValueError, match="Schema validation failed"):
                writer.write_table(table)
        finally:
            writer.close()

    def test_write_table_accepts_valid_table(self, tmp_path):
        """write_table accepts a correctly typed Arrow table built from records."""
        path = tmp_path / "test.feature.parquet"
        writer = FeatureWriter(path)
        schema = writer.arrow_schema
        # Use from_pylist to build a valid record batch (handles nullability correctly)
        record = make_feature_record()
        batch = pa.RecordBatch.from_pylist([record], schema=schema)
        table = pa.Table.from_batches([batch])

        with writer:
            writer.write_table(table)

        assert parquet_row_count(path) == 1


# ---------------------------------------------------------------------------
# Writer batching tests
# ---------------------------------------------------------------------------


class TestWriterBatching:
    """Test writer batching behavior (buffer fills and flushes)."""

    def test_batch_size_flush(self, tmp_path):
        """Records flush when batch_size is reached."""
        path = tmp_path / "test.feature.parquet"
        # Use small batch size to force multiple flushes
        with FeatureWriter(path, batch_size=2) as w:
            records = [make_feature_record(sequence=f"SEQ{i}") for i in range(5)]
            w.write_batch(records)

        assert parquet_row_count(path) == 5

    def test_remaining_buffer_flushed_on_close(self, tmp_path):
        """Remaining buffer records are flushed when writer is closed."""
        path = tmp_path / "test.feature.parquet"
        with FeatureWriter(path, batch_size=100) as w:
            # Write fewer than batch_size records
            records = [make_feature_record(sequence=f"SEQ{i}") for i in range(3)]
            w.write_batch(records)
            # At this point, records are still in buffer (3 < 100)
            # They should be flushed on close

        assert parquet_row_count(path) == 3

    def test_multiple_write_batch_calls(self, tmp_path):
        """Multiple write_batch calls accumulate correctly."""
        path = tmp_path / "test.feature.parquet"
        with FeatureWriter(path) as w:
            w.write_batch([make_feature_record(sequence="SEQ1")])
            w.write_batch([make_feature_record(sequence="SEQ2")])
            w.write_batch([make_feature_record(sequence="SEQ3")])

        assert parquet_row_count(path) == 3


# ---------------------------------------------------------------------------
# Writer context manager tests
# ---------------------------------------------------------------------------


class TestWriterContextManager:
    """Test writer context manager behavior."""

    def test_context_manager_writes_file(self, tmp_path):
        """Using writer as context manager produces a valid file."""
        path = tmp_path / "test.feature.parquet"
        with FeatureWriter(path) as w:
            w.write_batch([make_feature_record()])

        assert path.exists()
        assert parquet_row_count(path) == 1

    def test_context_manager_closes_properly(self, tmp_path):
        """After context manager exit, writer is closed (no open writer)."""
        path = tmp_path / "test.feature.parquet"
        writer = FeatureWriter(path)
        with writer:
            writer.write_batch([make_feature_record()])

        # After close, internal writer should be None
        assert writer._writer is None

    def test_manual_close(self, tmp_path):
        """Manual close() works correctly."""
        path = tmp_path / "test.feature.parquet"
        writer = FeatureWriter(path)
        writer.write_batch([make_feature_record()])
        writer.close()

        assert path.exists()
        assert parquet_row_count(path) == 1
        assert writer._writer is None


# ---------------------------------------------------------------------------
# Writer compression tests
# ---------------------------------------------------------------------------


class TestWriterCompression:
    """Test writer compression options."""

    def test_default_compression_zstd(self, tmp_path):
        """Default compression is zstd."""
        path = tmp_path / "test.feature.parquet"
        with FeatureWriter(path) as w:
            w.write_batch([make_feature_record()])

        meta = read_parquet_metadata(path)
        assert meta["compression_format"] == "zstd"

    def test_custom_compression(self, tmp_path):
        """Custom compression is recorded in footer."""
        path = tmp_path / "test.feature.parquet"
        with FeatureWriter(path, compression="snappy") as w:
            w.write_batch([make_feature_record()])

        meta = read_parquet_metadata(path)
        assert meta["compression_format"] == "snappy"
