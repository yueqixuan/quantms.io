"""Tests for QPX core layer: DuckDBEngine, LazyQuery, QueryResult, Parquet I/O."""

import pandas as pd
import pyarrow as pa
import pytest

from qpx.core.engine import DuckDBEngine, create_engine
from qpx.core.query import LazyQuery
from qpx.core.convert import QueryResult
from qpx.core.parquet_io import (
    read_parquet_metadata,
    read_parquet_schema,
    parquet_row_count,
)

# ---------------------------------------------------------------------------
# DuckDBEngine tests
# ---------------------------------------------------------------------------


class TestDuckDBEngine:
    """Test DuckDB connection management and configuration."""

    def test_creates_connection(self):
        """Engine creates a working DuckDB connection."""
        engine = DuckDBEngine()
        try:
            result = engine.execute("SELECT 42 AS answer")
            row = result.fetchone()
            assert row[0] == 42
        finally:
            engine.close()

    def test_custom_memory_and_threads(self):
        """Engine accepts custom memory_limit and threads."""
        engine = DuckDBEngine(memory_limit="1GB", threads=2)
        try:
            # Just verify it works — DuckDB accepts the settings
            result = engine.execute("SELECT current_setting('threads') AS t")
            row = result.fetchone()
            # DuckDB returns thread count as string
            assert int(row[0]) == 2
        finally:
            engine.close()

    def test_register_parquet_creates_view(self, feature_parquet):
        """register_parquet creates a queryable DuckDB view."""
        engine = DuckDBEngine()
        try:
            engine.register_parquet("feature", feature_parquet)
            result = engine.execute("SELECT COUNT(*) FROM feature")
            count = result.fetchone()[0]
            assert count == 3
        finally:
            engine.close()

    def test_register_parquet_replaces_existing(self, feature_parquet, tmp_path):
        """register_parquet with same name replaces the existing view."""
        from tests.conftest import make_feature_record
        from qpx.writers import FeatureWriter

        # Create a second file with 1 row
        path2 = tmp_path / "other.feature.parquet"
        with FeatureWriter(path2) as w:
            w.write_batch([make_feature_record()])

        engine = DuckDBEngine()
        try:
            engine.register_parquet("feature", feature_parquet)
            count1 = engine.execute("SELECT COUNT(*) FROM feature").fetchone()[0]
            assert count1 == 3

            engine.register_parquet("feature", path2)
            count2 = engine.execute("SELECT COUNT(*) FROM feature").fetchone()[0]
            assert count2 == 1
        finally:
            engine.close()

    def test_context_manager(self):
        """Engine works as a context manager."""
        with DuckDBEngine() as engine:
            result = engine.execute("SELECT 1")
            assert result.fetchone()[0] == 1
        # After exiting, connection should be closed

    def test_create_engine_factory(self):
        """create_engine() factory returns a DuckDBEngine."""
        engine = create_engine(memory_limit="2GB", threads=1)
        try:
            assert isinstance(engine, DuckDBEngine)
            result = engine.execute("SELECT 'ok'")
            assert result.fetchone()[0] == "ok"
        finally:
            engine.close()

    def test_execute_with_params(self):
        """Engine.execute supports parameterized queries."""
        engine = DuckDBEngine()
        try:
            result = engine.execute("SELECT ? + ? AS sum", [10, 20])
            assert result.fetchone()[0] == 30
        finally:
            engine.close()


# ---------------------------------------------------------------------------
# LazyQuery tests
# ---------------------------------------------------------------------------


class TestLazyQuery:
    """Test the immutable, chainable LazyQuery SQL builder."""

    @pytest.fixture
    def engine_with_feature(self, feature_parquet):
        """Provide an engine with a feature parquet registered."""
        engine = DuckDBEngine()
        engine.register_parquet("feature", feature_parquet)
        yield engine
        engine.close()

    def test_build_sql_default(self, engine_with_feature):
        """Default LazyQuery builds 'SELECT * FROM table_name'."""
        q = LazyQuery(engine_with_feature, "feature")
        assert q.build_sql() == "SELECT * FROM feature"

    def test_filter_builds_where_clause(self, engine_with_feature):
        """filter() builds correct SQL with WHERE clause."""
        q = LazyQuery(engine_with_feature, "feature")
        filtered = q.filter("charge > 2")
        sql = filtered.build_sql()
        assert "WHERE" in sql
        assert "charge > 2" in sql

    def test_select_builds_column_list(self, engine_with_feature):
        """select() builds correct SQL with column list."""
        q = LazyQuery(engine_with_feature, "feature")
        selected = q.select("sequence", "charge")
        sql = selected.build_sql()
        assert "sequence" in sql
        assert "charge" in sql

    def test_join_builds_join_clause(self, engine_with_feature, run_parquet):
        """join() builds correct SQL with JOIN clause."""
        engine_with_feature.register_parquet("run", run_parquet)
        q1 = LazyQuery(engine_with_feature, "feature")
        q2 = LazyQuery(engine_with_feature, "run")
        joined = q1.join(q2, on="run_file_name")
        sql = joined.build_sql()
        assert "JOIN" in sql
        assert "run_file_name" in sql

    def test_limit_builds_limit_clause(self, engine_with_feature):
        """limit() builds correct SQL with LIMIT clause."""
        q = LazyQuery(engine_with_feature, "feature")
        limited = q.limit(5)
        sql = limited.build_sql()
        assert "LIMIT 5" in sql

    def test_order_by_builds_order_clause(self, engine_with_feature):
        """order_by() builds correct SQL with ORDER BY clause."""
        q = LazyQuery(engine_with_feature, "feature")
        ordered = q.order_by("charge")
        sql = ordered.build_sql()
        assert "ORDER BY" in sql
        assert "charge ASC" in sql

    def test_order_by_desc(self, engine_with_feature):
        """order_by(desc=True) uses DESC direction."""
        q = LazyQuery(engine_with_feature, "feature")
        ordered = q.order_by("charge", desc=True)
        sql = ordered.build_sql()
        assert "charge DESC" in sql

    def test_count_returns_correct_count(self, engine_with_feature):
        """count() returns the correct row count."""
        q = LazyQuery(engine_with_feature, "feature")
        assert q.count() == 3

    def test_count_after_filter(self, engine_with_feature):
        """count() after filter returns filtered count."""
        q = LazyQuery(engine_with_feature, "feature")
        filtered = q.filter("is_decoy = false")
        assert filtered.count() == 2

    def test_chaining_is_immutable(self, engine_with_feature):
        """Chaining returns new instances; original is unchanged."""
        q = LazyQuery(engine_with_feature, "feature")
        original_sql = q.build_sql()

        filtered = q.filter("charge > 2")
        selected = q.select("sequence")
        limited = q.limit(1)

        # Original should be unchanged
        assert q.build_sql() == original_sql
        # Each derived query should be different
        assert filtered.build_sql() != original_sql
        assert selected.build_sql() != original_sql
        assert limited.build_sql() != original_sql
        # Derived queries should be different from each other
        assert filtered.build_sql() != selected.build_sql()

    def test_chained_filter_then_select(self, engine_with_feature):
        """Chaining filter().select() produces nested SQL."""
        q = LazyQuery(engine_with_feature, "feature")
        result = q.filter("charge = 2").select("sequence", "charge")
        sql = result.build_sql()
        assert "WHERE" in sql
        assert "charge = 2" in sql
        assert "sequence" in sql

    def test_execute_returns_results(self, engine_with_feature):
        """execute() returns a DuckDB result object."""
        q = LazyQuery(engine_with_feature, "feature")
        result = q.execute()
        rows = result.fetchall()
        assert len(rows) == 3

    def test_distinct_values(self, engine_with_feature):
        """distinct_values returns unique values for a column."""
        q = LazyQuery(engine_with_feature, "feature")
        runs = q.distinct_values("run_file_name")
        assert set(runs) == {"run_01", "run_02"}

    def test_group_by(self, engine_with_feature):
        """group_by builds correct SQL."""
        q = LazyQuery(engine_with_feature, "feature")
        grouped = q.group_by("run_file_name")
        sql = grouped.build_sql()
        assert "GROUP BY" in sql
        assert "COUNT(*)" in sql

    def test_source_property_without_sql(self, engine_with_feature):
        """source property returns table_name when no SQL is set."""
        q = LazyQuery(engine_with_feature, "feature")
        assert q.source == "feature"

    def test_source_property_with_sql(self, engine_with_feature):
        """source property wraps SQL in parens when SQL is set."""
        q = LazyQuery(engine_with_feature, "feature")
        filtered = q.filter("charge > 2")
        assert filtered.source.startswith("(")
        assert filtered.source.endswith(")")


# ---------------------------------------------------------------------------
# QueryResult tests
# ---------------------------------------------------------------------------


class TestQueryResult:
    """Test QueryResult materialization to different formats."""

    def test_to_df_returns_dataframe(self, feature_parquet):
        """to_df() returns a pandas DataFrame."""
        engine = DuckDBEngine()
        engine.register_parquet("feature", feature_parquet)
        raw = engine.execute("SELECT sequence, charge FROM feature")
        result = QueryResult(raw)
        df = result.to_df()
        assert isinstance(df, pd.DataFrame)
        assert len(df) == 3
        assert "sequence" in df.columns
        assert "charge" in df.columns
        engine.close()

    def test_to_arrow_returns_table(self, feature_parquet):
        """to_arrow() returns a PyArrow Table."""
        engine = DuckDBEngine()
        engine.register_parquet("feature", feature_parquet)
        raw = engine.execute("SELECT sequence, charge FROM feature")
        result = QueryResult(raw)
        table = result.to_arrow()
        assert isinstance(table, pa.Table)
        assert table.num_rows == 3
        engine.close()

    def test_fetchone(self, feature_parquet):
        """fetchone() returns a single row tuple."""
        engine = DuckDBEngine()
        engine.register_parquet("feature", feature_parquet)
        raw = engine.execute("SELECT COUNT(*) FROM feature")
        result = QueryResult(raw)
        row = result.fetchone()
        assert row[0] == 3
        engine.close()

    def test_fetchall(self, feature_parquet):
        """fetchall() returns list of row tuples."""
        engine = DuckDBEngine()
        engine.register_parquet("feature", feature_parquet)
        raw = engine.execute("SELECT sequence FROM feature ORDER BY sequence")
        result = QueryResult(raw)
        rows = result.fetchall()
        assert len(rows) == 3
        engine.close()

    def test_iter(self, feature_parquet):
        """QueryResult is iterable."""
        engine = DuckDBEngine()
        engine.register_parquet("feature", feature_parquet)
        raw = engine.execute("SELECT sequence FROM feature")
        result = QueryResult(raw)
        rows = list(result)
        assert len(rows) == 3
        engine.close()


# ---------------------------------------------------------------------------
# Parquet I/O utility tests
# ---------------------------------------------------------------------------


class TestParquetIO:
    """Test Parquet I/O utility functions."""

    def test_read_parquet_metadata(self, feature_parquet):
        """read_parquet_metadata returns footer key-value metadata."""
        meta = read_parquet_metadata(feature_parquet)
        assert isinstance(meta, dict)
        assert "qpx_version" in meta
        assert "file_type" in meta
        assert meta["file_type"] == "feature_file"

    def test_read_parquet_metadata_has_creator(self, feature_parquet):
        """Footer metadata should contain creator."""
        meta = read_parquet_metadata(feature_parquet)
        assert "creator" in meta

    def test_read_parquet_metadata_has_uuid(self, feature_parquet):
        """Footer metadata should contain a UUID."""
        meta = read_parquet_metadata(feature_parquet)
        assert "uuid" in meta
        assert len(meta["uuid"]) > 0

    def test_read_parquet_schema(self, feature_parquet):
        """read_parquet_schema returns a PyArrow schema."""
        schema = read_parquet_schema(feature_parquet)
        assert isinstance(schema, pa.Schema)
        assert "sequence" in schema.names

    def test_parquet_row_count(self, feature_parquet):
        """parquet_row_count returns correct count without reading data."""
        count = parquet_row_count(feature_parquet)
        assert count == 3

    def test_parquet_row_count_single_row(self, dataset_parquet):
        """parquet_row_count works for single-row files."""
        count = parquet_row_count(dataset_parquet)
        assert count == 1
