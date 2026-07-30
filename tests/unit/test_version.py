"""Tests for QPX spec-version stamping and the pre-1.1 PG load guard."""

import pyarrow as pa
import pyarrow.parquet as pq
import pytest

import qpx
from qpx.core.parquet_io import read_parquet_metadata
from qpx.version import (
    QPX_SPEC_VERSION,
    QpxVersionError,
    check_pg_file_compatible,
    parse_spec_version,
)
from qpx.writers import PgWriter
from tests.conftest import make_pg_record


def test_spec_version_is_the_format_version():
    """The spec version is the two-component format version, distinct from the build."""
    assert QPX_SPEC_VERSION == "1.1"
    assert parse_spec_version(QPX_SPEC_VERSION) == (1, 1)
    assert parse_spec_version("1.0") < parse_spec_version("1.1")


def test_written_pg_footer_carries_spec_version(tmp_path):
    """A freshly written PG file stamps qpx_version == the spec version (1.1)."""
    path = tmp_path / "test.pg.parquet"
    with PgWriter(path) as w:
        w.write_batch([make_pg_record()])

    meta = read_parquet_metadata(path)
    assert meta["qpx_version"] == "1.1"
    # writer_version (package build) is a separate provenance key.
    assert "writer_version" in meta
    # A well-formed 1.1 file passes the guard unchanged.
    check_pg_file_compatible(path)


def _write_old_pg_file(path, qpx_version="1.0"):
    """Synthesize a pre-1.1 PG file: scalar run_file_name, no grouped_runs."""
    table = pa.table({"anchor_protein": ["P12345"], "run_file_name": ["run_01"], "global_qvalue": [0.005]})
    table = table.replace_schema_metadata({b"qpx_version": qpx_version.encode()})
    pq.write_table(table, str(path))


def test_guard_raises_on_old_schema_pg_file(tmp_path):
    """A pre-1.1 PG file (scalar run_file_name) raises a clear version error."""
    path = tmp_path / "old.pg.parquet"
    _write_old_pg_file(path, qpx_version="1.0")

    with pytest.raises(QpxVersionError) as exc:
        check_pg_file_compatible(path)

    msg = str(exc.value)
    assert "grouped_runs" in msg
    assert "run_file_name" in msg
    assert "1.0" in msg  # names the offending file version
    assert "1.1" in msg  # names the required version


def test_read_pg_raises_on_old_schema_file(tmp_path):
    """The public read_pg() load path surfaces the clear error, not a raw DuckDB one."""
    path = tmp_path / "old.pg.parquet"
    _write_old_pg_file(path, qpx_version="1.0")

    with pytest.raises(QpxVersionError):
        qpx.read_pg(str(path))


def test_guard_noop_on_unstamped_but_grouped_file(tmp_path):
    """A file with grouped_runs but no qpx_version stamp is allowed (forward-only)."""
    path = tmp_path / "grouped.pg.parquet"
    table = pa.table({"anchor_protein": ["P1"], "grouped_runs": [["run_01"]]})
    pq.write_table(table, str(path))
    check_pg_file_compatible(path)  # must not raise
