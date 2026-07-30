"""Regression tests for qpx.mudata._attach_uns_metadata.

Previously _attach_uns_metadata only filtered None and float NaN; it missed
pd.NA (from pandas nullable dtypes when reading parquet). That caused
mdata.write() to fail on real data with:
    IORegistryError: No method registered for writing <class 'pandas._libs.missing.NAType'>

These tests ensure all scalar NA flavors are skipped while normal values
(including non-scalar dicts/lists) are preserved.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

mudata = pytest.importorskip("mudata")
anndata = pytest.importorskip("anndata")

from qpx.converters.orchestrator import BaseOrchestrator  # noqa: E402
from qpx.dataset import Dataset  # noqa: E402
from qpx.mudata import (  # noqa: E402
    _LABEL_FIELD_QUERIES,
    _attach_uns_metadata,
    _detect_intensity_label,
    build_mudata,
)
from qpx.writers import FeatureWriter, PgWriter, RunWriter  # noqa: E402
from tests.conftest import make_feature_record, make_pg_record, make_run_record  # noqa: E402


class _FakeResult:
    def __init__(self, df: pd.DataFrame) -> None:
        self._df = df

    def fetchdf(self) -> pd.DataFrame:
        return self._df


class _FakeEngine:
    """Minimal stand-in for qpx.core.engine.DuckDBEngine used by _attach_uns_metadata.

    _attach_uns_metadata only calls engine.execute(sql).fetchdf(), so this shim
    is sufficient and avoids spinning up a real DuckDB + parquet fixture.
    """

    def __init__(self, df: pd.DataFrame) -> None:
        self._df = df

    def execute(self, sql: str) -> _FakeResult:  # noqa: ARG002
        return _FakeResult(self._df)


class _FakeRowResult:
    """Returns a single-row fetchone() result."""

    def __init__(self, value) -> None:
        self._value = value

    def fetchone(self):
        return (self._value,) if self._value is not None else None

    def fetchdf(self) -> pd.DataFrame:
        return pd.DataFrame()


class _TableAwareEngine:
    """Fake engine that dispatches based on SQL table references."""

    def __init__(self, table_labels: dict[str, str]) -> None:
        self._table_labels = table_labels

    def execute(self, sql: str, params=None) -> _FakeRowResult:  # noqa: ARG002
        for table, label in self._table_labels.items():
            if table in sql:
                if "typeof" in sql.lower():
                    return _FakeRowResult('STRUCT("label" VARCHAR, intensity FLOAT)[]')
                return _FakeRowResult(label)
        return _FakeRowResult(None)


def _make_empty_mdata() -> "mudata.MuData":
    return mudata.MuData({"precursors": anndata.AnnData(X=np.zeros((1, 1)))})


def test_attach_uns_metadata_skips_all_na_flavors():
    """All scalar NA types must be filtered out of mdata.uns.

    Regression for IORegistryError on pd.NA when writing .h5mu.
    """
    df = pd.DataFrame(
        {
            "software_name": ["qpx"],
            "software_version": ["0.1.0"],
            "file_checksums": pd.array([pd.NA], dtype="object"),
            "file_row_counts": [None],
            "file_sizes_bytes": [np.nan],
            "creation_date": [pd.NaT],
            "total_structures": [5],
        }
    )
    mdata = _make_empty_mdata()

    _attach_uns_metadata(_FakeEngine(df), mdata)

    assert "file_checksums" not in mdata.uns
    assert "file_row_counts" not in mdata.uns
    assert "file_sizes_bytes" not in mdata.uns
    assert "creation_date" not in mdata.uns

    assert mdata.uns["software_name"] == "qpx"
    assert mdata.uns["software_version"] == "0.1.0"
    assert mdata.uns["total_structures"] == 5


def test_attach_uns_metadata_preserves_non_scalar_values():
    """Dicts and lists must not be misclassified as NA by pd.isna()."""
    df = pd.DataFrame(
        {
            "tags": [["a", "b"]],
            "config": [{"key": "value"}],
        }
    )
    mdata = _make_empty_mdata()

    _attach_uns_metadata(_FakeEngine(df), mdata)

    assert mdata.uns["tags"] == ["a", "b"]
    assert mdata.uns["config"] == {"key": "value"}


def test_attach_uns_metadata_allows_hdf5_write(tmp_path):
    """End-to-end: mdata.write() must succeed when dataset row contains pd.NA.

    This is the symptom we originally observed on real DIA-NN conversion output.
    """
    df = pd.DataFrame(
        {
            "software_name": ["qpx"],
            "file_checksums": pd.array([pd.NA], dtype="object"),
        }
    )
    mdata = _make_empty_mdata()
    _attach_uns_metadata(_FakeEngine(df), mdata)

    out = tmp_path / "roundtrip.h5mu"
    mdata.write(out)
    assert out.exists()
    assert out.stat().st_size > 0


def test_label_field_query_generates_table_specific_sql():
    """_LABEL_FIELD_QUERIES must embed the table name, not hardcode 'feature'."""
    sql_feat = _LABEL_FIELD_QUERIES[("label", "feature")]
    sql_pg = _LABEL_FIELD_QUERIES[("label", "pg")]
    assert "feature" in sql_feat and "pg" not in sql_feat
    assert "pg" in sql_pg and "feature" not in sql_pg


def test_detect_intensity_label_per_table():
    """DIA-NN feature uses 'raw', pg uses 'LFQ' — each must be detected independently.

    Regression: previously _detect_intensity_label only inspected the feature table,
    causing _build_protein_adata to query pg WHERE label='raw' and get zero rows.
    """
    engine = _TableAwareEngine({"feature": "raw", "pg": "LFQ"})
    assert _detect_intensity_label(engine, "feature") == "raw"
    assert _detect_intensity_label(engine, "pg") == "LFQ"


def _write_quant_bundle(tmp_path, prefix, feature_intensities, protein_intensities):
    with FeatureWriter(tmp_path / f"{prefix}.feature.parquet") as writer:
        writer.write_batch([make_feature_record(intensities=feature_intensities)])
    with PgWriter(tmp_path / f"{prefix}.pg.parquet") as writer:
        writer.write_batch([make_pg_record(intensities=protein_intensities)])

    run = make_run_record()
    run["samples"] = [
        {
            "sample_accession": f"{prefix}_{entry['label']}",
            "label": entry["label"],
            "biological_replicate": index,
            "technical_replicate": 1,
        }
        for index, entry in enumerate(feature_intensities, start=1)
    ]
    with RunWriter(tmp_path / f"{prefix}.run.parquet") as writer:
        writer.write_batch([run])


def test_orchestrator_mudata_exports_all_channels_from_requested_prefix(tmp_path):
    """Automatic MuData export must include every channel and never mix prefixes."""
    _write_quant_bundle(
        tmp_path,
        "a",
        [{"label": "TMT126", "intensity": 900.0}],
        [{"label": "TMT126", "intensity": 9000.0}],
    )
    _write_quant_bundle(
        tmp_path,
        "z",
        [
            {"label": "TMT126", "intensity": 100.0},
            {"label": "TMT127N", "intensity": 200.0},
        ],
        [
            {"label": "TMT126", "intensity": 1000.0},
            {"label": "TMT127N", "intensity": 2000.0},
        ],
    )

    output = BaseOrchestrator()._write_mudata(tmp_path, "z")
    assert output == tmp_path / "z.h5mu"

    mdata = mudata.read_h5mu(output)
    precursor = mdata.mod["precursors"]
    protein = mdata.mod["proteins"]
    assert list(precursor.obs_names) == ["run_01|TMT126", "run_01|TMT127N"]
    assert list(precursor.obs["sample_accession"]) == ["z_TMT126", "z_TMT127N"]
    np.testing.assert_allclose(precursor.X.toarray()[:, 0], [100.0, 200.0])
    np.testing.assert_allclose(protein.X.toarray()[:, 0], [1000.0, 2000.0])


def test_orchestrator_mudata_expected_failure_is_best_effort(tmp_path, monkeypatch):
    """Expected assembly failures must not fail the primary conversion."""
    closed = []

    class FakeDataset:
        def __init__(self, path, file_prefix=None):  # noqa: ARG002
            pass

        def close(self):
            closed.append(True)

    def fail_build(dataset, all_intensity_labels=False):  # noqa: ARG001
        raise ValueError("invalid MuData input")

    monkeypatch.setattr("qpx.dataset.Dataset", FakeDataset)
    monkeypatch.setattr("qpx.mudata.build_mudata", fail_build)

    output = BaseOrchestrator()._write_mudata(tmp_path, "broken")

    assert output is None
    assert closed == [True]


def test_build_mudata_explicit_label_keeps_run_level_contract(tmp_path):
    """Explicit single-label export keeps run observations and matches sample metadata."""
    _write_quant_bundle(
        tmp_path,
        "z",
        [
            {"label": "TMT126", "intensity": 100.0},
            {"label": "TMT127N", "intensity": 200.0},
        ],
        [
            {"label": "TMT126", "intensity": 1000.0},
            {"label": "TMT127N", "intensity": 2000.0},
        ],
    )

    dataset = Dataset(tmp_path, file_prefix="z")
    try:
        mdata = build_mudata(
            dataset,
            intensity_label="TMT127N",
            modalities=["precursors"],
        )
    finally:
        dataset.close()

    precursor = mdata.mod["precursors"]
    assert list(precursor.obs_names) == ["run_01"]
    assert list(precursor.obs["sample_accession"]) == ["z_TMT127N"]
    np.testing.assert_allclose(precursor.X.toarray()[:, 0], [200.0])
