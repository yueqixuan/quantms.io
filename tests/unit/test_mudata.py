"""Tests for MuData export (qpx.mudata)."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from scipy import sparse

ad = pytest.importorskip("anndata")
mu = pytest.importorskip("mudata")

from qpx.dataset import Dataset  # noqa: E402, I001
from qpx.mudata import (  # noqa: E402
    _build_differential_adata,
    _build_expression_adata,
    _build_feature_mapping,
    _build_precursor_adata,
    _build_protein_adata,
    _detect_intensity_label,
    _pivot_to_sparse,
    build_mudata,
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def ds(dataset_dir):
    """Open a Dataset from the shared dataset_dir fixture."""
    return Dataset(dataset_dir)


@pytest.fixture
def engine(ds):
    """Expose the DuckDB engine from the dataset."""
    return ds._engine


# ---------------------------------------------------------------------------
# TestPivotToSparse
# ---------------------------------------------------------------------------


class TestPivotToSparse:
    def test_basic_pivot(self):
        df = pd.DataFrame({
            "row": ["r1", "r1", "r2"],
            "col": ["c1", "c2", "c1"],
            "val": [1.0, 2.0, 3.0],
        })
        rows = pd.Index(["r1", "r2"])
        cols = pd.Index(["c1", "c2"])
        mat = _pivot_to_sparse(df, "row", "col", "val", rows, cols)

        assert sparse.issparse(mat)
        assert mat.shape == (2, 2)
        assert mat[0, 0] == pytest.approx(1.0)
        assert mat[0, 1] == pytest.approx(2.0)
        assert mat[1, 0] == pytest.approx(3.0)
        assert mat[1, 1] == 0.0


# ---------------------------------------------------------------------------
# TestBuildPrecursorAdata
# ---------------------------------------------------------------------------


class TestBuildPrecursorAdata:
    def test_shape(self, engine):
        adata = _build_precursor_adata(engine, "TMT126")
        # dataset_dir: PEPTIDEK|2 and ANOTHERK|3 in run_01 with TMT126
        # THIRDPEP|2 in run_02 with TMT127 -- NOT included for TMT126
        assert adata.n_obs == 1  # only run_01 has TMT126
        assert adata.n_vars == 2  # PEPTIDEK|2, ANOTHERK|3

    def test_index_format(self, engine):
        adata = _build_precursor_adata(engine, "TMT126")
        for name in adata.var_names:
            assert "|" in name, f"Expected pipe-separated precursor ID, got {name}"

    def test_sparse_x(self, engine):
        adata = _build_precursor_adata(engine, "TMT126")
        assert sparse.issparse(adata.X)

    def test_values(self, engine):
        adata = _build_precursor_adata(engine, "TMT126")
        # Both PEPTIDEK|2 and ANOTHERK|3 have intensity 1000.0 in TMT126
        dense = adata.X.toarray()
        assert np.any(dense > 0)


# ---------------------------------------------------------------------------
# TestBuildProteinAdata
# ---------------------------------------------------------------------------


class TestBuildProteinAdata:
    def test_shape(self, engine):
        adata = _build_protein_adata(engine, "TMT126")
        # pg: P12345 in run_01 with TMT126
        assert adata.n_obs == 1
        assert adata.n_vars == 1  # P12345

    def test_protein_in_var(self, engine):
        adata = _build_protein_adata(engine, "TMT126")
        assert "P12345" in adata.var_names

    def test_gene_name_in_var(self, engine):
        adata = _build_protein_adata(engine, "TMT126")
        assert "gene_name" in adata.var.columns
        # The pg fixture has gg_names=["GENE1"]
        assert adata.var.loc["P12345", "gene_name"] == "GENE1"

    def test_sparse_x(self, engine):
        adata = _build_protein_adata(engine, "TMT126")
        assert sparse.issparse(adata.X)


# ---------------------------------------------------------------------------
# TestBuildFeatureMapping
# ---------------------------------------------------------------------------


class TestBuildFeatureMapping:
    def test_square_symmetric(self, engine):
        # Build precursor and protein modalities first
        prec = _build_precursor_adata(engine, "TMT126")
        prot = _build_protein_adata(engine, "TMT126")
        mdata = mu.MuData({"precursors": prec, "proteins": prot})

        mat = _build_feature_mapping(engine, mdata)
        assert mat.shape[0] == mat.shape[1]
        # Symmetric: mat == mat.T
        diff = mat - mat.T
        assert diff.nnz == 0

    def test_has_edges(self, engine):
        prec = _build_precursor_adata(engine, "TMT126")
        prot = _build_protein_adata(engine, "TMT126")
        mdata = mu.MuData({"precursors": prec, "proteins": prot})

        mat = _build_feature_mapping(engine, mdata)
        assert mat.nnz > 0

    def test_correct_shape(self, engine):
        prec = _build_precursor_adata(engine, "TMT126")
        prot = _build_protein_adata(engine, "TMT126")
        mdata = mu.MuData({"precursors": prec, "proteins": prot})

        mat = _build_feature_mapping(engine, mdata)
        expected_size = prec.n_vars + prot.n_vars
        assert mat.shape == (expected_size, expected_size)


# ---------------------------------------------------------------------------
# TestBuildExpressionAdata
# ---------------------------------------------------------------------------


class TestBuildExpressionAdata:
    def test_returns_none_when_missing(self, dataset_dir):
        result = _build_expression_adata(dataset_dir)
        assert result is None

    def test_loads_h5ad(self, dataset_dir):
        # Create a dummy .ae.h5ad
        dummy = ad.AnnData(
            X=np.array([[1.0, 2.0], [3.0, 4.0]]),
            obs=pd.DataFrame(index=["s1", "s2"]),
            var=pd.DataFrame(index=["g1", "g2"]),
        )
        dummy.write_h5ad(dataset_dir / "test.ae.h5ad")

        result = _build_expression_adata(dataset_dir)
        assert result is not None
        assert result.n_obs == 2
        assert result.n_vars == 2


# ---------------------------------------------------------------------------
# TestBuildDifferentialAdata
# ---------------------------------------------------------------------------


class TestBuildDifferentialAdata:
    def test_returns_none_when_missing(self, dataset_dir):
        result = _build_differential_adata(dataset_dir)
        assert result is None

    def test_loads_and_restructures(self, dataset_dir):
        # Create a DE h5ad with uns["de_results"]
        de_df = pd.DataFrame({
            "contrast": ["A_vs_B", "A_vs_B", "C_vs_D", "C_vs_D"],
            "protein": ["P1", "P2", "P1", "P2"],
            "log2FC": [1.5, -0.5, 2.0, -1.0],
            "pvals_adj": [0.01, 0.05, 0.001, 0.1],
        })
        dummy = ad.AnnData(
            X=np.zeros((1, 1)),
            obs=pd.DataFrame(index=["dummy"]),
            var=pd.DataFrame(index=["dummy_var"]),
        )
        dummy.uns["de_results"] = de_df
        dummy.write_h5ad(dataset_dir / "test.de.h5ad")

        result = _build_differential_adata(dataset_dir)
        assert result is not None
        # Should be restructured: obs=contrasts, var=proteins
        assert result.n_obs == 2  # A_vs_B, C_vs_D
        assert result.n_vars == 2  # P1, P2
        assert "pvals_adj" in result.layers


# ---------------------------------------------------------------------------
# TestDetectIntensityLabel
# ---------------------------------------------------------------------------


class TestDetectIntensityLabel:
    def test_detects_label(self, engine):
        label = _detect_intensity_label(engine)
        assert label in ("TMT126", "TMT127")


# ---------------------------------------------------------------------------
# TestBuildMudata
# ---------------------------------------------------------------------------


class TestBuildMudata:
    def test_core_modalities_present(self, ds):
        mdata = build_mudata(ds, intensity_label="TMT126")
        assert "precursors" in mdata.mod
        assert "proteins" in mdata.mod

    def test_no_expression_without_file(self, ds):
        mdata = build_mudata(ds, intensity_label="TMT126")
        assert "expression" not in mdata.mod

    def test_with_expression(self, ds, dataset_dir):
        # Create an .ae.h5ad file
        dummy = ad.AnnData(
            X=np.array([[1.0]]),
            obs=pd.DataFrame(index=["s1"]),
            var=pd.DataFrame(index=["g1"]),
        )
        dummy.write_h5ad(dataset_dir / "test.ae.h5ad")

        mdata = build_mudata(ds, intensity_label="TMT126")
        assert "expression" in mdata.mod

    def test_modalities_filter(self, ds):
        mdata = build_mudata(ds, intensity_label="TMT126", modalities=["proteins"])
        assert "proteins" in mdata.mod
        assert "precursors" not in mdata.mod

    def test_uns_metadata(self, ds):
        mdata = build_mudata(ds, intensity_label="TMT126")
        assert "project_accession" in mdata.uns
        assert mdata.uns["project_accession"] == "PXD000001"

    def test_invalid_modality(self, ds):
        with pytest.raises(ValueError, match="Unknown modalities"):
            build_mudata(ds, intensity_label="TMT126", modalities=["invalid"])

    def test_auto_detect_label(self, ds):
        mdata = build_mudata(ds)
        # Should auto-detect and produce at least one modality
        assert len(mdata.mod) >= 1


# ---------------------------------------------------------------------------
# TestDatasetToMudata
# ---------------------------------------------------------------------------


class TestDatasetToMudata:
    def test_method_works(self, ds):
        mdata = ds.to_mudata(intensity_label="TMT126")
        assert isinstance(mdata, mu.MuData)
        assert "precursors" in mdata.mod

    def test_import_error(self, ds, monkeypatch):
        import builtins

        real_import = builtins.__import__

        def mock_import(name, *args, **kwargs):
            if name == "mudata":
                raise ImportError("no mudata")
            return real_import(name, *args, **kwargs)

        monkeypatch.setattr(builtins, "__import__", mock_import)
        with pytest.raises(ImportError, match="mudata is required"):
            ds.to_mudata(intensity_label="TMT126")
