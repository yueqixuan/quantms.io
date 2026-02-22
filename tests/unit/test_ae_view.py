"""Tests for AbsoluteExpressionView (PCA from .ae.h5ad)."""

import pytest

plotly = pytest.importorskip("plotly", reason="plotly required for plotting tests")

import numpy as np


@pytest.fixture
def ae_dataset_dir(dataset_dir):
    """Extend dataset_dir with a minimal .ae.h5ad file."""
    import anndata as ad
    import pandas as pd

    n_samples = 2
    n_proteins = 3
    X = np.array([[100.0, 200.0, 300.0], [400.0, 500.0, 600.0]], dtype=np.float32)

    obs = pd.DataFrame(
        {
            "organism": ["Homo sapiens", "Homo sapiens"],
            "disease": ["healthy", "cancer"],
        },
        index=["SAMPLE_01", "SAMPLE_02"],
    )
    var = pd.DataFrame(index=["P12345", "P67890", "P11111"])

    adata = ad.AnnData(X=X, obs=obs, var=var)
    adata.layers["ibaq"] = X * 2
    adata.layers["ibaq_log2"] = np.log2(np.where(X > 0, X, np.nan)).astype(np.float32)

    prefix = dataset_dir.name
    adata.write_h5ad(dataset_dir / f"{prefix}.ae.h5ad")
    return dataset_dir


class TestAbsoluteExpressionView:
    @pytest.mark.timeout(120)
    def test_pca_basic(self, ae_dataset_dir):
        from qpx import Dataset

        ds = Dataset(ae_dataset_dir)
        fig = ds.ae_view.pca()
        assert fig is not None
        assert hasattr(fig, "data")
        # Should have at least one trace
        assert len(fig.data) >= 1

    def test_pca_with_color(self, ae_dataset_dir):
        from qpx import Dataset

        ds = Dataset(ae_dataset_dir)
        fig = ds.ae_view.pca(color_by="disease")
        # Should have two traces (healthy, cancer)
        assert len(fig.data) == 2

    def test_pca_with_layer(self, ae_dataset_dir):
        from qpx import Dataset

        ds = Dataset(ae_dataset_dir)
        fig = ds.ae_view.pca(layer="ibaq")
        assert fig is not None

    def test_pca_no_h5ad_raises(self, dataset_dir):
        from qpx import Dataset

        ds = Dataset(dataset_dir)
        with pytest.raises(FileNotFoundError, match="No .ae.h5ad file found"):
            ds.ae_view.pca()

    def test_pca_fillna_median(self, ae_dataset_dir):
        from qpx import Dataset

        ds = Dataset(ae_dataset_dir)
        # ibaq_log2 layer contains NaN values — median fill should work
        fig = ds.ae_view.pca(layer="ibaq_log2", fillna="median")
        assert fig is not None
