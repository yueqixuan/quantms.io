"""Tests for SampleSummaryView."""

import pytest
from qpx import Dataset

_has_plotly = pytest.importorskip is not None  # always True; actual check below
try:
    import plotly  # noqa: F401
    _has_plotly = True
except ImportError:
    _has_plotly = False

requires_plotly = pytest.mark.skipif(not _has_plotly, reason="plotly not installed")


class TestSampleSummaryView:
    """Tests for per-run and per-sample peptide/protein counts."""

    def test_peptides_per_run(self, dataset_dir):
        ds = Dataset(dataset_dir)
        result = ds.sample_summary.peptides_per_run()
        df = result.to_df()
        assert "run_file_name" in df.columns
        assert "n_peptides" in df.columns
        # dataset_dir fixture: run_01 has PEPTIDEK + ANOTHERK, run_02 has THIRDPEP
        assert len(df) == 2
        run01 = df[df["run_file_name"] == "run_01"]
        assert run01["n_peptides"].iloc[0] == 2
        run02 = df[df["run_file_name"] == "run_02"]
        assert run02["n_peptides"].iloc[0] == 1

    def test_peptides_per_sample(self, dataset_dir):
        ds = Dataset(dataset_dir)
        result = ds.sample_summary.peptides_per_sample()
        df = result.to_df()
        assert "sample_accession" in df.columns
        assert "n_peptides" in df.columns
        assert len(df) == 2

    def test_proteins_per_sample(self, dataset_dir):
        ds = Dataset(dataset_dir)
        result = ds.sample_summary.proteins_per_sample()
        df = result.to_df()
        assert "sample_accession" in df.columns
        assert "n_proteins" in df.columns
        assert len(df) == 2

    @requires_plotly
    def test_plot_peptides_per_run(self, dataset_dir):
        ds = Dataset(dataset_dir)
        fig = ds.sample_summary.plot("peptides_per_run")
        assert fig is not None
        assert hasattr(fig, "data")

    @requires_plotly
    def test_plot_peptides_per_sample(self, dataset_dir):
        ds = Dataset(dataset_dir)
        fig = ds.sample_summary.plot("peptides_per_sample")
        assert fig is not None

    @requires_plotly
    def test_plot_proteins_per_sample(self, dataset_dir):
        ds = Dataset(dataset_dir)
        fig = ds.sample_summary.plot("proteins_per_sample")
        assert fig is not None

    def test_plot_unknown_metric_raises(self, dataset_dir):
        ds = Dataset(dataset_dir)
        with pytest.raises(ValueError, match="Unknown metric"):
            ds.sample_summary.plot("invalid_metric")
