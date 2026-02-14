"""View plotting tests."""

import pytest


@pytest.fixture
def ds(dataset_dir):
    import qpx

    d = qpx.open(str(dataset_dir))
    yield d
    d.close()


class TestIdentificationSummaryPlot:
    def test_plot_returns_figure(self, ds):
        import matplotlib.figure

        fig = ds.identification_summary.plot()
        assert isinstance(fig, matplotlib.figure.Figure)


class TestRunSummaryPlot:
    def test_plot_returns_figure(self, ds):
        import matplotlib.figure

        fig = ds.run_summary.plot()
        assert isinstance(fig, matplotlib.figure.Figure)


class TestQCViewPlot:
    def test_plot_returns_figure(self, ds):
        import matplotlib.figure

        fig = ds.qc_view.plot()
        assert isinstance(fig, matplotlib.figure.Figure)


class TestPlotSave:
    def test_save_to_png(self, ds, tmp_path):
        fig = ds.identification_summary.plot()
        out = tmp_path / "plot.png"
        fig.savefig(out)
        assert out.exists()
        assert out.stat().st_size > 0
