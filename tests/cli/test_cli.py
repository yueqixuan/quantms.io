"""CLI command tests — convert, transform, query, info, validate, ontology."""

from click.testing import CliRunner

from qpx.cli.main import qpx_main

# ---------------------------------------------------------------------------
# Convert
# ---------------------------------------------------------------------------


def _assert_help(result, *options):
    """Verify CLI help renders and contains expected options."""
    if result.exit_code != 0:
        raise AssertionError(f"exit_code={result.exit_code}, output={result.output}")
    for opt in options:
        if opt not in result.output:
            raise AssertionError(f"Missing option {opt} in help output")


class TestQuantMSConvertCLI:
    def test_quantms_help_renders(self):
        runner = CliRunner()
        result = runner.invoke(qpx_main, ["convert", "quantms", "--help"])
        _assert_help(result, "--mztab-path")


class TestDiaNNConvertCLI:
    def test_diann_help_renders(self):
        runner = CliRunner()
        result = runner.invoke(qpx_main, ["convert", "diann", "--help"])
        _assert_help(result, "--report-path")


class TestMaxQuantConvertCLI:
    def test_maxquant_help_renders(self):
        runner = CliRunner()
        result = runner.invoke(qpx_main, ["convert", "maxquant", "--help"])
        _assert_help(result, "--msms-file")


class TestFragPipeConvertCLI:
    def test_fragpipe_help_renders(self):
        runner = CliRunner()
        result = runner.invoke(qpx_main, ["convert", "fragpipe", "--help"])
        _assert_help(result, "--psm-file")

    def test_fragpipe_help_shows_new_options(self):
        runner = CliRunner()
        result = runner.invoke(qpx_main, ["convert", "fragpipe", "--help"])
        _assert_help(result, "--ion-file", "--peptide-file", "--pg-file")


class TestMzIdentMLConvertCLI:
    def test_mzidentml_help_renders(self):
        runner = CliRunner()
        result = runner.invoke(qpx_main, ["convert", "mzidentml", "--help"])
        _assert_help(result, "--mzid-path")

    def test_mzidentml_help_shows_new_options(self):
        runner = CliRunner()
        result = runner.invoke(qpx_main, ["convert", "mzidentml", "--help"])
        _assert_help(result, "--mgf-path", "--include-spectra", "--project-accession")

    def test_mzidentml_help_shows_enrich_pride(self):
        runner = CliRunner()
        result = runner.invoke(qpx_main, ["convert", "mzidentml", "--help"])
        _assert_help(result, "--enrich-pride")


class TestSdrfConvertCLI:
    def test_sdrf_help_renders(self):
        runner = CliRunner()
        result = runner.invoke(qpx_main, ["convert", "sdrf", "--help"])
        _assert_help(result, "--sdrf-file")


# ---------------------------------------------------------------------------
# Transform
# ---------------------------------------------------------------------------


class TestTransformGeneMapCLI:
    def test_genemap_help_renders(self):
        runner = CliRunner()
        result = runner.invoke(qpx_main, ["transform", "gene-map", "--help"])
        _assert_help(result, "--fasta")


class TestTransformQuantifyCLI:
    def test_quantify_help_renders(self):
        runner = CliRunner()
        result = runner.invoke(qpx_main, ["transform", "quantify", "--help"])
        _assert_help(result, "--feature-path", "--method", "--output")

    def test_quantify_help_shows_ibaq_options(self):
        runner = CliRunner()
        result = runner.invoke(qpx_main, ["transform", "quantify", "--help"])
        _assert_help(result, "--organism", "--ploidy", "--min-aa")
