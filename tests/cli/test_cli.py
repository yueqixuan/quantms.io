"""CLI command tests — convert, transform, query, info, validate, ontology."""

from click.testing import CliRunner

from qpx.cli.main import qpx_main

# ---------------------------------------------------------------------------
# Convert
# ---------------------------------------------------------------------------


class TestQuantMSConvertCLI:
    def test_quantms_help_renders(self):
        runner = CliRunner()
        result = runner.invoke(qpx_main, ["convert", "quantms", "--help"])
        assert result.exit_code == 0
        assert "--mztab-path" in result.output


class TestDiaNNConvertCLI:
    def test_diann_help_renders(self):
        runner = CliRunner()
        result = runner.invoke(qpx_main, ["convert", "diann", "--help"])
        assert result.exit_code == 0
        assert "--report-path" in result.output


class TestMaxQuantConvertCLI:
    def test_maxquant_help_renders(self):
        runner = CliRunner()
        result = runner.invoke(qpx_main, ["convert", "maxquant", "--help"])
        assert result.exit_code == 0
        assert "--msms-file" in result.output


class TestFragPipeConvertCLI:
    def test_fragpipe_help_renders(self):
        runner = CliRunner()
        result = runner.invoke(qpx_main, ["convert", "fragpipe", "--help"])
        assert result.exit_code == 0
        assert "--psm-file" in result.output

    def test_fragpipe_help_shows_new_options(self):
        runner = CliRunner()
        result = runner.invoke(qpx_main, ["convert", "fragpipe", "--help"])
        assert "--ion-file" in result.output
        assert "--peptide-file" in result.output
        assert "--pg-file" in result.output


class TestMzIdentMLConvertCLI:
    def test_mzidentml_help_renders(self):
        runner = CliRunner()
        result = runner.invoke(qpx_main, ["convert", "mzidentml", "--help"])
        assert result.exit_code == 0
        assert "--mzid-path" in result.output

    def test_mzidentml_help_shows_new_options(self):
        runner = CliRunner()
        result = runner.invoke(qpx_main, ["convert", "mzidentml", "--help"])
        assert "--mgf-path" in result.output
        assert "--include-spectra" in result.output
        assert "--project-accession" in result.output

    def test_mzidentml_help_shows_enrich_pride(self):
        runner = CliRunner()
        result = runner.invoke(qpx_main, ["convert", "mzidentml", "--help"])
        assert "--enrich-pride" in result.output


class TestSdrfConvertCLI:
    def test_sdrf_help_renders(self):
        runner = CliRunner()
        result = runner.invoke(qpx_main, ["convert", "sdrf", "--help"])
        assert result.exit_code == 0
        assert "--sdrf-file" in result.output


# ---------------------------------------------------------------------------
# Transform
# ---------------------------------------------------------------------------


class TestTransformGeneMapCLI:
    def test_genemap_help_renders(self):
        runner = CliRunner()
        result = runner.invoke(qpx_main, ["transform", "gene-map", "--help"])
        assert result.exit_code == 0
        assert "--fasta" in result.output


class TestTransformQuantifyCLI:
    def test_quantify_help_renders(self):
        runner = CliRunner()
        result = runner.invoke(qpx_main, ["transform", "quantify", "--help"])
        assert result.exit_code == 0
        assert "--feature-path" in result.output
        assert "--method" in result.output
        assert "--output" in result.output

    def test_quantify_help_shows_ibaq_options(self):
        runner = CliRunner()
        result = runner.invoke(qpx_main, ["transform", "quantify", "--help"])
        assert "--organism" in result.output
        assert "--ploidy" in result.output
        assert "--min-aa" in result.output
