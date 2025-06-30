"""
Commandline interface for quantmsio package allows generating the quantms.io file format from different sources and
 steps. The quantms.io specification is available in the docs folder of this repository.
"""

import logging

import click

from quantmsio import __version__ as __version__

# Convert commands
from quantmsio.commands.convert.diann import convert_diann_cmd, convert_diann_pg_cmd
from quantmsio.commands.convert.fragpipe import convert_fragpipe_psm_cmd
from quantmsio.commands.convert.maxquant import (
    convert_maxquant_feature_cmd,
    convert_maxquant_pg_cmd,
    convert_maxquant_psm_cmd,
)
from quantmsio.commands.convert.quantms import (
    convert_quantms_feature_cmd,
    convert_quantms_pg_cmd,
    convert_quantms_psm_cmd,
)

# Transform commands
from quantmsio.commands.transform.ae import convert_ibaq_absolute_cmd
from quantmsio.commands.transform.anndata import merge_ae_files_cmd
from quantmsio.commands.transform.de import convert_msstats_differential_cmd
from quantmsio.commands.transform.gene import map_gene_message_cmd
from quantmsio.commands.transform.ibaq import convert_ibaq_file_cmd
from quantmsio.commands.transform.spectra import map_spectrum_message_cmd
from quantmsio.commands.transform.uniprot import map_latest_uniprot_cmd
from quantmsio.commands.utils.attach import attach_file_to_json_cmd
from quantmsio.commands.utils.plot import plot_cmd

# Utility commands
from quantmsio.commands.utils.project import generate_pride_project_json_cmd
from quantmsio.commands.utils.stats import statistics_cmd

CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])


@click.version_option(
    version=__version__, package_name="quantmsio", message="%(package)s %(version)s"
)
@click.group(context_settings=CONTEXT_SETTINGS)
def cli() -> None:
    """
    quantmsio - A tool for converting and analyzing mass spectrometry proteomics data
    """
    logging.basicConfig(
        level=logging.INFO,
        datefmt="%H:%M:%S",
        format="[%(asctime)s] %(levelname).1s | %(name)s | %(message)s",
    )


@cli.group()
def convert():
    """Convert external formats to quantms.io format."""
    pass


@cli.group()
def transform():
    """Transform quantms.io data into different representations."""
    pass


@cli.group()
def visualize():
    """Visualize quantms.io data."""
    pass


@cli.group()
def stats():
    """Statistical analysis of quantms.io data."""
    pass


@cli.group()
def project():
    """Project management commands."""
    pass


# Convert commands
convert.add_command(convert_diann_cmd, name="diann")
convert.add_command(convert_diann_pg_cmd, name="diann-pg")
convert.add_command(convert_maxquant_psm_cmd, name="maxquant-psm")
convert.add_command(convert_maxquant_feature_cmd, name="maxquant-feature")
convert.add_command(convert_maxquant_pg_cmd, name="maxquant-pg")
convert.add_command(convert_fragpipe_psm_cmd, name="fragpipe")
convert.add_command(convert_quantms_psm_cmd, name="quantms-psm")
convert.add_command(convert_quantms_feature_cmd, name="quantms-feature")
convert.add_command(convert_quantms_pg_cmd, name="quantms-pg")


# Transform commands
transform.add_command(convert_ibaq_absolute_cmd, name="ae")
transform.add_command(convert_msstats_differential_cmd, name="de")
transform.add_command(map_spectrum_message_cmd, name="spectra")
transform.add_command(map_gene_message_cmd, name="gene")
transform.add_command(convert_ibaq_file_cmd, name="ibaq")
transform.add_command(map_latest_uniprot_cmd, name="uniprot")
transform.add_command(merge_ae_files_cmd, name="anndata")

# Visualization commands
visualize.add_command(plot_cmd, name="plot")

# Statistics commands
stats.add_command(statistics_cmd, name="analyze")

# Project commands
project.add_command(generate_pride_project_json_cmd, name="create")
project.add_command(attach_file_to_json_cmd, name="attach")


def quantms_io_main() -> None:
    """
    Main function to run the quantmsio command line interface
    :return: none
    """
    cli()


if __name__ == "__main__":
    quantms_io_main()
