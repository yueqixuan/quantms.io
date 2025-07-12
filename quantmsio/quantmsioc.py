"""
Commandline interface for quantmsio package allows generating the quantms.io file format from different sources and
 steps. The quantms.io specification is available in the docs folder of this repository.
"""

import logging

import click

from quantmsio import __version__ as __version__

# Convert commands
from quantmsio.commands.convert.quantms import convert as quantms_convert
from quantmsio.commands.convert.quantms import (
    convert_quantms_psm_cmd as quantms_psm_convert,
)
from quantmsio.commands.convert.quantms import create_duckdb_cmd as duckdb_creator
from quantmsio.commands.convert.quantms_project import (
    convert_quantms_project_cmd as quantms_project_convert,
)

# Transform commands
from quantmsio.commands.transform.ae import (
    convert_ibaq_absolute_cmd as ae_transform,
)
from quantmsio.commands.transform.anndata import (
    merge_ae_files_cmd as anndata_transform,
)
from quantmsio.commands.transform.de import (
    convert_msstats_differential_cmd as de_transform,
)
from quantmsio.commands.transform.gene import map_gene_message_cmd as gene_transform
from quantmsio.commands.transform.ibaq import convert_ibaq_file_cmd as ibaq_transform
from quantmsio.commands.transform.spectra import (
    map_spectrum_message_cmd as spectra_transform,
)
from quantmsio.commands.transform.uniprot import (
    map_latest_uniprot_cmd as uniprot_transform,
)

# Utility commands
from quantmsio.commands.utils.attach import attach_file_to_json_cmd as attach_utils
from quantmsio.commands.utils.plot import plot_cmd as plot_utils
from quantmsio.commands.utils.project import (
    generate_pride_project_json_cmd as project_utils,
)
from quantmsio.commands.utils.stats import statistics_cmd as stats_utils

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
convert.add_command(quantms_convert)
convert.add_command(quantms_psm_convert)
convert.add_command(duckdb_creator)
convert.add_command(quantms_project_convert)


# Transform commands
transform.add_command(ae_transform)
transform.add_command(de_transform)
transform.add_command(gene_transform)
transform.add_command(ibaq_transform)
transform.add_command(spectra_transform)
transform.add_command(uniprot_transform)
transform.add_command(anndata_transform)

# Visualization commands
visualize.add_command(plot_utils, name="plot")

# Statistics commands
stats.add_command(stats_utils, name="analyze")

# Project commands
project.add_command(project_utils, name="create")
project.add_command(attach_utils, name="attach")


def quantms_io_main() -> None:
    """
    Main function to run the quantmsio command line interface
    :return: none
    """
    cli()


if __name__ == "__main__":
    quantms_io_main()
