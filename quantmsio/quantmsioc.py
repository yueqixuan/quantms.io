"""
Commandline interface for quantmsio package allows generating the quantms.io file format from different sources and
 steps. The quantms.io specification is available in the docs folder of this repository.
"""

import logging
import click

from quantmsio import __version__ as __version__

# Convert commands
from quantmsio.commands.convert.feature import convert_feature_file
from quantmsio.commands.convert.psm import convert_psm_file, compare_set_of_psms
from quantmsio.commands.convert.diann import (
    diann_convert_to_parquet,
    diann_pg_convert_to_parquet,
)
from quantmsio.commands.convert.maxquant import (
    convert_maxquant_psm,
    convert_maxquant_feature,
)
from quantmsio.commands.convert.fragpipe import convert_fragpipe_psm

# Transform commands
from quantmsio.commands.transform.ae import convert_ibaq_absolute
from quantmsio.commands.transform.de import convert_msstats_differential
from quantmsio.commands.transform.spectra import map_spectrum_message_to_parquet
from quantmsio.commands.transform.gene import map_gene_message_to_parquet
from quantmsio.commands.transform.ibaq import convert_ibaq_file
from quantmsio.commands.transform.uniprot import map_latest_uniport
from quantmsio.commands.transform.anndata import merge_ae_files

# Utility commands
from quantmsio.commands.utils.project import generate_pride_project_json
from quantmsio.commands.utils.attach import attach_file_to_json
from quantmsio.commands.utils.plot import plot
from quantmsio.commands.utils.stats import statistics

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
convert.add_command(convert_feature_file, name="feature")
convert.add_command(convert_psm_file, name="psm")
convert.add_command(compare_set_of_psms, name="compare-psms")
convert.add_command(diann_convert_to_parquet, name="diann")
convert.add_command(diann_pg_convert_to_parquet, name="diann-pg")
convert.add_command(convert_maxquant_psm, name="maxquant-psm")
convert.add_command(convert_maxquant_feature, name="maxquant-feature")
convert.add_command(convert_fragpipe_psm, name="fragpipe")

# Transform commands
transform.add_command(convert_ibaq_absolute, name="ae")
transform.add_command(convert_msstats_differential, name="de")
transform.add_command(map_spectrum_message_to_parquet, name="spectra")
transform.add_command(map_gene_message_to_parquet, name="gene")
transform.add_command(convert_ibaq_file, name="ibaq")
transform.add_command(map_latest_uniport, name="uniprot")
transform.add_command(merge_ae_files, name="anndata")

# Visualization commands
visualize.add_command(plot, name="plot")

# Statistics commands
stats.add_command(statistics, name="analyze")

# Project commands
project.add_command(generate_pride_project_json, name="create")
project.add_command(attach_file_to_json, name="attach")


def quantms_io_main() -> None:
    """
    Main function to run the quantmsio command line interface
    :return: none
    """
    cli()


if __name__ == "__main__":
    quantms_io_main()
