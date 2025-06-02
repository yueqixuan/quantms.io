"""
quantmsio CLI main entry point
"""

import click
import logging
from quantmsio.utils.logger import get_logger
from quantmsio import __version__

# Import commands
from quantmsio.commands.feature_command import convert_feature_file
from quantmsio.commands.psm_command import convert_psm_file, compare_set_of_psms
from quantmsio.commands.diann_command import diann_convert_to_parquet, diann_pg_convert_to_parquet
from quantmsio.commands.ae_command import convert_ibaq_absolute
from quantmsio.commands.de_command import convert_msstats_differential
from quantmsio.commands.attach_file_command import attach_file_to_json
from quantmsio.commands.generate_spectra_message_command import map_spectrum_message_to_parquet
from quantmsio.commands.generate_gene_message_command import map_gene_message_to_parquet
from quantmsio.commands.plot_command import plot
from quantmsio.commands.statistic_command import statistics
from quantmsio.commands.maxquant_command import convert_maxquant_psm, convert_maxquant_feature
from quantmsio.commands.ibaq_command import convert_ibaq_file
from quantmsio.commands.fragpipe_command import convert_fragpipe_psm
from quantmsio.commands.map_latest_uniport_command import map_latest_uniport
from quantmsio.commands.anndata_command import merge_ae_files
from quantmsio.commands.project_command import generate_pride_project_json

logger = get_logger("quantmsio.cli")

@click.group()
@click.version_option(version=__version__)
def cli():
    """
    quantmsio - A tool for working with mass spectrometry data in the quantms.io format.

    This CLI provides commands for:
    - Converting data from various formats to quantms.io
    - Analyzing and processing quantms.io data
    - Visualizing results
    - Managing project metadata
    """
    pass


# Conversion commands group
@cli.group()
def convert():
    """Commands for converting data to quantms.io format."""
    pass

convert.add_command(convert_feature_file, name="features")
convert.add_command(convert_psm_file, name="psm")
convert.add_command(convert_maxquant_psm, name="maxquant-psm")
convert.add_command(convert_maxquant_feature, name="maxquant-features")
convert.add_command(convert_fragpipe_psm, name="fragpipe-psm")
convert.add_command(convert_msstats_differential, name="differential")
convert.add_command(convert_ibaq_absolute, name="absolute")
convert.add_command(diann_convert_to_parquet, name="diann")
convert.add_command(diann_pg_convert_to_parquet, name="diann-pg")


# Analysis commands group
@cli.group()
def analyze():
    """Commands for analyzing quantms.io data."""
    pass

analyze.add_command(compare_set_of_psms, name="compare-psms")
analyze.add_command(statistics, name="stats")


# Visualization commands group
@cli.group()
def visualize():
    """Commands for visualizing quantms.io data."""
    pass

visualize.add_command(plot)


# Project commands group
@cli.group()
def project():
    """Commands for managing project metadata."""
    pass

project.add_command(generate_pride_project_json, name="pride")
project.add_command(attach_file_to_json, name="attach")


# Transform commands group
@cli.group()
def transform():
    """Commands for data transformation operations."""
    pass

transform.add_command(map_spectrum_message_to_parquet, name="spectrum-parquet")
transform.add_command(map_gene_message_to_parquet, name="gene-parquet")
transform.add_command(map_latest_uniport, name="update-uniprot")
transform.add_command(merge_ae_files, name="merge-ae")
transform.add_command(convert_ibaq_file, name="ibaq")


def quantms_io_main():
    """Main entry point for the quantmsio CLI."""
    try:
        logger.info("Starting quantmsio CLI")
        cli()
    except Exception as e:
        logger.exception(str(e))
        raise
