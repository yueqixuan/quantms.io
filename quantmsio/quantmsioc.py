"""
Commandline interface for quantmsio package that provides tools for working with the quantms.io file format.
The quantms.io specification is available in the docs folder of this repository.
"""

import click
from typing import Optional, Dict, Any

from quantmsio import __version__
from quantmsio.utils.logger import setup_logging, get_logger
from quantmsio.commands.project_command import generate_pride_project_json
from quantmsio.commands.feature_command import convert_feature_file
from quantmsio.commands.psm_command import convert_psm_file, compare_set_of_psms
from quantmsio.commands.diann_command import (
    diann_convert_to_parquet,
    diann_pg_convert_to_parquet,
)
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

CONTEXT_SETTINGS: Dict[str, Any] = dict(help_option_names=["-h", "--help"])

def init_logging() -> None:
    """Initialize logging with default configuration."""
    setup_logging()
    logger = get_logger(__name__)
    logger.info("Starting quantmsio CLI")

@click.version_option(
    version=__version__,
    package_name="quantmsio",
    message="%(package)s %(version)s"
)
@click.group(context_settings=CONTEXT_SETTINGS)
def cli() -> None:
    """
    quantmsio - A tool for working with mass spectrometry data in the quantms.io format.
    
    This CLI provides commands for:
    - Converting data from various formats to quantms.io
    - Analyzing and processing quantms.io data
    - Visualizing results
    - Managing project metadata
    """
    init_logging()

# Conversion commands group
@cli.group()
def convert() -> None:
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
def analyze() -> None:
    """Commands for analyzing quantms.io data."""
    pass

analyze.add_command(compare_set_of_psms, name="compare-psms")
analyze.add_command(statistics, name="stats")

# Visualization commands group
@cli.group()
def visualize() -> None:
    """Commands for visualizing quantms.io data."""
    pass

visualize.add_command(plot)

# Project management commands group
@cli.group()
def project() -> None:
    """Commands for managing project metadata."""
    pass

project.add_command(generate_pride_project_json, name="pride")
project.add_command(attach_file_to_json, name="attach")

# Data transformation commands group
@cli.group()
def transform() -> None:
    """Commands for data transformation operations."""
    pass

transform.add_command(map_spectrum_message_to_parquet, name="spectrum-parquet")
transform.add_command(map_gene_message_to_parquet, name="gene-parquet")
transform.add_command(map_latest_uniport, name="update-uniprot")
transform.add_command(merge_ae_files, name="merge-ae")
transform.add_command(convert_ibaq_file, name="ibaq")

def quantms_io_main() -> None:
    """
    Main entry point for the quantmsio command line interface.
    Initializes logging and runs the CLI application.
    """
    cli()

if __name__ == "__main__":
    quantms_io_main()
