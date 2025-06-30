"""
Core converters for quantms.io formats (PSM, Feature, mzTab Protein Groups).
"""

import logging
from pathlib import Path
from typing import Optional

import click
import pyarrow.parquet as pq

from quantmsio.core.project import create_uuid_filename
from quantmsio.core.quantms.feature import Feature
from quantmsio.core.quantms.pg import MzTabProteinGroups
from quantmsio.core.quantms.psm import Psm


@click.group()
def convert():
    """Convert various formats to quantms.io format."""
    pass


# Feature conversion
@convert.command("quantms-feature")
@click.option(
    "--input-file",
    help="Input mzTab file path",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "--output-folder",
    help="Output directory for generated files",
    required=True,
    type=click.Path(file_okay=False, path_type=Path),
)
@click.option(
    "--output-prefix",
    help="Prefix for output files (final name will be {prefix}-{uuid}.feature.parquet)",
    default="feature",
)
@click.option(
    "--sdrf-file",
    help="SDRF file path",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "--msstats-file",
    help="MSstats input file path",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option("--verbose", help="Enable verbose logging", is_flag=True)
def convert_quantms_feature_cmd(
    input_file: Path,
    output_folder: Path,
    output_prefix: str,
    sdrf_file: Optional[Path] = None,
    msstats_file: Optional[Path] = None,
    verbose: bool = False,
):
    """Convert feature data from mzTab to quantms.io format."""
    logger = logging.getLogger("quantmsio.commands.convert.feature")
    if verbose:
        logger.setLevel(logging.DEBUG)

    try:
        # Create output directory if it doesn't exist
        output_folder.mkdir(parents=True, exist_ok=True)

        # Generate output filename with UUID
        filename = create_uuid_filename(output_prefix, ".feature.parquet")
        output_file = output_folder / filename

        feature = Feature(
            mztab_path=str(input_file),
            sdrf_path=str(sdrf_file) if sdrf_file else None,
            msstats_in_path=str(msstats_file) if msstats_file else None,
        )

        feature.write_feature_to_file(output_path=str(output_file))

    except Exception as e:
        logger.exception(f"Error in feature conversion: {str(e)}")
        raise click.ClickException(f"Error: {str(e)}\nCheck the logs for more details.")


# PSM conversion
@convert.command("quantms-psm")
@click.option(
    "--input-file",
    help="Input mzTab file path",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "--output-folder",
    help="Output directory for generated files",
    required=True,
    type=click.Path(file_okay=False, path_type=Path),
)
@click.option(
    "--output-prefix",
    help="Prefix for output files (final name will be {prefix}-{uuid}.psm.parquet)",
    default="psm",
)
@click.option(
    "--sdrf-file",
    help="SDRF file path",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option("--verbose", help="Enable verbose logging", is_flag=True)
def convert_quantms_psm_cmd(
    input_file: Path,
    output_folder: Path,
    output_prefix: str,
    sdrf_file: Optional[Path] = None,
    verbose: bool = False,
):
    """Convert PSM data from mzTab to quantms.io format."""
    logger = logging.getLogger("quantmsio.commands.convert.psm")
    if verbose:
        logger.setLevel(logging.DEBUG)

    try:
        # Create output directory if it doesn't exist
        output_folder.mkdir(parents=True, exist_ok=True)

        # Generate output filename with UUID
        filename = create_uuid_filename(output_prefix, ".psm.parquet")
        output_file = output_folder / filename

        psm = Psm(mztab_path=str(input_file))

        psm.write_psm_to_file(output_path=str(output_file))

    except Exception as e:
        logger.exception(f"Error in PSM conversion: {str(e)}")
        raise click.ClickException(f"Error: {str(e)}\nCheck the logs for more details.")


# mzTab Protein Group conversion
@convert.command("quantms-pg")
@click.option(
    "--input-file",
    help="Input mzTab file path",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "--msstats-file",
    help="Input msstats_in.csv file path for quantification",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "--sdrf-file",
    help="SDRF file path",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "--output-folder",
    help="Output directory for generated files",
    required=True,
    type=click.Path(file_okay=False, path_type=Path),
)
@click.option(
    "--output-prefix",
    help="Prefix for output files (final name will be {prefix}-{uuid}.pg.parquet)",
    default="pg",
)
@click.option(
    "--compute-topn/--no-compute-topn",
    help="Whether to compute TopN intensity",
    default=True,
)
@click.option(
    "--compute-ibaq/--no-compute-ibaq",
    help="Whether to compute iBAQ intensity",
    default=True,
)
@click.option(
    "--topn",
    help="Number of peptides to use for TopN intensity",
    default=3,
    type=int,
)
@click.option("--verbose", help="Enable verbose logging", is_flag=True)
def convert_quantms_pg_cmd(
    input_file: Path,
    msstats_file: Path,
    sdrf_file: Path,
    output_folder: Path,
    output_prefix: str,
    compute_topn: bool = True,
    compute_ibaq: bool = True,
    topn: int = 3,
    verbose: bool = False,
):
    """Convert protein groups from mzTab quantms TMT and LFQ data to quantms.io format using msstats for quantification.

    This command combines protein group definitions from mzTab with complete quantification
    data from msstats_in.csv. For LFQ data, it computes:

    - Default intensity: Sum of all peptidoform intensities
    - Additional intensities (optional):
      - TopN: Mean of top N peptide intensities (default N=3)
      - iBAQ: Intensity-based absolute quantification

    The output file will be named as {prefix}-{uuid}.pg.parquet.
    """
    logger = logging.getLogger("quantmsio.commands.convert.mztab")
    if verbose:
        logger.setLevel(logging.DEBUG)
        logger.info(
            "Converting mzTab protein groups to quantms.io format using msstats quantification"
        )

    try:
        # Create output directory if it doesn't exist
        output_folder.mkdir(parents=True, exist_ok=True)
        logger.info(f"Using output directory: {output_folder}")

        # Generate output filename with UUID
        filename = create_uuid_filename(output_prefix, ".pg.parquet")
        output_file = output_folder / filename
        logger.info(f"Will save protein groups file as: {filename}")

        # Perform conversion using OPTIMIZED msstats quantification (4-6x faster!)
        # Use context manager to ensure cleanup of any temporary files
        with MzTabProteinGroups(input_file) as mztab_pg:
            result_df = mztab_pg.quantify_from_msstats_optimized(
                str(msstats_file),
                str(sdrf_file),
                compute_topn=compute_topn,
                topn=topn,
                compute_ibaq=compute_ibaq,
            )

            # Convert to parquet and write
            table = mztab_pg._convert_to_parquet_format(result_df)
            pq.write_table(table, str(output_file))
            logger.info("Successfully wrote protein groups to parquet file")

    except Exception as e:
        logger.exception(f"Error in mzTab protein group conversion: {str(e)}")
        raise click.ClickException(f"Error: {str(e)}\nCheck the logs for more details.")
