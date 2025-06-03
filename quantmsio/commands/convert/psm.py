"""
Command to convert and compare PSM data in quantms.io format.
"""

from pathlib import Path
from typing import Optional, Tuple
import logging

import click
from quantmsio.core.project import create_uuid_filename
from quantmsio.core.psm import Psm
from quantmsio.operate.plots import plot_peptidoform_charge_venn, plot_sequence_venn
from quantmsio.utils.logger import get_logger


def convert_psm(
    mztab_file: Path,
    output_folder: Path,
    chunksize: int = 1000000,
    protein_file: Optional[Path] = None,
    output_prefix: Optional[str] = None,
    verbose: bool = False,
) -> None:
    """
    Convert PSM data from mzTab to parquet format.

    Args:
        mztab_file: mzTab file to extract protein information
        output_folder: Output directory for generated files
        chunksize: Read batch size (default: 1000000)
        protein_file: Optional protein file with specific requirements
        output_prefix: Optional prefix for output files
        verbose: Enable verbose logging
    """
    logger = get_logger("quantmsio.commands.psm")
    if verbose:
        logger.setLevel(logging.DEBUG)

    try:
        if not all([mztab_file, output_folder]):
            raise click.UsageError("Please provide all required parameters")

        # Ensure output directory exists
        output_folder = Path(output_folder)
        output_folder.mkdir(parents=True, exist_ok=True)

        # Set default prefix if not provided
        prefix = output_prefix or "psm"
        filename = create_uuid_filename(prefix, ".psm.parquet")
        output_path = output_folder / filename
        logger.info(f"Will save PSM file to: {output_path}")

        psm_manager = Psm(mztab_path=mztab_file)
        psm_manager.write_psm_to_file(
            output_path=str(output_path), chunksize=chunksize, protein_file=protein_file
        )
        logger.info(f"PSM file saved to: {output_path}")

    except Exception as e:
        logger.exception(f"Error in PSM conversion: {str(e)}")
        raise click.ClickException(f"Error: {str(e)}\nCheck the logs for more details.")


def compare_psm_sets(
    parquets: Tuple[Path, ...],
    tags: Tuple[str, ...],
    verbose: bool = False,
) -> None:
    """
    Compare a set of PSM parquet files.

    Args:
        parquets: Set of PSM parquet file paths
        tags: Set of labels for the parquet files
        verbose: Enable verbose logging
    """
    logger = get_logger("quantmsio.commands.psm")
    if verbose:
        logger.setLevel(logging.DEBUG)

    try:
        if len(parquets) != len(tags):
            raise click.UsageError(
                "Please provide same number of parquet files and labels"
            )

        plot_peptidoform_charge_venn(parquets, tags)
        plot_sequence_venn(parquets, tags)

    except Exception as e:
        logger.exception(f"Error in PSM comparison: {str(e)}")
        raise click.ClickException(f"Error: {str(e)}\nCheck the logs for more details.")


@click.command(
    "psm",
    short_help="Convert PSM data from mzTab to parquet format",
)
@click.option(
    "--mztab-file",
    help="mzTab file to extract protein information",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "--output-folder",
    help="Output directory for generated files",
    required=True,
    type=click.Path(file_okay=False, path_type=Path),
)
@click.option("--chunksize", help="Read batch size", default=1000000, type=int)
@click.option(
    "--protein-file",
    help="Protein file with specific requirements",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option("--output-prefix", help="Prefix for output files")
@click.option("--verbose", help="Enable verbose logging", is_flag=True)
def convert_psm_cmd(**kwargs):
    """
    Convert PSM data from mzTab to parquet format.
    
    This command takes an mzTab file and converts its PSM data to the quantms.io
    parquet format. The mzTab file is used to extract protein information.
    
    Example:
        quantmsio convert psm \\
            --mztab-file data.mztab \\
            --output-folder ./output \\
            --chunksize 1000000
    """
    convert_psm(**kwargs)


@click.command(
    "compare-psm",
    short_help="Compare PSM data from multiple parquet files",
)
@click.option(
    "-p",
    "--parquet",
    "parquets",
    help="PSM parquet file path",
    multiple=True,
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "-t",
    "--tag",
    "tags",
    help="Label for the parquet file",
    multiple=True,
    required=True,
)
@click.option("--verbose", help="Enable verbose logging", is_flag=True)
def compare_psm_cmd(**kwargs):
    """
    Compare PSM data from multiple parquet files.
    
    This command takes multiple PSM parquet files and generates comparison plots
    including peptidoform charge Venn diagrams and sequence Venn diagrams.
    
    Example:
        quantmsio compare psm \\
            -p file1.psm.parquet -t "Sample 1" \\
            -p file2.psm.parquet -t "Sample 2"
    """
    compare_psm_sets(**kwargs)
