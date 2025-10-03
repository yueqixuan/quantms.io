"""
Command to convert differential expression data to quantms.io format.
"""

import logging
from pathlib import Path
from typing import Optional

import click

from quantmsio.core.de import DifferentialExpressionHandler
from quantmsio.utils.file_utils import extract_protein_list
from quantmsio.utils.logger import get_logger


def convert_msstats_differential(
    msstats_file: Path,
    sdrf_file: Path,
    project_file: Optional[Path] = None,
    protein_file: Optional[Path] = None,
    fdr_threshold: float = 0.05,
    output_folder: Path = None,
    output_prefix: Optional[str] = None,
    delete_existing: bool = True,
    verbose: bool = False,
) -> None:
    """
    Convert MSstats differential file into quantms.io format.

    Args:
        msstats_file: MSstats differential file
        sdrf_file: SDRF file needed to extract metadata
        project_file: Optional quantms.io project file
        protein_file: Optional protein file with specific requirements
        fdr_threshold: FDR threshold to filter results (default: 0.05)
        output_folder: Output directory for generated files
        output_prefix: Optional prefix for output files
        delete_existing: Whether to delete existing files in output folder
        verbose: Enable verbose logging
    """
    logger = get_logger("quantmsio.commands.de")
    if verbose:
        logger.setLevel(logging.DEBUG)

    try:
        if not all([msstats_file, sdrf_file, output_folder]):
            raise click.UsageError("Please provide all required parameters")

        protein_list = extract_protein_list(protein_file) if protein_file else None
        protein_str = "|".join(protein_list) if protein_list else None

        de_handler = DifferentialExpressionHandler()
        if project_file:
            de_handler.load_project_file(project_file)

        de_handler.load_msstats_file(msstats_file, protein_str)
        de_handler.load_sdrf_file(sdrf_file)
        de_handler.set_fdr_threshold(fdr_threshold=fdr_threshold)

        de_handler.convert_msstats_to_quantms(
            output_folder=output_folder,
            output_file_prefix=output_prefix,
            delete_existing=delete_existing,
        )

    except Exception as e:
        logger.exception(f"Error in differential expression conversion: {str(e)}")
        raise click.ClickException(f"Error: {str(e)}\nCheck the logs for more details.")


@click.command(
    "differential",
    short_help="Convert a MSstats differential file into a quantms.io file format",
)
@click.option(
    "--msstats-file",
    help="MSstats differential file",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "--sdrf-file",
    help="SDRF file needed to extract metadata",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "--project-file",
    help="quantms.io project file",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "--protein-file",
    help="Protein file that meets specific requirements",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "--fdr-threshold",
    help="FDR threshold to filter results",
    default="0.05",
    type=float,
)
@click.option(
    "--output-folder",
    help="Output directory for generated files",
    required=True,
    type=click.Path(file_okay=False, path_type=Path),
)
@click.option("--output-prefix", help="Prefix for output files")
@click.option(
    "--delete-existing", help="Delete existing files in output folder", is_flag=True
)
@click.option("--verbose", help="Enable verbose logging", is_flag=True)
def convert_msstats_differential_cmd(**kwargs):
    """
    Convert MSstats differential expression data to quantms.io format.
    
    This command takes a differential expression file from MSstats and converts it
    to the quantms.io format, using metadata from an SDRF file. The results can be
    filtered by FDR threshold and specific proteins.
    
    Example:
        quantmsio convert differential \\
            --msstats-file data.msstats.txt \\
            --sdrf-file data.sdrf.tsv \\
            --output-folder ./output \\
            --fdr-threshold 0.05
    """
    convert_msstats_differential(**kwargs)
