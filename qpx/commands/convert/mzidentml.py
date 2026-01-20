"""mzIdentML-specific converters for QPX formats."""

import logging
from pathlib import Path
from typing import Optional

import click

from qpx.core.mzidentml import MzIdentML
from qpx.core.project import create_uuid_filename
from qpx.utils.logger import get_logger


@click.command(
    "convert-mzidentml",
    short_help="Convert mzIdentML to PSM parquet file in QPX",
)
@click.option(
    "--mzid-file",
    help="mzIdentML file (.mzid or .mzid.gz)",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "--output-folder",
    help="Output folder",
    required=True,
    type=click.Path(file_okay=False, path_type=Path),
)
@click.option(
    "--mzml-file",
    help="Optional single mzML file to attach spectra by scan",
    required=False,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "--mzml-folder",
    help="Optional folder containing mzML files (matches by reference_file_name)",
    required=False,
    type=click.Path(exists=True, file_okay=False, path_type=Path),
)
@click.option(
    "--output-prefix",
    help="Output file prefix",
    default=None,
)
@click.option(
    "--spectral-data",
    help="Include spectral data fields",
    is_flag=True,
)
@click.option(
    "--verbose",
    help="Enable verbose logging",
    is_flag=True,
)
def convert_mzidentml_file(
    mzid_file: Path,
    output_folder: Path,
    mzml_file: Optional[Path] = None,
    mzml_folder: Optional[Path] = None,
    output_prefix: Optional[str] = None,
    spectral_data: bool = False,
    verbose: bool = False,
) -> None:
    """Convert mzIdentML file to QPX PSM format.

    This command takes a mzIdentML file and converts it to the QPX parquet
    format for PSM data. Supports both .mzid and .mzid.gz compressed files.

    Example:
        qpxc convert mzidentml \\
            --mzid-file data.mzid \\
            --output-folder ./output \\
            --output-prefix pride-psm
    """
    logger = get_logger("qpx.commands.mzidentml")
    if verbose:
        logger.setLevel(logging.DEBUG)
        logger.debug("Verbose logging enabled")

    try:
        if not all([mzid_file, output_folder]):
            raise click.UsageError("ERROR: Please provide all required parameters")

        output_folder = Path(output_folder)
        output_folder.mkdir(parents=True, exist_ok=True)
        logger.info(f"Using output directory: {output_folder}")

        prefix = output_prefix or "psm"
        filename = create_uuid_filename(prefix, ".psm.parquet")
        output_path = output_folder / filename
        logger.info(f"Will save PSM file as: {filename}")

        # Validate mzML options
        if mzml_file and mzml_folder:
            raise click.UsageError("Cannot specify both --mzml-file and --mzml-folder")

        logger.info("Initializing mzIdentML parser...")
        parser = MzIdentML(
            mzid_path=mzid_file,
            mzml_path=mzml_file,
            mzml_folder=mzml_folder,
            spectral_data=spectral_data,
        )

        parser.to_parquet(output_path)
        logger.info(f"PSM file successfully saved to: {output_path}")
        logger.info(f"Total PSMs converted: {parser.get_psm_count():,}")

    except Exception as e:
        logger.error(f"Error in mzIdentML conversion: {str(e)}", exc_info=True)
        raise click.ClickException(f"Error: {str(e)}\nCheck the logs for more details.")
