"""
MaxQuant-specific converters for quantms.io formats.
"""

import logging
from pathlib import Path
from typing import Optional

import click

from quantmsio.core.maxquant.maxquant import MaxQuant
from quantmsio.core.project import create_uuid_filename
from quantmsio.utils.logger import get_logger


@click.group()
def convert():
    """Convert MaxQuant formats to quantms.io format."""
    pass


@convert.command(
    "maxquant-psm",
    short_help="Convert PSM data from MaxQuant msms.txt to parquet format",
)
@click.option(
    "--msms-file",
    help="MaxQuant msms.txt file to extract peptide information",
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
    "--batch-size",
    help="Read batch size",
    default=1000000,
    type=int,
)
@click.option(
    "--output-prefix",
    help="Prefix for output files",
)
@click.option("--verbose", help="Enable verbose logging", is_flag=True)
def convert_maxquant_psm_cmd(
    msms_file: Path,
    output_folder: Path,
    batch_size: int,
    output_prefix: Optional[str],
    verbose: bool = False,
):
    """
    Convert MaxQuant PSM data from msms.txt to parquet format.
    
    This command takes a MaxQuant msms.txt file and converts it to the quantms.io
    parquet format for PSM data.
    
    Example:
        quantmsio convert maxquant-psm \\
            --msms-file msms.txt \\
            --output-folder ./output \\
            --batch-size 1000000
    """
    logger = get_logger("quantmsio.commands.maxquant")
    if verbose:
        logger.setLevel(logging.DEBUG)
        logger.debug("Verbose logging enabled")

    try:
        if not all([msms_file, output_folder]):
            raise click.UsageError("ERROR: Please provide all required parameters")

        # Ensure output directory exists
        output_folder = Path(output_folder)
        output_folder.mkdir(parents=True, exist_ok=True)
        logger.info(f"Using output directory: {output_folder}")

        # Set default prefix if not provided
        prefix = output_prefix or "psm"
        filename = create_uuid_filename(prefix, ".psm.parquet")
        output_path = output_folder / filename
        logger.info(f"Will save PSM file as: {filename}")

        logger.info("Initializing MaxQuant PSM converter...")
        mq = MaxQuant()

        logger.info(f"Starting PSM conversion (batch size: {batch_size:,})...")
        mq.write_psm_to_file(
            msms_path=str(msms_file), output_path=str(output_path), chunksize=batch_size
        )
        logger.info(f"PSM file successfully saved to: {output_path}")

    except Exception as e:
        logger.error(f"Error in MaxQuant PSM conversion: {str(e)}", exc_info=True)
        raise click.ClickException(f"Error: {str(e)}\nCheck the logs for more details.")


@convert.command(
    "maxquant-feature",
    short_help="Convert feature data from MaxQuant evidence.txt to parquet format",
)
@click.option(
    "--evidence-file",
    help="MaxQuant evidence.txt file to extract peptide information",
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
    "--output-folder",
    help="Output directory for generated files",
    required=True,
    type=click.Path(file_okay=False, path_type=Path),
)
@click.option(
    "--protein-file",
    help="Protein file with specific requirements",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "--partitions",
    help="Field(s) used for splitting files (comma-separated)",
)
@click.option(
    "--batch-size",
    help="Read batch size",
    default=1000000,
    type=int,
)
@click.option(
    "--output-prefix",
    help="Prefix for output files",
)
@click.option("--verbose", help="Enable verbose logging", is_flag=True)
def convert_maxquant_feature_cmd(
    evidence_file: Path,
    sdrf_file: Path,
    output_folder: Path,
    protein_file: Optional[Path],
    partitions: Optional[str],
    batch_size: int,
    output_prefix: Optional[str],
    verbose: bool = False,
):
    """
    Convert MaxQuant feature data from evidence.txt to parquet format.

    This command takes a MaxQuant evidence.txt file and converts it to the quantms.io
    parquet format for feature data, using metadata from an SDRF file.

    Example:
        quantmsio convert maxquant-feature \\
            --evidence-file evidence.txt \\
            --sdrf-file data.sdrf.tsv \\
            --output-folder ./output \\
            --batch-size 1000000
    """
    logger = get_logger("quantmsio.commands.maxquant")
    if verbose:
        logger.setLevel(logging.DEBUG)
        logger.debug("Verbose logging enabled")

    try:
        if not all([evidence_file, sdrf_file, output_folder]):
            raise click.UsageError("ERROR: Please provide all required parameters")

        # Ensure output directory exists
        output_folder = Path(output_folder)
        output_folder.mkdir(parents=True, exist_ok=True)
        logger.info(f"Using output directory: {output_folder}")

        # Set default prefix if not provided
        prefix = output_prefix or "feature"
        filename = create_uuid_filename(prefix, ".feature.parquet")
        output_path = output_folder / filename
        logger.info(f"Will save feature file as: {filename}")

        logger.info("Initializing MaxQuant feature converter...")
        mq = MaxQuant()

        if not partitions:
            logger.info(f"Starting feature conversion (batch size: {batch_size:,})...")
            mq.write_feature_to_file(
                evidence_path=str(evidence_file),
                sdrf_path=str(sdrf_file),
                output_path=str(output_path),
                chunksize=batch_size,
                protein_file=str(protein_file) if protein_file else None,
            )
            logger.info(f"Feature file successfully saved to: {output_path}")
        else:
            logger.info(f"Starting partitioned feature conversion using: {partitions}")
            partition_list = partitions.split(",")
            mq.write_features_to_file(
                evidence_path=str(evidence_file),
                sdrf_path=str(sdrf_file),
                output_folder=str(output_folder),
                filename=filename,
                partitions=partition_list,
                chunksize=batch_size,
                protein_file=str(protein_file) if protein_file else None,
            )
            logger.info(
                f"Partitioned feature files successfully saved to: {output_folder}"
            )

    except Exception as e:
        logger.error(f"Error in MaxQuant feature conversion: {str(e)}", exc_info=True)
        raise click.ClickException(f"Error: {str(e)}\nCheck the logs for more details.")


@convert.command(
    "maxquant-pg",
    short_help="Convert MaxQuant proteinGroups.txt to quantms.io protein group format",
)
@click.option(
    "--protein-groups-file",
    help="MaxQuant proteinGroups.txt file",
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
    "--output-folder",
    help="Output directory for generated files",
    required=True,
    type=click.Path(file_okay=False, path_type=Path),
)
@click.option(
    "--protein-file",
    help="Protein file with specific requirements",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "--batch-size",
    help="Read batch size",
    default=1000000,
    type=int,
)
@click.option(
    "--output-prefix",
    help="Prefix for output files (will be appended with .pg.parquet)",
)
@click.option("--verbose", help="Enable verbose logging", is_flag=True)
def convert_maxquant_pg_cmd(
    protein_groups_file: Path,
    sdrf_file: Path,
    output_folder: Path,
    protein_file: Optional[Path],
    batch_size: int,
    output_prefix: Optional[str],
    verbose: bool = False,
):
    """
    Convert MaxQuant protein groups from proteinGroups.txt to parquet format.

    This command takes a MaxQuant proteinGroups.txt file and converts it to the quantms.io
    parquet format for protein groups, using metadata from an SDRF file.

    Example:
        quantmsio convert maxquant-pg \\
            --protein-groups-file proteinGroups.txt \\
            --sdrf-file data.sdrf.tsv \\
            --output-folder ./output \\
            --batch-size 1000000
    """
    logger = get_logger("quantmsio.commands.maxquant")
    if verbose:
        logger.setLevel(logging.DEBUG)
        logger.debug("Verbose logging enabled")

    try:
        if not all([protein_groups_file, sdrf_file, output_folder]):
            raise click.UsageError("ERROR: Please provide all required parameters")

        # Ensure output directory exists
        output_folder = Path(output_folder)
        output_folder.mkdir(parents=True, exist_ok=True)
        logger.info(f"Using output directory: {output_folder}")

        # Set default prefix if not provided
        prefix = output_prefix or "pg"
        filename = create_uuid_filename(prefix, ".pg.parquet")
        output_path = output_folder / filename
        logger.info(f"Will save protein groups file as: {filename}")

        logger.info("Initializing MaxQuant protein groups converter...")
        mq = MaxQuant()

        logger.info(
            f"Starting protein groups conversion (batch size: {batch_size:,})..."
        )
        mq.write_protein_groups_to_file(
            protein_groups_path=str(protein_groups_file),
            sdrf_path=str(sdrf_file),
            output_path=str(output_path),
            chunksize=batch_size,
            protein_file=str(protein_file) if protein_file else None,
        )
        logger.info(f"Protein groups file successfully saved to: {output_path}")

    except Exception as e:
        logger.error(
            f"Error in MaxQuant protein groups conversion: {str(e)}", exc_info=True
        )
        raise click.ClickException(f"Error: {str(e)}\nCheck the logs for more details.")
