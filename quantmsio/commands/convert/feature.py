"""
Command to convert msstats/mztab data to parquet format.
"""

import click
import logging
from pathlib import Path
from typing import Optional

from quantmsio.core.feature import Feature
from quantmsio.core.project import create_uuid_filename
from quantmsio.utils.logger import get_logger


def convert_feature(
    sdrf_file: Path,
    msstats_file: Path,
    mztab_file: Path,
    output_folder: Path,
    batch_size: int = 50,
    protein_file: Optional[Path] = None,
    partitions: Optional[str] = None,
    output_prefix: Optional[str] = None,
    duckdb_max_memory: Optional[str] = None,
    duckdb_threads: Optional[int] = None,
    verbose: bool = False,
) -> Path:
    """
    Convert msstats/mztab data to parquet format.

    Args:
        sdrf_file: SDRF file needed to extract metadata
        msstats_file: MSstats input file (main format to convert)
        mztab_file: mzTab file (used to extract protein information)
        output_folder: Output directory for generated files
        batch_size: Read batch size (default: 50)
        protein_file: Optional protein file with specific requirements
        partitions: Optional field(s) used for splitting files (comma-separated)
        output_prefix: Optional prefix for output files
        duckdb_max_memory: Optional maximum memory for DuckDB (e.g., "4GB")
        duckdb_threads: Optional number of threads for DuckDB
        verbose: Enable verbose logging

    Returns:
        Path: The path to the generated feature file
    """
    logger = get_logger("quantmsio.commands.feature")
    if verbose:
        logger.setLevel(logging.DEBUG)
        logger.debug("🔍 Verbose logging enabled")

    try:
        if not all([sdrf_file, msstats_file, mztab_file, output_folder]):
            raise click.UsageError("❌ Please provide all required parameters")

        # Ensure output directory exists
        output_folder = Path(output_folder)
        output_folder.mkdir(parents=True, exist_ok=True)
        logger.info(f"📂 Using output directory: {output_folder}")

        logger.info("🔄 Initializing feature manager...")
        feature_manager = Feature(
            mztab_path=mztab_file, sdrf_path=sdrf_file, msstats_in_path=msstats_file
        )

        # Set default prefix if not provided
        prefix = output_prefix or "feature"
        filename = create_uuid_filename(prefix, ".feature.parquet")
        output_path = output_folder / filename
        logger.info(f"📄 Will save feature file as: {filename}")
        logger.info(f"📍 Full output path: {output_path}")

        if not partitions:
            logger.info("🔄 Starting feature conversion (no partitions)...")
            feature_manager.write_feature_to_file(
                output_path=str(output_path),
                file_num=batch_size,
                protein_file=protein_file,
                duckdb_max_memory=duckdb_max_memory,
                duckdb_threads=duckdb_threads,
            )
            logger.info(f"✅ Feature file successfully saved to: {output_path}")
            return output_path
        else:
            logger.info(
                f"🔄 Starting partitioned feature conversion using: {partitions}"
            )
            partition_list = partitions.split(",")
            feature_manager.write_features_to_file(
                output_folder=str(output_folder),
                filename=filename,
                partitions=partition_list,
                file_num=batch_size,
                protein_file=protein_file,
                duckdb_max_memory=duckdb_max_memory,
                duckdb_threads=duckdb_threads,
            )
            logger.info(
                f"✅ Partitioned feature files successfully saved to: {output_folder}"
            )
            return output_path

    except Exception as e:
        logger.error(f"❌ Error in feature conversion: {str(e)}", exc_info=True)
        raise click.ClickException(
            f"❌ Error: {str(e)}\nCheck the logs for more details."
        )


@click.command(
    "feature",
    short_help="Convert msstats/mztab data to parquet format",
)
@click.option(
    "--sdrf-file",
    help="SDRF file needed to extract metadata",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "--msstats-file",
    help="MSstats input file (main format to convert)",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "--mztab-file",
    help="mzTab file (used to extract protein information)",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option("--batch-size", help="Read batch size", default=50, type=int)
@click.option(
    "--protein-file",
    help="Protein file with specific requirements",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "--output-folder",
    help="Output directory for generated files",
    required=True,
    type=click.Path(file_okay=False, path_type=Path),
)
@click.option(
    "--partitions", help="Field(s) used for splitting files (comma-separated)"
)
@click.option("--output-prefix", help="Prefix for output files")
@click.option("--duckdb-max-memory", help="Maximum memory for DuckDB (e.g., '4GB')")
@click.option("--duckdb-threads", help="Number of threads for DuckDB", type=int)
@click.option("--verbose", help="Enable verbose logging", is_flag=True)
def convert_feature_cmd(**kwargs):
    """
    Convert msstats/mztab data to parquet format.
    
    This command takes MSstats and mzTab files and converts them to the quantms.io
    parquet format. The mzTab file is used to extract protein information while
    the MSstats file provides the main data to convert.
    
    Example:
        quantmsio convert feature \\
            --msstats-file data.msstats.txt \\
            --mztab-file data.mztab \\
            --sdrf-file data.sdrf.tsv \\
            --output-folder ./output \\
            --file-num 50
    """
    convert_feature(**kwargs)
