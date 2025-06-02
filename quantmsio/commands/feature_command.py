"""
Command to convert feature data from various formats to quantms.io format.
"""

import click
from pathlib import Path
from typing import Optional, List
import logging
import time

from quantmsio.core.feature import Feature
from quantmsio.core.project import create_uuid_filename
from quantmsio.utils.logger import get_logger


def convert_feature(
    sdrf_file: Path,
    msstats_file: Path,
    mztab_file: Path,
    output_folder: Path,
    file_num: int = 50,
    protein_file: Optional[Path] = None,
    partitions: Optional[List[str]] = None,
    output_prefix: str = "feature",
    duckdb_max_memory: Optional[str] = None,
    duckdb_threads: Optional[int] = None,
    verbose: bool = False,
) -> None:
    """
    Convert feature data from MSstats/mzTab to parquet format.
    
    Args:
        sdrf_file: SDRF file path
        msstats_file: MSstats input file path
        mztab_file: mzTab file path
        output_folder: Output directory
        file_num: Read batch size
        protein_file: Optional protein file path
        partitions: Optional list of fields to partition by
        output_prefix: Prefix for output files
        duckdb_max_memory: Maximum memory for DuckDB
        duckdb_threads: Number of threads for DuckDB
        verbose: Enable verbose logging
    """
    logger = get_logger("quantmsio.commands.feature")
    if verbose:
        logger.setLevel(logging.DEBUG)
    
    start_time = time.time()
    logger.info("Starting feature conversion")
    
    try:
        # Log input parameters
        logger.debug(f"Input parameters:")
        logger.debug(f"  SDRF file: {sdrf_file}")
        logger.debug(f"  MSstats file: {msstats_file}")
        logger.debug(f"  mzTab file: {mztab_file}")
        logger.debug(f"  Output folder: {output_folder}")
        logger.debug(f"  Protein file: {protein_file}")
        logger.debug(f"  Batch size: {file_num}")
        logger.debug(f"  Partitions: {partitions}")
        logger.debug(f"  DuckDB max memory: {duckdb_max_memory}")
        logger.debug(f"  DuckDB threads: {duckdb_threads}")
        
        # Initialize feature manager
        logger.info("Initializing feature manager...")
        feature_manager = Feature(
            mztab_path=str(mztab_file),
            sdrf_path=str(sdrf_file),
            msstats_in_path=str(msstats_file)
        )
        logger.debug("Feature manager initialized successfully")
        
        # Generate output filename
        filename = create_uuid_filename(output_prefix, ".feature.parquet")
        output_path = output_folder / filename
        logger.debug(f"Generated output filename: {filename}")
        
        # Convert features
        if not partitions:
            logger.info(f"Writing features to single file: {output_path}")
            logger.debug("Starting feature conversion to single file...")
            
            # Add progress callback for feature manager
            def progress_callback(current: int, total: int, message: str):
                if verbose:
                    percent = (current / total) * 100 if total > 0 else 0
                    logger.debug(f"Progress: {percent:.1f}% - {message}")
            
            feature_manager.write_feature_to_file(
                output_path=str(output_path),
                file_num=file_num,
                protein_file=str(protein_file) if protein_file else None,
                duckdb_max_memory=duckdb_max_memory,
                duckdb_threads=duckdb_threads,
                progress_callback=progress_callback,
            )
            logger.info("Feature conversion to single file completed")
        else:
            logger.info(f"Writing features to partitioned files in: {output_folder}")
            logger.debug("Starting feature conversion to partitioned files...")
            feature_manager.write_features_to_file(
                output_folder=str(output_folder),
                filename=filename,
                partitions=partitions,
                file_num=file_num,
                protein_file=str(protein_file) if protein_file else None,
                duckdb_max_memory=duckdb_max_memory,
                duckdb_threads=duckdb_threads,
                progress_callback=progress_callback,
            )
            logger.info("Feature conversion to partitioned files completed")
        
        elapsed = time.time() - start_time
        logger.info(f"Completed feature conversion in {elapsed:.2f} seconds")
        
    except Exception as e:
        logger.exception(f"Error in feature conversion: {str(e)}")
        raise click.ClickException(f"Error: {str(e)}\nCheck the logs for more details.")


@click.command("features", short_help="Convert msstats/mztab to parquet file")
@click.option(
    "--sdrf-file",
    help="SDRF file",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "--msstats-file",
    help="MSstats input file",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "--mztab-file",
    help="mzTab file",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "--output-folder",
    help="Output directory",
    required=True,
    type=click.Path(file_okay=False, path_type=Path),
)
@click.option(
    "--file-num",
    help="Read batch size",
    type=click.INT,
    default=50,
)
@click.option(
    "--protein-file",
    help="Protein file",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "--partitions",
    help="Fields to partition by",
    multiple=True,
)
@click.option(
    "--output-prefix",
    help="Prefix for output files",
    default="feature",
)
@click.option("--duckdb-max-memory", help="Maximum memory for DuckDB")
@click.option("--duckdb-threads", help="Number of threads for DuckDB", type=click.INT)
@click.option("--verbose", help="Enable verbose logging", is_flag=True)
def convert_feature_file(**kwargs):
    """
    Convert feature data from MSstats/mzTab format to quantms.io parquet format.
    
    This command combines data from MSstats and mzTab files, using metadata from
    an SDRF file to create a standardized feature file in parquet format.
    
    The output can be optionally partitioned by specified fields and includes
    comprehensive protein information.
    
    Example:
        quantmsio features \\
            --sdrf-file data.sdrf.tsv \\
            --msstats-file data.msstats.txt \\
            --mztab-file data.mztab \\
            --output-folder ./output \\
            --partitions "sample,fraction"
    """
    convert_feature(**kwargs)
