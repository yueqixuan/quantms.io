"""
Command to convert feature data from various formats to quantms.io format.
"""

import click
from pathlib import Path
from typing import Optional, List
import logging

from quantmsio.commands.base_command import BaseCommand
from quantmsio.core.feature import Feature
from quantmsio.core.project import create_uuid_filename
from quantmsio.utils.logger import with_request_tracking


class FeatureCommand(BaseCommand):
    """
    Command to convert feature data from various formats to quantms.io format.
    Supports MSstats and mzTab input formats.
    """

    def __init__(self):
        super().__init__("feature")

    @with_request_tracking
    @BaseCommand.with_error_handling
    def convert_feature(
        self,
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
        """
        self.start_progress("feature conversion")
        
        # Initialize feature manager
        feature_manager = Feature(
            mztab_path=str(mztab_file),
            sdrf_path=str(sdrf_file),
            msstats_in_path=str(msstats_file)
        )
        
        # Generate output filename
        filename = create_uuid_filename(output_prefix, ".feature.parquet")
        output_path = output_folder / filename
        
        # Convert features
        if not partitions:
            self.logger.info(f"Writing features to single file: {output_path}")
            feature_manager.write_feature_to_file(
                output_path=str(output_path),
                file_num=file_num,
                protein_file=str(protein_file) if protein_file else None,
                duckdb_max_memory=duckdb_max_memory,
                duckdb_threads=duckdb_threads,
            )
        else:
            self.logger.info(f"Writing features to partitioned files in: {output_folder}")
            feature_manager.write_features_to_file(
                output_folder=str(output_folder),
                filename=filename,
                partitions=partitions,
                file_num=file_num,
                protein_file=str(protein_file) if protein_file else None,
                duckdb_max_memory=duckdb_max_memory,
                duckdb_threads=duckdb_threads,
            )
        
        self.end_progress("feature conversion")


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
@BaseCommand.common_options
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
    cmd = FeatureCommand()
    cmd.convert_feature(**kwargs)
