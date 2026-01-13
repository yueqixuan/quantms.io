from pathlib import Path

import click

from qpx.core.project import create_uuid_filename
from qpx.operate.tools import write_ibaq_feature


@click.command(
    "ibaq",
    short_help="Convert feature data to IBAQ format",
)
@click.option(
    "--feature-file",
    help="Feature file path",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "--sdrf-file",
    help="SDRF file for metadata extraction",
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
    help="Prefix for output files",
    default="ibaq",
)
def convert_ibaq_file_cmd(
    feature_file: Path,
    sdrf_file: Path,
    output_folder: Path,
    output_prefix: str,
):
    """Convert feature data to iBAQ format.
    
    Transforms feature-level quantification data into iBAQ (intensity-Based Absolute Quantification) 
    format. This command integrates feature quantification with sample metadata from SDRF files to 
    generate iBAQ values.
    
    Example:
        qpxc transform ibaq \\
            --feature-file ./output/feature.parquet \\
            --sdrf-file ./metadata.sdrf.tsv \\
            --output-folder ./output
    """
    if not all([feature_file, sdrf_file, output_folder]):
        raise click.UsageError("Please provide all required parameters")

    # Create output directory if it doesn't exist
    output_folder.mkdir(parents=True, exist_ok=True)

    output_path = output_folder / create_uuid_filename(output_prefix, ".ibaq.parquet")
    write_ibaq_feature(str(sdrf_file), str(feature_file), str(output_path))
