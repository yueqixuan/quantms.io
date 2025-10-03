from pathlib import Path
from typing import Optional

import click

from quantmsio.core.fragpipe.fragpipe import FragPipe


@click.command(
    "fragpipe-psm",
    short_help="Convert FragPipe PSMs from psm.tsv to parquet file in quantms.io",
)
@click.option(
    "--msms-file",
    type=click.Path(dir_okay=False, path_type=Path, file_okay=True, exists=True),
    help="the psm.tsv file, this will be used to extract the peptide information",
    required=True,
)
@click.option(
    "-o",
    "--output-folder",
    type=click.Path(dir_okay=True, path_type=Path, file_okay=False),
    help="Folder where the parquet file will be generated",
    required=True,
)
@click.option("-b", "--batch-size", help="Read batch size", default=1000000, type=int)
@click.option(
    "--output-prefix",
    help="Prefix of the parquet file needed to generate the file name",
    required=False,
)
def convert_fragpipe_psm_cmd(
    msms_file: Path,
    output_folder: Path,
    batch_size: int,
    output_prefix: Optional[str] = None,
):
    """Convert FragPipe PSMs to parquet format.

    Args:
        msms_file: PSM TSV file to extract peptide information
        output_folder: Output directory for generated files
        batch_size: Read batch size
        output_prefix: Optional prefix for output files
    """
    if not output_folder.exists():
        output_folder.mkdir(parents=True, exist_ok=True)
    converter = FragPipe(output_directory=output_folder)
    converter.write_psms_to_parquet(
        msms_file, batch_size=batch_size, output_prefix=output_prefix
    )
