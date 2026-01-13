import click

from qpx.operate.tools import generate_psms_of_spectrum


@click.command(
    "spectra",
    short_help="Map spectrum information from mzML to parquet format",
)
@click.option(
    "--parquet-path",
    help="PSM or feature parquet file path",
    required=True,
    type=click.Path(exists=True, dir_okay=False),
)
@click.option(
    "--mzml-directory",
    help="Directory containing mzML files",
    required=True,
    type=click.Path(exists=True, file_okay=False),
)
@click.option(
    "--output-folder",
    help="Output directory for generated files",
    required=True,
    type=click.Path(file_okay=False),
)
@click.option(
    "--file-num",
    help="Number of rows to read in each batch",
    default=10,
    type=int,
)
@click.option(
    "--partitions",
    help="Fields for splitting files (comma-separated)",
    required=False,
)
def map_spectrum_message_cmd(
    parquet_path: str,
    mzml_directory: str,
    output_folder: str,
    file_num: int,
    partitions: str = None,
):
    """Map spectrum information from mzML to parquet format.
    
    Enriches PSM or feature data with additional spectral information extracted from mzML files. 
    This command maps spectrum metadata and peak information to the corresponding peptide-spectrum matches.
    
    Example:
        qpxc transform spectra \\
            --parquet-path ./output/psm.parquet \\
            --mzml-directory ./mzml_files \\
            --output-folder ./output
    """
    if partitions:
        partitions = partitions.split(",")
    generate_psms_of_spectrum(
        parquet_path, mzml_directory, output_folder, file_num, partitions
    )
