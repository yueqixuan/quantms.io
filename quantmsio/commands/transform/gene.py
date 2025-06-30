import click

from quantmsio.operate.tools import generate_feature_of_gene


@click.command(
    "gene",
    short_help="Map gene information from FASTA to parquet format",
)
@click.option(
    "--parquet-path",
    help="PSM or feature parquet file path",
    required=True,
    type=click.Path(exists=True, dir_okay=False),
)
@click.option(
    "--fasta",
    help="FASTA file path",
    required=True,
    type=click.Path(exists=True, dir_okay=False),
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
@click.option(
    "--species",
    help="Species name (default: human)",
    default="human",
    required=False,
)
def map_gene_message_cmd(
    parquet_path: str,
    fasta: str,
    output_folder: str,
    file_num: int,
    partitions: str = None,
    species: str = "human",
):
    """Map gene information from FASTA file to parquet format.

    Args:
        parquet_path: PSM or feature parquet file path
        fasta: FASTA file path
        output_folder: Output directory for generated files
        file_num: Number of rows to read in each batch
        partitions: Optional fields for splitting files (comma-separated)
        species: Species name (default: human)
    """
    if partitions:
        partitions = partitions.split(",")
    generate_feature_of_gene(
        parquet_path, fasta, output_folder, file_num, partitions, species
    )
