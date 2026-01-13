import click

from qpx.operate.tools import generate_feature_of_gene


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
    """Map gene information from FASTA to parquet format.
    
    Maps gene names and information from a FASTA file to protein identifications in QPX 
    PSM or feature files. This command enriches protein data with gene-level metadata extracted 
    from FASTA headers.
    
    Example:
        qpxc transform gene \\
            --parquet-path ./output/psm.parquet \\
            --fasta proteins.fasta \\
            --output-folder ./output
    """
    if partitions:
        partitions = partitions.split(",")
    generate_feature_of_gene(
        parquet_path, fasta, output_folder, file_num, partitions, species
    )
