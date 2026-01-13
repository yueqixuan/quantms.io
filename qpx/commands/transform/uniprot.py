from pathlib import Path

import click

from qpx.core.project import create_uuid_filename
from qpx.operate.tools import map_peptide_to_protein


@click.command(
    "uniprot",
    short_help="Map feature data to latest UniProt version",
)
@click.option(
    "--feature-file",
    help="Feature file path",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "--fasta",
    help="UniProt FASTA file path",
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
    required=False,
)
def map_latest_uniprot_cmd(
    feature_file: Path,
    fasta: Path,
    output_folder: Path,
    output_prefix: str,
):
    """Map feature data to latest UniProt version.
    
    Maps peptides and features to the latest UniProt protein database using a FASTA file. 
    This command updates protein identifications to match current UniProt accessions and annotations.
    
    Example:
        qpxc transform uniprot \\
            --feature-file ./output/feature.parquet \\
            --fasta uniprot_human.fasta \\
            --output-folder ./output
    """
    if not all([feature_file, fasta, output_folder]):
        raise click.UsageError("Please provide all required parameters")

    if not output_prefix:
        output_prefix = "feature"

    filename = create_uuid_filename(output_prefix, ".feature.parquet")
    output_path = output_folder / filename
    map_peptide_to_protein(str(feature_file), str(fasta), str(output_folder), filename)
