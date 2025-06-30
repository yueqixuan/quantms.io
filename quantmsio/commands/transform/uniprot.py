from pathlib import Path

import click

from quantmsio.core.project import create_uuid_filename
from quantmsio.operate.tools import map_peptide_to_protein


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

    Args:
        feature_file: Feature file path
        fasta: UniProt FASTA file path
        output_folder: Output directory for generated files
        output_prefix: Optional prefix for output files
    """
    if not all([feature_file, fasta, output_folder]):
        raise click.UsageError("Please provide all required parameters")

    if not output_prefix:
        output_prefix = "feature"

    filename = create_uuid_filename(output_prefix, ".feature.parquet")
    output_path = output_folder / filename
    map_peptide_to_protein(str(feature_file), str(fasta), str(output_folder), filename)
