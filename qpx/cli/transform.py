"""
Transform subcommands — compute derived data from QPX datasets.

Subcommands:
    qpxc transform gene-map  — Map genes from FASTA
"""

import logging
from pathlib import Path

import click

logger = logging.getLogger("qpx.cli.transform")


# ---------------------------------------------------------------------------
# Transform group
# ---------------------------------------------------------------------------


@click.group()
def transform():
    """Transform QPX data into derived representations."""
    pass


# ---------------------------------------------------------------------------
# Gene Mapping
# ---------------------------------------------------------------------------


@transform.command("gene-map")
@click.option(
    "--parquet-path",
    help="QPX PSM or feature parquet file path",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "--fasta",
    help="FASTA database file path",
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
    "--species",
    help="Species name for gene mapping",
    default="human",
)
@click.option("--verbose", help="Enable verbose logging", is_flag=True)
def transform_gene_map_cmd(
    parquet_path: Path,
    fasta: Path,
    output_folder: Path,
    species: str,
    verbose: bool,
):
    """Map gene names from a FASTA file to QPX parquet data.

    Enriches protein identifications in QPX PSM or feature files with
    gene-level metadata extracted from FASTA database headers.

    \b
    Example:
        qpxc transform gene-map \\
            --parquet-path ./output/psm.parquet \\
            --fasta proteins.fasta \\
            --output-folder ./output \\
            --species human
    """
    if verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    output_folder = Path(output_folder)
    output_folder.mkdir(parents=True, exist_ok=True)

    import pandas as pd
    from qpx.transforms.gene_mapping import GeneMappingTransform

    mapping = GeneMappingTransform(fasta_path=str(fasta), species=species)
    df = pd.read_parquet(str(parquet_path))
    protein_col = (
        "pg_accessions" if "pg_accessions" in df.columns else "protein_accessions"
    )
    annotated = mapping.annotate_dataframe(df, protein_col=protein_col)
    output_path = output_folder / parquet_path.name
    annotated.to_parquet(str(output_path), index=False)

    click.echo(f"Gene mapping complete. Output: {output_path}")
