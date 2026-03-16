"""
Transform subcommands — compute derived data from QPX datasets.

Subcommands:
    qpxc transform gene-map   — Map genes from FASTA
"""

import logging
from pathlib import Path
from typing import Optional

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
    protein_col = "pg_accessions" if "pg_accessions" in df.columns else "protein_accessions"
    annotated = mapping.annotate_dataframe(df, protein_col=protein_col)
    output_path = output_folder / parquet_path.name
    annotated.to_parquet(str(output_path), index=False)

    click.echo(f"Gene mapping complete. Output: {output_path}")


# ---------------------------------------------------------------------------
# Protein Quantification via mokume
# ---------------------------------------------------------------------------

QUANT_METHODS = ["directlfq", "maxlfq", "topn", "top3", "ibaq", "sum"]
_SUPPORTED_SUFFIXES = {".parquet", ".tsv", ".csv"}


def _validate_quantify_inputs(method_lower: str, fasta: Optional[Path]) -> None:
    """Validate CLI inputs before running quantification."""
    if method_lower == "ibaq" and not fasta:
        raise click.UsageError("The --fasta option is required for the ibaq method")
    try:
        from mokume.quantification import is_directlfq_available
    except ImportError:
        raise click.UsageError("mokume is not installed. Install with: pip install mokume")
    if method_lower == "directlfq" and not is_directlfq_available():
        raise click.UsageError("DirectLFQ is not installed. Install with: pip install mokume[directlfq]")


def _run_ibaq(peptide_df, fasta, enzyme, normalize, min_aa, max_aa, ploidy, cpc, organism, output, verbose) -> None:
    """Run iBAQ quantification via mokume."""
    import tempfile

    from mokume.quantification import peptides_to_protein

    click.echo("Running iBAQ quantification ...")
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp_path = str(Path(tmpdir) / "peptides.parquet")
        peptide_df.to_parquet(tmp_path, index=False)
        peptides_to_protein(
            fasta=str(fasta),
            peptides=tmp_path,
            enzyme=enzyme,
            normalize=normalize,
            min_aa=min_aa,
            max_aa=max_aa,
            tpa=False,
            ruler=False,
            ploidy=ploidy,
            cpc=cpc,
            organism=organism,
            output=str(output),
            verbose=verbose,
            qc_report=str(output.with_suffix(".qc.pdf")),
        )
    click.echo(f"iBAQ results saved to {output}")


def _run_generic_quantify(peptide_df, method_lower, method, topn_n, threads, normalize):
    """Run non-iBAQ quantification and return result DataFrame."""
    from mokume.quantification import get_quantification_method

    click.echo(f"Using {method} quantification method")

    method_kwargs = {"topn": {"n": topn_n}, "top3": {"n": 3}, "maxlfq": {"threads": threads}}
    resolved = "topn" if method_lower == "top3" else method_lower
    quant_method = get_quantification_method(resolved, **method_kwargs.get(method_lower, {}))

    click.echo(
        f"Quantifying {peptide_df['ProteinName'].nunique()} proteins across {peptide_df['SampleID'].nunique()} samples ..."
    )
    result_df = quant_method.quantify(
        peptide_df,
        protein_column="ProteinName",
        peptide_column="PeptideCanonical",
        intensity_column="NormIntensity",
        sample_column="SampleID",
    )

    if normalize:
        intensity_cols = [c for c in result_df.columns if "Intensity" in c]
        if intensity_cols:
            int_col = intensity_cols[-1]
            result_df[f"{int_col}Norm"] = result_df.groupby("SampleID")[int_col].transform(lambda x: x / x.sum())
    return result_df


def _save_result(result_df, output: Path) -> None:
    """Save quantification results to the specified output format."""
    suffix = output.suffix.lower()
    if suffix not in _SUPPORTED_SUFFIXES:
        raise click.UsageError(f"Unsupported output format '{suffix}'. Use .parquet, .tsv, or .csv")
    output.parent.mkdir(parents=True, exist_ok=True)
    if suffix == ".parquet":
        result_df.to_parquet(str(output), index=False)
    elif suffix == ".tsv":
        result_df.to_csv(str(output), sep="\t", index=False)
    else:
        result_df.to_csv(str(output), index=False)


def _qpx_feature_to_peptide_df(feature_path: Path):
    """Read QPX feature.parquet and flatten to mokume peptide DataFrame.

    Mokume expects columns: ``ProteinName``, ``PeptideCanonical``,
    ``NormIntensity``, ``SampleID``.

    QPX stores intensities as a list of ``{label, intensity}`` structs.
    For label-free data each feature typically has a single intensity entry;
    for TMT/iTRAQ data the list contains one entry per channel.  This
    function explodes **all** intensity entries into separate rows so that
    multiplexed channels are preserved.  ``SampleID`` is set to
    ``<run_file_name>::<label>`` when multiple labels exist, or just
    ``<run_file_name>`` for single-label (LFQ) data.
    """
    import pandas as pd

    df = pd.read_parquet(str(feature_path))
    logger.info("Loaded %d QPX feature rows from %s", len(df), feature_path)

    # Filter decoys early
    if "is_decoy" in df.columns:
        df = df[~df["is_decoy"].astype(bool)]

    # Map protein column
    protein_col = "anchor_protein" if "anchor_protein" in df.columns else "pg_accessions"
    df["ProteinName"] = df[protein_col].apply(
        lambda x: x if isinstance(x, str) else (";".join(x) if isinstance(x, list) else str(x))
    )
    df["PeptideCanonical"] = df.get("sequence", df.get("peptidoform", ""))

    # Drop rows without intensities, then explode vectorised
    df = df[df["intensities"].apply(lambda x: x is not None and len(x) > 0)]
    if df.empty:
        return pd.DataFrame(columns=["ProteinName", "PeptideCanonical", "NormIntensity", "SampleID"])

    exploded = df[["ProteinName", "PeptideCanonical", "run_file_name", "intensities"]].explode("intensities")
    # Extract struct fields
    exploded["NormIntensity"] = exploded["intensities"].apply(lambda e: e.get("intensity", 0.0) if isinstance(e, dict) else 0.0)
    exploded["_label"] = exploded["intensities"].apply(lambda e: e.get("label", "") if isinstance(e, dict) else "")
    # Build SampleID: run::label when label present, else run
    has_label = exploded["_label"].astype(bool)
    exploded["SampleID"] = exploded["run_file_name"]
    exploded.loc[has_label, "SampleID"] = exploded.loc[has_label, "run_file_name"] + "::" + exploded.loc[has_label, "_label"]
    # Filter invalid intensities
    result = exploded.loc[
        exploded["NormIntensity"].notna() & (exploded["NormIntensity"] > 0),
        ["ProteinName", "PeptideCanonical", "NormIntensity", "SampleID"],
    ].copy()
    result = result[result["ProteinName"] != ""]

    logger.info(
        "Prepared %d peptide rows: %d proteins, %d samples",
        len(result),
        result["ProteinName"].nunique(),
        result["SampleID"].nunique(),
    )
    return result


@transform.command("quantify")
@click.option(
    "--feature-path",
    help="QPX feature.parquet file path",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "--method",
    help="Quantification method (directlfq, maxlfq, topn, top3, ibaq, sum)",
    type=click.Choice(QUANT_METHODS, case_sensitive=False),
    default="directlfq",
)
@click.option(
    "--fasta",
    help="FASTA database (required for ibaq method)",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=None,
)
@click.option(
    "--enzyme",
    help="Enzyme for iBAQ digestion (default: Trypsin)",
    default="Trypsin",
)
@click.option(
    "--topn-n",
    help="N for TopN method (default: 3)",
    default=3,
    type=int,
)
@click.option(
    "--threads",
    help="Parallel threads for MaxLFQ (-1 = all cores)",
    default=-1,
    type=int,
)
@click.option(
    "--output",
    "-o",
    help="Output file path (.parquet, .tsv, or .csv)",
    required=True,
    type=click.Path(path_type=Path),
)
@click.option("--normalize", help="Normalize quantification values", is_flag=True)
@click.option("--organism", help="Organism for iBAQ (default: human)", default="human")
@click.option("--ploidy", help="Ploidy for iBAQ ruler (default: 2)", default=2, type=int)
@click.option("--cpc", help="Cell copies per cell for iBAQ ruler (default: 200)", default=200, type=float)
@click.option("--min-aa", help="Min peptide length for iBAQ (default: 7)", default=7, type=int)
@click.option("--max-aa", help="Max peptide length for iBAQ (default: 30)", default=30, type=int)
@click.option("--verbose", help="Enable verbose logging", is_flag=True)
def transform_quantify_cmd(
    feature_path: Path,
    method: str,
    fasta: Optional[Path],
    enzyme: str,
    topn_n: int,
    threads: int,
    output: Path,
    normalize: bool,
    organism: str,
    ploidy: int,
    cpc: float,
    min_aa: int,
    max_aa: int,
    verbose: bool,
):
    """Compute protein quantification from QPX feature data using mokume.

    Reads a QPX feature.parquet file, extracts peptide-level intensities,
    and computes protein-level quantification using the selected method.

    \b
    Supported methods:
      directlfq  — DirectLFQ intensity traces (default)
      maxlfq     — MaxLFQ delayed normalization
      topn       — Average of N most intense peptides
      top3       — Average of 3 most intense peptides
      ibaq       — Intensity-Based Absolute Quantification (requires --fasta)
      sum        — Sum of all peptide intensities

    \b
    Examples:
        # DirectLFQ from QPX feature
        qpxc transform quantify \\
            --feature-path ./qpx_output/feature.parquet \\
            --method directlfq \\
            -o proteins_directlfq.parquet

        # iBAQ (requires FASTA)
        qpxc transform quantify \\
            --feature-path ./qpx_output/feature.parquet \\
            --method ibaq --fasta proteome.fasta \\
            -o proteins_ibaq.tsv

        # MaxLFQ with 8 threads
        qpxc transform quantify \\
            --feature-path ./qpx_output/feature.parquet \\
            --method maxlfq --threads 8 \\
            -o proteins_maxlfq.parquet
    """
    if verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    method_lower = method.lower()
    _validate_quantify_inputs(method_lower, fasta)

    # Step 1: Read QPX feature and flatten to mokume format
    click.echo(f"Reading QPX feature data from {feature_path}")
    peptide_df = _qpx_feature_to_peptide_df(feature_path)
    if peptide_df.empty:
        raise click.UsageError("No valid peptide intensities found in the feature file")

    output.parent.mkdir(parents=True, exist_ok=True)

    # Step 2: Run quantification
    if method_lower == "ibaq":
        _run_ibaq(peptide_df, fasta, enzyme, normalize, min_aa, max_aa, ploidy, cpc, organism, output, verbose)
        return

    result_df = _run_generic_quantify(peptide_df, method_lower, method, topn_n, threads, normalize)

    # Step 3: Save output
    _save_result(result_df, output)
    click.echo(f"Quantification complete: {result_df['ProteinName'].nunique()} proteins -> {output}")
