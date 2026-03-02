"""
Transform subcommands — compute derived data from QPX datasets.

Subcommands:
    qpxc transform gene-map   — Map genes from FASTA
    qpxc transform quantify   — Protein quantification via mokume
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
    protein_col = (
        "pg_accessions" if "pg_accessions" in df.columns else "protein_accessions"
    )
    annotated = mapping.annotate_dataframe(df, protein_col=protein_col)
    output_path = output_folder / parquet_path.name
    annotated.to_parquet(str(output_path), index=False)

    click.echo(f"Gene mapping complete. Output: {output_path}")


# ---------------------------------------------------------------------------
# Protein Quantification via mokume
# ---------------------------------------------------------------------------

QUANT_METHODS = ["directlfq", "maxlfq", "topn", "top3", "ibaq", "sum"]


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
    protein_col = (
        "anchor_protein" if "anchor_protein" in df.columns else "pg_accessions"
    )
    df["ProteinName"] = df[protein_col].apply(
        lambda x: (
            x
            if isinstance(x, str)
            else (";".join(x) if isinstance(x, list) else str(x))
        )
    )
    df["PeptideCanonical"] = df.get("sequence", df.get("peptidoform", ""))

    # Explode intensities list into one row per label/channel
    rows = []
    for _, row in df.iterrows():
        intensities = row.get("intensities")
        if intensities is None:
            continue
        run = row.get("run_file_name", "")
        protein = row["ProteinName"]
        peptide = row["PeptideCanonical"]
        for entry in intensities:
            try:
                intensity = entry["intensity"] if isinstance(entry, dict) else 0.0
            except (TypeError, KeyError):
                continue
            if not intensity or intensity <= 0:
                continue
            label = entry.get("label", "") if isinstance(entry, dict) else ""
            sample_id = f"{run}::{label}" if label else run
            rows.append(
                {
                    "ProteinName": protein,
                    "PeptideCanonical": peptide,
                    "NormIntensity": intensity,
                    "SampleID": sample_id,
                }
            )

    result = pd.DataFrame(rows)
    result = result.dropna(subset=["ProteinName", "NormIntensity"])
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
    help="Output file path (.parquet or .tsv)",
    required=True,
    type=click.Path(path_type=Path),
)
@click.option("--normalize", help="Normalize quantification values", is_flag=True)
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

    # Validate iBAQ requirement
    if method_lower == "ibaq" and not fasta:
        raise click.UsageError("The --fasta option is required for the ibaq method")

    # Import mokume
    try:
        from mokume.quantification import (
            get_quantification_method,
            is_directlfq_available,
        )
    except ImportError:
        raise click.UsageError(
            "mokume is not installed. Install with: pip install mokume"
        )

    if method_lower == "directlfq" and not is_directlfq_available():
        raise click.UsageError(
            "DirectLFQ is not installed. Install with: pip install mokume[directlfq]"
        )

    # Step 1: Read QPX feature and flatten to mokume format
    click.echo(f"Reading QPX feature data from {feature_path}")
    peptide_df = _qpx_feature_to_peptide_df(feature_path)

    if peptide_df.empty:
        raise click.UsageError("No valid peptide intensities found in the feature file")

    # Step 2: Run quantification
    if method_lower == "ibaq":
        # iBAQ uses peptides_to_protein directly (file-based API)
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
                min_aa=7,
                max_aa=30,
                tpa=False,
                ruler=False,
                ploidy=2,
                cpc=200,
                organism="human",
                output=str(output),
                verbose=verbose,
                qc_report=str(output.with_suffix(".qc.pdf")),
            )
        click.echo(f"iBAQ results saved to {output}")
        return

    # Non-iBAQ methods use the generic quantification API
    click.echo(f"Using {method} quantification method")

    if method_lower == "topn":
        quant_method = get_quantification_method(method, n=topn_n)
    elif method_lower == "top3":
        quant_method = get_quantification_method("topn", n=3)
    elif method_lower == "maxlfq":
        quant_method = get_quantification_method(method, threads=threads)
    elif method_lower == "directlfq":
        quant_method = get_quantification_method(method)
    else:
        quant_method = get_quantification_method(method)

    click.echo(
        f"Quantifying {peptide_df['ProteinName'].nunique()} proteins "
        f"across {peptide_df['SampleID'].nunique()} samples ..."
    )

    result_df = quant_method.quantify(
        peptide_df,
        protein_column="ProteinName",
        peptide_column="PeptideCanonical",
        intensity_column="NormIntensity",
        sample_column="SampleID",
    )

    # Normalize if requested
    if normalize:
        intensity_cols = [c for c in result_df.columns if "Intensity" in c]
        if intensity_cols:
            int_col = intensity_cols[-1]
            result_df[f"{int_col}Norm"] = result_df.groupby("SampleID")[
                int_col
            ].transform(lambda x: x / x.sum())

    # Step 3: Save output
    output.parent.mkdir(parents=True, exist_ok=True)
    if str(output).endswith(".parquet"):
        result_df.to_parquet(str(output), index=False)
    else:
        result_df.to_csv(str(output), sep="\t", index=False)

    click.echo(
        f"Quantification complete: {result_df['ProteinName'].nunique()} proteins -> {output}"
    )
