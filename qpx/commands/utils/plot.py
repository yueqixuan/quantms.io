from pathlib import Path
from typing import Optional

import click

from qpx.operate.plots import (
    plot_distribution_of_ibaq,
    plot_intensity_box_of_samples,
    plot_intensity_distribution_of_samples,
    plot_peptide_distribution_of_protein,
    plot_peptides_of_lfq_condition,
)

CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])


@click.group(name="plot", context_settings=CONTEXT_SETTINGS)
def plot_cmd() -> None:
    """Visualization commands for QPX data"""
    pass


@plot_cmd.command(
    "psm-peptides",
    short_help="Plot peptides by condition in LFQ",
)
@click.option(
    "--psm-parquet-path",
    help="PSM parquet file path",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "--sdrf-path",
    help="SDRF file path",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "--save-path",
    help="Output image path (e.g. plot.svg)",
    required=True,
    type=click.Path(dir_okay=False, path_type=Path),
)
def plot_peptides_cmd(
    psm_parquet_path: Path,
    sdrf_path: Path,
    save_path: Path,
) -> None:
    """Plot peptides by condition in label-free quantification (LFQ) experiments.
    
    Creates a visualization showing the distribution of identified peptides across different 
    experimental conditions. This plot helps assess data quality and completeness across samples.
    
    Example:
        qpxc visualize plot psm-peptides \\
            --psm-parquet-path ./output/psm.parquet \\
            --sdrf-path ./metadata.sdrf.tsv \\
            --save-path ./plots/peptides_by_condition.svg
    """
    plot_peptides_of_lfq_condition(
        psm_parquet_path=str(psm_parquet_path),
        sdrf_path=str(sdrf_path),
        save_path=str(save_path),
    )


@plot_cmd.command(
    "ibaq-distribution",
    short_help="Plot iBAQ distribution",
)
@click.option(
    "--ibaq-path",
    help="iBAQ file path",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "--save-path",
    help="Output image path (e.g. plot.svg)",
    required=True,
    type=click.Path(dir_okay=False, path_type=Path),
)
@click.option(
    "--select-column",
    help="Selected column in iBAQ file",
    required=False,
)
def plot_ibaq_distribution_cmd(
    ibaq_path: Path,
    save_path: Path,
    select_column: Optional[str],
) -> None:
    """Plot the distribution of iBAQ (intensity-Based Absolute Quantification) values.
    
    Creates a histogram or density plot showing the distribution of iBAQ values across proteins. 
    Useful for quality control and understanding the dynamic range of protein quantification.
    
    Example:
        qpxc visualize plot ibaq-distribution \\
            --ibaq-path ibaq_data.tsv \\
            --save-path ./plots/ibaq_distribution.svg
    """
    plot_distribution_of_ibaq(str(ibaq_path), str(save_path), select_column)


@plot_cmd.command(
    "kde-intensity",
    short_help="Plot KDE intensity distribution",
)
@click.option(
    "--feature-path",
    help="Feature file path",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "--save-path",
    help="Output image path (e.g. plot.svg)",
    required=True,
    type=click.Path(dir_okay=False, path_type=Path),
)
@click.option(
    "--num-samples",
    help="Number of samples to plot",
    default=10,
    type=int,
)
def plot_kde_intensity_distribution_cmd(
    feature_path: Path,
    save_path: Path,
    num_samples: int,
) -> None:
    """Plot Kernel Density Estimation (KDE) of intensity distributions across samples.
    
    Creates overlaid KDE plots showing intensity distributions for multiple samples, enabling 
    visual comparison of sample-to-sample variability and batch effects.
    
    Example:
        qpxc visualize plot kde-intensity \\
            --feature-path ./output/feature.parquet \\
            --save-path ./plots/intensity_kde.svg
    """
    plot_intensity_distribution_of_samples(
        str(feature_path), str(save_path), num_samples
    )


@plot_cmd.command(
    "peptide-distribution",
    short_help="Plot peptide distribution",
)
@click.option(
    "--feature-path",
    help="Feature file path",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "--save-path",
    help="Output image path (e.g. plot.svg)",
    required=True,
    type=click.Path(dir_okay=False, path_type=Path),
)
@click.option(
    "--num-samples",
    help="Number of samples to plot",
    default=20,
    type=int,
)
def plot_peptide_distribution_cmd(
    feature_path: Path,
    save_path: Path,
    num_samples: int,
) -> None:
    """Plot the distribution of peptides across proteins.
    
    Visualizes how many peptides are identified for each protein, providing insights into 
    protein coverage and identification confidence.
    
    Example:
        qpxc visualize plot peptide-distribution \\
            --feature-path ./output/feature.parquet \\
            --save-path ./plots/peptide_per_protein.svg
    """
    plot_peptide_distribution_of_protein(str(feature_path), str(save_path), num_samples)


@plot_cmd.command(
    "box-intensity",
    short_help="Plot intensity box plot",
)
@click.option(
    "--feature-path",
    help="Feature file path",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "--save-path",
    help="Output image path (e.g. plot.svg)",
    required=True,
    type=click.Path(dir_okay=False, path_type=Path),
)
@click.option(
    "--num-samples",
    help="Number of samples to plot",
    default=10,
    type=int,
)
def plot_box_intensity_distribution_cmd(
    feature_path: Path,
    save_path: Path,
    num_samples: int,
) -> None:
    """Plot box plots of intensity distributions across samples.
    
    Creates box plots showing the distribution of feature intensities for each sample, ideal 
    for identifying outliers and assessing normalization quality.
    
    Example:
        qpxc visualize plot box-intensity \\
            --feature-path ./output/feature.parquet \\
            --save-path ./plots/intensity_boxplot.svg
    """
    plot_intensity_box_of_samples(str(feature_path), str(save_path), num_samples)
