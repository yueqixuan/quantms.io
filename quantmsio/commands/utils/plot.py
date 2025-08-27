from pathlib import Path
from typing import Optional

import click

from quantmsio.operate.plots import (
    plot_distribution_of_ibaq,
    plot_intensity_box_of_samples,
    plot_intensity_distribution_of_samples,
    plot_peptide_distribution_of_protein,
    plot_peptides_of_lfq_condition,
)

CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])


@click.group(name="plot", context_settings=CONTEXT_SETTINGS)
def plot_cmd() -> None:
    """Visualization commands for quantms.io data"""
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
    """Plot peptides by condition in LFQ analysis.

    Args:
        psm_parquet_path: PSM parquet file path
        sdrf_path: SDRF file path
        save_path: Output image path
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
    """Plot iBAQ value distribution.

    Args:
        ibaq_path: iBAQ file path
        save_path: Output image path
        select_column: Selected column in iBAQ file
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
    """Plot KDE intensity distribution.

    Args:
        feature_path: Feature file path
        save_path: Output image path
        num_samples: Number of samples to plot
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
    """Plot peptide distribution.

    Args:
        feature_path: Feature file path
        save_path: Output image path
        num_samples: Number of samples to plot
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
    """Plot intensity box plot.

    Args:
        feature_path: Feature file path
        save_path: Output image path
        num_samples: Number of samples to plot
    """
    plot_intensity_box_of_samples(str(feature_path), str(save_path), num_samples)
