import sys
from pathlib import Path
from typing import Optional, TextIO

import click

from qpx.operate.statistics import IbaqStatistics, ParquetStatistics

CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])


@click.group(name="stats", context_settings=CONTEXT_SETTINGS)
def statistics_cmd() -> None:
    """Statistical analysis commands for QPX data"""
    pass


@statistics_cmd.command(
    "project-ae",
    short_help="Generate statistics for a project's AE data",
)
@click.option(
    "--absolute-path",
    help="Absolute expression file path",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "--parquet-path",
    help="PSM parquet file path",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "--save-path",
    help="Output statistics file path (e.g. stats.txt)",
    type=click.Path(dir_okay=False, path_type=Path),
)
def project_ae_statistics_cmd(
    absolute_path: Path,
    parquet_path: Path,
    save_path: Optional[Path],
) -> None:
    """Generate comprehensive statistics for a project's absolute expression data.
    
    Analyzes both absolute expression (AE) and PSM data to generate a complete statistical summary. 
    This command is useful for quality control and generating summary statistics for publications or reports.
    
    Example:
        qpxc stats analyze project-ae \\
            --absolute-path ./output/ae.parquet \\
            --parquet-path ./output/psm.parquet \\
            --save-path ./reports/project_statistics.txt
    """
    feature_statistics = ParquetStatistics(str(parquet_path))
    absolute_stats = IbaqStatistics(ibaq_path=str(absolute_path))

    def write_stats(file: TextIO, stats: ParquetStatistics) -> None:
        file.write(f"Number of proteins: {stats.get_number_of_proteins()}\n")
        file.write(f"Number of peptides: {stats.get_number_of_peptides()}\n")
        file.write(f"Number of samples: {stats.get_number_of_samples()}\n")
        file.write(f"Number of peptidoforms: {stats.get_number_of_peptidoforms()}\n")
        file.write(f"Number of msruns: {stats.get_number_msruns()}\n")

    def write_absolute_stats(file: TextIO, stats: IbaqStatistics) -> None:
        file.write(f"iBAQ Number of proteins: {stats.get_number_of_proteins()}\n")
        file.write(f"iBAQ Number of samples: {stats.get_number_of_samples()}\n")

    if save_path:
        with open(save_path, "w") as f:
            write_stats(f, feature_statistics)
            write_absolute_stats(f, absolute_stats)
    else:
        write_stats(sys.stdout, feature_statistics)
        write_absolute_stats(sys.stdout, absolute_stats)


@statistics_cmd.command(
    "psm",
    short_help="Generate statistics for a PSM parquet file",
)
@click.option(
    "--parquet-path",
    help="PSM parquet file path",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "--save-path",
    help="Output statistics file path (e.g. stats.txt)",
    type=click.Path(dir_okay=False, path_type=Path),
)
def psm_statistics_cmd(
    parquet_path: Path,
    save_path: Optional[Path],
) -> None:
    """Generate statistics for PSM (Peptide-Spectrum Match) data.
    
    Analyzes PSM data to generate detailed statistics about identifications, including protein, 
    peptide, and PSM counts.
    
    Example:
        qpxc stats analyze psm \\
            --parquet-path ./output/psm.parquet \\
            --save-path ./reports/psm_statistics.txt
    """

    def write_stats(file: TextIO, stats: ParquetStatistics) -> None:
        file.write(f"Number of proteins: {stats.get_number_of_proteins()}\n")
        file.write(f"Number of peptides: {stats.get_number_of_peptides()}\n")
        file.write(f"Number of peptidoforms: {stats.get_number_of_peptidoforms()}\n")
        file.write(f"Number of PSMs: {stats.get_number_of_psms()}\n")
        file.write(f"Number of msruns: {stats.get_number_msruns()}\n")

    feature_statistics = ParquetStatistics(str(parquet_path))
    if save_path:
        with open(save_path, "w") as f:
            write_stats(f, feature_statistics)
    else:
        write_stats(sys.stdout, feature_statistics)
