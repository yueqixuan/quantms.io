from pathlib import Path
import os
import subprocess
from typing import Optional

import click
from quantmsio.core.project import create_uuid_filename


def find_file(directory: str, pattern: str) -> Optional[Path]:
    """Find first file matching pattern in directory."""
    path = Path(directory)
    files = list(path.rglob(pattern))
    return files[0] if files else None


def get_project_prefix(sdrf_file: Path) -> str:
    """Extract project prefix from SDRF filename (e.g. 'PXD000865' from 'PXD000865.sdrf.tsv')."""
    filename = sdrf_file.name
    # Remove .sdrf.tsv and any variations like _openms_design.sdrf.tsv
    prefix = filename.split(".sdrf")[0].split("_openms")[0]
    return prefix


def check_dir(folder_path: str) -> None:
    """Create directory if it doesn't exist."""
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)


def run_task(command: list) -> bool:
    """Run command and return success status."""
    try:
        result = subprocess.run(command, capture_output=True, text=True, check=True)
        return True
    except subprocess.CalledProcessError as e:
        print(f"Error running command: {e.stderr}")
        return False


def quantmsio_workflow(
    base_folder: str, output_folder: str, prefix: Optional[str] = None
) -> None:
    """Convert quantms output to quantms.io format.

    Expected structure:
    base_folder/
        quant_tables/
            *.mzTab
            *msstats_in.csv
        sdrf/
            *.sdrf.tsv
        spectra/
            mzml_statistics/
    """
    # Setup paths
    quant_tables = Path(base_folder) / "quant_tables"
    sdrf_dir = Path(base_folder) / "sdrf"
    spectra_dir = Path(base_folder) / "spectra"

    # Find required files
    mztab_file = find_file(quant_tables, "*.mzTab")
    msstats_file = find_file(quant_tables, "*msstats_in.csv")
    sdrf_file = find_file(sdrf_dir, "*.sdrf.tsv")
    mzml_stats = spectra_dir / "mzml_statistics"

    if not all([mztab_file, msstats_file, sdrf_file, mzml_stats.exists()]):
        missing = []
        if not mztab_file:
            missing.append("mzTab file")
        if not msstats_file:
            missing.append("MSstats input file")
        if not sdrf_file:
            missing.append("SDRF file")
        if not mzml_stats.exists():
            missing.append("mzML statistics")
        raise click.UsageError(f"Missing required files: {', '.join(missing)}")

    # Use provided prefix or auto-detect from SDRF
    if prefix is None:
        prefix = get_project_prefix(sdrf_file)
        print(f"No prefix provided, auto-detected: {prefix}")

    # Create output directory
    check_dir(output_folder)

    # Convert features
    command_feature = [
        "quantmsioc",
        "convert",
        "feature",
        "--sdrf-file",
        str(sdrf_file),
        "--msstats-file",
        str(msstats_file),
        "--mztab-file",
        str(mztab_file),
        "--file-num",
        "30",
        "--output-folder",
        output_folder,
        "--duckdb-max-memory",
        "64GB",
        "--output-prefix",
        prefix,
    ]
    if not run_task(command_feature):
        print("Warning: Feature conversion failed")

    # Convert PSMs
    command_psm = [
        "quantmsioc",
        "convert",
        "psm",
        "--mztab-file",
        str(mztab_file),
        "--output-folder",
        output_folder,
        "--output-prefix",
        prefix,
    ]
    if not run_task(command_psm):
        print("Warning: PSM conversion failed")


@click.command(
    "quantms",
    short_help="Convert quantms project output to quantms.io format",
)
@click.option(
    "--quantms-dir",
    help="The quantms project directory containing quant_tables, sdrf, and spectra subdirectories",
    required=True,
    type=click.Path(exists=True, file_okay=False, path_type=Path),
)
@click.option(
    "--output-dir",
    help="Output directory for quantms.io files (defaults to 'quantms.io' in parent directory)",
    required=False,
    type=click.Path(file_okay=False, path_type=Path),
)
@click.option(
    "--prefix",
    help="Prefix for output files (e.g. 'PXD000865'). If not provided, will be auto-detected from SDRF filename",
    required=False,
    type=str,
)
def convert_quantms_project_cmd(
    quantms_dir: Path,
    output_dir: Optional[Path] = None,
    prefix: Optional[str] = None,
) -> None:
    """Convert a quantms project output to quantms.io format.

    The script expects a quantms output directory with:
    - quant_tables/ containing mzTab and MSstats files
    - sdrf/ containing SDRF files
    - spectra/ containing mzML statistics
    """
    # Default output to sibling quantms.io directory
    if not output_dir:
        output_dir = str(quantms_dir.parent / "quantms.io")

    quantmsio_workflow(str(quantms_dir), output_dir, prefix)


# @click.command(
#     "psm",
#     short_help="Convert quantms PSMs from psm.tsv to parquet file in quantms.io",
# )
# @click.option(
#     "--psm-file",
#     help="the psm.tsv file, this will be used to extract the peptide information",
#     required=True,
# )
# @click.option(
#     "--output-folder",
#     help="Folder where the parquet file will be generated",
#     required=True,
# )
# @click.option(
#     "--output-prefix",
#     help="Prefix of the parquet file needed to generate the file name",
#     required=False,
# )
# def convert_quantms_psm(
#     psm_file: str,
#     output_folder: str,
#     output_prefix: str,
# ):
#     """
#     :param psm_file: the psm.tsv file, this will be used to extract the peptide information
#     :param output_folder: Folder where the parquet file will be generated
#     :param output_prefix: Prefix of the Json file needed to generate the file name
#     """

#     if psm_file is None or output_folder is None:
#         raise click.UsageError("Please provide all the required parameters")

#     if not output_prefix:
#         output_prefix = "psm"

#     output_path = (
#         output_folder + "/" + create_uuid_filename(output_prefix, ".psm.parquet")
#     )
#     quantms_psm = QuantmsPSM(psm_file)
#     quantms_psm.write_psm_to_file(output_path)


# @click.command(
#     "feature",
#     short_help="Convert quantms feature from evidence.txt to parquet file in quantms.io",
# )
# @click.option(
#     "--feature-file",
#     help="the feature.tsv file, this will be used to extract the peptide information",
#     required=True,
# )
# @click.option(
#     "--output-folder",
#     help="Folder where the parquet file will be generated",
#     required=True,
# )
# @click.option(
#     "--output-prefix",
#     help="Prefix of the parquet file needed to generate the file name",
#     required=False,
# )
# def convert_quantms_feature(
#     feature_file: str,
#     output_folder: str,
#     output_prefix: str,
# ):
#     if feature_file is None or output_folder is None:
#         raise click.UsageError("Please provide all the required parameters")

#     if not output_prefix:
#         output_prefix = "feature"

#     filename = create_uuid_filename(output_prefix, ".feature.parquet")
#     output_path = output_folder + "/" + filename
#     quantms_feature = QuantmsFeature(feature_file)
#     quantms_feature.write_feature_to_file(output_path)
