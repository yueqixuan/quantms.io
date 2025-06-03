from pathlib import Path
import os
import subprocess
from typing import Optional

import click

from quantmsio.core.quantms import QuantmsPSM, QuantmsFeature
from quantmsio.core.project import create_uuid_filename


def is_diann(directory: str) -> bool:
    dirs = [
        diann_dir
        for diann_dir in os.listdir(directory)
        if os.path.isdir(os.path.join(directory, diann_dir))
    ]
    if "diannsummary" in dirs:
        return True
    else:
        return False


def check_dir(folder_path: str) -> None:
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)


def find_file(directory: str, kind: str) -> list[Path]:
    path = Path(directory)
    ae_files = list(path.rglob(f"*{kind}"))
    return ae_files


def run_task(command: list) -> bool:
    try:
        result = subprocess.run(command, capture_output=True, text=True, check=True)
        print(result.returncode)
        return True
    except subprocess.CalledProcessError as e:
        print(e.returncode)
        print(e.stderr)
        return False


def quantmsio_workflow(directory: str, output_root_folder: str) -> None:
    dirs = [
        os.path.join(directory, diann_dir)
        for diann_dir in os.listdir(directory)
        if os.path.isdir(os.path.join(directory, diann_dir))
    ]
    print(dirs)
    for diann_dir in dirs:
        if is_diann(diann_dir):
            report_file = find_file(diann_dir, "diann_report.tsv")[0]
            mzmlstatistics = os.path.join(diann_dir, "mzmlstatistics")
            sdrf_file = find_file(diann_dir, ".sdrf.tsv")[0]
            filename = os.path.basename(diann_dir)
            output_folder = os.path.join(directory, output_root_folder, filename)
            check_dir(output_folder)
            command = [
                "quantmsioc",
                "convert-diann",
                "--report_path",
                report_file,
                "--qvalue_threshold",
                "0.05",
                "--mzml_info_folder",
                mzmlstatistics,
                "--sdrf_path",
                sdrf_file,
                "--output_folder",
                output_folder,
                "--duckdb_max_memory",
                "64GB",
                "--output_prefix_file",
                filename,
            ]
            run_task(command)
        else:
            mztab_file = find_file(diann_dir, ".mzTab")
            if len(mztab_file) > 0:
                mztab_file = mztab_file[0]
                msstatsin_file = find_file(diann_dir, "msstats_in.csv")[0]
                sdrf_file = find_file(diann_dir, ".sdrf.tsv")[0]
                filename = os.path.basename(diann_dir)
                output_folder = os.path.join(directory, output_root_folder, filename)
                check_dir(output_folder)
                command_feature = [
                    "quantmsioc",
                    "convert-feature",
                    "--sdrf_file",
                    sdrf_file,
                    "--msstats_file",
                    msstatsin_file,
                    "--mztab_file",
                    mztab_file,
                    "--file_num",
                    "30",
                    "--output_folder",
                    output_folder,
                    "--duckdb_max_memory",
                    "64GB",
                    "--output_prefix_file",
                    filename,
                ]
                run_task(command_feature)
                command_psm = [
                    "quantmsioc",
                    "convert-psm",
                    "--mztab_file",
                    mztab_file,
                    "--output_folder",
                    output_folder,
                    "--output_prefix_file",
                    filename,
                ]
                run_task(command_psm)
            else:
                continue


@click.command(
    "convert-quantms—project",
    short_help="Generate quantmsio files from multiple quantms reports." "format",
)
@click.option(
    "--base_folder",
    help="The directory for storing project files",
    required=True,
)
@click.option(
    "--save_folder",
    help="The directory for saving the results",
    required=False,
)
def convert_quantms_project(
    base_folder: str,
    save_folder: Optional[str],
) -> None:
    if not save_folder:
        save_folder = "quantms_data"
    quantmsio_workflow(base_folder, save_folder)


@click.command(
    "convert-quantms-psm",
    short_help="Convert quantms PSMs from psm.tsv to parquet file in quantms.io",
)
@click.option(
    "--psm-file",
    help="the psm.tsv file, this will be used to extract the peptide information",
    required=True,
)
@click.option(
    "--output-folder",
    help="Folder where the parquet file will be generated",
    required=True,
)
@click.option(
    "--output-prefix",
    help="Prefix of the parquet file needed to generate the file name",
    required=False,
)
def convert_quantms_psm(
    psm_file: str,
    output_folder: str,
    output_prefix: str,
):
    """
    :param psm_file: the psm.tsv file, this will be used to extract the peptide information
    :param output_folder: Folder where the parquet file will be generated
    :param output_prefix: Prefix of the Json file needed to generate the file name
    """

    if psm_file is None or output_folder is None:
        raise click.UsageError("Please provide all the required parameters")

    if not output_prefix:
        output_prefix = "psm"

    output_path = (
        output_folder + "/" + create_uuid_filename(output_prefix, ".psm.parquet")
    )
    quantms_psm = QuantmsPSM(psm_file)
    quantms_psm.write_psm_to_file(output_path)


@click.command(
    "convert-quantms-feature",
    short_help="Convert quantms feature from evidence.txt to parquet file in quantms.io",
)
@click.option(
    "--feature-file",
    help="the feature.tsv file, this will be used to extract the peptide information",
    required=True,
)
@click.option(
    "--output-folder",
    help="Folder where the parquet file will be generated",
    required=True,
)
@click.option(
    "--output-prefix",
    help="Prefix of the parquet file needed to generate the file name",
    required=False,
)
def convert_quantms_feature(
    feature_file: str,
    output_folder: str,
    output_prefix: str,
):
    if feature_file is None or output_folder is None:
        raise click.UsageError("Please provide all the required parameters")

    if not output_prefix:
        output_prefix = "feature"

    filename = create_uuid_filename(output_prefix, ".feature.parquet")
    output_path = output_folder + "/" + filename
    quantms_feature = QuantmsFeature(feature_file)
    quantms_feature.write_feature_to_file(output_path)
