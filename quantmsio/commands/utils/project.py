from pathlib import Path

import click

from quantmsio.core.common import QUANTMSIO_VERSION
from quantmsio.core.project import check_directory


@click.command(
    "create",
    short_help="Generate a project file from original PRIDE accession",
)
@click.option(
    "--project-accession",
    help="PRIDE project accession",
    required=True,
)
@click.option(
    "--sdrf-file",
    help="SDRF file path for metadata extraction",
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
    "--software-name",
    help="Software name used to generate the data",
    required=False,
)
@click.option(
    "--software-version",
    help="Software version used to generate the data",
    required=False,
)
@click.option(
    "--delete-existing",
    help="Delete existing files in the output folder",
    is_flag=True,
)
def generate_pride_project_json_cmd(
    project_accession: str,
    sdrf_file: Path,
    software_name: str,
    software_version: str,
    output_folder: Path,
    delete_existing: bool,
):
    """Generate a project file from the original PRIDE accession and SDRF file.

    The project file definition is available in the docs folder of this repository:
    https://github.com/bigbio/quantms.io/blob/main/docs/PROJECT.md

    This command will generate the project file and attach the provided SDRF file.

    Args:
        project_accession: PRIDE project accession
        sdrf_file: SDRF file path for metadata extraction
        software_name: Software name used to generate the data
        software_version: Software version used to generate the data
        output_folder: Output directory for generated files
        delete_existing: Delete existing files in the output folder
    """
    if not all([project_accession, sdrf_file, output_folder]):
        raise click.UsageError("Please provide all required parameters")

    project_handler = check_directory(str(output_folder), project_accession)

    # Populate the project handler with metadata from PRIDE Archive and SDRF file
    project_handler.populate_from_pride_archive()
    project_handler.populate_from_sdrf(str(sdrf_file))
    project_handler.add_quantms_version(quantmsio_version=QUANTMSIO_VERSION)
    project_handler.add_software_provider(
        sortware_name=software_name, sortware_version=software_version
    )
    project_path = str(output_folder / "project.json")
    project_handler.add_sdrf_file(
        sdrf_file_path=str(sdrf_file),
        output_folder=str(output_folder),
        delete_existing=delete_existing,
    )
    project_handler.save_updated_project_info(output_file_name=project_path)
