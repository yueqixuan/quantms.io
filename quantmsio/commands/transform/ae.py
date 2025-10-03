import click

from quantmsio.core.ae import AbsoluteExpressionHander
from quantmsio.utils.file_utils import extract_protein_list


@click.command(
    "ae",
    short_help="Convert IBAQ absolute file into quantms.io format",
)
@click.option(
    "--ibaq-file",
    help="IBAQ file path",
    required=True,
    type=click.Path(exists=True, dir_okay=False),
)
@click.option(
    "--sdrf-file",
    help="SDRF file path",
    required=True,
    type=click.Path(exists=True, dir_okay=False),
)
@click.option(
    "--protein-file",
    help="Protein file that meets specific requirements",
    required=False,
    type=click.Path(exists=True, dir_okay=False),
)
@click.option(
    "--project-file",
    help="quantms.io project file",
    required=False,
    type=click.Path(exists=True, dir_okay=False),
)
@click.option(
    "--output-folder",
    help="Output directory for generated files",
    required=True,
    type=click.Path(file_okay=False),
)
@click.option(
    "--output-prefix",
    help="Prefix for output files",
    required=False,
)
@click.option(
    "--delete-existing",
    help="Delete existing files in output folder",
    is_flag=True,
)
def convert_ibaq_absolute_cmd(
    ibaq_file: str,
    sdrf_file: str,
    project_file: str,
    protein_file: str,
    output_folder: str,
    output_prefix: str,
    delete_existing: bool = True,
):
    """Convert IBAQ absolute file into quantms.io format.

    Args:
        ibaq_file: IBAQ file path
        sdrf_file: SDRF file path
        project_file: Optional quantms.io project file
        protein_file: Optional protein file with requirements
        output_folder: Output directory for generated files
        output_prefix: Optional prefix for output files
        delete_existing: Whether to delete existing files
    """
    protein_list = extract_protein_list(protein_file) if protein_file else None
    protein_str = "|".join(protein_list) if protein_list else None
    ae_handler = AbsoluteExpressionHander()
    if project_file:
        ae_handler.load_project_file(project_file)
    ae_handler.load_ibaq_file(ibaq_file, protein_str)
    ae_handler.load_sdrf_file(sdrf_file)
    ae_handler.convert_ibaq_to_quantms(
        output_folder=output_folder,
        output_file_prefix=output_prefix,
        delete_existing=delete_existing,
    )
