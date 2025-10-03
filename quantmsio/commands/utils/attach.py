import click

from quantmsio.core.project import ProjectHandler


@click.command(
    "attach-file",
    short_help="Register the file to project.json.",
)
@click.option("--project-file", help="the project.json file", required=True)
@click.option(
    "--attach-file", help="The path of the file that will be registered", required=True
)
@click.option(
    "--category",
    type=click.Choice(
        ["sdrf-file", "feature-file", "psm-file", "differential-file", "absolute-file"],
        case_sensitive=False,
    ),
    help="The type of file that will be registered",
    required=True,
)
@click.option(
    "--is-folder",
    help="A boolean value that indicates if the file is a folder or not",
    is_flag=True,
)
@click.option(
    "--partitions",
    help="The field used for splitting files, multiple fields are separated by ,",
    required=False,
)
@click.option("--replace-existing", help="Whether to delete old files", is_flag=True)
def attach_file_to_json_cmd(
    project_file, attach_file, category, is_folder, partitions, replace_existing
):
    if partitions:
        partitions = partitions.split(",")
    register = ProjectHandler(project_json_file=project_file)
    register.register_file(
        attach_file,
        category,
        is_folder=is_folder,
        partition_fields=partitions,
        replace_existing=replace_existing,
    )
    register.save_updated_project_info(output_file_name=project_file)
