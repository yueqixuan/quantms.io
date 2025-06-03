import click
from quantmsio.core.project import create_uuid_filename
from quantmsio.operate.tools import write_ibaq_feature


@click.command(
    "map-latest-uniprot",
    short_help="Map feature file to latest uniprot version",
)
@click.option(
    "--feature-file",
    help="feature file",
    required=True,
)
@click.option(
    "--output-folder",
    help="Folder where the Json file will be generated",
    required=True,
)
@click.option(
    "--output-prefix",
    help="Prefix of the parquet file needed to generate the file name",
    required=False,
)
def map_latest_uniprot(
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
    write_ibaq_feature(feature_file, output_path)
