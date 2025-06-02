import click
from quantmsio.core.project import create_uuid_filename
from quantmsio.operate.tools import write_ibaq_feature


@click.command(
    "ibaq",
    short_help="Create an iBAQ view from quantms.io feature data",
)
@click.option(
    "--feature_file",
    help="quantms.io feature file to process",
    required=True,
)
@click.option(
    "--sdrf_file",
    help="the SDRF file needed to extract metadata",
    required=True,
)
@click.option(
    "--output_folder",
    help="Folder where the iBAQ view will be generated",
    required=True,
)
@click.option(
    "--output_prefix_file",
    help="Prefix of the output iBAQ file",
    required=False,
)
def convert_ibaq_file(
    feature_file: str,
    sdrf_file: str,
    output_folder: str,
    output_prefix_file: str,
):
    """
    Create an iBAQ (intensity-based absolute quantification) view from quantms.io feature data.
    This command transforms existing quantms.io feature data into an iBAQ representation
    that can be consumed by ibaqpy.

    :param feature_file: quantms.io feature file to process
    :param sdrf_file: the SDRF file needed to extract metadata
    :param output_folder: Folder where the iBAQ view will be generated
    :param output_prefix_file: Prefix of the output iBAQ file
    """

    if feature_file is None or sdrf_file is None or output_folder is None:
        raise click.UsageError("Please provide all the required parameters")

    if not output_prefix_file:
        output_prefix_file = ""

    output_path = (
        output_folder + "/" + create_uuid_filename(output_prefix_file, ".ibaq.parquet")
    )
    write_ibaq_feature(sdrf_file, feature_file, output_path)
