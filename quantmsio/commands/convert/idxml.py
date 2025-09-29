import click
from pathlib import Path
from typing import Optional, Union

from quantmsio.core.project import create_uuid_filename
from quantmsio.core.idxml import IdXML


@click.command(
    "convert-idxml",
    short_help="Convert IdXML to PSM parquet file in quantms io",
)
@click.option(
    "--idxml-file",
    help="the IdXML file containing identifications",
    required=True,
)
@click.option(
    "--output-folder",
    help="Folder where the parquet file will be generated",
    required=True,
)
@click.option(
    "--mzml-file",
    help="Optional mzML to attach spectra by scan",
    required=False,
)
@click.option(
    "--output-prefix-file",
    help="Prefix of the parquet file needed to generate the file name",
    required=False,
)
@click.option(
    "--spectral-data",
    help="Spectral data fields (optional)",
    is_flag=True,
)
def convert_idxml_file(
    idxml_file: Union[Path, str],
    output_folder: str,
    mzml_file: Optional[Union[Path, str]],
    output_prefix_file: Optional[str],
    spectral_data: bool = False,
) -> None:

    if idxml_file is None or output_folder is None:
        raise click.UsageError("Please provide all the required parameters")

    if not output_prefix_file:
        output_prefix_file = "psm"

    parser = IdXML(
        idxml_path=idxml_file, mzml_path=mzml_file, spectral_data=spectral_data
    )
    output_path = (
        f"{output_folder}/{create_uuid_filename(output_prefix_file, '.psm.parquet')}"
    )
    parser.to_parquet(output_path)
