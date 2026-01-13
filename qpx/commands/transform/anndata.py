import click
import pandas as pd

from qpx.core.combiner import Combiner
from qpx.core.project import create_uuid_filename
from qpx.utils.file_utils import find_ae_files


@click.command(
    "anndata",
    short_help="Merge multiple AE files into a file in AnnData format.",
)
@click.option(
    "--directory",
    help="The directory for storing AE files",
    required=True,
    type=click.Path(exists=True, file_okay=False),
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
def merge_ae_files_cmd(
    directory: str,
    output_folder: str,
    output_prefix: str,
):
    """Merge multiple AE files into AnnData format.
    
    Combines multiple absolute expression (AE) files from a directory into a single AnnData object 
    (H5AD format). This command is useful for integrating data from multiple experiments for 
    downstream analysis with scanpy or other Python-based tools.
    
    Example:
        qpxc transform anndata \\
            --directory ./ae_files \\
            --output-folder ./output
    """
    ae_files = find_ae_files(directory)
    output_path = output_folder + "/" + create_uuid_filename(output_prefix, ".h5ad")
    ae_combiner = Combiner()
    if len(ae_files) == 0:
        raise click.UsageError("No AE files were found.")
    else:
        for file in ae_files:
            df = pd.read_csv(file, comment="#", sep="\t")
            adata = ae_combiner.transform_to_adata(df)
            ae_combiner.combine_adata(adata)
        ae_combiner.save_adata(output_path)
