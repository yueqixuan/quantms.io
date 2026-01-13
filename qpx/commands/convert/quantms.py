"""
Core converters for QPX formats (PSM, Feature, mzTab Protein Groups).
"""

import logging
import os
from pathlib import Path
import tempfile
from typing import Optional

import click
import pyarrow.parquet as pq

from qpx.core.idxml_utils.idxml import IdXmlPsm
from qpx.core.quantms.mztab import MzTabIndexer
from qpx.core.project import create_uuid_filename
from qpx.core.quantms.feature import Feature
from qpx.core.quantms.pg import MzTabProteinGroups
from qpx.core.quantms.psm import Psm


@click.group()
def convert():
    """Convert various formats to QPX format."""
    pass


# Feature conversion
@convert.command("quantms-feature")
@click.option(
    "--mztab-path",
    help="Input mzTab file path (required if creating a new indexer)",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "--database-path",
    help="DuckDB database file path (if exists, will be opened; if not, will be created if mztab-path is provided)",
    type=click.Path(dir_okay=False, path_type=Path),
)
@click.option(
    "--output-folder",
    help="Output directory for generated files",
    required=True,
    type=click.Path(file_okay=False, path_type=Path),
)
@click.option(
    "--output-prefix",
    help="Prefix for output files (final name will be {prefix}-{uuid}.feature.parquet)",
    default="feature",
)
@click.option(
    "--sdrf-file",
    help="SDRF file path",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "--msstats-file",
    help="MSstats input file path",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option("--verbose", help="Enable verbose logging", is_flag=True)
def convert_quantms_feature_cmd(
    mztab_path: Path,
    database_path: Path,
    output_folder: Path,
    output_prefix: str,
    sdrf_file: Path,
    msstats_file: Path,
    verbose: bool = False,
):
    """Convert feature data from mzTab to QPX format.
    
    Converts feature-level quantification data from mzTab format to the QPX standardized 
    format, including MSstats quantification data.
    
    Example:
        qpxc convert quantms-feature \\
            --mztab-path /path/to/data.mzTab \\
            --sdrf-file /path/to/metadata.sdrf.tsv \\
            --msstats-file /path/to/msstats_in.csv \\
            --output-folder ./output
    """
    logger = logging.getLogger("qpx.commands.convert.feature")
    if verbose:
        logger.setLevel(logging.DEBUG)
        # Set root logger to DEBUG as well for comprehensive logging
        logging.getLogger().setLevel(logging.DEBUG)

    try:
        # Create output directory if it doesn't exist
        output_folder.mkdir(parents=True, exist_ok=True)

        # Generate output filename with UUID
        filename = create_uuid_filename(output_prefix, ".feature.parquet")
        output_file = output_folder / filename

        # Determine how to open or create the indexer
        indexer = None
        if database_path and Path(database_path).exists():
            logger.info(f"Opening existing MzTabIndexer at {database_path}")
            indexer = MzTabIndexer.open(
                database_path=str(database_path),
                sdrf_path=sdrf_file,
            )
        elif database_path and mztab_path:
            logger.info(
                f"Creating new MzTabIndexer at {database_path} from {mztab_path}"
            )
            indexer = MzTabIndexer.create(
                mztab_path=str(mztab_path),
                msstats_path=msstats_file,
                sdrf_path=sdrf_file,
                database_path=str(database_path),
            )
        elif mztab_path:
            # No database_path provided, create a temporary one
            temp_db_file = tempfile.NamedTemporaryFile(suffix=".duckdb", delete=False)
            temp_db_path = temp_db_file.name
            temp_db_file.close()
            if os.path.exists(temp_db_path):
                os.unlink(temp_db_path)
            logger.info(
                f"Creating temporary MzTabIndexer at {temp_db_path} from {mztab_path}"
            )
            indexer = MzTabIndexer.create(
                mztab_path=str(mztab_path),
                msstats_path=msstats_file,
                sdrf_path=sdrf_file,
                database_path=temp_db_path,
            )
        else:
            raise click.ClickException(
                "You must provide either --database-path (existing) or --mztab-path (to create a new indexer)."
            )

        # Use new composition pattern
        feature = Feature(
            mztab_indexer=indexer,
        )

        feature.write_feature_to_file(output_path=str(output_file))

    except Exception as e:
        logger.exception(f"Error in feature conversion: {str(e)}")
        raise click.ClickException(f"Error: {str(e)}\nCheck the logs for more details.")


# PSM conversion
@convert.command("quantms-psm")
@click.option(
    "--mztab-path",
    help="Input mzTab file path (required if creating a new indexer)",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "--database-path",
    help="DuckDB database file path (if exists, will be opened; if not, will be created if mztab-path is provided)",
    type=click.Path(dir_okay=False, path_type=Path),
)
@click.option(
    "--output-folder",
    help="Output directory for generated files",
    required=True,
    type=click.Path(file_okay=False, path_type=Path),
)
@click.option(
    "--output-prefix",
    help="Prefix for output files (final name will be {prefix}-{uuid}.psm.parquet)",
    default="psm",
)
@click.option(
    "--spectral-data",
    help="Spectral data fields (optional)",
    is_flag=True,
)
@click.option("--verbose", help="Enable verbose logging", is_flag=True)
def convert_quantms_psm_cmd(
    mztab_path: Path,
    database_path: Path,
    output_folder: Path,
    output_prefix: str,
    spectral_data: bool = False,
    verbose: bool = False,
):
    """Convert PSM data from mzTab to QPX format.
    
    Converts PSM data from mzTab format to the QPX standardized parquet format. 
    Can work with existing DuckDB indexes or create new ones from mzTab files.
    
    Example:
        qpxc convert quantms-psm \\
            --mztab-path /path/to/data.mzTab \\
            --output-folder ./output \\
            --verbose
    """
    logger = logging.getLogger("qpx.commands.convert.psm")
    if verbose:
        logger.setLevel(logging.DEBUG)
        # Set root logger to DEBUG as well for comprehensive logging
        logging.getLogger().setLevel(logging.DEBUG)

    try:
        # Create output directory if it doesn't exist
        output_folder.mkdir(parents=True, exist_ok=True)

        # Generate output filename with UUID
        filename = create_uuid_filename(output_prefix, ".psm.parquet")
        output_file = output_folder / filename

        # Determine how to open or create the indexer
        if database_path and Path(database_path).exists():
            logger.info(f"Opening existing MzTabIndexer at {database_path}")
            indexer = MzTabIndexer.open(
                database_path=str(database_path),
            )
        elif database_path and mztab_path:
            logger.info(
                f"Creating new MzTabIndexer at {database_path} from {mztab_path}"
            )
            indexer = MzTabIndexer.create(
                mztab_path=str(mztab_path),
                database_path=str(database_path),
            )
        elif mztab_path:
            # No database_path provided, create a temporary one
            temp_db_file = tempfile.NamedTemporaryFile(suffix=".duckdb", delete=False)
            temp_db_path = temp_db_file.name
            temp_db_file.close()
            if os.path.exists(temp_db_path):
                os.unlink(temp_db_path)
            logger.info(
                f"Creating temporary MzTabIndexer at {temp_db_path} from {mztab_path}"
            )
            indexer = MzTabIndexer.create(
                mztab_path=str(mztab_path),
                database_path=temp_db_path,
            )
        else:
            raise click.ClickException(
                "You must provide either --database-path (existing) or --mztab-path (to create a new indexer)."
            )
        psm = Psm(indexer, spectral_data)
        psm.convert_to_parquet(output_path=str(output_file), chunksize=1000000)

    except Exception as e:
        logger.exception(f"Error in PSM conversion: {str(e)}")
        raise click.ClickException(f"Error: {str(e)}\nCheck the logs for more details.")


# idXML PSM conversion
@convert.command("idxml-psm")
@click.option(
    "--idxml-path",
    help="Input idXML file path",
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
    "--output-prefix",
    help="Prefix for output files (final name will be {prefix}-{uuid}.psm.parquet)",
    default="idxml-psm",
)
@click.option("--verbose", help="Enable verbose logging", is_flag=True)
@click.option(
    "--spectral-data",
    help="Spectral data fields (optional)",
    is_flag=True,
)
def convert_idxml_psm_cmd(
    idxml_path: Path,
    output_folder: Path,
    output_prefix: str = "idxml-psm",
    verbose: bool = False,
    spectral_data: bool = False,
):
    """Convert PSM data from idXML to QPX format."""
    logger = logging.getLogger("qpx.commands.convert.idxml")
    if verbose:
        logger.setLevel(logging.DEBUG)
        # Set root logger to DEBUG as well for comprehensive logging
        logging.getLogger().setLevel(logging.DEBUG)

    try:
        # Create output directory if it doesn't exist
        output_folder.mkdir(parents=True, exist_ok=True)

        # Generate output filename with UUID
        filename = create_uuid_filename(output_prefix, ".psm.parquet")
        output_file = output_folder / filename

        # Create IdXML PSM processor and convert to parquet
        logger.info(f"Processing idXML file: {idxml_path}")
        idxml_psm = IdXmlPsm(str(idxml_path), spectral_data)
        idxml_psm.convert_to_parquet(output_path=str(output_file))

        logger.info(
            f"Successfully converted {idxml_psm.get_psm_count()} PSMs to parquet format"
        )
        logger.info(f"Output file: {output_file}")

    except Exception as e:
        logger.exception(f"Error in idXML PSM conversion: {str(e)}")
        raise click.ClickException(f"Error: {str(e)}\nCheck the logs for more details.")


# mzTab Protein Group conversion
@convert.command("quantms-pg")
@click.option(
    "--mztab-path",
    help="Input mzTab file path (required if creating a new indexer)",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "--database-path",
    help="DuckDB database file path (if exists, will be opened; if not, will be created if mztab-path is provided)",
    type=click.Path(dir_okay=False, path_type=Path),
)
@click.option(
    "--msstats-file",
    help="Input msstats_in.csv file path for quantification",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "--sdrf-file",
    help="SDRF file path",
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
    "--output-prefix",
    help="Prefix for output files (final name will be {prefix}-{uuid}.pg.parquet)",
    default="pg",
)
@click.option(
    "--compute-topn/--no-compute-topn",
    help="Whether to compute TopN intensity",
    default=True,
)
@click.option(
    "--compute-ibaq/--no-compute-ibaq",
    help="Whether to compute iBAQ intensity",
    default=True,
)
@click.option(
    "--topn",
    help="Number of peptides to use for TopN intensity",
    default=3,
    type=int,
)
@click.option("--verbose", help="Enable verbose logging", is_flag=True)
def convert_quantms_pg_cmd(
    mztab_path: Path = None,
    database_path: Path = None,
    msstats_file: Path = None,
    sdrf_file: Path = None,
    output_folder: Path = None,
    output_prefix: str = "pg",
    compute_topn: bool = True,
    compute_ibaq: bool = True,
    topn: int = 3,
    verbose: bool = False,
):
    """Convert protein groups from mzTab quantms TMT and LFQ data to QPX format using msstats for quantification.

    This command combines protein group definitions from mzTab with complete quantification
    data from msstats_in.csv. Supports both TMT and LFQ data, with optional TopN and iBAQ 
    intensity calculations.

    Example:
        qpxc convert quantms-pg \\
            --mztab-path /path/to/data.mzTab \\
            --msstats-file /path/to/msstats_in.csv \\
            --sdrf-file /path/to/metadata.sdrf.tsv \\
            --output-folder ./output \\
            --compute-topn \\
            --compute-ibaq \\
            --topn 3
    """
    logger = logging.getLogger("qpx.commands.convert.mztab")
    if verbose:
        logger.setLevel(logging.DEBUG)
        # Set root logger to DEBUG as well for comprehensive logging
        logging.getLogger().setLevel(logging.DEBUG)
        logger.info(
            "Converting mzTab protein groups to QPX format using msstats quantification"
        )

    try:
        # Create output directory if it doesn't exist
        output_folder.mkdir(parents=True, exist_ok=True)
        logger.info(f"Using output directory: {output_folder}")

        # Generate output filename with UUID
        filename = create_uuid_filename(output_prefix, ".pg.parquet")
        output_file = output_folder / filename
        logger.info(f"Will save protein groups file as: {filename}")

        # Determine how to open or create the indexer
        indexer = None
        if database_path and Path(database_path).exists():
            logger.info(f"Opening existing MzTabIndexer at {database_path}")
            indexer = MzTabIndexer.open(
                database_path=str(database_path),
                sdrf_path=sdrf_file,
            )
        elif database_path and mztab_path:
            logger.info(
                f"Creating new MzTabIndexer at {database_path} from {mztab_path}"
            )
            indexer = MzTabIndexer.create(
                mztab_path=str(mztab_path),
                msstats_path=str(msstats_file),
                sdrf_path=str(sdrf_file),
                database_path=str(database_path),
            )
        elif mztab_path:
            # No database_path provided, create a temporary one
            temp_db_file = tempfile.NamedTemporaryFile(suffix=".duckdb", delete=False)
            temp_db_path = temp_db_file.name
            temp_db_file.close()
            if os.path.exists(temp_db_path):
                os.unlink(temp_db_path)
            logger.info(
                f"Creating temporary MzTabIndexer at {temp_db_path} from {mztab_path}"
            )
            indexer = MzTabIndexer.create(
                mztab_path=str(mztab_path),
                msstats_path=str(msstats_file),
                sdrf_path=str(sdrf_file),
                database_path=temp_db_path,
            )
        else:
            raise click.ClickException(
                "You must provide either --database-path (existing) or --mztab-path (to create a new indexer)."
            )

        # Perform conversion using OPTIMIZED msstats quantification (4-6x faster!)
        # Use context manager to ensure cleanup of any temporary files
        with MzTabProteinGroups(mztab_indexer=indexer) as mztab_pg:
            result_df = mztab_pg.quantify_from_msstats_optimized(
                compute_topn=compute_topn,
                topn=topn,
                compute_ibaq=compute_ibaq,
            )

            # Convert to parquet and write
            table = mztab_pg._convert_to_parquet_format(result_df)
            pq.write_table(table, str(output_file))
            logger.info(f"[Writer] Successfully wrote protein groups to: {output_file}")

    except Exception as e:
        logger.exception(f"Error in mzTab protein group conversion: {str(e)}")
        raise click.ClickException(f"Error: {str(e)}\nCheck the logs for more details.")


@convert.command("create-duckdb")
@click.option(
    "--input-file",
    help="Input mzTab file path",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "--msstats-file",
    help="MSstats input file path",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "--backend",
    help="Storage backend ('duckdb' or 'parquet'). Defaults to 'duckdb'.",
    default="duckdb",
    type=click.Choice(["duckdb", "parquet"]),
)
@click.option(
    "--output-folder",
    help="Output folder for the generated data (a .duckdb file or a directory of Parquet files).",
    required=True,
    type=click.Path(file_okay=False, dir_okay=True, path_type=Path),
)
@click.option(
    "--output-prefix",
    help="Prefix for the output file or directory.",
    default="quantms-data",
)
@click.option(
    "--max-memory",
    help="Maximum memory to use for DuckDB (e.g., '16GB')",
    default="16GB",
    type=str,
)
@click.option(
    "--threads",
    help="Number of threads to use for DuckDB",
    default=4,
    type=int,
)
@click.option("--verbose", help="Enable verbose logging", is_flag=True)
def create_duckdb_cmd(
    input_file: Path,
    msstats_file: Optional[Path] = None,
    backend: str = "duckdb",
    output_folder: Path = None,
    output_prefix: str = "quantms-data",
    max_memory: str = "16GB",
    threads: int = 4,
    verbose: bool = False,
):
    """Create a DuckDB database or a Parquet store from mzTab and optionally MSstats files."""
    logger = logging.getLogger("qpx.commands.convert.quantms")
    if verbose:
        logger.setLevel(logging.DEBUG)

    try:
        logger.info(f"Using '{backend}' backend.")
        logger.info(f"Using mzTab file: {input_file}")
        if msstats_file:
            logger.info(f"Using MSstats file: {msstats_file}")

        # Ensure the base output folder exists
        output_folder.mkdir(parents=True, exist_ok=True)

        # Determine the database path based on backend
        if backend == "duckdb":
            database_path = output_folder / f"{output_prefix}.duckdb"
            logger.info(f"Output will be a DuckDB file at: {database_path}")
        else:  # parquet backend
            database_path = output_folder / output_prefix
            logger.info(f"Output will be a Parquet directory at: {database_path}")

        # Check if output already exists
        if database_path.exists():
            raise click.ClickException(
                f"Output path already exists: {database_path}. Please choose a different output prefix or remove the existing file/directory."
            )

        # Create database with mzTab data using factory method
        indexer = MzTabIndexer.create(
            mztab_path=str(input_file),
            msstats_path=str(msstats_file) if msstats_file else None,
            backend=backend,
            database_path=str(database_path),
            max_memory=max_memory,
            worker_threads=threads,
        )

        logger.info("Processing complete!")

        # Log some basic statistics
        try:
            metadata_count = len(indexer.get_metadata())
            protein_count = indexer.get_protein_count()
            protein_details_count = indexer.get_protein_details_count()
            psm_count = indexer.get_psm_count()

            logger.info(f"Successfully processed:")
            logger.info(f"  - Metadata entries: {metadata_count:,}")
            logger.info(f"  - Proteins: {protein_count:,}")
            logger.info(f"  - Protein details: {protein_details_count:,}")
            logger.info(f"  - PSMs: {psm_count:,}")

            if indexer._has_msstats_table():
                msstats_data = indexer.get_msstats()
                if msstats_data is not None:
                    logger.info(f"  - MSstats entries: {len(msstats_data):,}")

        except Exception as stats_error:
            logger.warning(f"Could not retrieve statistics: {stats_error}")

    except Exception as e:
        logger.exception(f"Error during processing: {str(e)}")
        raise click.ClickException(f"Error: {str(e)}\nCheck the logs for more details.")
