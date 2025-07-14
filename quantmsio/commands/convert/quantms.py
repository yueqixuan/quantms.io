"""
Core converters for quantms.io formats (PSM, Feature, mzTab Protein Groups).
"""

import logging
from pathlib import Path
import tempfile
from typing import Optional

import click
import pyarrow.parquet as pq

from quantmsio.core.idxml.idxml import IdXmlPsm
from quantmsio.core.quantms.mztab import MzTabIndexer
from quantmsio.core.project import create_uuid_filename
from quantmsio.core.quantms.feature import Feature
from quantmsio.core.quantms.pg import MzTabProteinGroups
from quantmsio.core.quantms.psm import Psm


@click.group()
def convert():
    """Convert various formats to quantms.io format."""
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
    "--backend",
    help="Storage backend ('duckdb' or 'parquet'). Defaults to 'duckdb'.",
    default="duckdb",
    show_default=True,
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
    mztab_path: Path = None,
    database_path: Path = None,
    backend: str = "duckdb",
    output_folder: Path = None,
    output_prefix: str = "feature",
    sdrf_file: Path = None,
    msstats_file: Path = None,
    verbose: bool = False,
):
    """Convert feature data from mzTab to quantms.io format."""
    logger = logging.getLogger("quantmsio.commands.convert.feature")
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
            indexer = MzTabIndexer.open(str(database_path))
        elif database_path and mztab_path:
            logger.info(
                f"Creating new MzTabIndexer at {database_path} from {mztab_path}"
            )
            indexer = MzTabIndexer.create(
                mztab_path=str(mztab_path),
                backend=backend,
                database_path=str(database_path),
            )
        elif mztab_path:
            # No database_path provided, create a temporary one
            import tempfile

            temp_db_path = tempfile.mktemp(suffix=".duckdb")
            logger.info(
                f"Creating temporary MzTabIndexer at {temp_db_path} from {mztab_path}"
            )
            indexer = MzTabIndexer.create(
                mztab_path=str(mztab_path),
                backend=backend,
                database_path=temp_db_path,
            )
        else:
            raise click.ClickException(
                "You must provide either --database-path (existing) or --mztab-path (to create a new indexer)."
            )

        # Use new composition pattern
        feature = Feature(
            mztab_indexer=indexer,
            sdrf_path=str(sdrf_file) if sdrf_file else None,
            msstats_in_path=str(msstats_file) if msstats_file else None,
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
    "--backend",
    help="Storage backend ('duckdb' or 'parquet'). Defaults to 'duckdb'.",
    default="duckdb",
    show_default=True,
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
    "--sdrf-file",
    help="SDRF file path",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "--chunksize",
    help="Chunk size for writing to parquet file",
    default=10000,
    type=int,
)
@click.option("--verbose", help="Enable verbose logging", is_flag=True)
def convert_quantms_psm_cmd(
    mztab_path: Path = None,
    database_path: Path = None,
    backend: str = "duckdb",
    output_folder: Path = None,
    output_prefix: str = "psm",
    sdrf_file: Path = None,
    chunksize: int = 10000,
    verbose: bool = False,
):
    """Convert PSM data from mzTab to quantms.io format."""
    logger = logging.getLogger("quantmsio.commands.convert.psm")
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
            indexer = MzTabIndexer.open(str(database_path))
        elif database_path and mztab_path:
            logger.info(
                f"Creating new MzTabIndexer at {database_path} from {mztab_path}"
            )
            indexer = MzTabIndexer.create(
                mztab_path=str(mztab_path),
                backend=backend,
                database_path=str(database_path),
                sdrf_path=str(sdrf_file) if sdrf_file else None,
            )
        elif mztab_path:
            # Create database in the same folder as mzTab file
            mztab_parent = Path(mztab_path).parent
            db_path = mztab_parent / f"{Path(mztab_path).stem}.duckdb"
            logger.info(f"Creating MzTabIndexer at {db_path} from {mztab_path}")
            indexer = MzTabIndexer.create(
                mztab_path=str(mztab_path),
                backend=backend,
                database_path=str(db_path),
                sdrf_path=str(sdrf_file) if sdrf_file else None,
            )
        else:
            raise click.ClickException(
                "You must provide either --database-path (existing) or --mztab-path (to create a new indexer)."
            )
        psm = Psm(indexer)
        psm.convert_to_parquet(output_path=str(output_file), chunksize=chunksize)

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
def convert_idxml_psm_cmd(
    idxml_path: Path,
    output_folder: Path,
    output_prefix: str = "idxml-psm",
    verbose: bool = False,
):
    """Convert PSM data from idXML to quantms.io format."""
    logger = logging.getLogger("quantmsio.commands.convert.idxml")
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
        idxml_psm = IdXmlPsm(str(idxml_path))
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
    "--backend",
    help="Storage backend ('duckdb' or 'parquet'). Defaults to 'duckdb'.",
    default="duckdb",
    show_default=True,
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
    backend: str = "duckdb",
    msstats_file: Path = None,
    sdrf_file: Path = None,
    output_folder: Path = None,
    output_prefix: str = "pg",
    compute_topn: bool = True,
    compute_ibaq: bool = True,
    topn: int = 3,
    verbose: bool = False,
):
    """Convert protein groups from mzTab quantms TMT and LFQ data to quantms.io format using msstats for quantification.

    This command combines protein group definitions from mzTab with complete quantification
    data from msstats_in.csv. For LFQ data, it computes:

    - Default intensity: Sum of all peptidoform intensities
    - Additional intensities (optional):
      - TopN: Mean of top N peptide intensities (default N=3)
      - iBAQ: Intensity-based absolute quantification

    The output file will be named as {prefix}-{uuid}.pg.parquet.
    """
    logger = logging.getLogger("quantmsio.commands.convert.mztab")
    if verbose:
        logger.setLevel(logging.DEBUG)
        # Set root logger to DEBUG as well for comprehensive logging
        logging.getLogger().setLevel(logging.DEBUG)
        logger.info(
            "Converting mzTab protein groups to quantms.io format using msstats quantification"
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
            indexer = MzTabIndexer.open(str(database_path))
        elif database_path and mztab_path:
            logger.info(
                f"Creating new MzTabIndexer at {database_path} from {mztab_path}"
            )
            indexer = MzTabIndexer.create(
                mztab_path=str(mztab_path),
                backend=backend,
                database_path=str(database_path),
            )
        elif mztab_path:
            # No database_path provided, create a temporary one
            import tempfile

            temp_db_path = tempfile.mktemp(suffix=".duckdb")
            logger.info(
                f"Creating temporary MzTabIndexer at {temp_db_path} from {mztab_path}"
            )
            indexer = MzTabIndexer.create(
                mztab_path=str(mztab_path),
                backend=backend,
                database_path=temp_db_path,
            )
        else:
            raise click.ClickException(
                "You must provide either --database-path (existing) or --mztab-path (to create a new indexer)."
            )

        # Perform conversion using OPTIMIZED msstats quantification (4-6x faster!)
        # Use context manager to ensure cleanup of any temporary files
        with MzTabProteinGroups(mztab_path) as mztab_pg:
            result_df = mztab_pg.quantify_from_msstats_optimized(
                str(msstats_file),
                str(sdrf_file),
                compute_topn=compute_topn,
                topn=topn,
                compute_ibaq=compute_ibaq,
            )

            # Convert to parquet and write
            table = mztab_pg._convert_to_parquet_format(result_df)
            pq.write_table(table, str(output_file))
            logger.info("Successfully wrote protein groups to parquet file")

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
    logger = logging.getLogger("quantmsio.commands.convert.quantms")
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
