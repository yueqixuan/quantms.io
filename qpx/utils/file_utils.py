"""
File utility functions for qpx.
This module provides functions for file operations, optimized for performance.
"""

import logging
import os
from pathlib import Path
from typing import Dict, Iterator, List, Tuple, Union

import pandas as pd
import psutil
import pyarrow as pa
import pyarrow.parquet as pq

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def validate_file(mztab_path: Union[str, Path]) -> bool:
    """Validate that the mzTab file exists and is not empty."""
    if not Path(mztab_path).exists():
        raise FileNotFoundError(f"mzTab file not found: {mztab_path}")
    if Path(mztab_path).stat().st_size == 0:
        raise ValueError("mzTab file is empty")
    return Path(mztab_path).stat().st_size > 0


def validate_extension(file_path: Union[str, Path], extension: str) -> bool:
    """Validate that the file has the correct extension."""
    if not Path(file_path).exists():
        raise FileNotFoundError(f"File not found: {file_path}")
    if not Path(file_path).suffix == extension:
        raise ValueError(
            f"File {file_path} does not have the correct extension {extension}"
        )
    return True


def validate_files_exists_in_path(
    file_path: Union[str, Path], file_names: List[str]
) -> bool:
    """Validate that the files exist in the path."""
    # Ensure file_path is always a Path object
    path_obj = Path(file_path)

    if not path_obj.exists():
        return False

    for file_name in file_names:
        # Now we have Path / str which works correctly on all platforms
        if not (path_obj / file_name).exists():
            return False
    return True


def extract_protein_list(file: str) -> List[str]:
    """
    Extract a list of proteins from a file.

    Parameters:
    -----------
    file : str
        Path to the file containing protein accessions (one per line)

    Returns:
    --------
    List[str]
        List of protein accessions
    """
    # Use context manager for proper file handling
    with open(file, encoding="utf-8") as f:
        protein_list = [line.strip() for line in f]
    return protein_list


def delete_files_extension(folder: str, extension: str) -> None:
    """
    Delete all files with the given extension in the given folder.

    Parameters:
    -----------
    folder : str
        Folder path
    extension : str
        File extension to match (e.g., '.txt')
    """
    # Use Path for more reliable file operations
    folder_path = Path(folder)
    for file_path in folder_path.glob(f"*{extension}"):
        file_path.unlink()
        logger.info(f"Deleted {file_path}")


def extract_len(file_path: str, header: str) -> Tuple[int, int]:
    """
    Extract the length and position of a section in a tab-delimited file.
    Optimized for performance with proper file handling.

    Parameters:
    -----------
    file_path : str
        Path to the file
    header : str
        Header to search for

    Returns:
    --------
    Tuple[int, int]
        Length of the section and position in the file
    """
    map_tag = {"PSH": "PSM", "PEH": "PEP", "PRH": "PRT"}

    # Check file size first
    if os.stat(file_path).st_size == 0:
        raise ValueError(f"File {file_path} is empty")

    # Use context manager for proper file handling
    with open(file_path, "r") as f:
        pos = 0
        # Find the header
        for line in f:
            if line.split("\t")[0] == header:
                break
            pos = f.tell()

        # Count lines in the section
        file_len = 0
        for line in f:
            if line.split("\t")[0] != map_tag[header]:
                break
            file_len += 1

    return file_len, pos


def load_de_or_ae(path: str) -> Tuple[pd.DataFrame, str]:
    """
    Load differential expression or absolute expression file.
    Optimized to handle comments at the beginning of the file.

    Parameters:
    -----------
    path : str
        Path to the DE or AE file

    Returns:
    --------
    Tuple[pd.DataFrame, str]
        DataFrame containing the data and string containing the comments
    """
    content = []

    # Use context manager for proper file handling
    with open(path, encoding="utf-8") as f:
        # Read comment lines
        for line in f:
            if not line.startswith("#"):
                break
            content.append(line)

        # Reset file position to beginning
        f.seek(0)

        # Skip comment lines when reading with pandas
        df = pd.read_csv(f, sep="\t", comment="#")

    return df, "".join(content)


def read_large_parquet(
    parquet_path: str, batch_size: int = 500000
) -> Iterator[pd.DataFrame]:
    """
    Read a large parquet file in batches to reduce memory usage.

    Parameters:
    -----------
    parquet_path : str
        Path to the parquet file
    batch_size : int, optional
        Number of rows to read in each batch, defaults to 500000

    Yields:
    -------
    pd.DataFrame
        Batch of data from the parquet file
    """
    parquet_file = pq.ParquetFile(parquet_path)

    # Use memory-efficient batch processing
    for batch in parquet_file.iter_batches(batch_size=batch_size):
        yield batch.to_pandas()


def calculate_buffer_size(file_path: str) -> int:
    # Get the total available system memory
    total_memory = psutil.virtual_memory().available

    # Get the size of the file
    file_size = os.path.getsize(file_path)

    # Set the buffer size based on a fraction of available memory or a maximum size
    max_buffer_size = 1 * 1024 * 1024 * 1024  # 1GB
    fraction_of_memory = 0.4  # Adjust as needed

    return min(int(total_memory * fraction_of_memory), max_buffer_size, file_size)


def save_slice_file(
    parquet_table: pa.Table,
    pqwriters: Dict,
    output_folder: str,
    partitions: tuple,
    filename: str,
) -> Dict:
    """
    Save a parquet table to a file with partitioning.

    Parameters:
    -----------
    parquet_table : pyarrow.Table
        Table to save
    pqwriters : Dict
        Dictionary of parquet writers
    output_folder : str
        Base output folder
    partitions : tuple
        Partition values
    filename : str
        Output filename

    Returns:
    --------
    Dict
        Updated dictionary of parquet writers
    """
    # Create folder path using Path for better cross-platform compatibility
    folder_path = Path(output_folder)
    for part in partitions:
        folder_path = folder_path / str(part)

    # Create directory if it doesn't exist
    folder_path.mkdir(parents=True, exist_ok=True)

    # Create file path
    save_path = folder_path / filename

    # Create writer if it doesn't exist
    if partitions not in pqwriters:
        pqwriters[partitions] = pq.ParquetWriter(str(save_path), parquet_table.schema)

    # Write table
    pqwriters[partitions].write_table(parquet_table)

    return pqwriters


def save_file(
    parquet_table: pa.Table,
    pqwriter: pq.ParquetWriter | None,
    output_folder: str,
    filename: str,
) -> pq.ParquetWriter:
    """
    Save a parquet table to a file.

    Parameters:
    -----------
    parquet_table : pyarrow.Table
        Table to save
    pqwriter : ParquetWriter or None
        Existing parquet writer or None
    output_folder : str
        Output folder
    filename : str
        Output filename

    Returns:
    --------
    ParquetWriter
        Parquet writer object
    """
    # Create directory if it doesn't exist
    folder_path = Path(output_folder)
    folder_path.mkdir(parents=True, exist_ok=True)

    # Create file path
    save_path = folder_path / filename

    # Create writer if it doesn't exist
    if not pqwriter:
        pqwriter = pq.ParquetWriter(str(save_path), parquet_table.schema)

    # Write table
    pqwriter.write_table(parquet_table)

    return pqwriter


def close_file(
    pqwriters: Dict | None = None, pqwriter: pq.ParquetWriter | None = None
) -> None:
    """
    Close parquet writers.

    Parameters:
    -----------
    pqwriters : Dict, optional
        Dictionary of parquet writers
    pqwriter : ParquetWriter, optional
        Single parquet writer
    """
    if pqwriter:
        pqwriter.close()
    elif pqwriters:
        for writer in pqwriters.values():
            writer.close()


def find_ae_files(directory: str) -> List[Path]:
    """
    Find absolute expression files in a directory.

    Parameters:
    -----------
    directory : str
        Directory to search

    Returns:
    --------
    List[Path]
        List of absolute expression file paths
    """
    path = Path(directory)
    ae_files = list(path.rglob("*.absolute.tsv"))
    return ae_files


class ParquetBatchWriter:
    """Efficient batch writer for PSM data using PyArrow."""

    def __init__(
        self,
        output_path: str,
        schema: pa.Schema,
        batch_size: int = 1000000,
        compression: str = "gzip",
        file_metadata: dict = None,
    ):
        """Initialize batch writer.

        Args:
            output_path: Path to output Parquet file
            schema: PyArrow schema for the data
            batch_size: Number of records to accumulate before writing
            compression: Compression method to use
            file_metadata: Optional file-level metadata to store in parquet metadata
        """
        self.output_path = output_path
        self.schema = schema
        self.batch_size = batch_size
        self.batch_data = []
        self.parquet_writer = None
        self.compression = compression
        self.file_metadata = file_metadata or {}
        self.logger = logging.getLogger(__name__)

    def write_batch(self, records: List[dict]) -> None:
        """Write a batch of records.

        Args:
            records: List of record dictionaries
        """
        # Add records to current batch
        self.batch_data.extend(records)

        # Write batch if we've accumulated enough records
        if len(self.batch_data) >= self.batch_size:
            self._write_batch()

    def _write_batch(self) -> None:
        """Write accumulated batch data efficiently using PyArrow's streaming writer."""
        try:
            if self.batch_data:
                # Initialize writer lazily if not already created
                if self.parquet_writer is None:
                    # Prepare schema metadata by combining schema metadata with file metadata
                    schema_metadata = self.schema.metadata or {}

                    # Convert file_metadata to string format for parquet metadata
                    if self.file_metadata:
                        for key, value in self.file_metadata.items():
                            schema_metadata[f"file_metadata.{key}"] = str(value)

                    # Update schema with metadata
                    schema_with_metadata = self.schema.with_metadata(schema_metadata)

                    self.parquet_writer = pq.ParquetWriter(
                        where=self.output_path,
                        schema=schema_with_metadata,
                        compression=self.compression,
                    )

                # Create a RecordBatch directly from the current batch
                batch = pa.RecordBatch.from_pylist(self.batch_data, schema=self.schema)

                # Write the batch directly
                self.parquet_writer.write_batch(batch)
                self.batch_data = []

        except Exception as e:
            self.logger.error(
                f"Error during batch writing: {e}, file path: {self.output_path}"
            )
            raise

    def close(self) -> None:
        """Close the writer and write any remaining data."""
        try:
            # Write any remaining data
            if self.batch_data:
                self._write_batch()

            # Close the writer
            if self.parquet_writer:
                self.parquet_writer.close()

        except Exception as e:
            self.logger.error(f"Error closing writer: {e}")
            raise
