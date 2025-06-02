"""
Data validation utilities for quantms.io.
"""

from typing import Dict, Any, List, Optional, Union
import pyarrow as pa
import pyarrow.parquet as pq
import pandas as pd
import numpy as np
from pathlib import Path

from quantmsio.core.format import PSM_FIELDS, FEATURE_FIELDS, IBAQ_FIELDS, PG_FIELDS
from quantmsio.utils.logger import get_logger

logger = get_logger(__name__)


class ValidationError(Exception):
    """Exception raised for data validation errors."""

    pass


def validate_schema(
    data: Union[pd.DataFrame, pa.Table], schema: pa.Schema
) -> List[str]:
    """
    Validate that data conforms to the expected schema.

    Args:
        data: Data to validate
        schema: Expected schema

    Returns:
        List of validation errors (empty if validation passed)
    """
    errors = []

    # Convert pandas DataFrame to Arrow Table if needed
    if isinstance(data, pd.DataFrame):
        try:
            data = pa.Table.from_pandas(data)
        except Exception as e:
            errors.append(f"Failed to convert DataFrame to Arrow Table: {str(e)}")
            return errors

    # Check required fields
    for field in schema:
        if field.name not in data.column_names:
            errors.append(f"Missing required field: {field.name}")
            continue

        # Check data type
        try:
            data_type = data.schema.field(field.name).type
            if not data_type.equals(field.type):
                errors.append(
                    f"Invalid type for field {field.name}: "
                    f"expected {field.type}, got {data_type}"
                )
        except Exception as e:
            errors.append(f"Error checking type for field {field.name}: {str(e)}")

    return errors


def validate_feature_data(data: Union[pd.DataFrame, pa.Table]) -> List[str]:
    """
    Validate feature data against the feature schema.

    Args:
        data: Feature data to validate

    Returns:
        List of validation errors (empty if validation passed)
    """
    schema = pa.schema(FEATURE_FIELDS)
    return validate_schema(data, schema)


def validate_psm_data(data: Union[pd.DataFrame, pa.Table]) -> List[str]:
    """
    Validate PSM data against the PSM schema.

    Args:
        data: PSM data to validate

    Returns:
        List of validation errors (empty if validation passed)
    """
    schema = pa.schema(PSM_FIELDS)
    return validate_schema(data, schema)


def validate_ibaq_data(data: Union[pd.DataFrame, pa.Table]) -> List[str]:
    """
    Validate IBAQ data against the IBAQ schema.

    Args:
        data: IBAQ data to validate

    Returns:
        List of validation errors (empty if validation passed)
    """
    schema = pa.schema(IBAQ_FIELDS)
    return validate_schema(data, schema)


def validate_pg_data(data: Union[pd.DataFrame, pa.Table]) -> List[str]:
    """
    Validate protein group data against the PG schema.

    Args:
        data: Protein group data to validate

    Returns:
        List of validation errors (empty if validation passed)
    """
    schema = pa.schema(PG_FIELDS)
    return validate_schema(data, schema)


def validate_numeric_range(
    data: Union[pd.DataFrame, pa.Table],
    column: str,
    min_value: Optional[float] = None,
    max_value: Optional[float] = None,
) -> List[str]:
    """
    Validate that numeric values in a column are within the specified range.

    Args:
        data: Data containing the column to validate
        column: Name of the column to validate
        min_value: Minimum allowed value (inclusive)
        max_value: Maximum allowed value (inclusive)

    Returns:
        List of validation errors (empty if validation passed)
    """
    errors = []

    # Convert Arrow Table to pandas DataFrame if needed
    if isinstance(data, pa.Table):
        data = data.to_pandas()

    if column not in data.columns:
        return [f"Column not found: {column}"]

    # Check for non-numeric values
    non_numeric = data[column].apply(
        lambda x: not (isinstance(x, (int, float)) or (isinstance(x, np.number)))
    )
    if non_numeric.any():
        errors.append(
            f"Non-numeric values found in column {column} "
            f"at indices: {list(data.index[non_numeric])}"
        )

    # Check minimum value
    if min_value is not None:
        below_min = data[column] < min_value
        if below_min.any():
            errors.append(
                f"Values below minimum {min_value} found in column {column} "
                f"at indices: {list(data.index[below_min])}"
            )

    # Check maximum value
    if max_value is not None:
        above_max = data[column] > max_value
        if above_max.any():
            errors.append(
                f"Values above maximum {max_value} found in column {column} "
                f"at indices: {list(data.index[above_max])}"
            )

    return errors


def validate_unique_values(
    data: Union[pd.DataFrame, pa.Table], column: str
) -> List[str]:
    """
    Validate that values in a column are unique.

    Args:
        data: Data containing the column to validate
        column: Name of the column to validate

    Returns:
        List of validation errors (empty if validation passed)
    """
    errors = []

    # Convert Arrow Table to pandas DataFrame if needed
    if isinstance(data, pa.Table):
        data = data.to_pandas()

    if column not in data.columns:
        return [f"Column not found: {column}"]

    # Check for duplicate values
    duplicates = data[column].duplicated()
    if duplicates.any():
        errors.append(
            f"Duplicate values found in column {column} "
            f"at indices: {list(data.index[duplicates])}"
        )

    return errors


def validate_file_format(file_path: Path, expected_format: str) -> List[str]:
    """
    Validate that a file has the expected format.

    Args:
        file_path: Path to the file to validate
        expected_format: Expected file format (e.g., 'parquet', 'csv')

    Returns:
        List of validation errors (empty if validation passed)
    """
    errors = []

    if not file_path.exists():
        return [f"File not found: {file_path}"]

    # Check file extension
    if file_path.suffix.lower() != f".{expected_format.lower()}":
        errors.append(
            f"Invalid file format: expected .{expected_format}, "
            f"got {file_path.suffix}"
        )
        return errors

    # Try to read the file based on format
    try:
        if expected_format.lower() == "parquet":
            pq.read_metadata(file_path)
        elif expected_format.lower() in ["csv", "tsv"]:
            pd.read_csv(file_path, nrows=1)
    except Exception as e:
        errors.append(f"Failed to read {expected_format} file: {str(e)}")

    return errors
