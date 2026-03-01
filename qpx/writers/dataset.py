"""DatasetWriter — writes dataset.parquet files."""

from qpx.writers.base import BaseWriter
from qpx.core.data import DatasetSchema


class DatasetWriter(BaseWriter):
    _schema_class = DatasetSchema
