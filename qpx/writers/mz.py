"""MzWriter — writes mz.parquet files."""

from qpx.writers.base import BaseWriter
from qpx.core.data import MzSchema


class MzWriter(BaseWriter):
    _schema_class = MzSchema
