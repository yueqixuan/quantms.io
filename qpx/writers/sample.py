"""SampleWriter — writes sample.parquet files."""

from qpx.writers.base import BaseWriter
from qpx.core.data import SampleSchema


class SampleWriter(BaseWriter):
    _schema_class = SampleSchema
