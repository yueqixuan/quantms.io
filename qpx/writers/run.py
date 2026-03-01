"""RunWriter — writes run.parquet files."""

from qpx.writers.base import BaseWriter
from qpx.core.data import RunSchema


class RunWriter(BaseWriter):
    _schema_class = RunSchema
