"""PsmWriter — writes psm.parquet files."""

from qpx.writers.base import BaseWriter
from qpx.core.data import PsmSchema


class PsmWriter(BaseWriter):
    _schema_class = PsmSchema
