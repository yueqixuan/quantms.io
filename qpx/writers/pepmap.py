"""PepMapWriter — writes pepmap.parquet files."""

from qpx.writers.base import BaseWriter
from qpx.core.data import PepMapSchema


class PepMapWriter(BaseWriter):
    _schema_class = PepMapSchema
