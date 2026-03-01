"""PgWriter — writes pg.parquet files."""

from qpx.writers.base import BaseWriter
from qpx.core.data import PgSchema


class PgWriter(BaseWriter):
    _schema_class = PgSchema
