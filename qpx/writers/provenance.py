"""ProvenanceWriter — writes provenance.parquet files."""

from qpx.writers.base import BaseWriter
from qpx.core.data import ProvenanceSchema


class ProvenanceWriter(BaseWriter):
    _schema_class = ProvenanceSchema
