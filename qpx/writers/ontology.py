"""OntologyWriter — writes ontology.parquet files."""

from qpx.writers.base import BaseWriter
from qpx.core.data import OntologySchema


class OntologyWriter(BaseWriter):
    _schema_class = OntologySchema
