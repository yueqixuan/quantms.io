"""FeatureWriter — writes feature.parquet files."""

from qpx.writers.base import BaseWriter
from qpx.core.data import FeatureSchema


class FeatureWriter(BaseWriter):
    _schema_class = FeatureSchema
