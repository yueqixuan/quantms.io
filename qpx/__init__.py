"""QPX — quantitative proteomics data format library."""

from __future__ import annotations

from qpx._version import __version__  # noqa: F401
from qpx.dataset import Dataset

# --- Writing ---
from qpx.writers.feature import FeatureWriter  # noqa: F401
from qpx.writers.psm import PsmWriter  # noqa: F401
from qpx.writers.pg import PgWriter  # noqa: F401
from qpx.writers.mz import MzWriter  # noqa: F401
from qpx.writers.sample import SampleWriter  # noqa: F401
from qpx.writers.run import RunWriter  # noqa: F401
from qpx.writers.dataset import DatasetWriter  # noqa: F401
from qpx.writers.ontology import OntologyWriter  # noqa: F401
from qpx.writers.provenance import ProvenanceWriter  # noqa: F401
from qpx.writers.pepmap import PepMapWriter  # noqa: F401

# --- Data structure classes ---
from qpx.core.data.feature import Feature  # noqa: F401
from qpx.core.data.psm import PSM  # noqa: F401
from qpx.core.data.pg import PG  # noqa: F401
from qpx.core.data.mz import MzSpectra  # noqa: F401
from qpx.core.data.sample import Sample  # noqa: F401
from qpx.core.data.run import Run  # noqa: F401
from qpx.core.data.pepmap import PepMap  # noqa: F401

# --- Validation ---
from qpx.core.data.schema import ValidationResult, ValidationIssue  # noqa: F401

# --- Collections ---
from qpx.collection import DatasetCollection  # noqa: F401

# --- Views ---
from qpx import views  # noqa: F401


# --- Reading (data structures) ---


def open(path: str, structures: list[str] | None = None, **kwargs) -> Dataset:
    """Open a QPX dataset directory."""
    return Dataset(path, structures=structures, **kwargs)


def read_feature(path: str, **kwargs):
    """Open a single feature.parquet file as a standalone data structure."""
    return Feature.from_file(path, **kwargs)


def read_psm(path: str, **kwargs):
    """Open a single psm.parquet file."""
    return PSM.from_file(path, **kwargs)


def read_pg(path: str, **kwargs):
    """Open a single pg.parquet file."""
    return PG.from_file(path, **kwargs)


def read_mz(path: str, **kwargs):
    """Open a single mz.parquet file."""
    return MzSpectra.from_file(path, **kwargs)


def read_sample(path: str, **kwargs):
    """Open a single sample.parquet file."""
    return Sample.from_file(path, **kwargs)


def read_run(path: str, **kwargs):
    """Open a single run.parquet file."""
    return Run.from_file(path, **kwargs)


def read_pepmap(path: str, **kwargs):
    """Open a single pepmap.parquet file."""
    return PepMap.from_file(path, **kwargs)
