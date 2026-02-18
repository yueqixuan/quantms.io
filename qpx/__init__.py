"""QPX — quantitative proteomics data format library."""

from qpx._version import __version__
from qpx.dataset import Dataset

# --- Reading (data structures) ---


def open(path: str, structures: list[str] | None = None, **kwargs) -> Dataset:
    """Open a QPX dataset directory."""
    return Dataset(path, structures=structures, **kwargs)


def read_feature(path: str, **kwargs):
    """Open a single feature.parquet file as a standalone data structure."""
    from qpx.core.data.feature import Feature
    return Feature.from_file(path, **kwargs)


def read_psm(path: str, **kwargs):
    """Open a single psm.parquet file."""
    from qpx.core.data.psm import PSM
    return PSM.from_file(path, **kwargs)


def read_pg(path: str, **kwargs):
    """Open a single pg.parquet file."""
    from qpx.core.data.pg import PG
    return PG.from_file(path, **kwargs)


def read_mz(path: str, **kwargs):
    """Open a single mz.parquet file."""
    from qpx.core.data.mz import MzSpectra
    return MzSpectra.from_file(path, **kwargs)


def read_sample(path: str, **kwargs):
    """Open a single sample.parquet file."""
    from qpx.core.data.sample import Sample
    return Sample.from_file(path, **kwargs)


def read_run(path: str, **kwargs):
    """Open a single run.parquet file."""
    from qpx.core.data.run import Run
    return Run.from_file(path, **kwargs)


def read_pepmap(path: str, **kwargs):
    """Open a single pepmap.parquet file."""
    from qpx.core.data.pepmap import PepMap
    return PepMap.from_file(path, **kwargs)


# --- Writing ---
from qpx.writers.feature import FeatureWriter
from qpx.writers.psm import PsmWriter
from qpx.writers.pg import PgWriter
from qpx.writers.mz import MzWriter
from qpx.writers.sample import SampleWriter
from qpx.writers.run import RunWriter
from qpx.writers.dataset import DatasetWriter
from qpx.writers.ontology import OntologyWriter
from qpx.writers.provenance import ProvenanceWriter
from qpx.writers.pepmap import PepMapWriter

# --- Data structure classes ---
from qpx.core.data.feature import Feature
from qpx.core.data.psm import PSM
from qpx.core.data.pg import PG
from qpx.core.data.mz import MzSpectra
from qpx.core.data.sample import Sample
from qpx.core.data.run import Run
from qpx.core.data.pepmap import PepMap

# --- Validation ---
from qpx.core.data.schema import ValidationResult, ValidationIssue

# --- Collections ---
from qpx.collection import DatasetCollection

# --- Views ---
from qpx import views
