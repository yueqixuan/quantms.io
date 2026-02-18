"""DIA-NN integration test — PXD036609 full dataset.

Runs PG conversion from the full DIA-NN output (compressed report, 195k rows,
27 runs), produces a QPX dataset, and validates structures/schemas/queries.

This dataset does NOT include mzml_info files, so feature conversion is
skipped. Only PG conversion, ontology, provenance, and dataset are tested.

Files:
    tests/examples/diann/full/diann_report.tsv.gz          (~35 MB, 195k rows)
    tests/examples/diann/full/PXD036609.sdrf.tsv           (27 runs)
    tests/examples/diann/full/diann_report.pg_matrix.tsv   (846 PGs)

Markers:
    @pytest.mark.large_data — skipped in CI with ``-m "not large_data"``
"""

import pyarrow.parquet as pq
import pytest
from pathlib import Path

from qpx.dataset import Dataset

EXAMPLES_DIR = Path(__file__).parent.parent / "examples" / "diann" / "full"

_REPORT_GZ = EXAMPLES_DIR / "diann_report.tsv.gz"
_SDRF = EXAMPLES_DIR / "PXD036609.sdrf.tsv"
_PG_MATRIX = EXAMPLES_DIR / "diann_report.pg_matrix.tsv"

_PREFIX = "diann_full"


# ---------------------------------------------------------------------------
# Module-scoped fixture: convert once, open as Dataset
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def dataset(tmp_path_factory):
    """Run DIA-NN PG conversion and open as Dataset."""
    if not _REPORT_GZ.exists():
        pytest.skip(f"Test data not found: {_REPORT_GZ}")

    from qpx.converters.diann.converter import DiaNNConverter

    output = tmp_path_factory.mktemp("diann_full")

    converter = DiaNNConverter(
        report_path=str(_REPORT_GZ),
        sdrf_path=str(_SDRF),
    )
    converter.convert_pg(
        pg_matrix_path=str(_PG_MATRIX),
        output_folder=str(output),
        output_prefix=_PREFIX,
    )
    converter.write_ontology(output, prefix=_PREFIX)
    converter.write_provenance(output, prefix=_PREFIX)
    converter.write_dataset(output, prefix=_PREFIX, project_accession="PXD036609")

    ds = Dataset(output)
    yield ds
    ds.close()


@pytest.fixture(scope="module")
def pg_table(dataset):
    """Raw PyArrow table for pg.parquet."""
    path = Path(dataset.path) / f"{_PREFIX}.pg.parquet"
    if not path.exists():
        pytest.skip("pg.parquet was not produced")
    return pq.read_table(str(path))


# ---------------------------------------------------------------------------
# PG conversion tests
# ---------------------------------------------------------------------------


@pytest.mark.large_data
class TestPgConversion:
    """Validate pg.parquet from full DIA-NN dataset."""

    def test_has_rows(self, pg_table):
        assert pg_table.num_rows > 100, (
            f"Expected >100 PG rows, got {pg_table.num_rows}"
        )

    def test_key_columns_present(self, pg_table):
        expected = {"pg_accessions", "anchor_protein", "run_file_name", "intensities"}
        missing = expected - set(pg_table.column_names)
        assert not missing, f"Missing columns: {missing}"

    def test_protein_accessions_nonempty(self, pg_table):
        for accs in pg_table.column("pg_accessions").to_pylist()[:100]:
            assert accs is not None and len(accs) > 0
            assert all(isinstance(a, str) and len(a) > 0 for a in accs)

    def test_anchor_protein_nonempty(self, pg_table):
        anchors = pg_table.column("anchor_protein").to_pylist()[:100]
        assert all(isinstance(a, str) and len(a) > 0 for a in anchors)

    def test_intensities_structure(self, pg_table):
        for row in pg_table.column("intensities").to_pylist()[:100]:
            assert row is not None and len(row) > 0
            for entry in row:
                assert "label" in entry
                assert "intensity" in entry

    def test_intensities_nonnegative(self, pg_table):
        for row in pg_table.column("intensities").to_pylist()[:100]:
            if row is None:
                continue
            for entry in row:
                assert entry["intensity"] >= 0

    def test_multiple_runs(self, pg_table):
        """Full dataset should have multiple run files."""
        runs = set(pg_table.column("run_file_name").to_pylist())
        assert len(runs) > 5, f"Expected >5 runs, got {len(runs)}"

    def test_schema_validation(self, pg_table):
        from qpx.core.data import PgSchema

        errors = PgSchema.validate(pg_table)
        assert not errors, f"Schema validation errors: {errors}"


# ---------------------------------------------------------------------------
# Dataset structure tests
# ---------------------------------------------------------------------------


@pytest.mark.large_data
class TestDatasetStructure:
    """Validate all expected structures."""

    def test_has_pg(self, dataset):
        assert dataset.pg is not None
        assert dataset.pg.count() > 100

    def test_has_ontology(self, dataset):
        assert dataset.ontology is not None
        assert dataset.ontology.count() > 0

    def test_has_provenance(self, dataset):
        assert dataset.provenance is not None
        assert dataset.provenance.count() > 0

    def test_has_dataset_meta(self, dataset):
        assert dataset.dataset_meta is not None

    def test_available_structures(self, dataset):
        available = set(dataset.available_structures)
        expected = {"pg", "ontology", "provenance", "dataset"}
        assert expected.issubset(available), f"Missing: {expected - available}"


# ---------------------------------------------------------------------------
# PG querying via Dataset API
# ---------------------------------------------------------------------------


@pytest.mark.large_data
class TestPgQuerying:
    """Query protein groups through the Dataset API."""

    def test_unique_anchor_proteins(self, dataset):
        result = dataset.sql(
            "SELECT COUNT(DISTINCT anchor_protein) AS n FROM pg"
        )
        n = result.to_df()["n"].iloc[0]
        assert n > 50, f"Expected >50 unique proteins, got {n}"

    def test_unique_runs(self, dataset):
        result = dataset.sql(
            "SELECT COUNT(DISTINCT run_file_name) AS n FROM pg"
        )
        n = result.to_df()["n"].iloc[0]
        assert n > 5, f"Expected >5 runs, got {n}"

    def test_select_columns(self, dataset):
        df = dataset.pg.select("anchor_protein", "run_file_name").limit(20).to_df()
        assert set(df.columns) == {"anchor_protein", "run_file_name"}
        assert len(df) == 20

    def test_pg_to_dataframe(self, dataset):
        df = dataset.pg.limit(50).to_df()
        assert len(df) == 50
        assert "anchor_protein" in df.columns
        assert "pg_accessions" in df.columns


# ---------------------------------------------------------------------------
# Provenance tests
# ---------------------------------------------------------------------------


@pytest.mark.large_data
class TestProvenance:
    """Validate provenance records."""

    def test_has_diann_step(self, dataset):
        df = dataset.provenance.to_df()
        tool_names = df["tool_name"].tolist()
        assert any("DIA-NN" in str(t) for t in tool_names if t)

    def test_has_qpx_step(self, dataset):
        df = dataset.provenance.to_df()
        tool_names = df["tool_name"].tolist()
        assert any("qpx" in str(t).lower() for t in tool_names if t)

    def test_step_order_positive(self, dataset):
        df = dataset.provenance.to_df()
        assert (df["step_order"] > 0).all()


# ---------------------------------------------------------------------------
# Ontology tests
# ---------------------------------------------------------------------------


@pytest.mark.large_data
class TestOntology:
    """Validate ontology entries."""

    def test_has_entries(self, dataset):
        df = dataset.ontology.to_df()
        assert len(df) > 0

    def test_has_field_name(self, dataset):
        df = dataset.ontology.to_df()
        assert "field_name" in df.columns


# ---------------------------------------------------------------------------
# Dataset metadata tests
# ---------------------------------------------------------------------------


@pytest.mark.large_data
class TestDatasetMetadata:
    """Validate dataset.parquet metadata."""

    def test_project_accession(self, dataset):
        df = dataset.dataset_meta.to_df()
        assert df["project_accession"].iloc[0] == "PXD036609"

    def test_software_name(self, dataset):
        df = dataset.dataset_meta.to_df()
        assert df["software_name"].iloc[0] == "DIA-NN"

    def test_creation_date(self, dataset):
        df = dataset.dataset_meta.to_df()
        assert df["creation_date"].iloc[0] is not None
