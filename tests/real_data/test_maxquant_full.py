"""MaxQuant integration test — PXD001819 full dataset.

Runs MaxQuant conversion on the full PXD001819 dataset (compressed evidence,
proteinGroups, proper SDRF with 27 runs), produces a QPX dataset, and
validates structures/schemas/queries.

Files:
    tests/examples/maxquant/maxquant_full/evidence.txt.gz     (~22 MB, 244k rows)
    tests/examples/maxquant/maxquant_full/proteinGroups.txt    (1 074 PGs)
    tests/examples/maxquant/maxquant_full/PXD001819.sdrf.tsv   (27 runs)

Markers:
    @pytest.mark.large_data — skipped in CI with ``-m "not large_data"``
"""

import pyarrow.parquet as pq
import pytest
from pathlib import Path

from qpx.dataset import Dataset

EXAMPLES_DIR = Path(__file__).parent.parent / "examples" / "maxquant" / "maxquant_full"

_EVIDENCE = EXAMPLES_DIR / "evidence.txt.gz"
_PG = EXAMPLES_DIR / "proteinGroups.txt"
_SDRF = EXAMPLES_DIR / "PXD001819.sdrf.tsv"

_PREFIX = "mq_full"


# ---------------------------------------------------------------------------
# Module-scoped fixture: convert once, open as Dataset
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def dataset(tmp_path_factory):
    """Run full MaxQuant conversion and open as Dataset."""
    if not _EVIDENCE.exists():
        pytest.skip(f"Test data not found: {_EVIDENCE}")

    from qpx.converters.maxquant.converter import MaxQuantConverter

    output = tmp_path_factory.mktemp("mq_full")

    converter = MaxQuantConverter()
    converter.convert(
        output_folder=output,
        evidence_file=_EVIDENCE,
        protein_groups_file=_PG if _PG.exists() else None,
        sdrf_file=_SDRF if _SDRF.exists() else None,
        output_prefix=_PREFIX,
        project_accession="PXD001819",
    )

    ds = Dataset(output)
    yield ds
    ds.close()


@pytest.fixture(scope="module")
def feature_table(dataset):
    """Raw PyArrow table for feature.parquet."""
    path = Path(dataset.path) / f"{_PREFIX}.feature.parquet"
    if not path.exists():
        pytest.skip("feature.parquet was not produced")
    return pq.read_table(str(path))


@pytest.fixture(scope="module")
def pg_table(dataset):
    """Raw PyArrow table for pg.parquet."""
    path = Path(dataset.path) / f"{_PREFIX}.pg.parquet"
    if not path.exists():
        pytest.skip("pg.parquet was not produced")
    return pq.read_table(str(path))


# ---------------------------------------------------------------------------
# Feature conversion tests
# ---------------------------------------------------------------------------


@pytest.mark.large_data
class TestFeatureConversion:
    """Validate feature.parquet from full MaxQuant dataset."""

    def test_has_rows(self, feature_table):
        assert feature_table.num_rows > 1000, (
            f"Expected >1000 feature rows, got {feature_table.num_rows}"
        )

    def test_key_columns_present(self, feature_table):
        expected = {"sequence", "charge", "calculated_mz", "intensities", "run_file_name"}
        missing = expected - set(feature_table.column_names)
        assert not missing, f"Missing columns: {missing}"

    def test_sequences_nonempty(self, feature_table):
        seqs = feature_table.column("sequence").to_pylist()[:200]
        assert all(isinstance(s, str) and len(s) > 0 for s in seqs)

    def test_charges_in_range(self, feature_table):
        charges = feature_table.column("charge").to_pylist()[:200]
        assert all(1 <= c <= 10 for c in charges)

    def test_intensities_nonnegative(self, feature_table):
        for row in feature_table.column("intensities").to_pylist()[:200]:
            if row is None:
                continue
            for entry in row:
                assert entry["intensity"] >= 0

    def test_multiple_runs(self, feature_table):
        """Full dataset should have multiple run files."""
        runs = set(feature_table.column("run_file_name").to_pylist())
        assert len(runs) > 5, f"Expected >5 runs, got {len(runs)}"

    def test_schema_validation(self, feature_table):
        from qpx.core.data import FeatureSchema

        errors = FeatureSchema.validate(feature_table)
        assert not errors, f"Schema validation errors: {errors}"


# ---------------------------------------------------------------------------
# PG conversion tests
# ---------------------------------------------------------------------------


@pytest.mark.large_data
class TestPgConversion:
    """Validate pg.parquet from full MaxQuant dataset."""

    def test_has_rows(self, pg_table):
        assert pg_table.num_rows > 100, (
            f"Expected >100 PG rows, got {pg_table.num_rows}"
        )

    def test_key_columns_present(self, pg_table):
        expected = {"pg_accessions", "anchor_protein", "intensities"}
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

    def test_has_feature(self, dataset):
        assert dataset.feature is not None
        assert dataset.feature.count() > 1000

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

    def test_has_sample(self, dataset):
        """SDRF should produce sample.parquet."""
        if dataset.sample is None:
            pytest.skip("sample.parquet not produced (SDRF incomplete)")
        assert dataset.sample.count() > 0

    def test_has_run(self, dataset):
        """SDRF should produce run.parquet."""
        if dataset.run is None:
            pytest.skip("run.parquet not produced (SDRF incomplete)")
        assert dataset.run.count() > 0

    def test_available_structures(self, dataset):
        available = set(dataset.available_structures)
        expected = {"feature", "pg", "ontology", "provenance", "dataset"}
        assert expected.issubset(available), f"Missing: {expected - available}"


# ---------------------------------------------------------------------------
# Feature querying via Dataset API
# ---------------------------------------------------------------------------


@pytest.mark.large_data
class TestFeatureQuerying:
    """Query features through the Dataset API."""

    def test_filter_by_charge(self, dataset):
        c2 = dataset.feature.filter("charge = 2").count()
        c3 = dataset.feature.filter("charge = 3").count()
        total = dataset.feature.count()
        assert c2 > 0
        assert c2 + c3 <= total

    def test_unique_sequences(self, dataset):
        result = dataset.sql("SELECT COUNT(DISTINCT sequence) AS n FROM feature")
        n = result.to_df()["n"].iloc[0]
        assert n > 100, f"Expected >100 unique sequences, got {n}"

    def test_unique_runs(self, dataset):
        result = dataset.sql(
            "SELECT COUNT(DISTINCT run_file_name) AS n FROM feature"
        )
        n = result.to_df()["n"].iloc[0]
        assert n > 5, f"Expected >5 runs, got {n}"

    def test_select_columns(self, dataset):
        df = dataset.feature.select("sequence", "charge").limit(20).to_df()
        assert set(df.columns) == {"sequence", "charge"}
        assert len(df) == 20


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

    def test_pg_to_dataframe(self, dataset):
        df = dataset.pg.limit(20).to_df()
        assert len(df) == 20
        assert "anchor_protein" in df.columns


# ---------------------------------------------------------------------------
# Provenance tests
# ---------------------------------------------------------------------------


@pytest.mark.large_data
class TestProvenance:
    """Validate provenance records."""

    def test_has_maxquant_step(self, dataset):
        df = dataset.provenance.to_df()
        tool_names = df["tool_name"].tolist()
        assert any("MaxQuant" in str(t) for t in tool_names if t)

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
        assert df["project_accession"].iloc[0] == "PXD001819"

    def test_software_name(self, dataset):
        df = dataset.dataset_meta.to_df()
        assert df["software_name"].iloc[0] == "MaxQuant"

    def test_creation_date(self, dataset):
        df = dataset.dataset_meta.to_df()
        assert df["creation_date"].iloc[0] is not None
