"""Tests for QPX data layer: BaseStructure, Feature, PSM, PG and domain methods."""

import pandas as pd
import pyarrow as pa

from qpx.core.data.feature import Feature
from qpx.core.data.psm import PSM
from qpx.core.data.pg import PG
from qpx.core.data.sample import Sample
from qpx.core.data.run import Run
from qpx.core.data.dataset import DatasetMeta
from qpx.core.data.ontology import Ontology
from qpx.core.data.provenance import Provenance

# ---------------------------------------------------------------------------
# Feature: from_file and basic operations
# ---------------------------------------------------------------------------


class TestFeatureFromFile:
    """Test Feature data structure opened from standalone Parquet."""

    def test_from_file_opens_parquet(self, feature_parquet):
        """Feature.from_file opens a standalone Parquet and connects to DuckDB."""
        feat = Feature.from_file(feature_parquet)
        assert isinstance(feat, Feature)
        assert feat.count() == 3

    def test_from_file_with_engine_kwargs(self, feature_parquet):
        """Feature.from_file accepts engine kwargs."""
        feat = Feature.from_file(feature_parquet, memory_limit="1GB", threads=1)
        assert feat.count() == 3

    def test_schema_classmethod(self):
        """Feature.schema() returns the FeatureSchema arrow schema."""
        schema = Feature.schema()
        assert isinstance(schema, pa.Schema)
        assert "sequence" in schema.names


class TestFeatureFilter:
    """Test Feature.filter and domain-specific filters."""

    def test_filter_returns_filtered_results(self, feature_parquet):
        """filter() with a condition returns only matching rows."""
        feat = Feature.from_file(feature_parquet)
        filtered = feat.filter("charge = 2")
        assert filtered.count() == 2  # Two rows have charge=2

    def test_filter_returns_feature_instance(self, feature_parquet):
        """filter() returns a Feature instance (not BaseStructure)."""
        feat = Feature.from_file(feature_parquet)
        filtered = feat.filter("charge = 2")
        assert isinstance(filtered, Feature)

    def test_by_protein(self, feature_parquet):
        """by_protein() filters by anchor protein."""
        feat = Feature.from_file(feature_parquet)
        p12345 = feat.by_protein("P12345")
        assert p12345.count() == 2  # Two rows with P12345

    def test_by_run(self, feature_parquet):
        """by_run() filters by run_file_name."""
        feat = Feature.from_file(feature_parquet)
        run01 = feat.by_run("run_01")
        assert run01.count() == 2  # Two rows from run_01

    def test_filter_chaining(self, feature_parquet):
        """Chaining multiple filters works correctly."""
        feat = Feature.from_file(feature_parquet)
        result = feat.filter("charge = 2").filter("is_decoy = false")
        # Only the first row matches: charge=2 AND not decoy
        assert result.count() == 1


class TestFeatureSelect:
    """Test Feature.select returns subset of columns."""

    def test_select_returns_subset(self, feature_parquet):
        """select() returns only the specified columns."""
        feat = Feature.from_file(feature_parquet)
        selected = feat.select("sequence", "charge")
        df = selected.to_df()
        assert set(df.columns) == {"sequence", "charge"}
        assert len(df) == 3


class TestFeatureLimit:
    """Test Feature.limit limits row count."""

    def test_limit(self, feature_parquet):
        """limit() restricts the number of returned rows."""
        feat = Feature.from_file(feature_parquet)
        limited = feat.limit(2)
        assert limited.count() == 2


class TestFeatureCount:
    """Test Feature.count returns correct count."""

    def test_count_all(self, feature_parquet):
        """count() on unfiltered data returns total rows."""
        feat = Feature.from_file(feature_parquet)
        assert feat.count() == 3

    def test_count_after_filter(self, feature_parquet):
        """count() after filter returns filtered count."""
        feat = Feature.from_file(feature_parquet)
        assert feat.filter("is_decoy = true").count() == 1

    def test_len(self, feature_parquet):
        """len() works on Feature (delegates to count)."""
        feat = Feature.from_file(feature_parquet)
        assert len(feat) == 3


class TestFeatureMaterialization:
    """Test Feature materialization to DataFrame and Arrow."""

    def test_to_df(self, feature_parquet):
        """to_df() returns a pandas DataFrame."""
        feat = Feature.from_file(feature_parquet)
        df = feat.to_df()
        assert isinstance(df, pd.DataFrame)
        assert len(df) == 3

    def test_to_arrow(self, feature_parquet):
        """to_arrow() returns a PyArrow Table."""
        feat = Feature.from_file(feature_parquet)
        table = feat.to_arrow()
        assert isinstance(table, pa.Table)
        assert table.num_rows == 3


class TestFeatureDomainMethods:
    """Test Feature domain-specific methods."""

    def test_unique_proteins(self, feature_parquet):
        """unique_proteins() returns distinct anchor proteins."""
        feat = Feature.from_file(feature_parquet)
        proteins = feat.unique_proteins()
        assert set(proteins) == {"P12345", "P67890"}

    def test_peptide_intensities(self, feature_parquet):
        """peptide_intensities() flattens intensity array to long format."""
        feat = Feature.from_file(feature_parquet)
        result = feat.peptide_intensities()
        df = result.to_df()
        # Should have one row per (feature, label) combo
        assert "label" in df.columns
        assert "intensity" in df.columns
        assert len(df) == 3  # 3 features, each with 1 intensity

    def test_file_metadata(self, feature_parquet):
        """file_metadata property returns Parquet footer metadata."""
        feat = Feature.from_file(feature_parquet)
        meta = feat.file_metadata
        assert isinstance(meta, dict)
        assert "file_type" in meta
        assert meta["file_type"] == "feature_file"

    def test_repr(self, feature_parquet):
        """Feature repr shows class name, path, and row count."""
        feat = Feature.from_file(feature_parquet)
        r = repr(feat)
        assert "Feature" in r
        assert "rows=3" in r


# ---------------------------------------------------------------------------
# BaseStructure.iter_batches
# ---------------------------------------------------------------------------


class TestIterBatches:
    """Test iter_batches partitions data by column."""

    def test_iter_batches_by_run(self, feature_parquet):
        """iter_batches partitions by run_file_name."""
        feat = Feature.from_file(feature_parquet)
        batches = list(feat.iter_batches(partition_by="run_file_name", batch_size=10))
        # Should get one batch containing all distinct run values
        assert len(batches) == 1
        keys, df = batches[0]
        assert set(keys) == {"run_01", "run_02"}
        assert len(df) == 3

    def test_iter_batches_small_batch_size(self, feature_parquet):
        """iter_batches with batch_size=1 yields one batch per distinct value."""
        feat = Feature.from_file(feature_parquet)
        batches = list(feat.iter_batches(partition_by="run_file_name", batch_size=1))
        assert len(batches) == 2  # Two distinct runs

        for keys, df in batches:
            assert len(keys) == 1
            assert isinstance(df, pd.DataFrame)
            assert len(df) > 0

    def test_iter_batches_yields_correct_data(self, feature_parquet):
        """Each batch contains only data for its partition values."""
        feat = Feature.from_file(feature_parquet)
        for keys, df in feat.iter_batches(partition_by="run_file_name", batch_size=1):
            run_values = df["run_file_name"].unique()
            assert set(run_values) == set(keys)


# ---------------------------------------------------------------------------
# BaseStructure immutability
# ---------------------------------------------------------------------------


class TestBaseStructureImmutability:
    """Test that query operations return new instances (immutable)."""

    def test_filter_returns_new_instance(self, feature_parquet):
        """filter() returns a new Feature, not mutating the original."""
        feat = Feature.from_file(feature_parquet)
        original_count = feat.count()
        filtered = feat.filter("charge > 10")  # No rows match
        assert filtered.count() == 0
        assert feat.count() == original_count  # Original unchanged

    def test_select_returns_new_instance(self, feature_parquet):
        """select() returns a new Feature, not mutating the original."""
        feat = Feature.from_file(feature_parquet)
        selected = feat.select("sequence")
        assert selected is not feat

    def test_limit_returns_new_instance(self, feature_parquet):
        """limit() returns a new Feature, not mutating the original."""
        feat = Feature.from_file(feature_parquet)
        limited = feat.limit(1)
        assert limited.count() == 1
        assert feat.count() == 3  # Original unchanged


# ---------------------------------------------------------------------------
# Feature.order_by
# ---------------------------------------------------------------------------


class TestFeatureOrderBy:
    """Test Feature.order_by for sorting."""

    def test_order_by_asc(self, feature_parquet):
        """order_by returns rows in ascending order."""
        feat = Feature.from_file(feature_parquet)
        ordered = feat.order_by("charge")
        df = ordered.to_df()
        charges = df["charge"].tolist()
        assert charges == sorted(charges)

    def test_order_by_desc(self, feature_parquet):
        """order_by(desc=True) returns rows in descending order."""
        feat = Feature.from_file(feature_parquet)
        ordered = feat.order_by("charge", desc=True)
        df = ordered.to_df()
        charges = df["charge"].tolist()
        assert charges == sorted(charges, reverse=True)


# ---------------------------------------------------------------------------
# PSM data structure tests
# ---------------------------------------------------------------------------


class TestPSM:
    """Test PSM data structure and domain methods."""

    def test_from_file(self, psm_parquet):
        """PSM.from_file opens a standalone PSM Parquet."""
        psm = PSM.from_file(psm_parquet)
        assert isinstance(psm, PSM)
        assert psm.count() == 3

    def test_targets_only(self, psm_parquet):
        """targets_only() excludes decoy PSMs."""
        psm = PSM.from_file(psm_parquet)
        targets = psm.targets_only()
        assert targets.count() == 2
        df = targets.to_df()
        assert all(df["is_decoy"] == False)  # noqa: E712

    def test_by_run(self, psm_parquet):
        """by_run() filters PSMs by run_file_name."""
        psm = PSM.from_file(psm_parquet)
        run01 = psm.by_run("run_01")
        assert run01.count() == 2

    def test_filter(self, psm_parquet):
        """filter() works on PSM."""
        psm = PSM.from_file(psm_parquet)
        filtered = psm.filter("charge = 2")
        assert filtered.count() == 2

    def test_to_df(self, psm_parquet):
        """to_df() works on PSM."""
        psm = PSM.from_file(psm_parquet)
        df = psm.to_df()
        assert isinstance(df, pd.DataFrame)
        assert len(df) == 3


# ---------------------------------------------------------------------------
# PG data structure tests
# ---------------------------------------------------------------------------


class TestPG:
    """Test PG data structure and domain methods."""

    def test_from_file(self, pg_parquet):
        """PG.from_file opens a standalone PG Parquet."""
        pg = PG.from_file(pg_parquet)
        assert isinstance(pg, PG)
        assert pg.count() == 3

    def test_targets_only(self, pg_parquet):
        """targets_only() excludes decoy proteins."""
        pg = PG.from_file(pg_parquet)
        targets = pg.targets_only()
        assert targets.count() == 2

    def test_by_protein(self, pg_parquet):
        """by_protein() filters by anchor protein."""
        pg = PG.from_file(pg_parquet)
        p12345 = pg.by_protein("P12345")
        assert p12345.count() == 1

    def test_by_run(self, pg_parquet):
        """by_run() filters by run_file_name."""
        pg = PG.from_file(pg_parquet)
        run01 = pg.by_run("run_01")
        assert run01.count() == 2  # Two PGs in run_01

    def test_protein_intensities(self, pg_parquet):
        """protein_intensities() flattens PG intensities to long format."""
        pg = PG.from_file(pg_parquet)
        result = pg.protein_intensities()
        df = result.to_df()
        assert "label" in df.columns
        assert "intensity" in df.columns
        assert "anchor_protein" in df.columns
        # 3 PG rows, each with 1 intensity label
        assert len(df) == 3

    def test_to_arrow(self, pg_parquet):
        """to_arrow() works on PG."""
        pg = PG.from_file(pg_parquet)
        table = pg.to_arrow()
        assert isinstance(table, pa.Table)
        assert table.num_rows == 3


# ---------------------------------------------------------------------------
# Sample data structure tests
# ---------------------------------------------------------------------------


class TestSample:
    """Test Sample data structure."""

    def test_from_file(self, sample_parquet):
        """Sample.from_file opens a standalone sample Parquet."""
        sample = Sample.from_file(sample_parquet)
        assert isinstance(sample, Sample)
        assert sample.count() == 2

    def test_by_organism(self, sample_parquet):
        """by_organism() filters samples by organism."""
        sample = Sample.from_file(sample_parquet)
        human = sample.by_organism("Homo sapiens")
        assert human.count() == 2  # Both samples are Homo sapiens

    def test_to_df(self, sample_parquet):
        """to_df() works on Sample."""
        sample = Sample.from_file(sample_parquet)
        df = sample.to_df()
        assert isinstance(df, pd.DataFrame)
        assert "sample_accession" in df.columns


# ---------------------------------------------------------------------------
# Run data structure tests
# ---------------------------------------------------------------------------


class TestRun:
    """Test Run data structure."""

    def test_from_file(self, run_parquet):
        """Run.from_file opens a standalone run Parquet."""
        run = Run.from_file(run_parquet)
        assert isinstance(run, Run)
        assert run.count() == 2

    def test_to_df(self, run_parquet):
        """to_df() works on Run."""
        run = Run.from_file(run_parquet)
        df = run.to_df()
        assert isinstance(df, pd.DataFrame)
        assert "run_file_name" in df.columns
        assert "run_accession" in df.columns


# ---------------------------------------------------------------------------
# DatasetMeta data structure tests
# ---------------------------------------------------------------------------


class TestDatasetMeta:
    """Test DatasetMeta data structure."""

    def test_from_file(self, dataset_parquet):
        """DatasetMeta.from_file opens a standalone dataset Parquet."""
        ds = DatasetMeta.from_file(dataset_parquet)
        assert isinstance(ds, DatasetMeta)
        assert ds.count() == 1

    def test_to_df(self, dataset_parquet):
        """to_df() works on DatasetMeta."""
        ds = DatasetMeta.from_file(dataset_parquet)
        df = ds.to_df()
        assert "project_accession" in df.columns
        assert df["project_accession"].iloc[0] == "PXD000001"


# ---------------------------------------------------------------------------
# Ontology data structure tests
# ---------------------------------------------------------------------------


class TestOntology:
    """Test Ontology data structure."""

    def test_from_file(self, ontology_parquet):
        """Ontology.from_file opens a standalone ontology Parquet."""
        ont = Ontology.from_file(ontology_parquet)
        assert isinstance(ont, Ontology)
        assert ont.count() == 2

    def test_by_view(self, ontology_parquet):
        """by_view() filters ontology entries by view name."""
        ont = Ontology.from_file(ontology_parquet)
        feat_entries = ont.by_view("feature")
        assert feat_entries.count() == 1


# ---------------------------------------------------------------------------
# Provenance data structure tests
# ---------------------------------------------------------------------------


class TestProvenance:
    """Test Provenance data structure."""

    def test_from_file(self, provenance_parquet):
        """Provenance.from_file opens a standalone provenance Parquet."""
        prov = Provenance.from_file(provenance_parquet)
        assert isinstance(prov, Provenance)
        assert prov.count() == 2

    def test_to_df(self, provenance_parquet):
        """to_df() works on Provenance."""
        prov = Provenance.from_file(provenance_parquet)
        df = prov.to_df()
        assert "step_order" in df.columns
        assert "tool_name" in df.columns


# ---------------------------------------------------------------------------
# Join tests
# ---------------------------------------------------------------------------


class TestBaseStructureJoin:
    """Test join operations between data structures."""

    def test_feature_join_run(self, feature_parquet, run_parquet):
        """Feature.join(Run) joins on run_file_name."""
        feat = Feature.from_file(feature_parquet)
        # Register run in the same engine
        feat._engine.register_parquet("run", run_parquet)
        run = Run(engine=feat._engine, table_name="run", file_path=run_parquet)

        joined = feat.join(run, on="run_file_name")
        df = joined.to_df()
        assert len(df) == 3
        assert "run_accession" in df.columns
