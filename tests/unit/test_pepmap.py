"""Peptide-protein mapping tests."""
import pyarrow.parquet as pq


def make_mapping_record(**overrides):
    defaults = {
        "sequence": "PEPTIDEK",
        "protein_accession": "P12345",
        "protein_name": "BRCA1_HUMAN",
        "gene_name": "BRCA1",
        "start_position": 100,
        "end_position": 107,
        "is_unique": True,
        "is_proteotypic": True,
    }
    defaults.update(overrides)
    return defaults


class TestPepMapSchema:
    def test_schema_loads(self):
        from qpx.core.data import PepMapSchema

        field_names = [f.name for f in PepMapSchema.get_arrow_schema()]
        assert "sequence" in field_names
        assert "protein_accession" in field_names
        assert "is_unique" in field_names
        assert "is_proteotypic" in field_names

    def test_round_trip(self, tmp_path):
        from qpx.writers import PepMapWriter

        records = [
            make_mapping_record(),
            make_mapping_record(
                sequence="ANOTHERK",
                protein_accession="P67890",
                gene_name="TP53",
                is_unique=False,
                is_proteotypic=False,
            ),
        ]
        out = tmp_path / "test.pepmap.parquet"
        with PepMapWriter(out) as w:
            w.write_batch(records)
        table = pq.read_table(out)
        assert table.num_rows == 2
        assert table.column("sequence")[0].as_py() == "PEPTIDEK"
        assert table.column("protein_accession")[1].as_py() == "P67890"

    def test_nullable_fields(self, tmp_path):
        from qpx.writers import PepMapWriter

        record = make_mapping_record(
            protein_name=None,
            gene_name=None,
            start_position=None,
            end_position=None,
            is_unique=None,
            is_proteotypic=None,
        )
        out = tmp_path / "test.pepmap.parquet"
        with PepMapWriter(out) as w:
            w.write_batch([record])
        table = pq.read_table(out)
        assert table.num_rows == 1


class TestPepMapDataset:
    def test_dataset_discovers_mapping(self, tmp_path):
        from qpx.writers import PepMapWriter

        records = [make_mapping_record()]
        out = tmp_path / "exp.pepmap.parquet"
        with PepMapWriter(out) as w:
            w.write_batch(records)

        import qpx

        ds = qpx.open(str(tmp_path))
        assert ds.pepmap is not None
        assert ds.pepmap.count() == 1
        ds.close()

    def test_filter_by_protein(self, tmp_path):
        from qpx.writers import PepMapWriter

        records = [
            make_mapping_record(sequence="PEP1K", protein_accession="P12345"),
            make_mapping_record(sequence="PEP2K", protein_accession="P67890"),
            make_mapping_record(sequence="PEP1K", protein_accession="P67890"),
        ]
        out = tmp_path / "exp.pepmap.parquet"
        with PepMapWriter(out) as w:
            w.write_batch(records)

        import qpx

        ds = qpx.open(str(tmp_path))
        filtered = ds.pepmap.by_protein("P12345")
        assert filtered.count() == 1
        ds.close()

    def test_unique_peptides_filter(self, tmp_path):
        from qpx.writers import PepMapWriter

        records = [
            make_mapping_record(sequence="PEP1K", is_unique=True),
            make_mapping_record(
                sequence="PEP2K", protein_accession="P67890", is_unique=False
            ),
        ]
        out = tmp_path / "exp.pepmap.parquet"
        with PepMapWriter(out) as w:
            w.write_batch(records)

        import qpx

        ds = qpx.open(str(tmp_path))
        unique = ds.pepmap.unique_peptides()
        assert unique.count() == 1
        ds.close()

    def test_save_structure_round_trip(self, dataset_dir):
        import qpx

        ds = qpx.open(str(dataset_dir))
        records = [
            make_mapping_record(),
            make_mapping_record(sequence="ANOTHERK", protein_accession="P67890"),
        ]
        path = ds.save_structure(records, "pepmap")
        assert path.exists()

        ds.refresh()
        assert ds.pepmap is not None
        assert ds.pepmap.count() == 2
        ds.close()
