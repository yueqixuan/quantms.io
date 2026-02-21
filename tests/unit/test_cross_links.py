"""Cross-linking schema and round-trip tests."""

import pytest
import pyarrow.parquet as pq

from tests.conftest import make_psm_record


class TestCrossLinkSchema:
    def test_psm_schema_includes_cross_links(self):
        from qpx.core.data import PsmSchema

        field_names = [f.name for f in PsmSchema.get_arrow_schema()]
        assert "cross_links" in field_names

    def test_cross_links_is_list_of_struct(self):
        import pyarrow as pa
        from qpx.core.data import PsmSchema

        field = PsmSchema.get_arrow_schema().field("cross_links")
        assert pa.types.is_list(field.type)
        assert pa.types.is_struct(field.type.value_type)

    def test_cross_link_struct_fields(self):
        from qpx.core.data import PsmSchema

        xl_type = PsmSchema.get_arrow_schema().field("cross_links").type.value_type
        xl_field_names = [xl_type.field(i).name for i in range(xl_type.num_fields)]
        assert "xl_type" in xl_field_names
        assert "partner_sequence" in xl_field_names
        assert "partner_peptidoform" in xl_field_names
        assert "donor_position" in xl_field_names
        assert "acceptor_position" in xl_field_names
        assert "linker_name" in xl_field_names
        assert "linker_accession" in xl_field_names
        assert "linker_mass" in xl_field_names


class TestCrossLinkRoundTrip:
    def test_write_and_read_xl_psm(self, tmp_path):
        from qpx.writers import PsmWriter

        record = make_psm_record()
        record["cross_links"] = [
            {
                "xl_type": "inter",
                "partner_sequence": "PEPTIDEK",
                "partner_peptidoform": "PEPTIDEK",
                "donor_position": 5,
                "acceptor_position": 3,
                "linker_name": "DSS",
                "linker_accession": "XLMOD:02001",
                "linker_mass": 138.0681,
            }
        ]
        out = tmp_path / "xl.psm.parquet"
        with PsmWriter(out) as w:
            w.write_batch([record])
        table = pq.read_table(out)
        assert table.num_rows == 1
        xl = table.column("cross_links")[0].as_py()
        assert len(xl) == 1
        assert xl[0]["xl_type"] == "inter"
        assert xl[0]["linker_mass"] == pytest.approx(138.0681)
        assert xl[0]["partner_sequence"] == "PEPTIDEK"
        assert xl[0]["donor_position"] == 5
        assert xl[0]["acceptor_position"] == 3

    def test_null_cross_links_for_regular_psm(self, tmp_path):
        from qpx.writers import PsmWriter

        record = make_psm_record()
        out = tmp_path / "regular.psm.parquet"
        with PsmWriter(out) as w:
            w.write_batch([record])
        table = pq.read_table(out)
        assert table.column("cross_links")[0].as_py() is None

    def test_dead_end_cross_link(self, tmp_path):
        from qpx.writers import PsmWriter

        record = make_psm_record()
        record["cross_links"] = [
            {
                "xl_type": "dead-end",
                "partner_sequence": None,
                "partner_peptidoform": None,
                "donor_position": 5,
                "acceptor_position": None,
                "linker_name": "DSS",
                "linker_accession": None,
                "linker_mass": 156.0786,
            }
        ]
        out = tmp_path / "de.psm.parquet"
        with PsmWriter(out) as w:
            w.write_batch([record])
        table = pq.read_table(out)
        xl = table.column("cross_links")[0].as_py()
        assert xl[0]["xl_type"] == "dead-end"
        assert xl[0]["partner_sequence"] is None
        assert xl[0]["acceptor_position"] is None

    def test_multiple_cross_links(self, tmp_path):
        from qpx.writers import PsmWriter

        record = make_psm_record()
        record["cross_links"] = [
            {
                "xl_type": "inter",
                "partner_sequence": "PEPTIDEK",
                "partner_peptidoform": "PEPTIDEK",
                "donor_position": 5,
                "acceptor_position": 3,
                "linker_name": "DSSO",
                "linker_accession": None,
                "linker_mass": 158.0038,
            },
            {
                "xl_type": "intra",
                "partner_sequence": None,
                "partner_peptidoform": None,
                "donor_position": 2,
                "acceptor_position": 7,
                "linker_name": "DSSO",
                "linker_accession": None,
                "linker_mass": 158.0038,
            },
        ]
        out = tmp_path / "multi_xl.psm.parquet"
        with PsmWriter(out) as w:
            w.write_batch([record])
        table = pq.read_table(out)
        xl = table.column("cross_links")[0].as_py()
        assert len(xl) == 2
        assert xl[0]["xl_type"] == "inter"
        assert xl[1]["xl_type"] == "intra"
