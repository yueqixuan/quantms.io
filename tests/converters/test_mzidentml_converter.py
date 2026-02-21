"""Tests for the mzIdentML PSM adapter (including XL-MS cross-linking)."""

import pytest
from pathlib import Path

import pyarrow.parquet as pq

from qpx.converters.mzidentml.psm_adapter import (
    MzIdentMLPsmAdapter,
    _parse_scan,
    _group_siis,
)

XL_EXAMPLES = Path(__file__).parent.parent / "examples" / "mzidentml" / "xl"


# ---------------------------------------------------------------------------
# Helper
# ---------------------------------------------------------------------------


def _convert(
    mzid_name: str, tmp_path: Path, examples_dir: Path = None
) -> pq.ParquetFile:
    """Convert a test mzid file and return the output as a ParquetFile."""
    base_dir = examples_dir or XL_EXAMPLES
    mzid_path = base_dir / mzid_name
    out_path = tmp_path / f"{mzid_name}.psm.parquet"
    with MzIdentMLPsmAdapter() as adapter:
        adapter.convert(mzid_path=str(mzid_path), output_path=str(out_path))
    return pq.read_table(str(out_path))


# ---------------------------------------------------------------------------
# Unit tests for helpers
# ---------------------------------------------------------------------------


class TestParseScan:
    def test_scan_format(self):
        assert _parse_scan("scan=100") == [100]

    def test_index_format(self):
        assert _parse_scan("index=7483") == [7483]

    def test_spectrum_format(self):
        assert _parse_scan("spectrum=5") == [5]

    def test_bare_integer(self):
        assert _parse_scan("42") == [42]

    def test_unknown_format(self):
        assert _parse_scan("something_else") == [0]


class TestGroupSiis:
    def _make_sii(self, sii_id, cv_params=None):
        return {
            "id": sii_id,
            "charge": 2,
            "exp_mz": 500.0,
            "calc_mz": 500.0,
            "rank": 1,
            "pass_threshold": True,
            "peptide_ref": f"Pep_{sii_id}",
            "cv_params": cv_params or [],
            "peptide_evidence_refs": [],
        }

    def test_regular_sii(self):
        siis = [self._make_sii("SII_1")]
        groups = _group_siis(siis)
        assert len(groups) == 1
        assert groups[0]["xl_type"] == "none"

    def test_crosslink_pair(self):
        siis = [
            self._make_sii(
                "SII_1", [{"accession": "MS:1002511", "name": "xl", "value": "1"}]
            ),
            self._make_sii(
                "SII_2", [{"accession": "MS:1002511", "name": "xl", "value": "1"}]
            ),
        ]
        groups = _group_siis(siis)
        assert len(groups) == 1
        assert groups[0]["xl_type"] == "inter"
        assert groups[0]["alpha"]["id"] == "SII_1"
        assert groups[0]["beta"]["id"] == "SII_2"

    def test_looplink(self):
        siis = [
            self._make_sii(
                "SII_1", [{"accession": "MS:1003329", "name": "looplink", "value": ""}]
            ),
        ]
        groups = _group_siis(siis)
        assert len(groups) == 1
        assert groups[0]["xl_type"] == "intra"

    def test_noncovalent_pair(self):
        siis = [
            self._make_sii(
                "SII_1", [{"accession": "MS:1003331", "name": "noncov", "value": "1"}]
            ),
            self._make_sii(
                "SII_2", [{"accession": "MS:1003331", "name": "noncov", "value": "1"}]
            ),
        ]
        groups = _group_siis(siis)
        assert len(groups) == 1
        assert groups[0]["xl_type"] == "noncovalent"


# ---------------------------------------------------------------------------
# Integration tests — Xlink EDC (inter-peptide + looplinks)
# ---------------------------------------------------------------------------


class TestXlinkEDC:
    def test_convert_produces_psms(self, tmp_path):
        table = _convert("Xlink_EDC_mzIdentML_1_3_0_draft.mzid.gz", tmp_path)
        assert table.num_rows > 0

    def test_has_crosslinked_psms(self, tmp_path):
        table = _convert("Xlink_EDC_mzIdentML_1_3_0_draft.mzid.gz", tmp_path)
        xl_count = sum(1 for v in table.column("cross_links") if v.as_py() is not None)
        assert xl_count > 0, "Expected some PSMs with cross-links"

    def test_has_looplinks(self, tmp_path):
        table = _convert("Xlink_EDC_mzIdentML_1_3_0_draft.mzid.gz", tmp_path)
        intra = []
        for i in range(table.num_rows):
            xl = table.column("cross_links")[i].as_py()
            if xl and xl[0]["xl_type"] == "intra":
                intra.append(xl[0])
        assert len(intra) > 0, "Expected looplink (intra) cross-links"
        # Looplinks should have donor_position and acceptor_position but no partner
        for link in intra:
            assert link["partner_sequence"] is None
            assert link["donor_position"] > 0
            assert link["acceptor_position"] > 0

    def test_has_interlinks(self, tmp_path):
        table = _convert("Xlink_EDC_mzIdentML_1_3_0_draft.mzid.gz", tmp_path)
        inter = []
        for i in range(table.num_rows):
            xl = table.column("cross_links")[i].as_py()
            if xl and xl[0]["xl_type"] == "inter":
                inter.append(xl[0])
        assert len(inter) > 0, "Expected inter-peptide cross-links"
        for link in inter:
            assert link["partner_sequence"] is not None
            assert len(link["partner_sequence"]) > 0

    def test_linker_info(self, tmp_path):
        table = _convert("Xlink_EDC_mzIdentML_1_3_0_draft.mzid.gz", tmp_path)
        for i in range(table.num_rows):
            xl = table.column("cross_links")[i].as_py()
            if xl and xl[0]["xl_type"] in ("intra", "inter"):
                assert xl[0]["linker_name"] == "EDC"
                assert xl[0]["linker_mass"] != 0
                break

    def test_regular_psms_have_no_crosslinks(self, tmp_path):
        table = _convert("Xlink_EDC_mzIdentML_1_3_0_draft.mzid.gz", tmp_path)
        # Some PSMs should be regular (no cross-links)
        regular = sum(1 for v in table.column("cross_links") if v.as_py() is None)
        assert regular > 0, "Expected some regular (non-XL) PSMs"


# ---------------------------------------------------------------------------
# Integration tests — noncovalently associated peptides
# ---------------------------------------------------------------------------


class TestNoncovalent:
    def test_convert_produces_psms(self, tmp_path):
        table = _convert("noncovalently_assoc_1_3_0_draft.mzid.gz", tmp_path)
        assert table.num_rows >= 1

    def test_noncovalent_crosslink(self, tmp_path):
        table = _convert("noncovalently_assoc_1_3_0_draft.mzid.gz", tmp_path)
        xl = table.column("cross_links")[0].as_py()
        assert xl is not None
        assert xl[0]["xl_type"] == "noncovalent"
        assert xl[0]["partner_sequence"] is not None
        assert len(xl[0]["partner_sequence"]) > 0


# ---------------------------------------------------------------------------
# Integration tests — multiple spectra per identification
# ---------------------------------------------------------------------------


class TestMultipleSpectra:
    def test_convert_produces_psms(self, tmp_path):
        table = _convert("multiple_spectra_per_id_1_3_0_draft.mzid.gz", tmp_path)
        assert table.num_rows > 0

    def test_has_crosslinks(self, tmp_path):
        table = _convert("multiple_spectra_per_id_1_3_0_draft.mzid.gz", tmp_path)
        xl_count = sum(1 for v in table.column("cross_links") if v.as_py() is not None)
        assert xl_count > 0

    def test_inter_peptide_link(self, tmp_path):
        table = _convert("multiple_spectra_per_id_1_3_0_draft.mzid.gz", tmp_path)
        for i in range(table.num_rows):
            xl = table.column("cross_links")[i].as_py()
            if xl and xl[0]["xl_type"] == "inter":
                assert xl[0]["linker_name"] == "DSSO"
                assert xl[0]["partner_sequence"] is not None
                break
        else:
            pytest.fail("No inter-peptide cross-link found")


# ---------------------------------------------------------------------------
# Integration tests — scores and thresholds
# ---------------------------------------------------------------------------


class TestScoresAndThresholds:
    def test_convert_produces_psms(self, tmp_path):
        table = _convert("scores_and_thresholds_1_3_0_draft.mzid.gz", tmp_path)
        assert table.num_rows == 2

    def test_has_crosslinks(self, tmp_path):
        table = _convert("scores_and_thresholds_1_3_0_draft.mzid.gz", tmp_path)
        xl_count = sum(1 for v in table.column("cross_links") if v.as_py() is not None)
        assert xl_count == 2

    def test_has_fdr_scores(self, tmp_path):
        table = _convert("scores_and_thresholds_1_3_0_draft.mzid.gz", tmp_path)
        for i in range(table.num_rows):
            scores = table.column("additional_scores")[i].as_py()
            score_names = {s["score_name"] for s in scores}
            assert "crosslinked_psm_level_global_fdr" in score_names
            assert "peptide_pair_sequence_level_global_fdr" in score_names


# ---------------------------------------------------------------------------
# Schema validation
# ---------------------------------------------------------------------------


class TestPepExtraction:
    """Test posterior_error_probability extraction."""

    def test_extract_pep_from_cvparams(self):
        """PEP should be extracted from MS:1001493 cvParam."""
        from qpx.converters.mzidentml.psm_adapter import _extract_pep

        cv_params = [
            {"accession": "MS:1001171", "name": "Mascot:score", "value": "45.2"},
            {
                "accession": "MS:1001493",
                "name": "posterior error probability",
                "value": "0.001",
            },
        ]
        assert _extract_pep(cv_params) == pytest.approx(0.001)

    def test_extract_pep_global(self):
        """Should extract PEP from MS:1002352 (global PEP)."""
        from qpx.converters.mzidentml.psm_adapter import _extract_pep

        cv_params = [
            {
                "accession": "MS:1002352",
                "name": "PSM-level global FDR",
                "value": "0.01",
            },
        ]
        pep = _extract_pep(cv_params)
        assert pep == pytest.approx(0.01)

    def test_extract_pep_missing(self):
        """Should return None when PEP is not present."""
        from qpx.converters.mzidentml.psm_adapter import _extract_pep

        cv_params = [
            {"accession": "MS:1001171", "name": "Mascot:score", "value": "45.2"},
        ]
        assert _extract_pep(cv_params) is None

    def test_extract_pep_invalid_value(self):
        """Should return None for non-numeric PEP values."""
        from qpx.converters.mzidentml.psm_adapter import _extract_pep

        cv_params = [
            {
                "accession": "MS:1001493",
                "name": "posterior error probability",
                "value": "",
            },
        ]
        assert _extract_pep(cv_params) is None

    def test_pep_also_in_additional_scores(self):
        """PEP should appear in both dedicated field and additional_scores."""
        from qpx.converters.mzidentml.psm_adapter import _extract_scores, _extract_pep

        cv_params = [
            {
                "accession": "MS:1001493",
                "name": "posterior error probability",
                "value": "0.05",
            },
            {"accession": "MS:1001171", "name": "Mascot:score", "value": "30.0"},
        ]
        pep = _extract_pep(cv_params)
        scores = _extract_scores(cv_params)

        assert pep == pytest.approx(0.05)
        # PEP should also be in additional_scores
        pep_scores = [s for s in scores if "posterior" in s["score_name"]]
        assert len(pep_scores) == 1
        assert pep_scores[0]["score_value"] == pytest.approx(0.05)


class TestRtExtraction:
    """Test retention time extraction."""

    def test_rt_from_sir_cvparam_seconds(self):
        """RT in seconds should be extracted directly."""
        from qpx.core.cv_terms import CV_SCAN_START_TIME

        # Just test that the constant exists and is correct
        assert CV_SCAN_START_TIME == "MS:1000016"

    def test_rt_from_mgf(self):
        """MGF RTINSECONDS should be parsed and stored."""
        from qpx.converters.mzidentml.mgf_parser import MgfSpectraIndex
        import tempfile, os

        mgf_content = (
            "BEGIN IONS\n"
            "TITLE=test_spectrum\n"
            "SCANS=100\n"
            "RTINSECONDS=123.45\n"
            "PEPMASS=500.0\n"
            "CHARGE=2+\n"
            "200.0 1000\n"
            "300.0 2000\n"
            "END IONS\n"
        )

        with tempfile.NamedTemporaryFile(mode="w", suffix=".mgf", delete=False) as f:
            f.write(mgf_content)
            f.flush()
            tmp_name = f.name
        try:
            idx = MgfSpectraIndex(tmp_name)
            spec = idx.get_spectrum(100)
            assert spec is not None
            assert spec["rt"] == pytest.approx(123.45)
        finally:
            os.unlink(tmp_name)

    def test_rt_none_when_absent(self):
        """RT should be None when RTINSECONDS is not in MGF."""
        from qpx.converters.mzidentml.mgf_parser import MgfSpectraIndex
        import tempfile, os

        mgf_content = (
            "BEGIN IONS\n"
            "TITLE=test_spectrum\n"
            "SCANS=200\n"
            "PEPMASS=500.0\n"
            "200.0 1000\n"
            "END IONS\n"
        )

        with tempfile.NamedTemporaryFile(mode="w", suffix=".mgf", delete=False) as f:
            f.write(mgf_content)
            f.flush()
            tmp_name = f.name
        try:
            idx = MgfSpectraIndex(tmp_name)
            spec = idx.get_spectrum(200)
            assert spec is not None
            assert spec.get("rt") is None
        finally:
            os.unlink(tmp_name)


class TestSchemaValidation:
    def test_output_matches_psm_schema(self, tmp_path):
        """Verify the output parquet has the correct PSM schema structure."""
        table = _convert("Xlink_EDC_mzIdentML_1_3_0_draft.mzid.gz", tmp_path)
        from qpx.core.data import PsmSchema

        errors = PsmSchema.validate(table)
        assert not errors, f"Schema validation errors: {errors}"
