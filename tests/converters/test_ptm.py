"""Test shared PTM utilities (flat ptm.py)."""

from qpx.converters.ptm import (
    UNIMOD_MASS,
    mass_to_unimod,
    build_proforma,
    from_proforma,
    compute_precursor_mz,
    parse_unimod_delta,
    AA_MASS,
)


class TestUnimodMassRegistry:
    def test_oxidation_mass(self):
        assert abs(UNIMOD_MASS[35] - 15.994915) < 0.0001

    def test_acetyl_mass(self):
        assert abs(UNIMOD_MASS[1] - 42.010565) < 0.0001

    def test_mass_to_unimod_oxidation(self):
        assert mass_to_unimod(15.9949) == 35

    def test_mass_to_unimod_unknown(self):
        assert mass_to_unimod(999.0) is None


class TestBuildProforma:
    def test_no_mods(self):
        assert build_proforma("PEPTIDEK", []) == "PEPTIDEK"

    def test_internal_mod(self):
        assert build_proforma("PEPTMIDEK", [(5, "UNIMOD:35")]) == "PEPTM[UNIMOD:35]IDEK"

    def test_nterm_mod(self):
        assert build_proforma("PEPTIDEK", [(0, "UNIMOD:1")]) == "[UNIMOD:1]-PEPTIDEK"

    def test_empty_sequence(self):
        assert build_proforma("", []) == ""


class TestFromProforma:
    def test_unimod_tag(self):
        meta = {"UNIMOD:35": ("Oxidation", ["M"], ["Anywhere"])}
        result = from_proforma("M[UNIMOD:35]PEPTIDEK", "MPEPTIDEK", meta=meta)
        assert result is not None
        assert len(result) == 1
        assert result[0]["name"] == "Oxidation"
        assert result[0]["accession"] == "UNIMOD:35"

    def test_no_mods(self):
        assert from_proforma("PEPTIDEK", "PEPTIDEK", meta=None) is None

    def test_nterm(self):
        meta = {"UNIMOD:1": ("Acetyl", ["X"], ["N-term"])}
        result = from_proforma("[UNIMOD:1]-PEPTIDEK", "PEPTIDEK", meta=meta)
        assert result is not None
        assert result[0]["positions"][0]["position"] == 0

    def test_site_scores_attached(self):
        """Site scores are attached to the correct modification position."""
        site_scores = {
            1: [
                {
                    "score_name": "phosphors_site_probability",
                    "score_value": 0.95,
                    "higher_better": True,
                }
            ],
        }
        result = from_proforma(
            "M[UNIMOD:35]PEPTIDEK",
            "MPEPTIDEK",
            meta=None,
            site_scores=site_scores,
        )
        assert result is not None
        pos = result[0]["positions"][0]
        assert pos["position"] == 1
        assert pos["scores"] is not None
        assert pos["scores"][0]["score_value"] == 0.95

    def test_site_scores_none_when_not_provided(self):
        """Without site_scores, positions still have scores=None."""
        result = from_proforma("M[UNIMOD:35]PEPTIDEK", "MPEPTIDEK", meta=None)
        assert result is not None
        assert result[0]["positions"][0]["scores"] is None

    def test_site_scores_partial_positions(self):
        """Only positions in the site_scores dict get scores attached."""
        site_scores = {
            5: [
                {
                    "score_name": "phospho_prob",
                    "score_value": 0.8,
                    "higher_better": True,
                }
            ],
        }
        # Two mods: position 1 and position 5; only position 5 has scores
        result = from_proforma(
            "M[UNIMOD:35]PEPTS[UNIMOD:21]IDEK",
            "MPEPTSIDEK",
            meta=None,
            site_scores=site_scores,
        )
        assert result is not None
        # Find the phospho mod (position 5)
        for mod in result:
            for pos in mod["positions"]:
                if pos["position"] == 5:
                    assert pos["scores"] is not None
                    assert pos["scores"][0]["score_value"] == 0.8
                elif pos["position"] == 1:
                    assert pos["scores"] is None


class TestComputePrecursorMz:
    def test_basic(self):
        mz = compute_precursor_mz("PEPTM(UniMod:35)IDE", 2)
        assert mz is not None
        assert mz > 0

    def test_empty(self):
        assert compute_precursor_mz("", 2) is None
        assert compute_precursor_mz("PEPTIDE", 0) is None
