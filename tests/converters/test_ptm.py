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


class TestComputePrecursorMz:
    def test_basic(self):
        mz = compute_precursor_mz("PEPTM(UniMod:35)IDE", 2)
        assert mz is not None
        assert mz > 0

    def test_empty(self):
        assert compute_precursor_mz("", 2) is None
        assert compute_precursor_mz("PEPTIDE", 0) is None
