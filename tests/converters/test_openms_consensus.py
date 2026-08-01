"""consensusXML -> QPX converter (interim path).

Builds a tiny ConsensusMap with pyopenms so the extraction runs in CI without a
large fixture file.
"""

import duckdb
import pytest

pytest.importorskip("pyopenms")
import pyopenms as oms  # noqa: E402

from qpx.converters.openms_consensus.converter import OpenMSConsensusConverter  # noqa: E402
from qpx.converters.openms_consensus.feature_adapter import to_proforma  # noqa: E402


@pytest.mark.parametrize(
    ("openms_seq", "expected"),
    [
        ("PEPTIDEK", "PEPTIDEK"),
        (".(TMT6plex)THSQEEM(Oxidation)QHMQR", "[UNIMOD:737]-THSQEEM[UNIMOD:35]QHMQR"),
        ("C(Carbamidomethyl)PEPTIDEK", "C[UNIMOD:4]PEPTIDEK"),
        ("PEPTIDER.(Amidated)", "PEPTIDER-[UNIMOD:2]"),
    ],
)
def test_to_proforma(openms_seq, expected):
    assert to_proforma(oms.AASequence.fromString(openms_seq)) == expected


def _write_tmt_consensusxml(path):
    """A 2-channel isobaric consensusXML: one peptide, 2 TMT channels, 1 protein."""
    cm = oms.ConsensusMap()
    cm.setExperimentType("labeled_MS2")
    headers = {}
    for idx, label in ((0, "tmt6plex_126"), (1, "tmt6plex_127")):
        h = oms.ColumnHeader()
        h.filename = "run_01.mzML"
        h.label = label
        h.size = 1
        h.unique_id = idx + 1
        headers[idx] = h
    cm.setColumnHeaders(headers)

    cf = oms.ConsensusFeature()
    cf.setCharge(2)
    cf.setRT(100.0)
    cf.setMZ(450.25)
    for idx, inten in ((0, 1000.0), (1, 2000.0)):
        peak = oms.Peak2D()
        peak.setRT(100.0)
        peak.setMZ(450.25)
        peak.setIntensity(inten)
        cf.insert(idx, peak, idx)

    pid = oms.PeptideIdentification()
    pid.setMZ(450.26)
    pid.setRT(100.0)
    pid.setMetaValue("spectrum_reference", "controllerType=0 controllerNumber=1 scan=42")
    hit = oms.PeptideHit()
    hit.setSequence(oms.AASequence.fromString("PEPTIDEK"))
    hit.setCharge(2)
    hit.setMetaValue("target_decoy", "target")
    hit.setMetaValue("Posterior Error Probability_score", 0.001)
    ev = oms.PeptideEvidence()
    ev.setProteinAccession("P12345")
    hit.setPeptideEvidences([ev])
    pid.setHits([hit])
    cf.setPeptideIdentifications([pid])
    cm.push_back(cf)

    prot = oms.ProteinIdentification()
    ph = oms.ProteinHit()
    ph.setAccession("P12345")
    prot.setHits([ph])
    cm.setProteinIdentifications([prot])

    oms.ConsensusXMLFile().store(str(path), cm)


def test_consensusxml_to_qpx_feature_has_channels_pg_is_identification_only(tmp_path):
    cx = tmp_path / "test.consensusXML"
    _write_tmt_consensusxml(cx)
    out = tmp_path / "out"
    written = OpenMSConsensusConverter().convert(str(cx), str(out), output_prefix="t", structures=("feature", "psm", "pg"))
    con = duckdb.connect()

    feat = con.execute(
        f"SELECT peptidoform, charge, run_file_name, intensities FROM read_parquet('{written['feature']}')"
    ).fetchall()
    assert len(feat) == 1
    pep, charge, run, intensities = feat[0]
    assert pep == "PEPTIDEK" and charge == 2 and run == "run_01"
    labels = {e["label"]: e["intensity"] for e in intensities}
    assert labels == {"TMT126": 1000.0, "TMT127": 2000.0}  # both channels, canonicalized, quant kept

    psm = con.execute(f"SELECT peptidoform, scan FROM read_parquet('{written['psm']}')").fetchall()
    assert psm and psm[0][0] == "PEPTIDEK" and list(psm[0][1]) == [42]

    pg = con.execute(f"SELECT anchor_protein, label, intensity FROM read_parquet('{written['pg']}')").fetchall()
    assert pg and pg[0][0] == "P12345"
    # interim: pg is identification-only — no protein intensity, no label
    assert all(label is None and intensity is None for _, label, intensity in pg)
