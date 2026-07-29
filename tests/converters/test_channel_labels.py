"""Plex-aware channel-label resolution for the QuantMS feature adapter.

Regression coverage for the qpx intensity-label bug: TMT/iTRAQ channels were
labeled from a plex-agnostic table (TMT10 channel 10 wrongly labeled TMT131N)
and iTRAQ was not mapped at all. Labels now come from the shared sdrf-pipelines
channel_map and are plex-aware.
"""

import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from qpx.converters.channel_labels import (
    experiment_type_from_labels,
    relabel_intensities_parquet,
    resolve_channel_labels,
)


def _resolve(experiment_type, channel_indices, sdrf_labels=None):
    return resolve_channel_labels(experiment_type, sdrf_labels, channel_indices)


def test_tmt10_channel10_is_tmt131_not_tmt131n():
    labels = _resolve("TMT", list(range(1, 11)))  # 10 channels -> tmt10plex
    assert labels[1] == "TMT126"
    assert labels[10] == "TMT131"  # NOT TMT131N (the historical bug)
    assert "TMT131N" not in labels.values()


def test_tmt11_uses_131n_131c():
    labels = _resolve("TMT", list(range(1, 12)))  # 11 channels -> tmt11plex
    assert labels[10] == "TMT131N"
    assert labels[11] == "TMT131C"


def test_experiment_type_from_labels():
    assert experiment_type_from_labels({"TMT126", "TMT131C"}) == "TMT"
    assert experiment_type_from_labels({"ITRAQ114"}) == "iTRAQ"
    assert experiment_type_from_labels({"label free sample"}) == "LFQ"
    assert experiment_type_from_labels(None) == "LFQ"


def _write_intensities_parquet(path, rows):
    itype = pa.list_(pa.struct([("label", pa.string()), ("intensity", pa.float32())]))
    pq.write_table(pa.table({"intensities": pa.array(rows, type=itype)}), path)


def test_relabel_openms_feature_filename_labels(tmp_path):
    """OpenMS feature -out_qpx: filename in every label -> reporter names by position."""
    src = tmp_path / "in.parquet"
    dst = tmp_path / "out.parquet"
    fname = "run01.mzML"
    _write_intensities_parquet(src, [[{"label": fname, "intensity": float(i)} for i in range(1, 11)]])
    labels = resolve_channel_labels("TMT", {f"TMT12{i}" for i in range(6)} | {"TMT131"}, range(1, 11))
    relabel_intensities_parquet(str(src), str(dst), labels, is_lfq=False)
    got = [e["label"] for e in pq.read_table(dst).column("intensities").to_pylist()[0]]
    assert got[0] == "TMT126"
    assert got[9] == "TMT131"  # tmt10plex, position 10
    assert fname not in got


def test_relabel_openms_pg_index_labels(tmp_path):
    """OpenMS pg -out_qpx: bare index labels '1','2',... -> reporter names."""
    src = tmp_path / "pg_in.parquet"
    dst = tmp_path / "pg_out.parquet"
    _write_intensities_parquet(src, [[{"label": str(i), "intensity": 1.0} for i in range(1, 12)]])
    labels = resolve_channel_labels("TMT", None, range(1, 12))  # tmt11plex by count
    relabel_intensities_parquet(str(src), str(dst), labels, is_lfq=False)
    got = [e["label"] for e in pq.read_table(dst).column("intensities").to_pylist()[0]]
    assert got == ["TMT126", "TMT127N", "TMT127C", "TMT128N", "TMT128C", "TMT129N",
                   "TMT129C", "TMT130N", "TMT130C", "TMT131N", "TMT131C"]


def test_run_label_is_canonical_join_key(tmp_path):
    """run.samples[].label must be the SAME canonical label that feature/pg carry:
    the SDRF ontology form 'label free sample' is normalized to 'LFQ'."""
    import pyarrow.parquet as pq

    from qpx.converters.sdrf import SdrfConverter

    sdrf = "tests/examples/quantms/dda-lfq-full/PXD007683-LFQ.sdrf.tsv"
    with SdrfConverter() as conv:
        conv.convert(
            sdrf_path=sdrf,
            sample_output=str(tmp_path / "s.parquet"),
            run_output=str(tmp_path / "r.parquet"),
        )
    run = pq.read_table(tmp_path / "r.parquet")
    labels = {s["label"] for row in run.column("samples").to_pylist() for s in (row or [])}
    assert labels == {"LFQ"}


def test_relabel_lfq_collapses_to_lfq(tmp_path):
    src = tmp_path / "lfq_in.parquet"
    dst = tmp_path / "lfq_out.parquet"
    _write_intensities_parquet(src, [[{"label": "run01.mzML", "intensity": 5.0}]])
    relabel_intensities_parquet(str(src), str(dst), {}, is_lfq=True)
    assert pq.read_table(dst).column("intensities").to_pylist()[0][0]["label"] == "LFQ"


@pytest.mark.parametrize(
    ("exp_type", "n", "last_label"),
    [
        ("TMT", 6, "TMT131"),
        ("TMT", 16, "TMT134N"),
        ("TMT", 18, "TMT135N"),
        ("iTRAQ", 4, "ITRAQ117"),
        ("iTRAQ", 8, "ITRAQ121"),
    ],
)
def test_plex_resolution_by_channel_count(exp_type, n, last_label):
    labels = _resolve(exp_type, list(range(1, n + 1)))
    assert len(labels) == n
    assert labels[n] == last_label


def test_itraq_is_now_mapped():
    labels = _resolve("iTRAQ", [1, 2, 3, 4])
    assert labels == {1: "ITRAQ114", 2: "ITRAQ115", 3: "ITRAQ116", 4: "ITRAQ117"}


def test_non_isobaric_returns_empty():
    assert _resolve("LFQ", [1]) == {}
    assert _resolve("TMT", []) == {}


def test_sdrf_labels_override_channel_count():
    """SDRF ground truth wins: a TMT11 SDRF with an empty 11th channel in the
    data still resolves to tmt11plex (ch10 = TMT131N), where count-inference
    from indices 1..10 alone would wrongly pick tmt10plex (ch10 = TMT131)."""
    tmt11_labels = {
        "TMT126", "TMT127N", "TMT127C", "TMT128N", "TMT128C", "TMT129N",
        "TMT129C", "TMT130N", "TMT130C", "TMT131N", "TMT131C",
    }
    labels = _resolve("TMT", list(range(1, 11)), sdrf_labels=tmt11_labels)
    assert labels[10] == "TMT131N"  # SDRF-declared tmt11plex, not tmt10plex
    # without the SDRF, count-inference falls back to tmt10plex
    assert _resolve("TMT", list(range(1, 11)))[10] == "TMT131"


def test_unresolvable_sdrf_labels_fall_back_to_count():
    labels = _resolve("TMT", list(range(1, 11)), sdrf_labels={"garbage", "nonsense"})
    assert labels[10] == "TMT131"  # fell back to count-based tmt10plex


def test_tmt_fixture_labels_are_reporter_names(tmp_path):
    """End-to-end: the TMT11 example converts to TMT126..TMT131C reporter labels."""
    from pathlib import Path

    import pyarrow.parquet as pq

    from qpx.converters.quantms.converter import QuantMSConverter

    base = Path("tests/examples/quantms/dda-plex-full")
    mztab = base / "PXD007683TMT.sdrf_openms_design_openms.mzTab.gz"
    msstats = base / "PXD007683TMT.sdrf_openms_design_msstats_in.csv.gz"
    sdrf = base / "PXD007683-TMT.sdrf.tsv"

    out = tmp_path / "out"
    out.mkdir()
    QuantMSConverter(
        mztab_path=str(mztab), sdrf_file=str(sdrf), msstats_file=str(msstats)
    ).convert(output_folder=out, output_prefix="tmt", structures=["feature"])
    feature_files = list(out.glob("*.feature.parquet"))
    assert feature_files, "no feature parquet produced"

    labels = set()
    for entry in pq.read_table(feature_files[0]).column("intensities").to_pylist():
        for ch in entry or []:
            labels.add(ch["label"])
    # Reporter names, never bare numeric indices.
    assert "TMT126" in labels
    assert labels <= {
        "TMT126", "TMT127N", "TMT127C", "TMT128N", "TMT128C",
        "TMT129N", "TMT129C", "TMT130N", "TMT130C", "TMT131N", "TMT131C",
    }
    assert not any(lbl.isdigit() for lbl in labels)
