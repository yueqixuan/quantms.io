"""FragPipe PG converter tests focused on grouped_runs expansion.

FragPipe ``combined_protein.tsv`` reports intensities per *experiment*. Each
experiment can aggregate several raw files, so ``grouped_runs`` must expand to
the experiment's member ``run_file_name`` values (via ``experiment_to_runs``);
absent that mapping it falls back to the single experiment token.
"""

from __future__ import annotations

import pyarrow.parquet as pq
import pytest

from qpx.converters.fragpipe.pg_adapter import FragPipePgAdapter

_COMBINED_PROTEIN_TSV = (
    "Protein\tGene\tDescription\t"
    "exp1 Total Intensity\texp1 MaxLFQ Intensity\texp1 Spectral Count\texp1 Unique Spectral Count\t"
    "exp2 Total Intensity\texp2 MaxLFQ Intensity\texp2 Spectral Count\texp2 Unique Spectral Count\t"
    "Combined Total Peptides\tCombined Unique Peptides\n"
    "P12345\tGENE1\tsome protein\t"
    "1000\t900\t10\t8\t"
    "2000\t1800\t20\t16\t"
    "5\t3\n"
)

_EXPERIMENT_ANNOTATION_TSV = "file\tsample\n/data/run_B.raw\texp1\n/data/run_A.mzML\texp1\nC:\\\\data\\\\run_C.RAW\texp2\n"


def _write_input(tmp_path):
    path = tmp_path / "combined_protein.tsv"
    path.write_text(_COMBINED_PROTEIN_TSV)
    return path


def _grouped_runs_by_intensity(table):
    """Return {total_intensity: grouped_runs} so records are identifiable per experiment."""
    df = table.to_pandas()
    out = {}
    for _, row in df.iterrows():
        total = float(row["intensities"][0]["intensity"])
        out[total] = list(row["grouped_runs"])
    return out


def test_grouped_runs_expand_to_member_run_files(tmp_path):
    """With an experiment->runs mapping, grouped_runs holds real run file names."""
    protein_path = _write_input(tmp_path)
    out_path = tmp_path / "fragpipe.pg.parquet"

    with FragPipePgAdapter() as adapter:
        adapter.convert(
            protein_path=str(protein_path),
            output_path=str(out_path),
            experiment_to_runs={"exp1": ["run_A", "run_B"], "exp2": ["run_C"]},
        )

    table = pq.read_table(str(out_path))
    grouped = _grouped_runs_by_intensity(table)
    assert grouped[1000.0] == ["run_A", "run_B"]  # exp1 expanded to its two fractions
    assert grouped[2000.0] == ["run_C"]


def test_grouped_runs_fallback_to_experiment_token(tmp_path):
    """Without a mapping, grouped_runs falls back to the single experiment token."""
    protein_path = _write_input(tmp_path)
    out_path = tmp_path / "fragpipe.pg.parquet"

    with FragPipePgAdapter() as adapter:
        adapter.convert(
            protein_path=str(protein_path),
            output_path=str(out_path),
        )

    table = pq.read_table(str(out_path))
    grouped = _grouped_runs_by_intensity(table)
    assert grouped[1000.0] == ["exp1"]
    assert grouped[2000.0] == ["exp2"]


def test_experiment_annotation_expands_and_normalizes_member_runs(tmp_path):
    """The official FragPipe annotation file drives grouped_runs directly."""
    protein_path = _write_input(tmp_path)
    annotation_path = tmp_path / "experiment_annotation.tsv"
    annotation_path.write_text(_EXPERIMENT_ANNOTATION_TSV)
    out_path = tmp_path / "fragpipe.pg.parquet"

    with FragPipePgAdapter() as adapter:
        adapter.convert(
            protein_path=str(protein_path),
            output_path=str(out_path),
            experiment_annotation_path=str(annotation_path),
        )

    grouped = _grouped_runs_by_intensity(pq.read_table(str(out_path)))
    assert grouped[1000.0] == ["run_A", "run_B"]
    assert grouped[2000.0] == ["run_C"]


def test_experiment_annotation_and_explicit_mapping_are_mutually_exclusive(
    tmp_path,
):
    """Ambiguous mapping inputs fail before writing output."""
    protein_path = _write_input(tmp_path)
    annotation_path = tmp_path / "experiment_annotation.tsv"
    annotation_path.write_text(_EXPERIMENT_ANNOTATION_TSV)

    with FragPipePgAdapter() as adapter:
        with pytest.raises(ValueError, match="either experiment_to_runs"):
            adapter.convert(
                protein_path=str(protein_path),
                output_path=str(tmp_path / "fragpipe.pg.parquet"),
                experiment_to_runs={"exp1": ["run_A"]},
                experiment_annotation_path=str(annotation_path),
            )
