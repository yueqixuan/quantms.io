"""Validate Spectronaut conversion fidelity — compare source TSV vs output parquet."""

import argparse
import sys

import duckdb
import numpy as np
import pyarrow.parquet as pq

sys.stdout.reconfigure(line_buffering=True)

_DEFAULT_REPORT = r"D:\Git-repository\Bigbio\QPX_data\test_Spectronaut\spectropiper_zenodo\extracted\HYE_Exploris480_SN19_Report_SpectroPipeR (Normal).tsv"
_DEFAULT_FEAT = r"D:\Git-repository\Bigbio\QPX_data\test_Spectronaut\spectropiper_zenodo\qpx_output_hye\hye.feature.parquet"
_DEFAULT_PG = r"D:\Git-repository\Bigbio\QPX_data\test_Spectronaut\spectropiper_zenodo\qpx_output_hye\hye.pg.parquet"


def main():
    parser = argparse.ArgumentParser(description="Validate Spectronaut conversion")
    parser.add_argument("--report", default=_DEFAULT_REPORT, help="Source Spectronaut TSV")
    parser.add_argument("--feature", default=_DEFAULT_FEAT, help="Output feature.parquet")
    parser.add_argument("--pg", default=_DEFAULT_PG, help="Output pg.parquet")
    args = parser.parse_args()
    REPORT = args.report
    FEAT_PQ = args.feature
    PG_PQ = args.pg
    print(f"Report: {REPORT}")
    print(f"Feature: {FEAT_PQ}")
    print(f"PG: {PG_PQ}")
    print()

    con = duckdb.connect()
    safe_path = REPORT.replace("'", "''")
    con.execute(
        f"CREATE VIEW report AS SELECT * FROM read_csv_auto('{safe_path}', "
        "delim='\\t', header=true, auto_detect=true, null_padding=true)"
    )

    feat = pq.read_table(FEAT_PQ)
    pg_tbl = pq.read_table(PG_PQ)

    print("=" * 60)
    print("1. ROW COUNT (precursor-level)")
    print("=" * 60)

    src_all = con.execute(
        "SELECT COUNT(*) FROM ("
        '  SELECT DISTINCT "R.FileName", "EG.ModifiedSequence", '
        '    "PEP.StrippedSequence", "FG.Charge" '
        '  FROM report WHERE "PG.ProteinGroups" IS NOT NULL'
        ")"
    ).fetchone()[0]

    src_qv = con.execute(
        "SELECT COUNT(*) FROM ("
        '  SELECT DISTINCT "R.FileName", "EG.ModifiedSequence", '
        '    "PEP.StrippedSequence", "FG.Charge" '
        '  FROM report WHERE "PG.ProteinGroups" IS NOT NULL '
        '    AND CAST("EG.Qvalue" AS DOUBLE) < 0.01'
        ")"
    ).fetchone()[0]

    out_rows = feat.num_rows
    print(f"  Source precursors (all):     {src_all:>12,}")
    print(f"  Source precursors (q<0.01):  {src_qv:>12,}")
    print(f"  Output feature rows:         {out_rows:>12,}")
    match = "YES" if out_rows == src_qv else f"NO (delta={out_rows - src_qv})"
    print(f"  Match (q<0.01): {match}")

    print()
    print("=" * 60)
    print("2. RUN COUNT")
    print("=" * 60)
    src_runs = con.execute('SELECT COUNT(DISTINCT "R.FileName") FROM report').fetchone()[0]
    out_runs = len(set(feat.column("run_file_name").to_pylist()))
    print(f"  Source runs:  {src_runs}")
    print(f"  Output runs:  {out_runs}")
    print(f"  Match: {'YES' if src_runs == out_runs else 'NO'}")

    print()
    print("=" * 60)
    print("3. PROTEIN GROUP COUNT")
    print("=" * 60)
    src_pg = con.execute('SELECT COUNT(DISTINCT "PG.ProteinGroups") FROM report WHERE "PG.ProteinGroups" IS NOT NULL').fetchone()[
        0
    ]
    out_pg_accs = set()
    for accs in pg_tbl.column("pg_accessions").to_pylist():
        if accs:
            out_pg_accs.add(";".join(accs))
    print(f"  Source unique PG:  {src_pg:>8,}")
    print(f"  Output unique PG:  {len(out_pg_accs):>8,}")
    delta_pg = src_pg - len(out_pg_accs)
    print(f"  Match: {'YES' if delta_pg == 0 else f'DELTA={delta_pg}'}")

    print()
    print("=" * 60)
    print("4. UNIQUE SEQUENCES")
    print("=" * 60)
    src_seqs = con.execute('SELECT COUNT(DISTINCT "PEP.StrippedSequence") FROM report').fetchone()[0]
    out_seqs = len(set(feat.column("sequence").to_pylist()))
    print(f"  Source:  {src_seqs:>8,}")
    print(f"  Output:  {out_seqs:>8,}")
    delta_seq = src_seqs - out_seqs
    print(f"  Match: {'YES' if delta_seq == 0 else f'DELTA={delta_seq}'}")

    print()
    print("=" * 60)
    print("5. INTENSITY SPOT CHECK (sample 1000 precursors)")
    print("=" * 60)

    # Compare FG.Quantity from source vs intensities[0].intensity in output
    src_sample = con.execute(
        'SELECT "R.FileName", "EG.ModifiedSequence", "PEP.StrippedSequence", '
        '  "FG.Charge", CAST("FG.Quantity" AS DOUBLE) AS qty '
        "FROM ("
        '  SELECT DISTINCT "R.FileName", "EG.ModifiedSequence", '
        '    "PEP.StrippedSequence", "FG.Charge", '
        '    FIRST("FG.Quantity") AS "FG.Quantity" '
        "  FROM report "
        '  WHERE "PG.ProteinGroups" IS NOT NULL '
        '    AND CAST("EG.Qvalue" AS DOUBLE) < 0.01 '
        '  GROUP BY "R.FileName", "EG.ModifiedSequence", '
        '    "PEP.StrippedSequence", "FG.Charge"'
        ") USING SAMPLE 1000"
    ).fetchdf()

    # Build lookup from output
    out_df = feat.select(["run_file_name", "peptidoform", "sequence", "charge", "intensities"]).to_pandas()
    out_df["raw_intensity"] = out_df["intensities"].apply(lambda x: x[0]["intensity"] if x else 0.0)

    matched = 0
    mismatched = 0
    rel_errors = []
    for _, row in src_sample.iterrows():
        run = row["R.FileName"]
        seq = row["PEP.StrippedSequence"]
        charge = int(row["FG.Charge"])
        src_qty = float(row["qty"]) if row["qty"] is not None else 0.0

        mask = (
            (out_df["sequence"] == seq)
            & (out_df["charge"] == charge)
            & (out_df["run_file_name"].str.contains(run.replace(".raw", "").replace(".mzML", "")))
        )
        matches = out_df.loc[mask]
        if len(matches) > 0:
            out_qty = float(matches.iloc[0]["raw_intensity"])
            if src_qty > 0:
                rel_err = abs(out_qty - src_qty) / src_qty
                rel_errors.append(rel_err)
            matched += 1
        else:
            mismatched += 1

    print(f"  Sampled precursors:  {len(src_sample)}")
    print(f"  Found in output:     {matched}")
    print(f"  Not found:           {mismatched}")
    if rel_errors:
        arr = np.array(rel_errors)
        print("  Relative error (intensity):")
        print(f"    median: {np.median(arr):.6f}")
        print(f"    mean:   {np.mean(arr):.6f}")
        print(f"    max:    {np.max(arr):.6f}")
        print(f"    exact (err<1e-6): {np.sum(arr < 1e-6)} / {len(arr)}")

    print()
    print("=" * 60)
    print("6. MODIFICATION COVERAGE")
    print("=" * 60)
    src_mod_count = con.execute(
        'SELECT COUNT(DISTINCT "EG.ModifiedSequence") FROM report WHERE "EG.ModifiedSequence" LIKE \'%[%\''
    ).fetchone()[0]
    src_unmod_count = con.execute(
        'SELECT COUNT(DISTINCT "EG.ModifiedSequence") FROM report WHERE "EG.ModifiedSequence" NOT LIKE \'%[%\''
    ).fetchone()[0]

    out_pf = feat.column("peptidoform").to_pylist()
    out_with_mod = sum(1 for p in set(out_pf) if "[" in str(p))
    out_without_mod = sum(1 for p in set(out_pf) if "[" not in str(p))

    print(f"  Source modified peptidoforms:    {src_mod_count:>8,}")
    print(f"  Output modified peptidoforms:    {out_with_mod:>8,}")
    print(f"  Source unmodified peptidoforms:  {src_unmod_count:>8,}")
    print(f"  Output unmodified peptidoforms:  {out_without_mod:>8,}")

    # Check unique mod types in source
    src_mods = con.execute(
        "SELECT DISTINCT regexp_extract_all(\"EG.ModifiedSequence\", '\\[([^\\]]+)\\]') AS mods "
        "FROM report WHERE \"EG.ModifiedSequence\" LIKE '%[%' LIMIT 5000"
    ).fetchdf()
    all_mods = set()
    for mlist in src_mods["mods"]:
        if len(mlist) > 0:
            for m in mlist:
                all_mods.add(m)
    print(f"  Source modification types: {sorted(all_mods)}")

    con.close()
    print("\n" + "=" * 60)
    print("VALIDATION COMPLETE")
    print("=" * 60)


if __name__ == "__main__":
    main()
