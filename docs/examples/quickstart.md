# Quick Start Examples

## Example 1: Convert MaxQuant Data

Convert MaxQuant output files to QPX format.

```bash
# Convert PSM data only
qpxc convert maxquant \
    --msms-file tests/examples/maxquant/maxquant_simple/msms.txt \
    --output-folder ./output \
    --structures psm

# Convert all structures (PSM + feature + PG)
qpxc convert maxquant \
    --msms-file tests/examples/maxquant/maxquant_full/msms.txt.gz \
    --evidence-file tests/examples/maxquant/maxquant_full/evidence.txt.gz \
    --sdrf-file tests/examples/maxquant/maxquant_full/PXD001819.sdrf.tsv \
    --protein-groups-file tests/examples/maxquant/maxquant_full/proteinGroups.txt \
    --output-folder ./output
```

**Expected Output**:

- `output/psm-{uuid}.psm.parquet`
- `output/feature-{uuid}.feature.parquet`
- `output/pg-{uuid}.pg.parquet`

## Example 2: Convert DIA-NN Data

Process DIA-NN report files for Data-Independent Acquisition (DIA) data.

```bash
# Convert to feature format
qpxc convert diann \
    --report-path tests/examples/diann/small/diann_report.tsv \
    --qvalue-threshold 0.05 \
    --mzml-info-folder tests/examples/diann/small/mzml \
    --sdrf-path tests/examples/diann/small/PXD019909-DIA.sdrf.tsv \
    --output-folder ./output

# Convert features + protein groups
qpxc convert diann \
    --report-path tests/examples/diann/full/diann_report.tsv.gz \
    --pg-matrix-path tests/examples/diann/full/diann_report.pg_matrix.tsv \
    --sdrf-file tests/examples/diann/full/PXD036609.sdrf.tsv \
    --output-folder ./output
```

**Expected Output**:

- `output/feature-{uuid}.feature.parquet`
- `output/pg-{uuid}.pg.parquet`

## Example 3: Absolute Expression Analysis

Calculate absolute protein expression from iBAQ data.

```bash
qpxc transform ae \
    --ibaq-file tests/examples/AE/PXD016999.1-ibaq.tsv \
    --sdrf-file tests/examples/AE/PXD016999-first-instrument.sdrf.tsv \
    --project-file tests/examples/AE/project.json \  # v1.0; will be dataset.parquet in v2.0
    --output-folder ./output \
    --output-prefix ae_analysis
```

**Expected Output**:

- `output/ae_analysis-{uuid}.absolute.parquet`

---

[← Back to Examples Overview](index.md) | [View Workflows →](workflows.md)
