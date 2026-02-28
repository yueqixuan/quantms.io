# Workflow Examples

Complete processing workflows for different analysis scenarios.

## Complete Label-Free Quantification (LFQ) Workflow

Process a complete LFQ experiment from raw MaxQuant output to publication-ready results.

```bash
#!/bin/bash

# Define paths
RAW_DIR="data/raw"
OUTPUT_DIR="output"
PLOTS_DIR="plots"
REPORTS_DIR="reports"
PROJECT_DIR="project"

# Create directories
mkdir -p $OUTPUT_DIR $PLOTS_DIR $REPORTS_DIR $PROJECT_DIR

# Step 1: Convert MaxQuant data (all structures in one command)
echo "Converting MaxQuant data..."
qpxc convert maxquant \
    --msms-file $RAW_DIR/msms.txt \
    --evidence-file $RAW_DIR/evidence.txt.gz \
    --sdrf-file $RAW_DIR/experiment.sdrf.tsv \
    --protein-groups-file $RAW_DIR/proteinGroups.txt \
    --output-folder $OUTPUT_DIR \
    --verbose

# Step 2: Calculate absolute expression
echo "Calculating absolute expression..."
qpxc transform ae \
    --ibaq-file $RAW_DIR/ibaq.tsv \
    --sdrf-file $RAW_DIR/experiment.sdrf.tsv \
    --output-folder $OUTPUT_DIR

# Step 3: Validate the QPX output
echo "Validating dataset..."
qpxc validate --dataset-folder $OUTPUT_DIR

# Step 4: Inspect dataset info
echo "Dataset info..."
qpxc info --dataset-folder $OUTPUT_DIR

echo "Workflow complete!"
echo "Results:"
echo "  - Converted data: $OUTPUT_DIR/"
```

## Differential Expression Analysis Workflow

Analyze differential expression from MSstats output.

```bash
#!/bin/bash

# Convert differential expression data
qpxc transform differential \
    --msstats-file tests/examples/DE/PXD033169.sdrf_openms_design_msstats_in_comparisons.csv \
    --sdrf-file tests/examples/DE/PXD033169.sdrf.tsv \
    --project-file tests/examples/DE/project.json \  # v1.0; dataset.parquet in v2.0
    --fdr-threshold 0.05 \
    --output-folder ./output \
    --verbose

echo "Differential expression analysis complete!"
```

---

[← Back to Examples Overview](index.md) | [View Integration Examples →](integration.md)
