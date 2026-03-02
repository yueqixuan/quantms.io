# CLI Reference

The qpx command-line tool provides a comprehensive set of commands for converting, transforming, querying, and validating mass spectrometry proteomics data.

## Overview

The qpx CLI is organized into five main command groups:

### [Convert Commands](convert.md)

Convert various mass spectrometry data formats to the QPX standard format:

- **QuantMS Conversion**: Convert QuantMS mzTab format to QPX data files
- **DIA-NN Conversion**: Convert DIA-NN reports to feature and protein group formats
- **MaxQuant Conversion**: Convert MaxQuant PSM, feature, and protein group data
- **FragPipe Conversion**: Convert FragPipe PSM data to QPX format
- **mzIdentML Conversion**: Convert mzIdentML format PSM data
- **SDRF Conversion**: Convert SDRF metadata to QPX sample/run format

[View detailed documentation →](convert.md)

### [Transform Commands](transform.md)

Transform and process data within the QPX ecosystem:

- **iBAQ Transformation**: Process iBAQ quantification files
- **Absolute Expression (AE)**: Convert iBAQ absolute expression data to QPX format ([format specification](../spec/absolute.md))
- **Differential Expression (DE)**: Convert MSstats differential expression analysis results
- **Gene Mapping**: Map gene information to protein data

[View detailed documentation →](transform.md)

### Query Commands *(coming soon)*

Query QPX datasets using SQL or structured commands:

- Execute SQL queries on QPX Parquet files
- Filter and aggregate data
- Join multiple datasets
- Export query results

### Info Commands *(coming soon)*

Inspect QPX datasets and metadata:

- Display dataset schema and statistics
- Show sample and run information
- List available columns and data types
- Validate data integrity

### Validate Commands *(coming soon)*

Validate QPX data against schemas:

- Check data conformance to QPX specifications
- Validate required fields and data types
- Verify referential integrity
- Generate validation reports

## Python API

Visualization, statistics, and project management functionality are available through the Python API:

```python
import qpx

# Load a dataset
dataset = qpx.Dataset("path/to/dataset")

# Generate visualizations
dataset.psm.plot.distribution()
dataset.feature.plot.intensity_heatmap()

# Compute statistics
stats = dataset.psm.stats()
summary = dataset.get_summary()

# Project management
project = qpx.Project.from_pride("PXD000001")
project.attach_files(["psm.parquet", "feature.parquet"])
```

See the API documentation for more details.

## Quick Start

### Installation

```bash
pip install qpx
```

### Basic Usage

View all available commands:

```bash
qpxc --help
```

View help for a specific command group:

```bash
qpxc convert --help
```

View detailed help for a specific command:

```bash
qpxc convert diann --help
```

## Common Options

Most commands support the following common options:

- `--verbose`: Enable verbose logging for debugging
- `--output-folder`: Specify the output directory
- `--output-prefix`: Specify the output file prefix

## Example Data Processing Workflow

A typical data processing workflow:

```bash
# 1. Convert raw data from MaxQuant
qpxc convert maxquant \
    --evidence-file evidence.txt \
    --msms-file msms.txt \
    --output-folder ./output

# 2. Transform to absolute expression data
qpxc transform ae \
    --ibaq-file ibaq.tsv \
    --sdrf-file metadata.sdrf.tsv \
    --output-folder ./output

# 3. Query the data
qpxc query \
    --dataset ./output \
    --sql "SELECT protein_accession, AVG(intensity) FROM psm GROUP BY protein_accession" \
    --output results.tsv

# 4. Validate the dataset
qpxc validate \
    --dataset ./output \
    --report validation_report.txt

# 5. Inspect dataset information
qpxc info \
    --dataset ./output \
    --show-schema \
    --show-stats
```

For visualization and statistical analysis, use the Python API:

```python
import qpx

# Load the dataset
dataset = qpx.Dataset("./output")

# Generate visualizations
dataset.psm.plot.distribution(save_path="./plots/psm_distribution.svg")
dataset.feature.plot.intensity_distribution(save_path="./plots/intensity_dist.svg")

# Generate statistics
stats = dataset.psm.stats()
print(stats.summary())

# Save statistical report
stats.save_report("./stats/report.txt")
```

## Getting Help

- Each command provides detailed help information using the `--help` parameter
- See [Format Specification](../spec/index.md) for output file formats and detailed schema information
- Check the [Troubleshooting Guide](../troubleshooting.md) for common issues
- Visit the [GitHub Repository](https://github.com/bigbio/qpx) to report issues
