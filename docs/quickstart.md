# Quick Start

Get started with qpx in minutes - from installation to your first data conversion.

## Quick Start Flow

![Quick Start Flow](images/workflow2.png)

## Prerequisites

Before installing qpx, ensure you have:

- **Python 3.10 or higher** - Check with `python --version`
- **pip** - Python package manager (included with Python)
- **Optional**: conda/mamba for environment management

## Installation

=== "pip (Recommended)"

    ```bash
    pip install qpx
    ```

=== "conda"

    ```bash
    # Create a new environment
    conda create -n qpx python=3.10
    conda activate qpx

    # Install qpx
    pip install qpx
    ```

=== "From Source"

    ```bash
    git clone https://github.com/bigbio/qpx.git
    cd qpx
    pip install .
    ```

## Verify Installation

After installation, verify qpx is working correctly:

```bash
# Check version
qpxc --version

# View available commands
qpxc --help
```

You should see output similar to:

```
Usage: qpxc [OPTIONS] COMMAND [ARGS]...

  qpx command line interface for proteomics data processing.

Options:
  --version  Show the version and exit.
  --help     Show this message and exit.

Commands:
  convert    Convert proteomics data formats to QPX format
  project    Project management commands
  stats      Statistical analysis commands
  transform  Data transformation commands
  visualize  Visualization commands
```

## Your First Conversion

Let's convert some sample MaxQuant data to QPX format.

### Step 1: Download Sample Data

```bash
# Create a working directory
mkdir qpx-tutorial && cd qpx-tutorial

# Download sample MaxQuant msms.txt file
curl -L -o msms.txt \
  "https://raw.githubusercontent.com/bigbio/qpx/main/tests/examples/maxquant/maxquant_simple/msms.txt"
```

### Step 2: Convert to QPX Format

```bash
# Convert MaxQuant PSM data to QPX parquet format
qpxc convert maxquant-psm \
    --msms-file msms.txt \
    --output-folder ./output \
    --verbose
```

### Step 3: Verify the Output

```bash
# List the generated files
ls -la output/

# You should see a file like: psm-{uuid}.psm.parquet
```

### Step 4: Inspect the Data (Optional)

```python
# Using Python to read the parquet file
import pyarrow.parquet as pq

table = pq.read_table("output/psm-*.psm.parquet")
df = table.to_pandas()

print(f"Total PSMs: {len(df)}")
print(f"Columns: {list(df.columns)}")
print(df.head())
```

## What's Next?

Now that you've completed your first conversion, explore more:

| Next Step | Description |
|-----------|-------------|
| [Examples Overview](examples-overview.md) | More conversion and analysis examples |
| [Convert Commands](cli-convert.md) | All available data converters |
| [Transform Commands](cli-transform.md) | Data transformation tools |
| [Format Specification](format-specification.md) | Understanding QPX data formats |

## Common Commands

```bash
# Convert DIA-NN data
qpxc convert diann --report-path report.tsv --output-folder ./output

# Convert FragPipe data
qpxc convert fragpipe --psm-file psm.tsv --output-folder ./output

# Generate statistics
qpxc stats analyze psm --parquet-path ./output/psm.parquet

# Create visualizations
qpxc visualize plot ibaq-distribution --ibaq-path ./output/ae.parquet
```

## Need Help?

- Run `qpxc <command> --help` for detailed command help
- Check the [Troubleshooting](troubleshooting.md) guide
- Visit our [GitHub Issues](https://github.com/bigbio/qpx/issues) for support
