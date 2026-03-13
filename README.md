# QPX

[![Python application](https://github.com/bigbio/qpx/actions/workflows/python-app.yml/badge.svg?branch=dev)](https://github.com/bigbio/qpx/actions/workflows/python-app.yml)
[![Upload Python Package](https://github.com/bigbio/qpx/actions/workflows/python-publish.yml/badge.svg)](https://github.com/bigbio/qpx/actions/workflows/python-publish.yml)
[![Codacy Badge](https://app.codacy.com/project/badge/Grade/e71a662e8d4f483094576c1d8f8888c3)](https://app.codacy.com/gh/bigbio/qpx/dashboard?utm_source=gh&utm_medium=referral&utm_content=&utm_campaign=Badge_grade)
[![Codacy Badge](https://app.codacy.com/project/badge/Coverage/e71a662e8d4f483094576c1d8f8888c3)](https://app.codacy.com/gh/bigbio/qpx/dashboard?utm_source=gh&utm_medium=referral&utm_content=&utm_campaign=Badge_Coverage)
[![PyPI version](https://badge.fury.io/py/qpx.svg)](https://badge.fury.io/py/qpx)

A Python package for working with mass spectrometry data in the QPX format.

## Features

- Convert data from various mass spectrometry formats to QPX format
- Analyze and process QPX data
- Visualize results
- Manage project metadata
- Transform data between different formats

## Installation

> **Note:** QPX is not yet available on PyPI. Please install directly from GitHub until the first official release.

### Install from GitHub (Recommended)

```bash
# Install the latest version directly from GitHub:
pip install git+https://github.com/bigbio/qpx.git
```

### Install from Source

```bash
# Clone the repository
git clone https://github.com/bigbio/qpx.git
cd qpx

# Install the package locally
pip install .
```

### Install and build with uv

[uv](https://docs.astral.sh/uv/) is a fast Python package installer and resolver. The project supports PEP 621 and can be installed, built, and published with uv.

**Prerequisites:** [Install uv](https://docs.astral.sh/uv/getting-started/installation/) (e.g. `curl -LsSf https://astral.sh/uv/install.sh | sh` or `pip install uv`).

```bash
# Install from GitHub
uv pip install "qpx @ git+https://github.com/bigbio/qpx.git"

# With optional extras (transforms, plotting)
uv pip install "qpx[transforms,plotting] @ git+https://github.com/bigbio/qpx.git"
```

**From a local clone:**

```bash
git clone https://github.com/bigbio/qpx.git
cd qpx

# Create a venv, install the project and its dependencies (recommended)
uv sync

# Or install in editable mode with optional dev dependencies
uv sync --extra dev

# Run the CLI without installing globally
uv run qpxc --help
```

**Build distributable packages** (sdist and wheel in `dist/`):

```bash
uv build
```

**Publish to PyPI** (after configuring credentials or trusted publishing):

```bash
uv build
uv publish
```

Both Poetry and uv can be used on this repo: the `pyproject.toml` includes a PEP 621 `[project]` section for uv/pip and `[tool.poetry]` for Poetry.

### Development Installation

For development with all dependencies:

```bash
# Using uv (recommended for fast installs)
uv sync --extra dev

# Using Poetry
poetry install

# Or using pip
pip install -e ".[dev]"
```

### System Dependencies

QPX depends on pyOpenMS, which requires certain system libraries. If you encounter errors related to missing shared libraries (e.g., `libglib-2.0.so.0`), install the required system dependencies:

**Ubuntu/Debian:**

```bash
sudo apt-get update
sudo apt-get install -y libglib2.0-0
```

**macOS:**

```bash
brew install glib
```

**Using Conda/Mamba (Recommended for pyOpenMS):**

Using mamba (faster dependency resolution):

```bash
mamba env create -f environment.yml
conda activate qpx
pip install git+https://github.com/bigbio/qpx.git
```

Or with conda:

```bash
conda env create -f environment.yml
conda activate qpx
pip install git+https://github.com/bigbio/qpx.git
```

## Usage

The package provides a command-line interface (CLI) with several command groups:

### Main CLI

```bash
Usage: cli [OPTIONS] COMMAND [ARGS]...

  qpx - A tool for converting and analyzing mass spectrometry proteomics
  data

Options:
  --version   Show the version and exit.
  -h, --help  Show this message and exit.

Commands:
  convert    Convert external formats to QPX format.
  project    Project management commands.
  stats      Statistical analysis of QPX data.
  transform  Transform QPX data into different representations.
  visualize  Visualize QPX data.
```

### Convert Commands

Convert data from various external formats to QPX:

```bash
Usage: convert [OPTIONS] COMMAND [ARGS]...

  Convert external formats to QPX format.

Options:
  --help  Show this message and exit.

Commands:
  diann             Convert DIA-NN report to QPX format
  diann-pg          Convert DIA-NN report to protein group format
  fragpipe          Convert FragPipe PSMs from psm.tsv to parquet file in
                    QPX
  idxml             Convert IdXML to PSM parquet file in QPX
  idxml-batch       Convert multiple IdXML files to a single merged PSM parquet
                    file
  maxquant-feature  Convert feature data from MaxQuant evidence.txt to parquet
                    format
  maxquant-pg       Convert MaxQuant proteinGroups.txt to QPX protein
                    group format
  maxquant-psm      Convert PSM data from MaxQuant msms.txt to parquet format
  quantms-feature   Convert feature data from mzTab to QPX format.
  quantms-pg        Convert protein groups from mzTab quantms TMT and LFQ...
  quantms-psm       Convert PSM data from mzTab to QPX format.
```

### Transform Commands

Transform data within the QPX ecosystem:

```bash
Usage: transform [OPTIONS] COMMAND [ARGS]...

  Transform QPX data into different representations.

Options:
  --help  Show this message and exit.

Commands:
  ae            Convert IBAQ absolute file into QPX format
  anndata       Merge multiple AE files into a file in AnnData format.
  differential  Convert a MSstats differential file into a QPX file
                format
  gene          Map gene information from FASTA to parquet format
  ibaq          Convert feature data to IBAQ format
  spectra       Map spectrum information from mzML to parquet format
  uniprot       Map feature data to latest UniProt version
```

### Visualization Commands

Visualize QPX data:

```bash
Usage: visualize [OPTIONS] COMMAND [ARGS]...

  Visualize QPX data.

Options:
  --help  Show this message and exit.

Commands:
  plot  Visualization commands for QPX data
```

### Statistics Commands

Analyze QPX data:

```bash
Usage: stats [OPTIONS] COMMAND [ARGS]...

  Statistical analysis of QPX data.

Options:
  --help  Show this message and exit.

Commands:
  analyze  Statistical analysis commands for QPX data
```

### Project Management Commands

Manage project metadata:

```bash
Usage: project [OPTIONS] COMMAND [ARGS]...

  Project management commands.

Options:
  --help  Show this message and exit.

Commands:
  attach  Register the file to project.json.
  create  Generate a project file from original PRIDE accession
```

## Configuration

Most commands support a `--verbose` flag that enables more detailed logging to stdout. The CLI uses standard logging configuration and does not require environment variables.

## Development

### Project Structure

```
qpx/
├── cli/                    # Click CLI (entry point: qpx.cli.main:main)
│   ├── main.py             # Top-level CLI group
│   └── convert.py          # convert subcommands (maxquant, diann, quantms, fragpipe, mzidentml, sdrf)
├── converters/             # Tool-specific converters
│   ├── quantms/            # QuantMS (mzTab) converter
│   ├── diann/              # DIA-NN converter
│   ├── maxquant/           # MaxQuant converter
│   ├── fragpipe/           # FragPipe converter
│   ├── mzidentml/          # mzIdentML converter
│   └── sdrf.py             # Shared SDRF converter
├── core/                   # Core logic & formats
│   ├── data/               # Schema definitions (YAML + Python)
│   │   └── schemas/        # YAML schema files for all structures
│   ├── engine.py           # DuckDB engine wrapper
│   ├── scores.py           # Score normalization & ontology
│   └── ontology/           # OBO ontology registry
├── writers/                # Parquet writers (one per structure)
├── views/                  # Analytical views (protein, peptide, QC)
└── dataset.py              # Main Dataset class entry point
```

### Contributing

1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Run tests
5. Submit a pull request

## License

This project is licensed under the Apache-2.0 License - see the LICENSE file for details.

## Core contributors and collaborators

The project is run by different groups:

- Yasset Perez-Riverol (PRIDE Team, European Bioinformatics Institute - EMBL-EBI, U.K.)
- Ping Zheng (Chongqing Key Laboratory of Big Data for Bio Intelligence, Chongqing University of Posts and Telecommunications, Chongqing, China)

IMPORTANT: If you contribute with the following specification, please make sure to add your name to the list of contributors.

## Code of Conduct

As part of our efforts toward delivering open and inclusive science, we follow the [Contributor Covenant Code of Conduct for Open Source Projects](https://www.contributor-covenant.org/version/2/0/code_of_conduct/).

## How to cite

## Copyright notice

    Copyright 2025 BigBio

    Licensed under the Apache License, Version 2.0.
    See the LICENSE file for details.
