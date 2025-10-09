# quantms.io

[![Python application](https://github.com/bigbio/quantms.io/actions/workflows/python-app.yml/badge.svg?branch=dev)](https://github.com/bigbio/quantms.io/actions/workflows/python-app.yml)
[![Upload Python Package](https://github.com/bigbio/quantms.io/actions/workflows/python-publish.yml/badge.svg)](https://github.com/bigbio/quantms.io/actions/workflows/python-publish.yml)
[![Codacy Badge](https://app.codacy.com/project/badge/Grade/e71a662e8d4f483094576c1d8f8888c3)](https://app.codacy.com/gh/bigbio/quantms.io/dashboard?utm_source=gh&utm_medium=referral&utm_content=&utm_campaign=Badge_grade)
[![Codacy Badge](https://app.codacy.com/project/badge/Coverage/e71a662e8d4f483094576c1d8f8888c3)](https://app.codacy.com/gh/bigbio/quantms.io/dashboard?utm_source=gh&utm_medium=referral&utm_content=&utm_campaign=Badge_Coverage)
[![PyPI version](https://badge.fury.io/py/quantmsio.svg)](https://badge.fury.io/py/quantmsio)

A Python package for working with mass spectrometry data in the quantms.io format.

## Features

- Convert data from various mass spectrometry formats to quantms.io format
- Analyze and process quantms.io data
- Visualize results
- Manage project metadata
- Transform data between different formats

## Installation

### Install from PyPI

```bash
# To install the stable release from PyPI:
pip install quantmsio
```

### Install from Source (Without PyPI)

```bash
# Fork the repository on GitHub

# Clone the repository
git clone https://github.com/your-username/quantms.io.git
cd quantms.io

# Install the package locally
pip install .
```

## Usage

The package provides a command-line interface (CLI) with several command groups:

### Convert Commands

Convert data from various external formats to quantms.io:

```bash
# Convert quantms files
quantmsioc convert quantms-psm [OPTIONS]
quantmsioc convert quantms-feature [OPTIONS]
quantmsioc convert quantms-pg [OPTIONS]

# Convert MaxQuant files
quantmsioc convert maxquant-psm [OPTIONS]
quantmsioc convert maxquant-feature [OPTIONS]
quantmsioc convert maxquant-pg [OPTIONS]

# Convert FragPipe files
quantmsioc convert fragpipe [OPTIONS]

# Convert DIA-NN files
quantmsioc convert diann [OPTIONS]
quantmsioc convert diann-pg [OPTIONS]

# Convert IdXML to PSM parquet
quantmsioc convert idxml [OPTIONS]
```

### Analysis Commands

Analyze quantms.io data:

```bash
# PSM statistics
quantmsioc stats analyze psm [OPTIONS]

# Project AE (iBAQ) and PSM parquet statistics
quantmsioc stats analyze project-ae [OPTIONS]
```

### Visualization Commands

Visualize quantms.io data:

```bash
# Plot peptides by condition in LFQ
quantmsioc visualize plot psm-peptides [OPTIONS]

# Plot iBAQ distribution
quantmsioc visualize plot ibaq-distribution [OPTIONS]

# Plot KDE intensity distribution
quantmsioc visualize plot kde-intensity [OPTIONS]

# Plot peptide distribution
quantmsioc visualize plot peptide-distribution [OPTIONS]

# Plot intensity box plot
quantmsioc visualize plot box-intensity [OPTIONS]
```

### Project Management Commands

Manage project metadata:

```bash
# Generate project.json from a PRIDE accession and SDRF
quantmsioc project create [OPTIONS]

# Attach files to project.json
quantmsioc project attach [OPTIONS]
```

### Data Transformation Commands

Transform data within the quantms.io ecosystem:

```bash
# Generate iBAQ feature file
quantmsioc transform ibaq [OPTIONS]

# Convert IBAQ absolute file
quantmsioc transform ae [OPTIONS]

# Merge AE files into AnnData (.h5ad)
quantmsioc transform anndata [OPTIONS]

# Convert MSstats differential file
quantmsioc transform differential [OPTIONS]

# Map gene information from FASTA
quantmsioc transform gene [OPTIONS]

# Map spectrum information from mzML
quantmsioc transform spectra [OPTIONS]

# Update UniProt mappings
quantmsioc transform uniprot [OPTIONS]
```

## Configuration

Most commands support a `--verbose` flag that enables more detailed logging to stdout. The CLI uses standard logging configuration and does not require environment variables.

## Development

### Project Structure

```
quantmsio/
├── __init__.py
├── quantmsioc.py           # CLI entry point (poetry script: quantmsioc)
├── commands/               # CLI command groups
│   ├── convert/            # Converters: quantms, maxquant, diann, idxml, fragpipe
│   ├── transform/          # Transforms: ibaq, ae, gene, spectra, anndata, differential, uniprot
│   └── utils/              # Utility CLIs: project(create/attach), stats(analyze), plot
├── core/                   # Core logic & formats
│   ├── quantms/            # quantms feature/psm/pg, mztab helpers
│   ├── diann/, maxquant/, fragpipe/, idxml_utils/ ...
│   └── project.py, duckdb.py, format.py, common.py
├── operate/                # High-level operations (stats, plotting, tools)
│   ├── plots.py, query.py, statistics.py, tools.py
│   └── ...
└── utils/                  # Utilities
    ├── logger.py           # Basic logger getter
    ├── file_utils.py       # File helpers (e.g., AE file discovery)
    ├── pride_utils.py      # PRIDE archive helpers
    ├── mztab_utils.py      # mzTab helpers
    ├── system.py           # System utilities
    └── constants.py        # Constants and configurations
```

### Recent Improvements

1. **CLI Structure**

   - Commands are organized into `convert`, `transform`, `visualize`, `stats`, and `project`
   - Improved help messages and documentation

2. **Code Organization**
   - Separation of concerns across `core`, `commands`, `operate`, and `utils`
   - Modular design with clearer entry points

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

    This information is free; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This information is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this work; if not, write to the Free Software
    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
