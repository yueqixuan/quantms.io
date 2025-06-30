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

```bash
pip install quantmsio
```

## Usage

The package provides a command-line interface (CLI) with several command groups:

### Convert Commands

Convert data from various external formats to quantms.io:

```bash
# Convert feature files
quantmsio convert features [OPTIONS]

# Convert PSM files
quantmsio convert psm [OPTIONS]

# Convert MaxQuant files
quantmsio convert maxquant-psm [OPTIONS]
quantmsio convert maxquant-feature [OPTIONS]
quantmsio convert maxquant-pg [OPTIONS]

# Convert FragPipe files
quantmsio convert fragpipe-psm [OPTIONS]

# Convert DiaNN files
quantmsio convert diann [OPTIONS]
quantmsio convert diann-pg [OPTIONS]

# Convert expression data
quantmsio convert differential [OPTIONS]
quantmsio convert absolute [OPTIONS]
```

### Analysis Commands

Analyze quantms.io data:

```bash
# Compare sets of PSMs
quantmsio analyze compare-psms [OPTIONS]

# Generate statistics
quantmsio analyze stats [OPTIONS]
```

### Visualization Commands

Visualize quantms.io data:

```bash
# Create plots
quantmsio visualize plot [OPTIONS]
```

### Project Management Commands

Manage project metadata:

```bash
# Generate PRIDE project JSON
quantmsio project pride [OPTIONS]

# Attach files to JSON
quantmsio project attach [OPTIONS]
```

### Data Transformation Commands

Transform data within the quantms.io ecosystem:

```bash
# Generate iBAQ view
quantmsio transform ibaq [OPTIONS]

# Convert spectrum data to Parquet
quantmsio transform spectrum-parquet [OPTIONS]

# Convert gene data to Parquet
quantmsio transform gene-parquet [OPTIONS]

# Update UniProt mappings
quantmsio transform update-uniprot [OPTIONS]

# Merge AE files
quantmsio transform merge-ae [OPTIONS]
```

## Configuration

The package can be configured using environment variables:

- `QUANTMSIO_LOG_LEVEL`: Set logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL)
- `QUANTMSIO_LOG_FILE`: Path to log file
- `QUANTMSIO_LOG_FORMAT`: Custom log format string
- `QUANTMSIO_LOG_DATE_FORMAT`: Custom date format for logs
- `QUANTMSIO_LOG_JSON`: Enable JSON-formatted logs if set to "true"

Example:

```bash
export QUANTMSIO_LOG_LEVEL=DEBUG
export QUANTMSIO_LOG_FILE=/path/to/logs/quantmsio.log
export QUANTMSIO_LOG_JSON=true
```

## Development

### Project Structure

```
quantmsio/
├── __init__.py
├── quantmsioc.py          # CLI entry point
├── commands/              # Command implementations
├── core/                  # Core functionality
├── operate/              # Operation-specific code
└── utils/                # Utility functions
    ├── logger.py         # Logging configuration
    ├── file_utils.py     # File handling utilities
    ├── pride_utils.py    # PRIDE-specific utilities
    └── constants.py      # Constants and configurations
```

### Recent Improvements

1. **Enhanced CLI Structure**
   - Organized commands into logical groups
   - Improved help messages and documentation
   - Better command naming and organization

2. **Improved Logging System**
   - Structured logging support
   - JSON log format option
   - Request tracking with unique IDs
   - Automatic log rotation
   - Environment-based configuration

3. **Code Organization**
   - Better separation of concerns
   - More modular design
   - Improved type hints
   - Enhanced documentation

### Contributing

1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Run tests
5. Submit a pull request

## License

This project is licensed under the MIT License - see the LICENSE file for details.

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

