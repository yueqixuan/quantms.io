# Installation

## Requirements

- **Python**: 3.10 or higher
- **Operating System**: Linux, macOS, or Windows
- **Memory**: Minimum 4GB RAM (8GB+ recommended for large datasets)

## Installation Methods

### Standard Installation (Recommended)

Install the latest stable version from PyPI:

```bash
pip install qpx
```

### Conda Environment

For isolated environment management using conda:

```bash
# Create a new environment with Python 3.10+
conda create -n qpx python=3.10
conda activate qpx

# Install qpx
pip install qpx
```

### Using environment.yml

Clone the repository and use the provided environment file:

```bash
git clone https://github.com/bigbio/qpx.git
cd qpx

# Create environment from file
conda env create -f environment.yml
conda activate qpx
```

### From Source

Install directly from the GitHub repository:

```bash
git clone https://github.com/bigbio/qpx.git
cd qpx
pip install .
```

### Development Installation

For contributing to the project:

```bash
# 1. Fork the repository on GitHub
# Visit https://github.com/bigbio/qpx and click "Fork"

# 2. Clone your fork
git clone https://github.com/YOUR-USERNAME/qpx
cd qpx

# 3. Install dependencies using Poetry (recommended)
pip install poetry
poetry install

# Or using pip with development dependencies
pip install -r requirements.txt
pip install -e ".[dev]"
```

## Verify Installation

After installation, verify everything is working:

```bash
# Check the installed version
qpxc --version

# View available commands
qpxc --help

# Test a simple command
qpxc convert --help
```

## Optional Dependencies

Some features require additional packages:

| Feature | Package | Installation |
|---------|---------|--------------|
| AnnData export | anndata | `pip install anndata` |
| Gene mapping | mygene | `pip install mygene` |
| Advanced plots | plotly | `pip install plotly` |

Install all optional dependencies:

```bash
pip install qpx[all]
```

## Upgrading

To upgrade to the latest version:

```bash
pip install --upgrade qpx
```

## Uninstalling

To remove qpx:

```bash
pip uninstall qpx
```

## Troubleshooting

If you encounter installation issues, see the [Troubleshooting Guide](troubleshooting.md) or:

- Check [GitHub Issues](https://github.com/bigbio/qpx/issues) for known problems
- Ensure Python 3.10+ is installed: `python --version`
- Try installing in a fresh virtual environment
