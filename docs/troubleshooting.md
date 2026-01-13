# Troubleshooting

Common issues and solutions when using qpx.

## Installation Issues

### ModuleNotFoundError: No module named 'qpx'

**Problem**: Python cannot find the qpx package after installation.

**Solutions**:

1. Ensure you're using the correct Python environment:

   ```bash
   which python
   which qpxc
   ```

2. Reinstall in the active environment:

   ```bash
   pip install --force-reinstall qpx
   ```

3. If using conda, ensure the environment is activated:

   ```bash
   conda activate qpx
   ```

### Python version incompatibility

**Problem**: Installation fails with Python version errors.

**Solution**: qpx requires Python 3.10 or higher. Check your version:

```bash
python --version
```

If needed, install a compatible Python version:

```bash
# Using conda
conda create -n qpx python=3.10
conda activate qpx
pip install qpx
```

### Missing dependencies

**Problem**: Import errors for packages like `venn`, `pyopenms`, or `anndata`.

**Solution**: Install the missing package:

```bash
pip install venn pyopenms anndata
```

Or install all optional dependencies:

```bash
pip install qpx[all]
```

## Conversion Issues

### File not found errors

**Problem**: `FileNotFoundError` when running convert commands.

**Solutions**:

1. Use absolute paths:

   ```bash
   qpxc convert maxquant-psm \
       --msms-file /full/path/to/msms.txt \
       --output-folder /full/path/to/output
   ```

2. Verify the file exists:
   ```bash
   ls -la path/to/file.txt
   ```

### Memory errors with large files

**Problem**: `MemoryError` or system becomes unresponsive with large datasets.

**Solutions**:

1. Process files in batches if supported by the command

2. Increase available memory or use a machine with more RAM

3. For DIA-NN reports, use the `--qvalue-threshold` to filter data:
   ```bash
   qpxc convert diann \
       --report-path report.tsv \
       --qvalue-threshold 0.01 \
       --output-folder ./output
   ```

### Invalid file format errors

**Problem**: `ValueError` or parsing errors when reading input files.

**Solutions**:

1. Verify the file format matches the expected format for the converter

2. Check for file corruption:

   ```bash
   head -20 input_file.txt
   ```

3. Ensure the file encoding is UTF-8:
   ```bash
   file input_file.txt
   ```

## Output Issues

### Empty output files

**Problem**: Parquet files are created but contain no data.

**Solutions**:

1. Check if input data passes quality filters (q-value, PEP thresholds)

2. Verify column names match expected format for the software

3. Use `--verbose` flag to see processing details:
   ```bash
   qpxc convert maxquant-psm --msms-file msms.txt --output-folder ./ --verbose
   ```

### Missing columns in output

**Problem**: Expected columns are not present in the output Parquet file.

**Solutions**:

1. Check if the input file contains the required source columns

2. For spectral data, ensure `--spectral-data` flag is used:

   ```bash
   qpxc convert maxquant-psm --msms-file msms.txt --output-folder ./ --spectral-data
   ```

3. Review the [Format Specification](format-specification.md) for required vs optional fields

## SDRF Issues

### Sample name mismatches

**Problem**: Samples in data files don't match SDRF sample names.

**Solutions**:

1. Ensure `source name` column in SDRF matches file names (without extension)

2. Check for whitespace or case sensitivity issues:

   ```python
   import pandas as pd
   sdrf = pd.read_csv('experiment.sdrf.tsv', sep='\t')
   print(sdrf['source name'].unique())
   ```

### Missing factor values

**Problem**: Factor values are not extracted from SDRF.

**Solution**: Ensure factor columns follow the format `factor value[factor_name]`:

```text
source name    factor value[disease]    factor value[organism part]
sample1        healthy                  liver
sample2        cancer                   liver
```

## Performance Issues

### Slow processing

**Solutions**:

1. Use SSD storage for input/output files

2. Increase available RAM

3. For large datasets, consider processing samples in parallel

4. Use compressed input files (.gz) to reduce I/O

### High memory usage

**Solutions**:

1. Close other applications to free memory

2. Process smaller batches of data

3. Use streaming/chunked processing where available

## Getting More Help

If your issue isn't listed here:

1. **Search existing issues**: [GitHub Issues](https://github.com/bigbio/qpx/issues)

2. **Enable verbose logging**: Add `--verbose` to any command for detailed output

3. **Create a new issue**: Include:
   - qpx version (`qpxc --version`)
   - Python version (`python --version`)
   - Operating system
   - Complete error message
   - Minimal reproducible example

4. **Community support**: Check the [Community](community.md) page for additional resources
