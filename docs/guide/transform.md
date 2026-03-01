# Transform Commands

!!! note "v2.0 Migration"
    The `--project-file` parameter in the examples below refers to `project.json` (v1.0 format). In QPX v2.0, this will be replaced by `dataset.parquet`. See [Dataset Metadata](../spec/dataset.md) for the new format.

Transform and process data within the QPX ecosystem.

```python exec="1" session="doc_utils" result="ansi"
import click
import textwrap

def get_click_type_display(param):
    param_type = param.type
    type_str = str(param_type)
    if 'Path' in type_str:
        if hasattr(param_type, 'dir_okay') and not param_type.dir_okay:
            return 'FILE'
        elif hasattr(param_type, 'file_okay') and not param_type.file_okay:
            return 'DIRECTORY'
        else:
            return 'PATH'
    elif isinstance(param_type, click.types.FloatParamType):
        return 'FLOAT'
    elif isinstance(param_type, click.types.IntParamType):
        return 'INTEGER'
    elif param.is_flag:
        return 'FLAG'
    else:
        return 'TEXT'

def generate_params_table(command):
    table = '<table>\n<thead>\n<tr>\n'
    table += '<th>Parameter</th><th>Type</th><th>Required</th><th>Default</th><th>Description</th>\n'
    table += '</tr>\n</thead>\n<tbody>\n'
    for param in command.params:
        if isinstance(param, click.Option) and param.name not in ['help']:
            param_names = param.opts
            param_name = param_names[0] if param_names else f"--{param.name}"
            param_type = get_click_type_display(param)
            required = 'Yes' if param.required else 'No'

            # Extract default value from help text if it contains "(default: ...)"
            description = param.help or ''
            default_from_help = None
            if '(default:' in description.lower():
                import re
                match = re.search(r'\(default:\s*([^)]+)\)', description, re.IGNORECASE)
                if match:
                    default_from_help = match.group(1).strip()

            # Determine default value display
            if default_from_help:
                default = default_from_help
            elif param.default is not None and str(param.default) != 'Sentinel.UNSET':
                if param.is_flag:
                    default = '-'
                elif isinstance(param.default, (int, float)):
                    default = str(param.default)
                elif isinstance(param.default, str):
                    default = f'<code>{param.default}</code>'
                else:
                    default = str(param.default)
            else:
                default = '-'

            table += f'<tr>\n<td><code>{param_name}</code></td>\n<td>{param_type}</td>\n<td>{required}</td>\n<td>{default}</td>\n<td>{description}</td>\n</tr>\n'
    table += '</tbody>\n</table>'
    return table

def generate_description(command):
    if command.help:
        help_text = command.help
        if 'Example' in help_text:
            description = help_text.split('Example')[0].strip()
        else:
            description = help_text.strip()
        lines = description.split('\n')
        if len(lines) > 1:
            description = '\n'.join(lines[1:]).strip()
            return f'<p>{description}</p>'
    return ''

def generate_example(command, default_text=''):
    if command.help and 'Example' in command.help:
        example_section = command.help.split('Example')[1]
        if ':' in example_section:
            example_section = example_section.split(':', 1)[1]
        example_section = textwrap.dedent(example_section).strip()
        output = ''
        if default_text:
            output += f'<p>{default_text}</p>\n'
        output += f'<pre><code class="language-bash">{example_section}</code></pre>'
        return output
    return ''
```

## Overview

The `transform` command group provides tools for processing and transforming QPX data into various downstream formats. These commands enable absolute and differential expression analysis, metadata mapping, and data format conversions.

## Available Commands

- [ibaq](#ibaq) - Compute iBAQ from features
- [ae](#ae) - Create absolute expression data
- [de](#de) - Create differential expression data
- [gene-map](#gene-map) - Map genes from FASTA

---

## ibaq

Compute iBAQ from feature data.

### Description {#ibaq-description}

```python exec="1" html="1" session="doc_utils"
from qpx.cli.transform import transform_ibaq_cmd
print(generate_description(transform_ibaq_cmd))
```

### Parameters {#ibaq-parameters}

```python exec="1" html="1" session="doc_utils"
from qpx.cli.transform import transform_ibaq_cmd
print(generate_params_table(transform_ibaq_cmd))
```

### Usage Examples {#ibaq-examples}

#### Basic Example {#ibaq-example-basic}

```python exec="1" html="1" session="doc_utils"
from qpx.cli.transform import transform_ibaq_cmd
print(generate_example(transform_ibaq_cmd, 'Compute iBAQ from feature data:'))
```

#### With Custom Prefix {#ibaq-example-prefix}

```bash
qpxc transform ibaq \
    --feature-file ./output/feature.parquet \
    --sdrf-file ./metadata.sdrf.tsv \
    --output-folder ./output \
    --output-prefix ibaq_quantification
```

### Output Files {#ibaq-output}

- **Output**: `{output-prefix}-{uuid}.ibaq.parquet`
- **Format**: Parquet file containing iBAQ quantification values
- **Content**: Protein-level iBAQ values per sample

### Best Practices {#ibaq-best-practices}

- Ensure feature file contains all necessary quantification data
- Verify SDRF metadata matches sample identifiers in feature file
- Use iBAQ output for absolute protein quantification analysis

---

## ae

Create absolute expression data from iBAQ.

### Description {#ae-description}

```python exec="1" html="1" session="doc_utils"
from qpx.cli.transform import transform_ae_cmd
print(generate_description(transform_ae_cmd))
```

**Format Specification**: For details about the AE format structure and fields, see the [Absolute Expression Format Specification](../spec/absolute.md).

### Parameters {#ae-parameters}

```python exec="1" html="1" session="doc_utils"
from qpx.cli.transform import transform_ae_cmd
print(generate_params_table(transform_ae_cmd))
```

### Usage Examples {#ae-examples}

#### Basic Example {#ae-example-basic}

```python exec="1" html="1" session="doc_utils"
from qpx.cli.transform import transform_ae_cmd
print(generate_example(transform_ae_cmd, 'Convert iBAQ data with default settings:'))
```

#### With Project Metadata {#ae-example-metadata}

```bash
qpxc transform ae \
    --ibaq-file tests/examples/AE/PXD016999.1-ibaq.tsv \
    --sdrf-file tests/examples/AE/PXD016999-first-instrument.sdrf.tsv \
    --project-file tests/examples/AE/project.json \
    --output-folder ./output \
    --output-prefix ae_with_metadata \
    --delete-existing
```

#### Filter Specific Proteins {#ae-example-filter}

```bash
qpxc transform ae \
    --ibaq-file tests/examples/AE/PXD016999.1-ibaq.tsv \
    --sdrf-file tests/examples/AE/PXD016999-first-instrument.sdrf.tsv \
    --protein-file tests/examples/fasta/Homo-sapiens.fasta \
    --output-folder ./output
```

### Input File Formats {#ae-input-formats}

**iBAQ File**: Tab-separated file with protein accessions and iBAQ intensities

```
ProteinID    Sample1    Sample2    Sample3
P12345       1000000    950000     1050000
Q67890       500000     480000     520000
```

**SDRF File**: Standard PRIDE SDRF format with sample metadata

### Output Files {#ae-output}

- **Output**: `{output-prefix}-{uuid}.absolute.parquet`
- **Format**: Parquet file containing absolute expression quantification
- **Schema**: Conforms to [QPX absolute expression specification](../spec/absolute.md)

### Common Issues {#ae-issues}

**Issue**: Mismatched sample names between iBAQ and SDRF

- **Solution**: Ensure column names in iBAQ file match sample identifiers in SDRF

**Issue**: Missing protein accessions

- **Solution**: Provide `--protein-file` to filter and validate protein IDs

### Best Practices {#ae-best-practices}

- Always provide project metadata file when available for better data provenance
- Use `--delete-existing` flag carefully to avoid accidental data loss
- Validate SDRF file format before processing
- Check sample name consistency across input files

---

## de

Create differential expression data from MSstats.

### Description {#de-description}

```python exec="1" html="1" session="doc_utils"
from qpx.cli.transform import transform_de_cmd
print(generate_description(transform_de_cmd))
```

### Parameters {#de-parameters}

```python exec="1" html="1" session="doc_utils"
from qpx.cli.transform import transform_de_cmd
print(generate_params_table(transform_de_cmd))
```

### Usage Examples {#de-examples}

#### Basic Example {#de-example-basic}

```python exec="1" html="1" session="doc_utils"
from qpx.cli.transform import transform_de_cmd
print(generate_example(transform_de_cmd, 'Convert MSstats differential expression data:'))
```

#### With Custom FDR Threshold {#de-example-fdr}

```bash
qpxc transform de \
    --msstats-file tests/examples/DE/PXD033169.sdrf_openms_design_msstats_in_comparisons.csv \
    --sdrf-file tests/examples/DE/PXD033169.sdrf.tsv \
    --fdr-threshold 0.01 \
    --output-folder ./output \
    --output-prefix de_stringent \
    --verbose
```

#### With Project Metadata {#de-example-metadata}

```bash
qpxc transform de \
    --msstats-file tests/examples/DE/PXD033169.sdrf_openms_design_msstats_in_comparisons.csv \
    --sdrf-file tests/examples/DE/PXD033169.sdrf.tsv \
    --project-file tests/examples/DE/project.json \
    --fdr-threshold 0.05 \
    --output-folder ./output \
    --delete-existing
```

### Input File Format {#de-input-format}

**MSstats File**: CSV file with comparison results

```
Protein,Label,log2FC,SE,Tvalue,DF,pvalue,adj.pvalue
P12345,Condition2-Condition1,2.5,0.3,8.33,10,0.0001,0.001
Q67890,Condition2-Condition1,-1.8,0.4,-4.5,10,0.002,0.01
```

### Output Files {#de-output}

- **Output**: `{output-prefix}-{uuid}.differential.parquet`
- **Format**: Parquet file containing differential expression results
- **Schema**: Conforms to QPX differential expression specification

### Common Issues {#de-issues}

**Issue**: No significant results after FDR filtering

- **Solution**: Increase `--fdr-threshold` or check input data quality

**Issue**: Memory errors with large comparison files

- **Solution**: Process comparisons in batches or increase available memory

### Best Practices {#de-best-practices}

- Use FDR threshold of 0.05 or lower for publication-quality results
- Enable verbose mode to monitor filtering statistics
- Validate comparison group names match SDRF metadata
- Include project file for complete data provenance

---

## gene-map

Map gene information from FASTA to parquet format.

### Description {#gene-map-description}

```python exec="1" html="1" session="doc_utils"
from qpx.cli.transform import transform_gene_map_cmd
print(generate_description(transform_gene_map_cmd))
```

### Parameters {#gene-map-parameters}

```python exec="1" html="1" session="doc_utils"
from qpx.cli.transform import transform_gene_map_cmd
print(generate_params_table(transform_gene_map_cmd))
```

### Usage Examples {#gene-map-examples}

#### Basic Example {#gene-map-example-basic}

```python exec="1" html="1" session="doc_utils"
from qpx.cli.transform import transform_gene_map_cmd
print(generate_example(transform_gene_map_cmd, 'Map gene information to parquet file:'))
```

#### With Partitioning {#gene-map-example-partition}

```bash
qpxc transform gene-map \
    --parquet-path ./output/feature.parquet \
    --fasta tests/examples/fasta/Homo-sapiens.fasta \
    --output-folder ./output \
    --file-num 20 \
    --partitions run_file_name \
    --species human
```

### Output Files {#gene-map-output}

- **Output**: Enhanced parquet file(s) with gene information
- **Format**: Parquet file in output folder
- **Added Fields**: Gene names and metadata from FASTA headers

### Best Practices {#gene-map-best-practices}

- Use species-specific FASTA files for accurate gene annotation
- Adjust `--file-num` based on available memory for large files
- Use partitioning for better file organization in large datasets

---

## Related Commands

- [Convert Commands](convert.md) - Convert raw data to QPX format
- [Visualization Commands](visualize.md) - Visualize transformed data
- [Statistics Commands](stats.md) - Analyze transformed data
