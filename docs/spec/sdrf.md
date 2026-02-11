# Sample Metadata (SDRF)

The SDRF (Sample and Data Relationship Format) describes the relationship between biological samples, data files, and experimental factors in a proteomics experiment. In QPX, SDRF serves as the **ingestion format** for sample metadata -- at conversion time, the original bracket-style column names are translated to clean snake_case field names and the file is stored as `sdrf.parquet`.

## Overview

SDRF is a tab-delimited format originally defined by the [proteomics-sample-metadata](https://github.com/bigbio/proteomics-sample-metadata) project. Each row represents a single data file (or channel within a multiplexed file), and columns describe the biological sample, technical parameters, and experimental factors associated with that file.

!!! info "SDRF as an ingestion format"
    QPX does **not** invent its own sample metadata model. It adopts the community SDRF standard and converts it to a query-friendly Parquet representation. The original `.sdrf.tsv` file is always preserved alongside the converted Parquet for provenance.

Key properties of SDRF in QPX:

- **One row per data-file-channel combination** -- for label-free experiments each row corresponds to one raw file; for multiplexed experiments (TMT/iTRAQ) each row corresponds to one channel within a raw file
- **Bracket columns are translated** -- `characteristics[organism]` becomes `organism`, `comment[data file]` becomes `reference_file_name`
- **Factor values are preserved** -- columns like `factor value[disease]` become `factor_disease` and indicate which sample properties are experimental variables

---

## Column Name Mapping

When SDRF is converted to QPX, bracket-style column names are mapped to snake_case field names according to the following registry. This mapping is deterministic and applied automatically during ingestion.

### Sample Columns

Sample columns describe biological properties of the specimen.

| SDRF Column | QPX Field Name | Category |
|---|---|---|
| `characteristics[organism]` | `organism` | sample |
| `characteristics[organism part]` | `organism_part` | sample |
| `characteristics[disease]` | `disease` | sample |
| `characteristics[cell line]` | `cell_line` | sample |
| `characteristics[cell type]` | `cell_type` | sample |
| `characteristics[sex]` | `sex` | sample |
| `characteristics[age]` | `age` | sample |
| `characteristics[individual]` | `individual` | sample |
| `characteristics[biological replicate]` | `biological_replicate` | sample |

### File Columns

File columns describe technical and instrument-related properties of the data acquisition.

| SDRF Column | QPX Field Name | Category |
|---|---|---|
| `comment[data file]` | `reference_file_name` | file |
| `comment[label]` | `labeling_channel` | file |
| `comment[fraction identifier]` | `fraction` | file |
| `comment[technical replicate]` | `technical_replicate` | file |
| `comment[instrument]` | `instrument` | file |
| `comment[cleavage agent details]` | `enzyme` | file |
| `comment[proteomics data acquisition method]` | `acquisition_method` | file |

### Factor Columns

Factor columns mark which sample properties are experimental variables in the study design.

| SDRF Column | QPX Field Name | Category |
|---|---|---|
| `factor value[disease]` | `factor_disease` | factor |
| `factor value[organism part]` | `factor_organism_part` | factor |

!!! note
    Factor columns are open-ended. Any `factor value[X]` column in the SDRF is mapped to `factor_X` (with spaces replaced by underscores). The table above shows the most common factor columns; additional factors are mapped using the same convention.

---

## Entity Model (v2.0)

In v2.0, every column in the SDRF is assigned to one of three **semantic levels**. The level determines how the column is used in downstream queries and which library-level views expose it.

### Semantic Levels

1. **Sample columns** (biological properties)
    - Describe the biological specimen: `organism`, `organism_part`, `disease`, `cell_line`, `cell_type`, `sex`, `age`, `individual`, `biological_replicate`
    - Used to define unique biological samples and to derive aggregated project-level lists (e.g. all organisms in the experiment)

2. **File columns** (technical properties)
    - Describe data acquisition: `reference_file_name`, `labeling_channel`, `fraction`, `technical_replicate`, `instrument`, `enzyme`, `acquisition_method`
    - Vary at the run/channel level and are not meaningful for defining biological samples

3. **Factor columns** (experimental design)
    - Mark which sample properties are independent variables: `factor_disease`, `factor_organism_part`, etc.
    - Used by differential expression analysis to define contrasts

The `column_semantics` map in [experiment.parquet](project.md) defines which semantic level each column belongs to:

```json
{
  "organism": "sample",
  "organism_part": "sample",
  "disease": "sample",
  "reference_file_name": "file",
  "labeling_channel": "file",
  "instrument": "file",
  "factor_disease": "factor"
}
```

### Library-Level Views

The QPX library provides two convenience views built on top of the SDRF:

| View | Access | Description |
|---|---|---|
| **Samples** | `project.samples` | Deduplicated biological samples -- only sample-level columns. Each row represents a unique biological specimen. |
| **Runs** | `project.runs` | All file + channel rows -- all columns. Each row represents a data-file-channel combination as written in the original SDRF. |

!!! example "Samples vs Runs"
    In a TMT-10plex experiment with 3 fractions per plex, `project.runs` contains 30 rows (10 channels x 3 fractions) while `project.samples` contains 10 rows (one per biological sample, with fraction and file columns removed).

    ```python
    # Get unique biological samples
    samples = project.samples
    print(len(samples))  # 10

    # Get all file-channel rows
    runs = project.runs
    print(len(runs))     # 30
    ```

---

## Multi-valued Metadata

In rare cases a single SDRF row may describe a sample derived from multiple organisms (e.g. a spike-in experiment) or multiple cell lines (pooled samples). QPX handles this with a **pipe-delimiter convention**.

### Convention

Multi-valued metadata is stored as a `pa.string()` column with individual values separated by the pipe character (`|`):

```
Homo sapiens|Saccharomyces cerevisiae
```

This keeps the Parquet schema flat (no nested lists for sample-level columns) while still supporting the uncommon multi-organism or multi-disease case.

### Parsing

The QPX library provides a helper to split pipe-delimited values:

```python
# Parse a multi-valued organism field
organisms = project.parse_multi("organism")
# Returns: ["Homo sapiens", "Saccharomyces cerevisiae"]
```

!!! warning "Pipe delimiter is for exceptional cases only"
    The vast majority of proteomics experiments have single-valued sample metadata. The pipe delimiter should only be used for genuine pooled or spike-in scenarios -- not as a general-purpose list encoding. If more than ~5% of rows in a column contain pipes, consider restructuring the experimental design.

### SQL Queries

When querying `sdrf.parquet` with DuckDB or similar engines, use `string_split` to unnest pipe-delimited values:

```sql
SELECT DISTINCT unnest(string_split(organism, '|')) AS organism
FROM sdrf
```

---

## Proposed sdrf.parquet Format

In v2.0, the SDRF TSV is converted to a Parquet file with clean snake_case column names. The conversion is lossless -- every column and every row from the original TSV is preserved.

### Conversion Rules

1. **Column renaming** -- bracket columns are mapped to snake_case names using the registry above
2. **Type inference** -- all columns are stored as `pa.string()` (sample metadata is inherently textual)
3. **Row preservation** -- one row per data-file-channel combination, matching the original TSV
4. **File naming** -- the output file is named `{accession}.sdrf.parquet` (e.g. `PXD014414.sdrf.parquet`)

### Provenance

The original SDRF TSV is always preserved alongside the Parquet conversion:

```
PXD014414/
  PXD014414.sdrf.tsv        # Original SDRF (provenance)
  PXD014414.sdrf.parquet     # Converted Parquet (query-optimized)
  PXD014414.experiment.parquet
  PXD014414.feature.parquet
  ...
```

!!! tip "Reading sdrf.parquet"
    ```python
    import pyarrow.parquet as pq

    sdrf = pq.read_table("PXD014414.sdrf.parquet")
    print(sdrf.column_names)
    # ['organism', 'organism_part', 'disease', 'reference_file_name',
    #  'labeling_channel', 'fraction', 'instrument', 'enzyme', ...]

    # Derive project-level organism list (replaces project.json 'organisms' field)
    organisms = sdrf.column("organism").unique().to_pylist()
    print(organisms)  # ['Homo sapiens']
    ```

### Example Rows

| organism | organism_part | disease | reference_file_name | labeling_channel | fraction | instrument |
|---|---|---|---|---|---|---|
| Homo sapiens | breast | breast cancer | 20190101_sample01 | LFQ | 1 | Q Exactive HF |
| Homo sapiens | breast | normal | 20190101_sample02 | LFQ | 1 | Q Exactive HF |
| Homo sapiens | breast | breast cancer | 20190101_sample03 | LFQ | 1 | Q Exactive HF |

---

## Related Pages

- [Project / Experiment Metadata](project.md)
- [Feature View](feature.md)
- [Modifications](modifications.md)
- [QPX Format Overview](index.md)
