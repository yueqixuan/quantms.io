# File Extensions & Naming

QPX defines a consistent file naming convention that encodes the project identity, file uniqueness, data view, and serialization format directly in the filename. This enables both humans and automated tools to identify file contents without reading metadata.

## Convention

All QPX files follow this naming pattern:

```
{PREFIX}-{UUID}.{view}.{format}
```

Where:

- **`{PREFIX}`** -- Usually the ProteomeXchange project accession (e.g., `PXD000000`). For non-ProteomeXchange datasets, any short identifier may be used.
- **`{UUID}`** -- A UUID v4 identifier ([RFC 4122](https://www.rfc-editor.org/rfc/rfc4122)). Optional but recommended. Ensures globally unique filenames even when the same project is processed multiple times.
- **`{view}`** -- One of the QPX spec-defined views (e.g., `psm`, `feature`, `pg`, `peptide`, `protein`, `mz`, `absolute`, `differential`, `sdrf`).
- **`{format}`** -- The serialization format extension (`parquet`, `tsv`, or `json`).

!!! note
    The UUID component is separated from the prefix by a hyphen (`-`), while the view and format are separated by dots (`.`). This makes it straightforward to parse filenames programmatically.

## Examples

| File Name | View | Format |
|-----------|------|--------|
| `PXD000000-943a8f02.psm.parquet` | PSM | Parquet |
| `PXD000000-943a8f02.feature.parquet` | Feature | Parquet |
| `PXD000000-943a8f02.pg.parquet` | Protein Group | Parquet |
| `PXD000000-943a8f02.peptide.parquet` | Peptide | Parquet |
| `PXD000000-943a8f02.protein.parquet` | Protein | Parquet |
| `PXD000000-943a8f02.mz.parquet` | Mass Spectra | Parquet |
| `PXD000000-943a8f02.absolute.tsv` | Absolute Expression | TSV |
| `PXD000000-943a8f02.differential.tsv` | Differential Expression | TSV |
| `PXD000000.sdrf.tsv` | Sample Metadata | TSV |

## Metadata Files (Fixed Names)

Metadata files use fixed names without a UUID component. These files are project-level artifacts that are not expected to have multiple versions within a single project:

| File Name Pattern | Description |
|-------------------|-------------|
| `{ACCESSION}.experiment.parquet` | Experiment-level metadata |
| `{ACCESSION}.sdrf.parquet` | Sample metadata in Parquet format |
| `{ACCESSION}.sdrf.tsv` | Original SDRF file, preserved for provenance |

!!! info "Why preserve the original SDRF?"
    The `.sdrf.tsv` file is kept alongside the `.sdrf.parquet` conversion to maintain a direct link to the original experimental metadata as submitted to ProteomeXchange. This ensures full provenance and allows downstream tools that expect TSV-formatted SDRF files to work without conversion.

## UUID Generation

UUIDs ensure globally unique filenames across projects, reprocessing runs, and institutions. QPX uses UUID version 4, which is generated from random numbers and does not contain any identifying information about the host or time of generation.

In Python, UUIDs can be generated using the standard library:

```python
import uuid

file_uuid = uuid.uuid4()
filename = f"PXD012345-{file_uuid}.psm.parquet"
# Example output: PXD012345-943a8f02-0527-4528-b1a3-b96de99ebe75.psm.parquet
```

!!! tip "Short UUIDs in examples"
    Throughout the QPX documentation, UUIDs are often shortened to the first 8 characters (e.g., `943a8f02`) for readability. In practice, the full 36-character UUID should be used to guarantee uniqueness.

## See Also

- [Serialization](serialization.md) -- serialization formats and Parquet details
- [QPX Format Overview](index.md) -- overview of all QPX views and their relationships
