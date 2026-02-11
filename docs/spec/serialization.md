# Serialization & Parquet Details

QPX uses a combination of serialization formats, each chosen for its suitability to the data it stores. The primary format is Apache Parquet for columnar proteomics data, with TSV and JSON used for expression matrices and project metadata, respectively.

## Supported Serialization Formats

| Format | Extension | Used By | Description |
|--------|-----------|---------|-------------|
| Apache Parquet | `.parquet` | PSM, Feature, PG, Peptide, Protein, MZ | Columnar storage with compression |
| TSV | `.tsv` | AE, DE (v1.0), SDRF | Tab-separated values |
| JSON | `.json` | Project (v1.0) | Key-value metadata |

!!! note
    The choice of format for each view is driven by its access patterns. Views that benefit from column-level filtering and large-scale analytical queries use Parquet. Expression matrices and sample metadata remain in TSV for broad tool compatibility.

## Parquet Format

Apache Parquet is a columnar storage format designed for efficient analytical processing. Unlike row-oriented formats (CSV, TSV, mzTab), Parquet stores data column by column, which provides significant advantages for proteomics workflows:

- **File metadata and column metadata**: Each Parquet file includes a footer that describes the schema, row groups, column statistics (min/max values), and custom key-value metadata. This allows query engines to skip irrelevant data without reading the full file.
- **Column-oriented design**: Analytical queries that access a subset of columns (e.g., retrieving only `sequence` and `precursor_charge` from a PSM file) read only the relevant column chunks, reducing I/O dramatically compared to row-based formats.
- **Broad ecosystem support**: Parquet is natively supported by PyArrow, DuckDB, Polars, Apache Spark, pandas, R (via the `arrow` package), and many other tools. This makes QPX files immediately usable across languages and frameworks without custom parsers.

## Compression

Parquet's columnar layout achieves substantial compression ratios because values within a column tend to have similar data types and distributions. The following benchmarks illustrate the storage savings and write performance for real-world proteomics datasets:

| Project | Type | Original Size (GB) | Parquet Size (MB) | PSM Write Time (s) | Feature Write Time (s) |
|---------|------|--------------------:|-------------------:|--------------------:|-----------------------:|
| PXD046440 | MaxQuant | 48 | 337/343 | 985.27 | 678.47 |
| PXD016999 | mzTab | 160 | 155/228 | 539.00 | 3554.53 |
| PXD019909 | DIA-NN | 1.9 | 195 | -- | 229.48 |

!!! tip "Compression algorithms"
    QPX supports the following compression algorithms for Parquet files:

    - **Snappy** -- Fast compression/decompression with moderate ratio (default)
    - **Gzip** -- Higher compression ratio, slower speed
    - **LZO** -- Optimized for decompression speed
    - **None** -- No compression

    In addition to block-level compression, Parquet applies encoding schemes at the column level:

    - **RLE (Run-Length Encoding)** -- Efficient for columns with many repeated values (e.g., `is_decoy`, `precursor_charge`)
    - **Dictionary encoding** -- Replaces repeated string values with integer codes (e.g., `reference_file_name`, `sequence`)

## Parquet Features

QPX leverages the following Parquet capabilities:

- **Columnar Storage**: Each column is stored independently, enabling query engines to read only the columns needed for a given analysis. This reduces I/O for analytical queries that access a subset of fields.
- **Efficient Compression**: Column-level compression with Snappy, Gzip, or LZO, combined with RLE and dictionary encoding, achieves compression ratios of 10x--100x on typical proteomics datasets.
- **Schema Evolution**: Columns can be added, removed, or modified across QPX versions without breaking compatibility with existing files. Readers that encounter unknown columns can safely ignore them.
- **Complex Data Types**: Parquet natively supports nested structs, arrays, and maps. QPX uses these for fields like `modifications` (array of structs), `additional_scores` (array of structs), and spectral arrays (`mz_array`, `intensity_array`).

## Parquet Slicing (Partitioning)

QPX supports partitioning Parquet files by any field, enabling directory-based data organization. This is particularly useful for large projects where loading all data into memory is impractical. Partitioning by `sample_accession` or `reference_file_name` allows tools to read only the data for specific samples or MS runs.

Example directory structure partitioned by sample accession:

```
PXD004683/
├── sample_accession_1/
│   ├── file1.parquet
│   └── file2.parquet
├── sample_accession_2/
│   └── file3.parquet
```

!!! warning "Practical limits of Parquet metadata"
    Parquet stores all file-level and column-level metadata in the file footer. When reading a partitioned dataset, the footer of every partition file must be loaded into memory to reconstruct the full schema and row group statistics. For datasets with thousands of partitions, this metadata overhead can become significant. Consider using fewer, larger partitions rather than many small ones.

## File-Level Metadata

Every QPX Parquet file includes metadata as key-value pairs stored in the Parquet file footer. These fields identify the file's origin, format version, and provenance:

| Key | Description |
|-----|-------------|
| `qpx_version` | Version of the QPX specification used to generate the file |
| `software_provider` | Name and version of the software that produced the data |
| `project_accession` | ProteomeXchange or project accession (e.g., `PXD012345`) |
| `project_title` | Human-readable title of the project |
| `scan_format` | Format of scan identifiers: `scan`, `index`, `nativeId`, or `multiple` |
| `creator` | Name of the tool or person who created the file |
| `file_type` | QPX view type (e.g., `psm_file`, `feature_file`, `pg_file`) |
| `creation_date` | ISO 8601 date when the file was created |
| `uuid` | UUID v4 identifier for the file |
| `compression_format` | Compression algorithm used: `snappy`, `gzip`, `lzo`, or `none` |

!!! tip "Writing Parquet with metadata in Python"
    ```python
    import pyarrow as pa
    import pyarrow.parquet as pq

    # Build your table
    table = pa.table({"sequence": ["PEPTIDER", "ANOTHERPEPTIDE"]})

    # Define file-level metadata
    file_metadata = {
        'qpx_version': '1.0',
        'software_provider': 'QuantMS 1.3.0',
        'project_accession': 'PXD012345',
        'file_type': 'psm_file',
        'creation_date': '2021-01-01',
        'uuid': '943a8f02-0527-4528-b1a3-b96de99ebe75'
    }

    # Merge with existing schema metadata and write
    existing_metadata = table.schema.metadata or {}
    merged_metadata = {
        **existing_metadata,
        **{k.encode(): v.encode() for k, v in file_metadata.items()}
    }
    table = table.replace_schema_metadata(merged_metadata)
    pq.write_table(table, 'psm_data.parquet')
    ```

!!! tip "Reading file metadata in Python"
    ```python
    import pyarrow.parquet as pq

    parquet_file = pq.ParquetFile("psm_data.parquet")
    metadata = parquet_file.schema_arrow.metadata
    for key, value in metadata.items():
        print(f"{key.decode()}: {value.decode()}")
    ```

## See Also

- [File Naming](file-naming.md) -- file extension and naming conventions for QPX files
- [Versioning](versioning.md) -- version format and backward compatibility rules
- [QPX Format Overview](index.md) -- overview of all QPX views and their relationships
