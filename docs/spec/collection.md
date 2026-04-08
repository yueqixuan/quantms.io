# Collections

A **collection** is a group of QPX datasets stored under a common path. Collections enable cross-dataset operations -- browsing available datasets, querying across all of them, and building materialized indexes for fast search.

Collections are **programmatic, not persisted**. There is no `collection.parquet` file. The collection is defined by a directory convention: any subfolder that contains a `*.dataset.parquet` file is recognized as a QPX dataset within the collection.

## Use Cases

- **Public repositories**: A large-scale reanalysis effort produces QPX datasets for dozens or hundreds of public proteomics projects. Users browse the collection to find datasets by organism or instrument, then download individual projects.
- **Lab-scale analysis**: A research group stores all their QPX-converted experiments under a shared directory. Cross-dataset queries find which experiments identified a protein of interest.
- **AI/ML training**: Researchers discover and download datasets matching specific criteria (organism, size, instrument) to build training corpora for machine learning models.
- **Multi-cohort studies**: Datasets from different cohorts or conditions are grouped into a collection for comparative analysis without merging into a single dataset.

## Directory Convention

A collection is any directory where immediate subfolders contain QPX datasets:

```text
my_collection/                          # Collection root
├── PXD000561/                          # QPX dataset
│   ├── PXD000561.dataset.parquet
│   ├── PXD000561.sample.parquet
│   ├── PXD000561.run.parquet
│   ├── PXD000561.sdrf.tsv
│   └── psm/                           # Hive-partitioned PSMs
│       ├── run_file_name=run_01/
│       │   └── part-0.parquet
│       └── run_file_name=run_02/
│           └── part-0.parquet
│
├── PXD002137/                          # Another QPX dataset
│   ├── PXD002137.dataset.parquet
│   └── ...
│
└── _index/                             # Collection indexes (optional)
    └── peptide/
        ├── sequence_prefix=AA/data.parquet
        └── ...
```

### Discovery Rules

1. Scan immediate subfolders of the collection root.
2. A subfolder is a **dataset** if it contains at least one `*.dataset.parquet` file.
3. Subfolders starting with `_` (e.g., `_index/`) are **reserved** for collection infrastructure and are not treated as datasets.
4. Nested collections (collections inside collections) are not supported.

## Opening a Collection

```python
import qpx

# Local directory
coll = qpx.open_collection("/data/my_collection/")

# S3
coll = qpx.open_collection(
    "s3://bucket/my_collection/",
    s3_config={"region": "us-east-1", "anonymous": True},
)

# List discovered datasets
print(coll.datasets)
# ['PXD000561', 'PXD002137', 'PXD012636', ...]
```

## Cross-Dataset Queries

Collections register all datasets in a single DuckDB engine, enabling SQL queries across datasets:

```python
# How many PSMs per project?
coll.sql("""
    SELECT project_accession, COUNT(*) AS num_psms
    FROM psm
    GROUP BY project_accession
    ORDER BY num_psms DESC
""")

# Find a peptide across all projects
coll.sql("""
    SELECT project_accession, peptidoform, charge, COUNT(*) AS spectra
    FROM psm
    WHERE sequence = 'PEPTIDEK'
    GROUP BY project_accession, peptidoform, charge
""")
```

### Derived Collection Summary

The collection summary is computed on demand from existing QPX structures, not stored:

```python
# Aggregated view — derived from dataset.parquet + sample.parquet + run.parquet
summary = coll.summary()

# Returns a DataFrame with one row per dataset:
#   project_accession | organisms          | instruments       | num_runs | num_psms
#   PXD000561         | [Homo sapiens]     | [Q Exactive HF]  | 24       | 4,200,000
#   PXD012636         | [Homo sapiens, ... | [Orbitrap Fusion] | 48       | 12,300,000
```

This is equivalent to:

```sql
SELECT d.project_accession,
       LIST(DISTINCT s.organism) AS organisms,
       COUNT(DISTINCT r.run_file_name) AS num_runs
FROM dataset d
JOIN sample s USING (project_accession)
JOIN run r USING (project_accession)
GROUP BY d.project_accession
```

### Species Search

Because organism information lives in `sample.parquet`, species search works without any additional infrastructure:

```python
# Find all mouse datasets
coll.sql("""
    SELECT DISTINCT d.project_accession, d.project_title
    FROM dataset d
    JOIN sample s USING (project_accession)
    WHERE list_contains(s.organism, 'Mus musculus')
""")
```

## Bulk Download

Users download individual datasets from a collection:

```bash
# Download one project
aws s3 sync s3://bucket/my_collection/PXD000561/ ./PXD000561/

# Download all human datasets (using the summary to filter first)
for project in $(python3 -c "
import qpx
coll = qpx.open_collection('s3://bucket/my_collection/', s3_config={'anonymous': True})
for row in coll.sql(\"\"\"
    SELECT d.project_accession FROM dataset d
    JOIN sample s USING (project_accession)
    WHERE list_contains(s.organism, 'Homo sapiens')
\"\"\").fetchall():
    print(row[0])
"); do
    aws s3 sync "s3://bucket/my_collection/$project/" "./$project/"
done
```

## Indexes

Collections can optionally include **materialized indexes** for fast search across datasets. See [Collection Indexes](collection-indexes.md) for details.

## Limitations

- **No manifest file**: The collection is discovered from the filesystem. If a dataset is added or removed, the collection reflects the change on next open.
- **No cross-dataset joins on partitioned columns**: Hive partition columns (e.g., `run_file_name`) are local to each dataset. Cross-dataset queries that filter on partition columns work, but the partition column values may overlap across datasets.
- **Memory**: Opening a large collection registers all datasets in DuckDB. For collections with hundreds of datasets, consider using `structures=["dataset", "sample"]` to limit which views are loaded.
