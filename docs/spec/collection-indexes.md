# Collection Indexes

An **index** is a materialized, pre-computed data structure that enables fast search across all datasets in a [collection](collection.md). Unlike the collection itself (which is programmatic), indexes are **persisted as Parquet files** because the alternative -- scanning terabytes of PSM data for every query -- is impractical.

## Why Indexes?

Consider a collection with 72 datasets totaling 1.18 TB of PSM data. Without an index, searching for a peptide means reading all 1.18 TB:

```python
# Without index: scans ALL parquet files (~minutes)
coll.sql("SELECT * FROM psm WHERE sequence = 'PEPTIDEK'")

# With index: reads only the relevant partition (~milliseconds)
coll.index("peptide").search("PEPTIDEK")
```

An index trades **one-time build cost** (hours) for **instant queries** (milliseconds). Think of it as a library catalog: the catalog is small and fast to search, and it tells you which book (dataset/file) to open.

## Index Types

QPX defines a standard structure for indexes. The initial spec includes one index type, with the framework designed to support future types.

### Peptide Index

Maps peptide sequences to the datasets and files that contain them.

**Stored at**: `{collection_root}/_index/peptide/`

**Partitioned by**: First two amino acids of the sequence (Hive-style), yielding ~400 partitions (`AA`, `AC`, `AD`, ..., `YW`, `YY`). Each partition is typically 1-5 MB compressed.

```text
_index/
└── peptide/
    ├── sequence_prefix=AA/data.parquet
    ├── sequence_prefix=AC/data.parquet
    ├── sequence_prefix=AD/data.parquet
    ├── ...
    ├── sequence_prefix=PE/data.parquet    ← "PEPTIDEK" lives here
    └── sequence_prefix=YY/data.parquet
```

#### Schema

Each row in the peptide index represents one (peptide, peptidoform, dataset) combination:

| Field | Type | Description |
|-------|------|-------------|
| `sequence` | string | Unmodified peptide amino acid sequence |
| `peptidoform` | string | Peptide with modifications in ProForma notation |
| `project_accession` | string | Dataset identifier (e.g., `PXD000561`) |
| `charge_states` | array[int16] | Distinct precursor charge states observed |
| `spectra_count` | int32 | Number of PSMs for this peptide in this dataset |
| `best_pep` | float64 | Best (lowest) posterior error probability across all PSMs |
| `protein_accessions` | array[string] | All proteins this peptide maps to in this dataset |

!!! note "What is NOT in the index"
    The index does not store spectral data (m/z arrays, intensities), per-PSM details, or run-level information. It is a lightweight pointer -- "this peptide exists in this dataset with this confidence." To get full PSM details, query the dataset's `psm.parquet` directly.

#### Querying

```python
# Using the collection API
coll = qpx.open_collection("/data/msnet/")

# Search by exact sequence
results = coll.index("peptide").search("PEPTIDEK")
# Returns DataFrame:
#   project_accession | peptidoform           | charge_states | spectra_count | best_pep
#   PXD000561         | PEPTIDEK              | [2, 3]        | 142           | 1.2e-15
#   PXD000561         | PEP(Phospho)TIDEK     | [2]           | 8             | 3.4e-06
#   PXD002137         | PEPTIDEK              | [2, 3, 4]     | 891           | 8.7e-22

# Search by prefix (wildcard)
results = coll.index("peptide").search_prefix("PEPTID")

# Search by peptidoform (with modifications)
results = coll.index("peptide").search_peptidoform("PEP(Phospho)TIDEK")
```

Under the hood, this translates to:

```sql
-- DuckDB reads only the 'PE' partition (~1-5 MB), not the full index
SELECT *
FROM read_parquet('_index/peptide/sequence_prefix=PE/data.parquet')
WHERE sequence = 'PEPTIDEK'
```

### Future Index Types

The index framework is designed to be extensible. Potential future indexes:

| Index | Partition key | Use case |
|-------|--------------|----------|
| **protein** | First 2 chars of accession | "Which datasets identified P04637 (TP53)?" |
| **spectrum** | Binned precursor m/z | "Find spectra matching this m/z and charge" |
| **modification** | Modification name | "Which datasets have phosphorylation data?" |

Each would follow the same pattern: Hive-partitioned Parquet under `_index/{type}/`.

## Building Indexes

Indexes are built from the collection's QPX datasets using a batch process:

```python
# Build the peptide index for the entire collection
coll = qpx.open_collection("/data/msnet/")
coll.build_index("peptide")
```

This:

1. Reads `psm.parquet` from every dataset in the collection.
2. Filters out decoy PSMs (`is_decoy = false`).
3. Aggregates per (sequence, peptidoform, project_accession).
4. Partitions by `sequence_prefix` (first 2 amino acids).
5. Writes Hive-partitioned Parquet files to `_index/peptide/`.

### Build Algorithm

```sql
-- Pseudocode for the peptide index build
INSERT INTO _index/peptide/
SELECT
    sequence,
    peptidoform,
    project_accession,
    LIST(DISTINCT charge)                    AS charge_states,
    COUNT(*)                                 AS spectra_count,
    MIN(posterior_error_probability)          AS best_pep,
    LIST(DISTINCT UNNEST(protein_accessions)) AS protein_accessions,
    LEFT(sequence, 2)                        AS sequence_prefix
FROM psm
WHERE is_decoy = false
GROUP BY sequence, peptidoform, project_accession
```

### Rebuilding

Indexes can go stale when datasets are added or removed. The rebuild process is idempotent -- it drops the existing index and recreates it from scratch:

```python
# Rebuild after adding new datasets
coll.build_index("peptide", force=True)

# Check if index is stale (compares dataset list vs index metadata)
coll.index("peptide").is_stale()
```

### Index Metadata

Each index stores a metadata file at `_index/{type}/_metadata.json`:

```json
{
    "index_type": "peptide",
    "qpx_version": "1.0",
    "built_at": "2026-04-08T14:30:00Z",
    "datasets_included": ["PXD000561", "PXD002137", "..."],
    "total_entries": 14400000,
    "partitions": 400,
    "compression": "zstd"
}
```

The `datasets_included` list enables staleness detection: if the current collection has datasets not in this list, the index is stale.

## Storage

Indexes are lightweight compared to the source data:

| Component | Size | Notes |
|-----------|------|-------|
| Source PSM data | ~1.18 TB | The actual data |
| Peptide index | ~0.6 GB | ~0.05% of source |
| Protein index (future) | ~0.1 GB | Estimated |

## Conventions

- Indexes live under `_index/` at the collection root.
- The `_` prefix signals "infrastructure, not a dataset" -- collection discovery skips it.
- Each index type gets its own subdirectory: `_index/peptide/`, `_index/protein/`, etc.
- Partition column names use the `{field}_prefix` naming pattern.
- Index Parquet files use ZSTD compression.
- Index metadata is stored as `_metadata.json` (not Parquet) for easy inspection.
