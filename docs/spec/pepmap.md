# Peptide-Protein Mapping

`pepmap.parquet` is a **deduplicated mapping** between peptide sequences and the proteins they belong to. Each row represents a unique (sequence, protein_accession) pair.

The YAML schema definition lives at `qpx/core/data/schemas/pepmap.yaml` in the source repository.

## Use Cases

- **Peptide uniqueness analysis**: Identify peptides that map to a single protein (proteotypic) vs. shared peptides.
- **Protein coverage**: Determine which regions of a protein are covered by identified peptides.
- **Gene-level rollup**: Map peptides to genes via the `gene_name` field.
- **Space-efficient storage**: At scale (10M+ PSMs), storing protein mappings as a separate table avoids duplicating protein accession strings across millions of PSM rows.

## Schema

| Field | Description | Type | Required |
|-------|-------------|------|----------|
| `sequence` | Unmodified peptide amino acid sequence | string | yes |
| `protein_accession` | UniProt protein accession | string | yes |
| `protein_name` | Protein entry name (e.g., `BRCA1_HUMAN`) | string | no |
| `gene_name` | Gene symbol (e.g., `BRCA1`) | string | no |
| `start_position` | Start position in protein sequence (1-indexed) | int32 | no |
| `end_position` | End position in protein sequence (1-indexed) | int32 | no |
| `is_unique` | `true` if peptide maps to exactly one protein | bool | no |
| `is_proteotypic` | `true` if peptide is proteotypic (unique to one gene) | bool | no |

!!! note "Primary key"
    The primary key is `(sequence, protein_accession)`. Each unique peptide-protein pair appears exactly once.

## Python API

```python
import qpx

ds = qpx.open("PXD014414/")

# Access the mapping
mapping = ds.pepmap

# Filter by protein
brca1_peptides = mapping.by_protein("P38398")
print(brca1_peptides.to_df())

# Filter by peptide
proteins = mapping.by_peptide("PEPTIDEK")
print(proteins.to_df())

# Get only unique (proteotypic) peptides
unique = mapping.unique_peptides()
print(f"Unique peptides: {unique.count()}")
```

## Relationship to Other Views

The `pepmap` complements the `protein_accessions` field on PSM and Feature views:

- **PSM/Feature `protein_accessions`**: A convenience list on each row. Simple but duplicates strings at scale.
- **`pepmap.parquet`**: A deduplicated lookup table with additional fields (gene name, positions, uniqueness flags). More space-efficient for large datasets.

Both approaches are valid. Use `protein_accessions` for quick lookups. Use `pepmap` when you need protein metadata or are working at scale.

## Writing a Peptide-Protein Map

```python
from qpx.writers import PepMapWriter

records = [
    {
        "sequence": "PEPTIDEK",
        "protein_accession": "P12345",
        "protein_name": "BRCA1_HUMAN",
        "gene_name": "BRCA1",
        "start_position": 100,
        "end_position": 107,
        "is_unique": True,
        "is_proteotypic": True,
    },
]

with PepMapWriter("exp.pepmap.parquet") as w:
    w.write_batch(records)
```

## Related Pages

- [PSM View](psm.md) -- carries optional `protein_accessions` list
- [Feature View](feature.md) -- `pg_accessions` and `anchor_protein` fields
- [Protein Group View](pg.md) -- `pg_accessions` and protein group inference
