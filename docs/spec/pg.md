# Protein Group View

The protein group (PG) view is a tabular Parquet file that contains the details of protein groups identified and quantified per raw file. It captures the relationship between protein groups and the raw files in which they were detected, including peptide counts, feature counts, quality metrics, and intensity-based quantification.

This view is analogous to outputs from tools such as MaxQuant (`proteinGroups.txt`), DIA-NN (`pg_matrix`), and FragPipe protein group reports.

## Use cases

- Retrieve all protein groups identified or quantified in a given raw file.
- Compute protein group abundance by file and condition.
- Store and query FDR q-values for protein groups at both the run and experiment level.
- Support downstream statistical analysis by providing per-file protein-level quantification.

## Schema

Fields marked with **(PK)** are primary keys and MUST NOT be null. Fields marked with **(nullable)** may have null values. See the full Avro schema in [`pg.avsc`](../pg.avsc).

### Identity

| Field | Description | Type | Required |
|-------|-------------|------|----------|
| `pg_accessions` | Protein accessions of all proteins within this group | `array[string]` | Yes (PK) |
| `pg_names` | Descriptive names for the proteins in the group | `array[string]` | No |
| `gg_accessions` | Gene group accessions as a string array | `array[string]` | No |
| `gg_names` | Gene names corresponding to the proteins in the group | `array[string]` | No |
| `anchor_protein` | Representative protein of the group (leading protein) | `string` | No |
| `reference_file_name` | The raw file containing the identified/quantified protein group | `string` | Yes (PK) |

### Counts

| Field | Description | Type | Required |
|-------|-------------|------|----------|
| `peptide_counts` | Peptide sequence counts for this protein group in this file | `struct` | No |
| `peptide_counts.unique_sequences` | Number of peptide sequences unique to this protein group within this file | `int` | -- |
| `peptide_counts.total_sequences` | Total number of peptide sequences identified for this protein group in this file | `int` | -- |
| `feature_counts` | Peptide feature counts (peptide-charge combinations) for this protein group | `struct` | No |
| `feature_counts.unique_features` | Number of unique peptide features specific to this protein group within this file | `int` | -- |
| `feature_counts.total_features` | Total number of peptide features identified for this protein group in this file | `int` | -- |

### Quality

| Field | Description | Type | Required |
|-------|-------------|------|----------|
| `global_qvalue` | Global q-value of the protein group at the experiment level | `float64` | No |
| `pg_qvalue` | Protein group q-value at the run level (DIA-NN specific) | `float64` | No |
| `is_decoy` | Decoy indicator: 1 if decoy, 0 if target | `int32` | Yes |
| `contaminant` | Contaminant indicator: 1 if contaminant, 0 otherwise | `int32` | No |
| `sequence_coverage` | Percentage of the protein sequence covered by identified peptides | `float32` | No |
| `molecular_weight` | Molecular weight of the protein in kDa | `float32` | No |

### Quantification

| Field | Description | Type | Required |
|-------|-------------|------|----------|
| `intensities` | Primary intensity-based abundance of the protein group across channels. See [Intensities](intensities.md) | `array[struct]` | No |
| `additional_intensities` | Derived/processed intensity values (normalized, LFQ, iBAQ, etc.). See [Intensities](intensities.md) | `array[struct]` | No |
| `additional_scores` | Additional scores and metrics (posterior error probability, confidence, etc.). See [Scores](scores.md) | `array[struct]` | No |

Each entry in `intensities` contains:

| Sub-field | Description | Type |
|-----------|-------------|------|
| `sample_accession` | Sample accession from SDRF | `string` |
| `channel` | Labeling channel identifier | `string` |
| `intensity` | Raw intensity value | `float32` |

Each entry in `additional_intensities` contains:

| Sub-field | Description | Type |
|-----------|-------------|------|
| `sample_accession` | Sample accession from SDRF | `string` |
| `channel` | Labeling channel identifier | `string` |
| `intensities` | Array of name-value pairs for derived intensities | `array[struct{intensity_name, intensity_value}]` |

### Peptide detail

| Field | Description | Type | Required |
|-------|-------------|------|----------|
| `peptides` | Number of peptides per individual protein in the protein group | `array[struct]` | Yes |

Each entry in `peptides` contains:

| Sub-field | Description | Type |
|-----------|-------------|------|
| `protein_name` | Protein accession | `string` |
| `peptide_count` | Number of peptides for this protein | `int` |

### CV params

| Field | Description | Type | Required |
|-------|-------------|------|----------|
| `cv_params` | Optional list of controlled vocabulary parameters for additional metadata | `array[struct{cv_name, cv_value}]` | No |

## Example

```json
{
  "pg_accessions": ["P04217", "A0A024R4E5"],
  "pg_names": ["Alpha-1B-glycoprotein", "Alpha-1B-glycoprotein variant"],
  "gg_accessions": ["A1BG"],
  "gg_names": ["A1BG"],
  "anchor_protein": "P04217",
  "reference_file_name": "20230101_sample_01",
  "peptide_counts": {
    "unique_sequences": 12,
    "total_sequences": 18
  },
  "feature_counts": {
    "unique_features": 24,
    "total_features": 36
  },
  "global_qvalue": 0.001,
  "pg_qvalue": 0.005,
  "is_decoy": 0,
  "contaminant": 0,
  "sequence_coverage": 45.2,
  "molecular_weight": 54.3,
  "intensities": [
    {
      "sample_accession": "Sample-1",
      "channel": "TMT126",
      "intensity": 1.5e8
    }
  ],
  "additional_intensities": [
    {
      "sample_accession": "Sample-1",
      "channel": "TMT126",
      "intensities": [
        {"intensity_name": "LFQ", "intensity_value": 1.2e8},
        {"intensity_name": "iBAQ", "intensity_value": 3.4e7}
      ]
    }
  ],
  "peptides": [
    {"protein_name": "P04217", "peptide_count": 15},
    {"protein_name": "A0A024R4E5", "peptide_count": 3}
  ],
  "additional_scores": [
    {"score_name": "Posterior Error Probability", "score_value": 0.0001, "higher_better": false}
  ],
  "cv_params": null
}
```

## Notes

!!! note "Relationship to other views"
    The PG view provides per-file protein group quantification. For a simpler per-sample protein summary without inference details, see the [Protein Summary View](protein.md). For downstream absolute or differential expression results, see the [Absolute Expression](absolute.md) and [Differential Expression](differential.md) views.

!!! tip "Tool equivalences"
    - **MaxQuant**: Corresponds to `proteinGroups.txt`. Fields like `pg_accessions` map to "Protein IDs", `is_decoy` to "Reverse", and `contaminant` to "Potential contaminant".
    - **DIA-NN**: Corresponds to the protein group matrix. `pg_accessions` maps to "Protein.Group", `global_qvalue` to "Global.PG.Q.Value", and `intensities` to "PG.Quantity".
    - **FragPipe**: `pg_accessions` maps to "Group" + "Indistinguishable Proteins".

!!! warning "Primary key constraints"
    The combination of `pg_accessions` and `reference_file_name` forms the composite primary key. Both fields MUST NOT be null. Each record represents a single protein group observed in a single raw file.
