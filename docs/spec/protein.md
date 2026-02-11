# Protein Summary View

The protein summary view is a simplified report of proteins identified and quantified in the experiment. Unlike the [Protein Group View](pg.md), it does not contain detailed inference information about protein groups. Instead, it provides a concise per-sample protein abundance summary designed for fast reporting and integration with downstream expression formats.

## Use cases

- Fast reports of proteins quantified/identified in an experiment for web interfaces and search engines.
- Simple protein abundance summaries across samples without protein group inference details.
- Connection to [Absolute Expression](absolute.md) and [Differential Expression](differential.md) formats to provide context on protein identification coverage.

## Schema

See the full Avro schema in [`protein.avsc`](../protein.avsc).

| Field | Description | Type | Required |
|-------|-------------|------|----------|
| `protein_accessions` | Protein accession identifiers (UniProt) | `array[string]` | Yes |
| `abundance` | Abundance of the given protein in the sample/experiment | `float32` | No |
| `sample_accession` | Sample accession in the SDRF (`source name` column) | `string` | Yes |
| `global_qvalue` | Global q-value for the protein or protein group | `float32` | No |
| `is_decoy` | Decoy indicator: 1 if decoy, 0 if target | `int32` | No |
| `best_id_score` | Best search engine score for the identification | `struct{name, value}` | No |
| `additional_scores` | List of subsidiary identification scores | `array[struct{name, value}]` | No |
| `gg_accessions` | Gene accessions corresponding to every protein | `array[string]` | No |
| `gg_names` | Gene names corresponding to every protein | `array[string]` | No |
| `number_peptides` | Total number of peptides for the protein | `int32` | No |
| `number_psms` | Total number of peptide spectrum matches | `int32` | No |
| `number_unique_peptides` | Total number of unique peptides | `int32` | No |

### Score structures

The `best_id_score` field is a struct with:

| Sub-field | Description | Type |
|-----------|-------------|------|
| `name` | Name of the score type (e.g., "Mascot score", "Andromeda score") | `string` |
| `value` | Score value | `float32` |

Each entry in `additional_scores` has the same structure.

## Example

```json
{
  "protein_accessions": ["P04217"],
  "abundance": 1.5e8,
  "sample_accession": "Sample-1",
  "global_qvalue": 0.001,
  "is_decoy": 0,
  "best_id_score": {
    "name": "Posterior Error Probability",
    "value": 0.0001
  },
  "additional_scores": [
    {"name": "Andromeda score", "value": 125.4}
  ],
  "gg_accessions": ["ENSG00000121410"],
  "gg_names": ["A1BG"],
  "number_peptides": 18,
  "number_psms": 42,
  "number_unique_peptides": 12
}
```

## Notes

!!! note "Simpler than the PG view"
    The protein summary view intentionally omits protein group inference details (e.g., `peptide_counts` struct, `feature_counts` struct, `anchor_protein`, `contaminant` flag). It is designed for fast lookups and reporting, not for detailed protein group analysis. Use the [Protein Group View](pg.md) when you need full inference context.

!!! tip "Connection to expression formats"
    The protein summary view serves as a bridge between identification-level data and expression-level data. The `protein_accessions` field links to the `protein` field in both the [Absolute Expression](absolute.md) and [Differential Expression](differential.md) views.
