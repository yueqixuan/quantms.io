# Peptide View

The peptide view is a summary representation of quantified peptides by sample. Each row represents a single peptide (modified or unmodified) with its abundance in a given sample, along with gene annotations and the best identification score observed across all features or PSMs.

## Use Cases

- **Quick reports**: Provides a concise summary of all peptides quantified in an experiment, suitable for dashboards and summary tables.
- **Peptide abundance summaries**: One row per peptide per sample, with a single abundance value, making it straightforward to aggregate and compare across conditions.
- **Web interfaces**: The simple flat structure (no nested arrays for intensities or spectral data) makes this view well-suited for search engines, REST APIs, and web-based data browsers.

## Schema

| Field | Description | Type | Required |
|-------|-------------|------|----------|
| `sequence` | Unmodified peptide amino acid sequence | string | yes |
| `peptidoform` | Peptide sequence with modifications in ProForma notation | string | yes |
| `modifications` | Structured list of modifications with name, accession, position, and localization scores | array[struct], null | no |
| `gg_accessions` | Gene group accessions | array[string], null | no |
| `gg_names` | Gene names | array[string], null | no |
| `best_id_score` | Best search engine score from all features/PSMs identified for this peptide | array[struct], null | no |
| `sample_accession` | Sample accession in the SDRF (the `source name` column) | string, null | no |
| `abundance` | Peptide abundance in the given sample | float32, null | no |

!!! note "Score structure"
    The `best_id_score` field uses the same score structure as `additional_scores` in other views: each entry contains `score_name` (string), `score_value` (float64), and `higher_better` (boolean, nullable).

## Shared Fields

Several fields in the peptide view use structures shared across other QPX views:

- For details on the `modifications` field structure, see [Modifications](modifications.md).
- For details on score structures used in `best_id_score`, see [Scores](scores.md).

## Example

```json
{
  "sequence": "AADLLTSFLGHK",
  "peptidoform": "AADLLTSFLGHK",
  "modifications": null,
  "gg_accessions": ["ENSG00000121410"],
  "gg_names": ["A1BG"],
  "best_id_score": [
    {
      "score_name": "global_qvalue",
      "score_value": 0.0005,
      "higher_better": false
    },
    {
      "score_name": "posterior_error_probability",
      "score_value": 1.23e-10,
      "higher_better": false
    }
  ],
  "sample_accession": "Sample-1",
  "abundance": 15234.5
}
```

### Peptide with Modifications

```json
{
  "sequence": "VLHPLEGAVVIIFK",
  "peptidoform": "[UniMod:1]-VLHPLEGAVVIIFK",
  "modifications": [
    {
      "name": "Acetyl",
      "accession": "UniMod:1",
      "positions": [
        {"position": "N-term.0", "scores": null}
      ]
    }
  ],
  "gg_accessions": ["ENSG00000244734"],
  "gg_names": ["HBB"],
  "best_id_score": [
    {
      "score_name": "andromeda_score",
      "score_value": 142.56,
      "higher_better": true
    }
  ],
  "sample_accession": "Sample-3",
  "abundance": 8923.1
}
```

## File Metadata

Peptide Parquet files store file-level metadata as key-value pairs in the Parquet footer. The following metadata fields are defined:

| Field | Description |
|-------|-------------|
| `qpx_version` | Version of the QPX format used to generate the file |
| `software_provider` | Name and version of the software that generated the data |
| `creator` | Name of the tool or person who created the file |
| `file_type` | Type of the file (value: `peptide_file`) |
| `creation_date` | Date when the file was created |
| `uuid` | Unique identifier for the file |

!!! note "No scan_format field"
    Unlike the PSM and feature views, the peptide view does not include a `scan_format` metadata field because it does not reference individual scans.

## Notes

- **Summary view**: The peptide view is intentionally simpler than the PSM or feature views. It does not contain spectral data, individual scan references, or per-channel intensity breakdowns. Its purpose is to provide a quick, flat summary of peptide quantification across samples.
- **Relationship to feature view**: The peptide view can be derived from the [Feature View](feature.md) by aggregating features across files and channels into a single abundance value per sample. The feature view retains the full detail (per-file, per-channel intensities, scan references, protein group q-values), while the peptide view provides the rolled-up summary.
- **Relationship to PSM view**: The `best_id_score` in the peptide view represents the best score observed across all PSMs or features for that peptide, collapsing detailed spectrum-level information into a single confidence metric.
- **Granularity**: A peptide row can represent either a modified peptide (peptidoform) or an unmodified sequence, depending on the use case and the level of granularity chosen during generation. The `peptidoform` field always carries the modification information regardless.
