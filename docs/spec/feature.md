# Feature View

The feature view captures quantified peptide information at the MS run level. Each row represents a peptide feature -- a quantified peptidoform in a specific reference file -- including its intensity across samples and channels, protein group mappings, and gene annotations.

## Use Cases

- **Quantified peptide information**: Stores peptide intensities linked to sample metadata, enabling downstream quantitative analysis and integration with SDRF annotations.
- **Peptide-level statistics**: Enables algorithms that operate at the peptide level (e.g., peptide-to-protein rollup, normalization, missing value imputation).
- **Integration with sample metadata**: Each feature carries sample accession and channel information, connecting quantification data to experimental design described in the SDRF.

## Schema

### Core Identification Fields

These fields are shared with the PSM view and describe the peptide identification associated with the feature.

| Field | Description | Type | Required |
|-------|-------------|------|----------|
| `sequence` | Unmodified peptide amino acid sequence | string | yes |
| `peptidoform` | Peptide sequence with modifications in ProForma notation | string | yes |
| `modifications` | Structured list of modifications with name, accession, position, and localization scores | array[struct], null | no |
| `precursor_charge` | Charge state of the precursor ion | int32 | yes |
| `posterior_error_probability` | Posterior error probability (PEP) for the peptide or PSM match | float64, null | no |
| `is_decoy` | Decoy indicator: 1 if the peptide is a decoy, 0 if target | int32 | yes |
| `calculated_mz` | Theoretical peptide mass-to-charge ratio based on identified sequence and modifications | float32 | yes |
| `observed_mz` | Experimental peptide mass-to-charge ratio (in Da) | float32 | yes |
| `rt` | Precursor retention time (in seconds) | float32, null | no |
| `rt_start` | Start of the retention time window for the feature | float32, null | no |
| `rt_stop` | End of the retention time window for the feature | float32, null | no |
| `predicted_rt` | Predicted retention time of the peptide (in seconds) | float32, null | no |
| `ion_mobility` | Ion mobility value for the precursor ion | float32, null | no |
| `start_ion_mobility` | Start ion mobility value for the precursor ion | float32, null | no |
| `stop_ion_mobility` | Stop ion mobility value for the precursor ion | float32, null | no |
| `additional_scores` | List of score structures with name, value, and direction indicator | array[struct], null | no |
| `cv_params` | Optional list of controlled vocabulary parameters for additional metadata | array[struct], null | no |

### Quantification Fields

| Field | Description | Type | Required |
|-------|-------------|------|----------|
| `intensities` | Primary intensity-based abundance of the feature across samples and channels | array[struct], null | no |
| `additional_intensities` | Derived/processed intensity values (e.g., normalized, LFQ, iBAQ) as named key-value pairs per sample and channel | array[struct], null | no |
| `reference_file_name` | The reference file name that contains the feature | string | yes |

!!! tip "Intensity structure"
    For details on the `intensities` and `additional_intensities` data structures, including examples for LFQ and TMT experiments, see [Intensities](intensities.md).

### Protein and Gene Fields

| Field | Description | Type | Required |
|-------|-------------|------|----------|
| `pg_accessions` | Protein group accessions of all proteins that the peptide maps to | array[string], null | no |
| `anchor_protein` | One protein accession that represents the protein group | string | yes |
| `unique` | Unique peptide indicator: 1 if the peptide maps to a single protein, 0 otherwise | int32, null | no |
| `pg_global_qvalue` | Global q-value of the protein group at the experiment level | float64, null | no |
| `gg_accessions` | Gene group accessions | array[string], null | no |
| `gg_names` | Gene names | array[string], null | no |

### Spectra Reference

| Field | Description | Type | Required |
|-------|-------------|------|----------|
| `scan_reference_file_name` | The reference file containing the best PSM that identified the feature (may differ from `reference_file_name`) | string, null | no |
| `scan` | Scan number or index of the spectrum used to identify the feature | string | yes |

## Shared Fields

Several fields in the feature view use structures shared across other QPX views:

- For details on the `modifications` field structure, see [Modifications](modifications.md).
- For details on `intensities` and `additional_intensities`, see [Intensities](intensities.md).
- For details on `additional_scores` and score semantics, see [Scores](scores.md).
- For details on `cv_params` usage and recommended terms, see [CV Parameters](cv-params.md).

## Example

### Feature with TMT Intensities

```json
{
  "sequence": "AADLLTSFLGHK",
  "peptidoform": "AADLLTSFLGHK",
  "modifications": null,
  "precursor_charge": 3,
  "posterior_error_probability": 1.23e-10,
  "is_decoy": 0,
  "calculated_mz": 424.2345,
  "observed_mz": 424.2350,
  "rt": 2345.67,
  "rt_start": 2340.12,
  "rt_stop": 2351.23,
  "predicted_rt": 2348.00,
  "ion_mobility": null,
  "start_ion_mobility": null,
  "stop_ion_mobility": null,
  "reference_file_name": "20200101_TMT_fraction01",
  "scan": "15234",
  "scan_reference_file_name": "20200101_TMT_fraction01",
  "pg_accessions": ["P04217", "P04217-2"],
  "anchor_protein": "P04217",
  "unique": 1,
  "pg_global_qvalue": 0.001,
  "gg_accessions": ["ENSG00000121410"],
  "gg_names": ["A1BG"],
  "intensities": [
    {"sample_accession": "Sample-1", "channel": "TMT126", "intensity": 15234.5},
    {"sample_accession": "Sample-2", "channel": "TMT127N", "intensity": 18456.7},
    {"sample_accession": "Sample-3", "channel": "TMT127C", "intensity": 12890.3},
    {"sample_accession": "Sample-4", "channel": "TMT128N", "intensity": 21045.8}
  ],
  "additional_intensities": [
    {
      "sample_accession": "Sample-1",
      "channel": "TMT126",
      "intensities": [
        {"intensity_name": "normalize_intensity", "intensity_value": 0.1523},
        {"intensity_name": "ibaq", "intensity_value": 4567.8}
      ]
    },
    {
      "sample_accession": "Sample-2",
      "channel": "TMT127N",
      "intensities": [
        {"intensity_name": "normalize_intensity", "intensity_value": 0.1846},
        {"intensity_name": "ibaq", "intensity_value": 5432.1}
      ]
    }
  ],
  "additional_scores": [
    {
      "score_name": "global_qvalue",
      "score_value": 0.0012,
      "higher_better": false
    }
  ],
  "cv_params": null
}
```

### Feature with LFQ Intensity

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
  "precursor_charge": 2,
  "is_decoy": 0,
  "calculated_mz": 782.4721,
  "observed_mz": 782.4725,
  "rt": 3567.89,
  "reference_file_name": "20200101_LFQ_rep1",
  "scan": "28901",
  "anchor_protein": "P68871",
  "pg_accessions": ["P68871"],
  "unique": 1,
  "gg_names": ["HBB"],
  "intensities": [
    {"sample_accession": "Sample-1", "channel": "LFQ", "intensity": 98765.4}
  ],
  "additional_intensities": null
}
```

## File Metadata

Feature Parquet files store file-level metadata as key-value pairs in the Parquet footer. The following metadata fields are defined:

| Field | Description |
|-------|-------------|
| `qpx_version` | Version of the QPX format used to generate the file |
| `software_provider` | Name and version of the software that generated the data |
| `scan_format` | Format of scan identifiers: `scan`, `index`, `nativeId`, or `multiple` |
| `creator` | Name of the tool or person who created the file |
| `file_type` | Type of the file (value: `feature_file`) |
| `creation_date` | Date when the file was created |
| `uuid` | Unique identifier for the file |
| `compression_format` | Compression algorithm used: `gzip`, `snappy`, `lzo`, or `none` |

!!! tip "Reading file metadata in Python"
    ```python
    import pyarrow.parquet as pq

    parquet_file = pq.ParquetFile("experiment.feature.parquet")
    metadata = parquet_file.schema_arrow.metadata
    for key, value in metadata.items():
        print(f"{key.decode()}: {value.decode()}")
    ```

## Notes

!!! note "Works for both DDA and DIA"
    Unlike the PSM view (which is DDA-specific), the feature view supports both DDA and DIA workflows. For DIA experiments, the feature view is the primary peptide-level output.

- **Relationship to PSM view**: A feature aggregates one or more PSMs into a single quantified peptide entry. The `scan` and `scan_reference_file_name` fields link back to the best PSM that identified the feature. In DDA-LFQ workflows, both views may exist; in DIA workflows, only the feature view is typically generated.
- **Relationship to protein group view**: The feature view contains protein group mappings (`pg_accessions`, `anchor_protein`, `pg_global_qvalue`) that connect each peptide to its inferred protein groups. The protein group (PG) view provides the complete protein-level perspective with aggregated intensities and peptide counts.
- **Intensities vs additional_intensities**: Use `intensities` for raw/primary measurements across experimental channels (TMT/iTRAQ tags) or samples (LFQ). Use `additional_intensities` for derived values such as normalized intensities, LFQ intensities, or iBAQ values. This separation keeps experimental design aspects distinct from data processing aspects.
- **Protein group accessions**: The `pg_accessions` field should contain all proteins within a protein group. The `anchor_protein` is the representative protein selected by the search engine or inference algorithm to represent the group.
