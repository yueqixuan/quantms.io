# PSM View

The PSM (Peptide Spectrum Match) view captures spectrum-level identification results. Each row represents a single match between a mass spectrum and a peptide sequence, including the identification scores, optional spectral data, and protein mappings.

## Use Cases

- **AI/ML training**: Provides peptide-spectrum pairs with optional spectral arrays (m/z, intensity, ion types) for training intensity prediction, de novo sequencing, and clustering models.
- **Spectrum-level analysis**: Enables detailed inspection of individual identifications, including retention time, charge state, and search engine scores.
- **DDA identification results**: Designed primarily for data-dependent acquisition (DDA) workflows where each spectrum yields one or more peptide identifications.

## Schema

### Core Identification Fields

| Field | Description | Type | Required |
|-------|-------------|------|----------|
| `sequence` | Unmodified peptide amino acid sequence | string | yes |
| `peptidoform` | Peptide sequence with modifications in ProForma notation | string | yes |
| `modifications` | Structured list of modifications with name, accession, position, and localization scores | array[struct], null | no |
| `precursor_charge` | Charge state of the precursor ion | int32 | yes |
| `posterior_error_probability` | Posterior error probability (PEP) for the peptide-spectrum match | float64, null | no |
| `is_decoy` | Decoy indicator: 1 if the PSM is a decoy, 0 if target | int32 | yes |
| `calculated_mz` | Theoretical peptide mass-to-charge ratio based on identified sequence and modifications | float32 | yes |
| `observed_mz` | Experimental peptide mass-to-charge ratio (in Da) | float32 | yes |
| `rt` | MS2 scan precursor retention time (in seconds) | float32, null | no |
| `predicted_rt` | Predicted retention time of the peptide (in seconds) | float32, null | no |
| `reference_file_name` | Spectrum file name without path or extension | string | yes |
| `scan` | Scan number or index of the identified spectrum | string | yes |
| `additional_scores` | List of score structures with name, value, and direction indicator | array[struct], null | no |
| `cv_params` | Optional list of controlled vocabulary parameters for additional metadata | array[struct], null | no |

### Protein Mapping Fields

| Field | Description | Type | Required |
|-------|-------------|------|----------|
| `protein_accessions` | Protein accessions of all proteins that the peptide maps to | array[string], null | no |

### Spectral Data Fields (optional)

These fields are optional and may not exist in the file at all. They are included when spectral data is requested during conversion (e.g., using the `--spectral-data` flag). When present, they enable AI/ML workflows that require direct access to spectral information alongside identifications.

| Field | Description | Type | Required |
|-------|-------------|------|----------|
| `ion_mobility` | Ion mobility value for the precursor ion | float32, null | no |
| `number_peaks` | Number of peaks in the spectrum | int32, null | no |
| `mz_array` | Array of m/z values for the spectrum | array[float32], null | no |
| `intensity_array` | Array of intensity values for the spectrum | array[float32], null | no |
| `charge_array` | Array of fragment ion charge values | array[int32], null | no |
| `ion_type_array` | Array of fragment ion type annotations (e.g., b, y, a) | array[string], null | no |
| `ion_mobility_array` | Array of fragment ion mobility values | array[float32], null | no |

!!! note "Nullable vs Optional"
    Core fields marked as "not required" are **nullable** -- the column always exists in the file but individual values may be null. Spectral data fields are **optional** -- the column itself may be absent from the file entirely, depending on conversion settings.

## Shared Fields

Several fields in the PSM view use structures shared across other QPX views:

- For details on the `modifications` field structure, see [Modifications](modifications.md).
- For details on `additional_scores` and score semantics, see [Scores](scores.md).
- For details on `cv_params` usage and recommended terms, see [CV Parameters](cv-params.md).

## Example

### Basic PSM Record

```json
{
  "sequence": "AAAAAAAAAAGAAGGR",
  "peptidoform": "_(Acetyl (Protein N-term))AAAAAAAAAAGAAGGR_",
  "precursor_charge": 2,
  "scan": "42164",
  "rt": 5140.98,
  "calculated_mz": 635.3311,
  "observed_mz": 635.3315,
  "is_decoy": 0,
  "posterior_error_probability": 5.58e-20,
  "predicted_rt": null,
  "reference_file_name": "20200101_sample_A",
  "protein_accessions": ["Q86U42-2", "Q86U42"],
  "modifications": [
    {
      "name": "Acetyl",
      "accession": "UniMod:1",
      "positions": [
        {
          "position": "N-term.0",
          "scores": []
        }
      ]
    }
  ],
  "additional_scores": [
    {
      "score_name": "andromeda_score",
      "score_value": 175.73,
      "higher_better": true
    },
    {
      "score_name": "andromeda_delta_score",
      "score_value": 160.47,
      "higher_better": true
    },
    {
      "score_name": "parent_ion_fraction",
      "score_value": 0.0,
      "higher_better": true
    }
  ],
  "cv_params": [
    {"cv_name": "dissociation method", "cv_value": "HCD"},
    {"cv_name": "normalized collision energy", "cv_value": "28"}
  ]
}
```

### PSM with Spectral Data

When spectral arrays are included, the record also contains peak-level data:

```json
{
  "sequence": "AAAAAAAAAAGAAGGR",
  "peptidoform": "_(Acetyl (Protein N-term))AAAAAAAAAAGAAGGR_",
  "precursor_charge": 2,
  "scan": "42164",
  "rt": 5140.98,
  "calculated_mz": 635.3311,
  "observed_mz": 635.3315,
  "is_decoy": 0,
  "reference_file_name": "20200101_sample_A",
  "number_peaks": 21,
  "mz_array": [175.119, 289.163, 360.200, 431.236, 488.258],
  "intensity_array": [1234.5, 5678.9, 3456.7, 2345.6, 1234.5],
  "charge_array": [1, 1, 1, 1, 1],
  "ion_type_array": ["y1", "y2", "b3", "b4", "b5"]
}
```

## File Metadata

PSM Parquet files store file-level metadata as key-value pairs in the Parquet footer. The following metadata fields are defined:

| Field | Description |
|-------|-------------|
| `qpx_version` | Version of the QPX format used to generate the file |
| `software_provider` | Name and version of the software that generated the data |
| `scan_format` | Format of scan identifiers: `scan`, `index`, `nativeId`, or `multiple` |
| `creator` | Name of the tool or person who created the file |
| `file_type` | Type of the file (value: `psm_file`) |
| `creation_date` | Date when the file was created |
| `uuid` | Unique identifier for the file |
| `compression_format` | Compression algorithm used: `gzip`, `snappy`, `lzo`, or `none` |

!!! tip "Reading file metadata in Python"
    ```python
    import pyarrow.parquet as pq

    parquet_file = pq.ParquetFile("experiment.psm.parquet")
    metadata = parquet_file.schema_arrow.metadata
    for key, value in metadata.items():
        print(f"{key.decode()}: {value.decode()}")
    ```

## Notes

!!! warning "DDA-specific view"
    The PSM view is designed primarily for **DDA** (data-dependent acquisition) methods. It is **not recommended** for DIA experiments, where the feature view should be used instead. Generating a PSM file for DIA data would produce duplicated information relative to the feature view.

- **Relationship to feature view**: The PSM view captures individual spectrum matches, while the [Feature View](feature.md) aggregates these into quantified peptide features with intensity data. A single feature may correspond to multiple PSMs across different scans.
- **Protein inference**: Protein inference results should not be the primary focus of the PSM view. However, `protein_accessions` is provided for use cases such as peptide filtering and protein-level browsing. For full protein group information, use the protein group (PG) view.
- **Spectral arrays**: The `mz_array` and `intensity_array` are parallel arrays of the same length, matching the `number_peaks` value. For large-scale spectral storage, the dedicated mass spectra (mz) view is recommended.
- **Recommended additional scores**: `global_qvalue` (experiment-level PSM q-value), `rank` (peptide rank in search results), and `pg_global_qvalue` (protein group q-value used for filtering).
