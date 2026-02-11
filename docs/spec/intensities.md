# Intensities

Intensity values in QPX are captured through two complementary fields: `intensities` for primary/raw measurements and `additional_intensities` for derived or processed values. This separation enforces a clean boundary between **experimental design** (which channels and samples were measured) and **data processing** (how raw values were normalized or transformed).

## Semantic separation

```mermaid
graph LR
    RAW["Instrument<br/>Raw signal"] --> INT["intensities<br/><i>Primary measurements</i>"]
    INT --> NORM["additional_intensities<br/><i>Normalized / LFQ / iBAQ</i>"]

    style RAW fill:#e1f5fe
    style INT fill:#e8f5e9
    style NORM fill:#fff3e0
```

| Field | Contains | Semantics |
| ----- | -------- | --------- |
| `intensities` | Raw intensity values from the instrument or quantification tool | Experimental design: channels, samples, raw measurements |
| `additional_intensities` | Normalized, LFQ, iBAQ, or other algorithm-derived values | Data processing: transformations applied to primary intensities |

This design means a data consumer can always find the raw measurement in `intensities` and any post-processed variants in `additional_intensities`, regardless of the tool that produced the file.

## Primary intensities

The `intensities` field is an array of structs. Each element represents one intensity measurement for a specific sample and channel within a single reference file.

### Struct definition

```
intensities: array[struct{
    sample_accession: string,  -- Sample identifier (typically the "source name" from SDRF)
    channel:          string,  -- Channel label (e.g. "LFQ", "TMT126", "iTRAQ114")
    intensity:        float    -- Raw intensity value
}]
```

### TMT example

In a TMT-10plex experiment, a single feature row contains one intensity per channel -- each channel maps to a different biological sample via the SDRF.

```json
[
  {"sample_accession": "Sample-1", "channel": "TMT126",  "intensity": 15234.7},
  {"sample_accession": "Sample-2", "channel": "TMT127N", "intensity": 18902.3},
  {"sample_accession": "Sample-3", "channel": "TMT127C", "intensity": 12045.1},
  {"sample_accession": "Sample-4", "channel": "TMT128N", "intensity": 22310.5},
  {"sample_accession": "Sample-5", "channel": "TMT128C", "intensity": 9871.6}
]
```

### LFQ example

In a label-free quantification experiment, each feature row has a single intensity value. The channel is set to `"LFQ"` by convention.

```json
[
  {"sample_accession": "Sample-1", "channel": "LFQ", "intensity": 98765.4}
]
```

!!! note
    For LFQ experiments, each feature row corresponds to one MS run and one sample, so the `intensities` array always has a single element. For multiplexed experiments (TMT, iTRAQ), the array contains one element per channel.

## Additional intensities

The `additional_intensities` field captures derived values that result from post-processing the primary measurements -- normalization, LFQ inference, iBAQ computation, and so on. Each element mirrors the sample/channel structure of primary intensities but adds a nested array of named intensity types.

### Struct definition

```
additional_intensities: array[struct{
    sample_accession: string,  -- Matches the corresponding primary intensity entry
    channel:          string,  -- Matches the corresponding primary intensity entry
    intensities:      array[struct{
        intensity_name:  string,  -- Algorithm or metric name (e.g. "normalize_intensity")
        intensity_value: float    -- Computed value
    }]
}]
```

### Example

```json
[
  {
    "sample_accession": "Sample-1",
    "channel": "LFQ",
    "intensities": [
      {"intensity_name": "normalize_intensity", "intensity_value": 0.1234},
      {"intensity_name": "lfq", "intensity_value": 23456.7},
      {"intensity_name": "ibaq", "intensity_value": 4567.8}
    ]
  },
  {
    "sample_accession": "Sample-2",
    "channel": "LFQ",
    "intensities": [
      {"intensity_name": "normalize_intensity", "intensity_value": 0.1456},
      {"intensity_name": "lfq", "intensity_value": 25890.2},
      {"intensity_name": "ibaq", "intensity_value": 5102.3}
    ]
  }
]
```

!!! tip
    Common `intensity_name` values include `normalize_intensity`, `lfq`, `ibaq`, and `maxlfq`. Tool-specific names are acceptable when no standard name exists.

## When to use which field

| Scenario | Use `intensities` | Use `additional_intensities` |
| -------- | ------------------ | ---------------------------- |
| Raw reporter ion signal from TMT channels | Yes | -- |
| Single MS1 precursor intensity in LFQ | Yes | -- |
| Median-normalized intensity | -- | Yes |
| MaxLFQ protein intensity | -- | Yes |
| iBAQ value | -- | Yes |
| Algorithm-specific rescaled intensity | -- | Yes |

## Where intensities are used

| View                          | `intensities` | `additional_intensities` | Notes |
| ----------------------------- | :-----------: | :----------------------: | ----- |
| Feature (`feature_file`)      | Yes           | Yes                      | Per-feature, per-run measurements |
| Protein Group (`pg_file`)     | Yes           | Yes                      | Aggregated at protein group level |
| Protein Summary (`protein_file`) | --         | --                       | Uses `abundance` field instead |

!!! note
    The Protein Summary view uses a simpler `abundance` field (a single float per sample) rather than the full intensities struct, since it serves as a lightweight report format.

## Further reading

- [Scores & CV Terms](scores.md) -- additional structured metadata for features and protein groups
- [QPX Format Overview](index.md) -- full list of views and concepts
