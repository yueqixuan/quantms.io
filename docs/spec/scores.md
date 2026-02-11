# Scores & CV Terms

QPX provides two structured fields for attaching metadata to records beyond the core columns: `additional_scores` for numeric quality metrics and `cv_params` for controlled vocabulary annotations. Both are stored as arrays of structs, making them extensible without schema changes.

## Additional scores

The `additional_scores` field captures search engine scores, quality metrics, and other numeric values associated with a PSM, feature, peptide, or protein group. Each entry is self-describing: it carries its name, value, and an indication of whether higher values are better.

### Struct definition

```
additional_scores: array[struct{
    score_name:    string,       -- Score identifier (e.g. "Comet:xcorr")
    score_value:   float,        -- Numeric value
    higher_better: bool, null    -- true = higher is better; false = lower is better; null = unknown
}]
```

The `higher_better` field makes the data self-describing so that downstream consumers can interpret scores without looking up external ontology definitions.

### Example

```json
[
  {"score_name": "Comet:xcorr",    "score_value": 3.42,   "higher_better": true},
  {"score_name": "global_qvalue",  "score_value": 0.0012, "higher_better": false},
  {"score_name": "rank",           "score_value": 1.0,    "higher_better": false}
]
```

### Recommended score names

The following score names are commonly used across QPX views. It is recommended to use HUPO-PSI MS ontology terms where possible.

| Score name | Typical view(s) | Direction | Description |
| ---------- | ---------------- | --------- | ----------- |
| `global_qvalue` | PSM, Feature | lower is better | Global q-value at the experiment level |
| `pg_global_qvalue` | PSM, Feature | lower is better | Protein group global q-value used to filter at the protein group level |
| `rank` | PSM | lower is better | Rank of the peptide in the search engine results (1 = best) |
| `Comet:xcorr` | PSM | higher is better | Cross-correlation score from the Comet search engine |
| `Comet:deltacn` | PSM | higher is better | Delta CN score from Comet |
| `Comet:expect` | PSM | lower is better | Expectation value from Comet |
| `MSGF:RawScore` | PSM | higher is better | Raw score from MS-GF+ |
| `MSGF:SpecEValue` | PSM | lower is better | Spectral E-value from MS-GF+ |
| `Sage:hyperscore` | PSM | higher is better | Hyperscore from Sage search engine |
| `DIA-NN:Q.Value` | Feature | lower is better | Run-level q-value from DIA-NN |
| `DIA-NN:Global.Q.Value` | Feature | lower is better | Global q-value from DIA-NN |
| `DIA-NN:CScore` | Feature | higher is better | Confidence score from DIA-NN |
| `consensus_support` | PSM | higher is better | Number of search engines supporting the identification |

!!! tip
    When using tool-specific scores, prefix the score name with the tool name and a colon (e.g., `Comet:xcorr`, `DIA-NN:Q.Value`). This avoids name collisions between different tools.

### Protein-level additional scores

At the protein group level, `additional_scores` follows the same struct definition. The `score_value` is a float that applies to the entire protein group entry. For per-protein scores within a group, the values array index corresponds to the `pg_accessions` array index.

## Controlled vocabulary terms

The `cv_params` field stores key-value annotations drawn from controlled vocabularies. Unlike `additional_scores` (which are always numeric), CV parameters can carry string values or be value-less (presence-only).

### Struct definition

```
cv_params: array[struct{
    cv_name:  string,        -- Term name from a controlled vocabulary
    cv_value: string, null   -- Term value; null if the term is a flag (presence-only)
}]
```

### Example

```json
[
  {"cv_name": "ms level",          "cv_value": "2"},
  {"cv_name": "deconvoluted data", "cv_value": null},
  {"cv_name": "prot:FDR threshold", "cv_value": "0.01"},
  {"cv_name": "number of unmatched peaks", "cv_value": "3"}
]
```

!!! note
    The `cv_name` is always required. The `cv_value` is optional -- a `null` value indicates a boolean/flag-style term where mere presence is meaningful (e.g., "deconvoluted data").

## Where these fields are used

| View | `additional_scores` | `cv_params` |
| ---- | :-----------------: | :---------: |
| PSM (`psm_file`) | Yes | Yes |
| Feature (`feature_file`) | Yes | Yes |
| Peptide (`peptide_file`) | Yes (as `best_id_score`) | -- |
| Protein Group (`pg_file`) | Yes | -- |
| MZ (`mz_file`) | -- | Yes |

!!! warning
    Do not store search engine scores in `cv_params`. Numeric scores belong in `additional_scores`, which provides the `higher_better` flag and enforces float typing. Reserve `cv_params` for non-numeric annotations and metadata.

## Further reading

- [Modifications](modifications.md) -- localization scores use a similar score struct pattern
- [Intensities](intensities.md) -- additional processed values for features and protein groups
- [QPX Format Overview](index.md) -- full list of views and concepts
