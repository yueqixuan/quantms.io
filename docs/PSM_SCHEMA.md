# PSM Schema Documentation and Examples

## Summary

This PR provides refined PSM (Peptide Spectrum Match) schema documentation and examples for multiple configurations, enabling team discussion and feedback before release.

## What's Included

### 1. Schema Review

- Verified consistency between `docs/psm.avsc` (Avro schema) and `qpx/core/format.py` (PyArrow schema)
- All 22 PSM fields have correct types, nullable settings, and descriptions
- Note: Minor field order difference in `protein_accessions` position (does not affect functionality)

### 2. Documentation

The PSM section in `docs/format-specification.md` (L696-851) is complete with:

- Complete field description tables
- Tool mapping (MaxQuant, DIA-NN, FragPipe, mzTab)
- Use case descriptions (DDA/DIA, with/without spectral data)
- File metadata documentation

### 3. Example Generation (Demonstration)

The following MaxQuant PSM examples can be generated to demonstrate the conversion workflow:

| Configuration | Description           | Command                                                                             |
| ------------- | --------------------- | ----------------------------------------------------------------------------------- |
| Basic PSM     | Without spectral data | `qpxc convert maxquant-psm --msms-file msms.txt --output-folder ./`                 |
| PSM + Spectra | With spectral arrays  | `qpxc convert maxquant-psm --msms-file msms.txt --output-folder ./ --spectral-data` |

See "How to Generate Examples" section below for complete commands.

### 4. MaxQuant to PSM Conversion Verification

Verified the conversion from MaxQuant `msms.txt` to PSM parquet format:

- **Field Mapping**: All 22 fields correctly mapped via `MAXQUANT_PSM_MAP`
- **Modifications**: Correctly parsed with name, accession, positions structure
- **Additional Scores**: Correctly collected andromeda_score, andromeda_delta_score, parent_ion_fraction
- **Spectral Arrays**: mz_array, intensity_array, charge_array, ion_type_array all correctly processed

## PSM Schema Fields (22 total)

### Field Classification

Fields are classified as:

- **PK** (Primary Key): Must not be null, required for data integrity
- **nullable**: Column always exists, but value can be null
- **optional**: Column may not exist in file at all

### Core Identification Fields (14)

| Field                       | Type         | Classification | Description                                    |
| --------------------------- | ------------ | -------------- | ---------------------------------------------- |
| sequence                    | string       | PK             | Unmodified peptide sequence                    |
| peptidoform                 | string       | PK             | Peptide sequence with modifications (ProForma) |
| modifications               | list[struct] | nullable       | Modification details with positions and scores |
| precursor_charge            | int32        | PK             | Precursor ion charge                           |
| posterior_error_probability | float32      | nullable       | PEP value                                      |
| is_decoy                    | int32        | required       | Decoy indicator (1=decoy, 0=target)            |
| calculated_mz               | float32      | required       | Theoretical m/z                                |
| observed_mz                 | float32      | required       | Experimental m/z                               |
| additional_scores           | list[struct] | nullable       | Search engine scores                           |
| predicted_rt                | float32      | nullable       | Predicted retention time (seconds)             |
| reference_file_name         | string       | PK             | Reference file name                            |
| cv_params                   | list[struct] | nullable       | CV parameters                                  |
| scan                        | string       | PK             | Scan number                                    |
| rt                          | float32      | nullable       | MS2 retention time (seconds)                   |

### Protein Mapping Fields (1)

| Field              | Type         | Classification | Description        |
| ------------------ | ------------ | -------------- | ------------------ |
| protein_accessions | list[string] | nullable       | Protein accessions |

### Spectral Data Fields (7) - Optional

**Note**: These fields are **optional** - they may not exist in the file at all. Use `--spectral-data` flag during conversion to include these columns.

| Field              | Type          | Classification | Description                          |
| ------------------ | ------------- | -------------- | ------------------------------------ |
| ion_mobility       | float32       | optional       | Ion mobility value                   |
| number_peaks       | int32         | optional       | Number of peaks                      |
| mz_array           | list[float32] | optional       | m/z values array                     |
| intensity_array    | list[float32] | optional       | Intensity values array               |
| charge_array       | list[int32]   | optional       | Fragment ion charge array            |
| ion_type_array     | list[string]  | optional       | Ion type annotations (b, y, a, etc.) |
| ion_mobility_array | list[float32] | optional       | Fragment ion mobility array          |

**Nullable vs Optional**: These fields are *optional* (column may be absent), not just *nullable* (column exists with null values). See [Field Classification](format-specification.md#field-classification) for details.

## How to Generate Examples

```bash
# Basic PSM (without spectral data)
qpxc convert maxquant-psm \
    --msms-file tests/examples/maxquant/maxquant_simple/msms.txt \
    --output-folder examples/psm/ \
    --output-prefix maxquant-dda \
    --verbose

# PSM with spectral data
qpxc convert maxquant-psm \
    --msms-file tests/examples/maxquant/maxquant_simple/msms.txt \
    --output-folder examples/psm/ \
    --output-prefix maxquant-dda-spectra \
    --spectral-data \
    --verbose
```

## Example Files

Pre-generated example files are available in `examples/psm/`:

| File                           | Source   | Description                     |
| ------------------------------ | -------- | ------------------------------- |
| `maxquant-basic.psm.parquet`   | MaxQuant | Basic PSM without spectral data |
| `maxquant-spectra.psm.parquet` | MaxQuant | PSM with spectral arrays        |
| `idxml-basic.psm.parquet`      | OpenMS   | IdXML converted PSM             |
| `quantms-lfq.psm.parquet`      | quantms  | mzTab LFQ converted PSM         |

## Example Data Snippets

### Basic PSM Record

```json
{
  "sequence": "AAAAAAAAAAGAAGGR",
  "peptidoform": "_(Acetyl (Protein N-term))AAAAAAAAAAGAAGGR_",
  "precursor_charge": 2,
  "scan": "42164",
  "rt": 5140.98,
  "posterior_error_probability": 5.58e-20,
  "protein_accessions": ["Q86U42-2", "Q86U42"],
  "modifications": [
    {
      "name": "Acetyl",
      "accession": "UniMod:1",
      "positions": [{ "position": "N-term.0", "scores": [] }]
    }
  ],
  "additional_scores": [
    { "score_name": "andromeda_score", "score_value": 175.73, "higher_better": true },
    { "score_name": "andromeda_delta_score", "score_value": 160.47, "higher_better": true },
    { "score_name": "parent_ion_fraction", "score_value": 0.0, "higher_better": true }
  ]
}
```

### PSM with Spectral Data

When `--spectral-data` is enabled, spectral arrays are populated:

```json
{
  "number_peaks": 21,
  "mz_array": [175.119, 289.163, 360.200, 431.236, 488.258, ...],
  "intensity_array": [1234.5, 5678.9, ...],
  "charge_array": [1, 1, 1, ...],
  "ion_type_array": ["y1", "y2", "b3", ...]
}
```

### PSM with Fragmentation Information (via cv_params)

Fragmentation method and collision energy should be stored as CV terms in the `cv_params` field:

```json
{
  "sequence": "PEPTIDESEQ",
  "cv_params": [
    { "cv_name": "dissociation method", "cv_value": "HCD" },
    { "cv_name": "normalized collision energy", "cv_value": "28" }
  ]
}
```

**Why human-readable names?** QPX uses readable term names (like `dissociation method`) instead of ontology accessions (like `MS:1000044`) to align with successful omics formats such as GTF and AnnData. This makes data self-documenting while the specification provides formal definitions. See [Design Philosophy](format-specification.md#psm-cv-params) for details.

**Common fragmentation CV terms:**

| CV Name | Example Values |
|---------|----------------|
| dissociation method | HCD, CID, ETD, ECD, UVPD |
| collision energy | 28, 35 (in eV) |
| normalized collision energy | 28, 30 (percentage) |

Full reference with PSI-MS accessions: [Common CV Terms](format-specification.md#common-cv-terms)

### De Novo / No Protein Association

For De Novo analysis where PSMs have no protein mapping, `protein_accessions` is nullable:

```json
{
  "sequence": "PEPTIDESEQ",
  "protein_accessions": null
}
```

## How to Read Examples

```python
import pyarrow.parquet as pq

# Read PSM file
table = pq.read_table("examples/psm/maxquant-basic.psm.parquet")
df = table.to_pandas()

# View sample data
print(df[['sequence', 'peptidoform', 'precursor_charge', 'scan']].head())

# View modifications
print(df['modifications'].head())

# View additional scores
print(df['additional_scores'].head())
```

## Related Files

```
docs/
|-- psm.avsc                    # Avro schema definition
+-- format-specification.md     # Format specification (PSM: L700-851)

qpx/core/
|-- format.py                   # PyArrow schema (PEPTIDE_FIELDS, PSM_UNIQUE_FIELDS)
|-- common.py                   # PSM_SCHEMA, MAXQUANT_PSM_MAP
+-- maxquant/maxquant.py        # MaxQuant converter (process_psm_file)

tests/examples/maxquant/
|-- maxquant_simple/msms.txt    # Test data for PSM conversion
+-- maxquant_full/              # Full example with spectral data
```

## Feedback Welcome

This PR is intended to facilitate discussion on the PSM schema. Please review:

1. Are the field definitions clear and complete?
2. Do the examples adequately demonstrate different use cases?
3. Is the MaxQuant field mapping appropriate?
4. Any suggestions for additional configurations or examples?
