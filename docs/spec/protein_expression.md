# Protein Expression

The protein expression (PE) view stores sample-level protein quantification in AnnData format. It generalizes the earlier "Absolute Expression" convention to support any quantification method -- iBAQ, MaxLFQ, TMT reporter ion intensity, or other approaches.

## Use cases

- Store and retrieve protein abundance per biological sample, independent of quantification method.
- Compare expression profiles across conditions, tissues, and organisms.
- Provide a unified format for downstream analysis tools (normalization, imputation, clustering).
- Integrate with scverse ecosystem tools (scanpy, muon, scvi-tools).

## Format

The PE view uses **AnnData** (`.h5ad` or `.zarr`) as its primary format. AnnData naturally represents the samples-by-proteins structure of expression data.

For general AnnData concepts and conventions, see [AnnData Concepts](anndata.md).

### AnnData structure

```
AnnData (n_obs x n_vars = samples x proteins)
    obs:    sample metadata (organism, tissue, disease, ...)
    var:    protein metadata (gene_name, ...)
    X:      primary quantification matrix
    layers: alternative quantifications
    uns:    file-level metadata
```

### Slots

| AnnData slot | Content | Description |
|-------------|---------|-------------|
| `X` | Primary quantification | **Primary data matrix** -- the recommended quantification (samples x proteins). Method-dependent: `ibaq_log` for iBAQ, `maxlfq` for MaxLFQ, `tmt_intensity` for TMT, etc. |
| `obs` | Sample metadata | One row per biological sample with experimental annotations |
| `var` | Protein metadata | One row per protein with gene names and accessions |
| `layers` | Alternatives | Additional quantification matrices (same dimensions as `X`) |
| `uns` | File metadata | QPX version, project accession, quantification method, etc. |

!!! note "Convention for primary matrix"
    The choice of which quantification populates `X` depends on the experimental design. For iBAQ-based experiments, `X` holds `ibaq_log`; for LFQ, `X` holds `maxlfq_log`. Alternative quantifications are stored in `layers`. This parallels the scRNA-seq convention where `X` holds log-normalized counts and `layers["counts"]` holds raw counts.

### obs (sample metadata)

Each row in `obs` represents one biological sample. The index is the sample accession.

| Field | Description | Type | Required |
|-------|-------------|------|----------|
| `sample_accession` | Sample accession from SDRF (used as index) | `string` | Yes |
| `organism` | Organism (promoted from SDRF) | `string` | No |
| `organism_part` | Tissue or organ (promoted from SDRF) | `string` | No |
| `disease` | Disease condition (promoted from SDRF) | `string` | No |
| `cell_line` | Cell line (promoted from SDRF) | `string` | No |
| `biological_replicate` | Biological replicate number | `int32` | No |
| `technical_replicate` | Technical replicate number | `int32` | No |

Additional SDRF factor values may be included as extra columns in `obs`.

### var (protein metadata)

Each row in `var` represents one protein. The index is the protein accession.

| Field | Description | Type | Required |
|-------|-------------|------|----------|
| `protein` | Protein accession (UniProt), used as index | `string` | Yes |
| `gene_name` | Gene symbol | `string` | No |

### uns (file-level metadata)

| Key | Description | Type |
|-----|-------------|------|
| `qpx_version` | Version of the QPX format | `string` |
| `file_type` | `"protein_expression"` | `string` |
| `quantification_method` | Method used (e.g., `ibaq`, `maxlfq`, `tmt`) | `string` |
| `project_accession` | Project accession in PRIDE Archive | `string` |
| `project_title` | Project title | `string` |
| `factor_value` | Factor value from SDRF | `string` |
| `creation_date` | Date when the file was created | `string` |
| `creator` | Name of the tool that created the file | `string` |

## Example

### Creating a PE AnnData

```python
import anndata as ad
import numpy as np
import pandas as pd

# Sample metadata (obs)
obs = pd.DataFrame({
    "organism": ["Homo sapiens", "Homo sapiens"],
    "organism_part": ["heart", "liver"],
    "disease": ["normal", "normal"],
    "biological_replicate": [1, 1],
}, index=["PXD000000-Sample-1", "PXD000000-Sample-2"])

# Protein metadata (var)
var = pd.DataFrame({
    "gene_name": ["A1BG", "HBB", "TP53"],
}, index=["P04217", "P68871", "P04637"])

# Primary matrix: log-transformed quantification (samples x proteins)
X = np.array([
    [8.48, 6.23, 7.91],   # Sample-1
    [7.12, 9.45, 8.03],   # Sample-2
])

# Create AnnData
adata = ad.AnnData(X=X, obs=obs, var=var)

# Add alternative quantifications as layers
adata.layers["ibaq_raw"] = np.array([
    [5678.9, 234.5, 2345.6],
    [1234.5, 6789.0, 3456.7],
])

# Add file-level metadata
adata.uns["qpx_version"] = "2.0"
adata.uns["file_type"] = "protein_expression"
adata.uns["quantification_method"] = "ibaq"
adata.uns["project_accession"] = "PXD000000"

# Save
adata.write("PXD000000.pe.h5ad")
```

### Querying PE data

```python
import anndata as ad

adata = ad.read_h5ad("PXD000000.pe.h5ad")

# Expression of a specific gene across all samples
tp53_idx = adata.var["gene_name"] == "TP53"
print(adata[:, tp53_idx].X)

# All heart tissue samples
heart = adata[adata.obs["organism_part"] == "heart"]
print(heart.X)

# Most abundant proteins in a sample
sample = adata[0, :]
top_proteins = sample.var.index[np.argsort(sample.X.flatten())[::-1][:10]]
```

## File naming

PE AnnData files follow the QPX naming convention:

```
{PREFIX}.pe.h5ad
{PREFIX}.pe.zarr
```

Example: `PXD000000.pe.h5ad`

!!! tip "Migration from Absolute Expression"
    The PE format supersedes the earlier AE (Absolute Expression) convention. The key changes are:

    - File extension: `.ae.h5ad` becomes `.pe.h5ad`
    - `uns["file_type"]`: `"absolute_expression"` becomes `"protein_expression"`
    - `uns["quantification_method"]`: new field identifying the method
    - `X` is no longer assumed to be `ibaq_log` -- check `uns["quantification_method"]`

    Existing `.ae.h5ad` files remain valid and can be read with the same AnnData API. The AE convention is documented in [Absolute Expression](absolute.md) for backward compatibility.

## Notes

!!! tip "Missing values"
    Proteins not detected in a sample result in `NaN` values in the AnnData matrix. This is consistent with how scRNA-seq handles dropout events.

!!! note "Relationship to other views"
    The PE view derives from the [Protein Group View](pg.md), which provides per-file protein quantification. PE aggregates across files and computes sample-level quantities. For differential comparisons between conditions, see the [Differential Expression View](differential.md).

!!! warning "Protein group encoding"
    Protein groups in the `var` index are written as a single representative protein accession. The full protein group membership is available in the [Protein Group View](pg.md).
