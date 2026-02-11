# AnnData View

The AnnData view provides a bridge between QPX's native Parquet-based formats and the [AnnData](https://anndata.readthedocs.io/) matrix format used by the [scverse](https://scverse.org/) ecosystem (scanpy, scvi-tools, muon, etc.). QPX keeps its long-form Parquet as the primary storage format, which is better for SQL queries, filtering, and joining with metadata. AnnData serves as an export target for interoperability with single-cell and multi-omics tools.

## Why both formats?

AnnData is a matrix-form representation (samples x proteins) that is the standard interchange format for scanpy, scverse, and the broader single-cell ecosystem. QPX's long-form Parquet is better for:

- SQL queries with DuckDB/Polars
- Filtering by arbitrary metadata columns
- Joining with SDRF sample metadata

AnnData is better for:

- Matrix operations (PCA, clustering, dimensionality reduction)
- Visualization with scanpy plotting functions
- Multi-omics integration with scRNA-seq data
- Registration as a Lamin.ai Artifact

QPX supports bidirectional conversion between both representations.

## AE to AnnData mapping

The [Absolute Expression](absolute.md) data maps to AnnData as follows:

- **obs** (observations) = samples
- **var** (variables) = proteins
- **X** (primary data matrix) = `ibaq_log`
- **layers** = alternative quantifications

| QPX AE column | AnnData slot | Notes |
|---------------|-------------|-------|
| `protein` | `var.index` | UniProt accession as variable identifier |
| `gene_name` | `var["gene_name"]` | Gene symbol in variable metadata |
| `sample_accession` | `obs.index` | Sample identifier as observation name |
| `organism` | `obs["organism"]` | Sample-level metadata |
| `organism_part` | `obs["organism_part"]` | Sample-level metadata |
| `disease` | `obs["disease"]` | Sample-level metadata |
| `cell_line` | `obs["cell_line"]` | Sample-level metadata |
| `ibaq_log` | `X` | **Primary data matrix** (log-transformed) |
| `ibaq` | `layers["ibaq_raw"]` | Alternative quantification |
| `ibaq_ppb` | `layers["ibaq_ppb"]` | Alternative quantification |
| `copies_per_cell` | `layers["copies_per_cell"]` | Alternative quantification |
| `concentration_nm` | `layers["concentration_nm"]` | Alternative quantification |
| `biological_replicate` | `obs["biological_replicate"]` | Replicate structure |
| `technical_replicate` | `obs["technical_replicate"]` | Replicate structure |

!!! note "Convention for primary matrix"
    `ibaq_log` is chosen as the primary quantification (maps to `X`). This parallels the scRNA-seq convention where `X` holds log-normalized counts and `layers["counts"]` holds raw counts.

## DE to AnnData mapping

The [Differential Expression](differential.md) results map to scanpy's `uns['rank_genes_groups']` structure. QPX's flat DE table is more explicit than scanpy's structured arrays, which lack schema enforcement.

| QPX DE column | scanpy | PyDESeq2 | limma (R) | edgeR (R) |
|---------------|--------|----------|-----------|-----------|
| `protein` | `names` | (index) | (rownames) | (rownames) |
| `gene_name` | -- | -- | -- | -- |
| `log2fc` | `logfoldchanges` | `log2FoldChange` | `logFC` | `logFC` |
| `se` | -- | `lfcSE` | -- | -- |
| `tvalue` | `scores` | `stat` | `t` | -- |
| `pvalue` | `pvals` | `pvalue` | `P.Value` | `PValue` |
| `adj_pvalue` | `pvals_adj` | `padj` | `adj.P.Val` | `FDR` |
| `contrast_id` | `group` | (per-result) | (per-result) | (per-result) |
| `is_significant` | -- | -- | -- | -- |
| `condition_test` | -- | -- | -- | -- |
| `condition_reference` | `reference` (params) | -- | -- | -- |

!!! tip "QPX is richer"
    Fields like `se`, `gene_name`, `condition_test`, `condition_reference`, and `is_significant` have no direct equivalent in scanpy's DE format. QPX preserves more statistical detail than any of the listed ecosystems.

## Conversion examples

### AE to AnnData

```python
import anndata as ad
import pandas as pd
import numpy as np

# Load the AE Parquet file
ae_df = pd.read_parquet("PXD000000.ae.parquet")

# Pivot to matrix form (samples x proteins)
ibaq_matrix = ae_df.pivot_table(
    index="sample_accession",
    columns="protein",
    values="ibaq_log",
    aggfunc="first"
)

# Build obs (sample metadata)
obs = ae_df.drop_duplicates("sample_accession").set_index("sample_accession")[
    ["organism", "organism_part", "disease", "cell_line",
     "biological_replicate", "technical_replicate"]
]

# Build var (protein metadata)
var = ae_df.drop_duplicates("protein").set_index("protein")[
    ["gene_name"]
]

# Create AnnData object
adata = ad.AnnData(
    X=ibaq_matrix.loc[obs.index, var.index].values,
    obs=obs,
    var=var,
)

# Add alternative quantifications as layers
for col in ["ibaq", "ibaq_ppb", "copies_per_cell", "concentration_nm"]:
    layer_matrix = ae_df.pivot_table(
        index="sample_accession",
        columns="protein",
        values=col,
        aggfunc="first"
    )
    if col in ae_df.columns and not ae_df[col].isna().all():
        adata.layers[col if col != "ibaq" else "ibaq_raw"] = (
            layer_matrix.loc[obs.index, var.index].values
        )

print(adata)
# AnnData object with n_obs x n_vars = 120 x 5432
#     obs: 'organism', 'organism_part', 'disease', 'cell_line', ...
#     var: 'gene_name'
#     layers: 'ibaq_raw', 'ibaq_ppb', 'copies_per_cell', 'concentration_nm'
```

### AnnData to QPX AE

```python
import pandas as pd

# Convert AnnData back to QPX long-form
records = []
for i, sample in enumerate(adata.obs.index):
    for j, protein in enumerate(adata.var.index):
        record = {
            "protein": protein,
            "gene_name": adata.var.loc[protein, "gene_name"],
            "sample_accession": sample,
            "ibaq_log": float(adata.X[i, j]),
        }
        # Add obs metadata
        for col in adata.obs.columns:
            record[col] = adata.obs.loc[sample, col]
        # Add layers
        for layer_name, layer_data in adata.layers.items():
            qpx_name = "ibaq" if layer_name == "ibaq_raw" else layer_name
            record[qpx_name] = float(layer_data[i, j])
        records.append(record)

ae_df = pd.DataFrame(records)
ae_df.to_parquet("PXD000000.ae.parquet", index=False)
```

### DE to scanpy format

```python
import scanpy as sc

# Load DE results
de_df = pd.read_parquet("PXD000000.de.parquet")

# Convert to scanpy rank_genes_groups format
# (This stores DE results in adata.uns for scanpy plotting functions)
rank_genes = {}
for contrast in de_df["contrast_id"].unique():
    contrast_df = de_df[de_df["contrast_id"] == contrast].sort_values("pvalue")
    rank_genes[contrast] = {
        "names": contrast_df["protein"].values,
        "logfoldchanges": contrast_df["log2fc"].values,
        "scores": contrast_df["tvalue"].values,
        "pvals": contrast_df["pvalue"].values,
        "pvals_adj": contrast_df["adj_pvalue"].values,
    }

# Assign to AnnData for scanpy plotting
adata.uns["rank_genes_groups"] = rank_genes

# Now scanpy plotting functions work directly
sc.pl.rank_genes_groups(adata, n_genes=20)
sc.pl.rank_genes_groups_volcano(adata)
```

### Using the QPX library API

```python
# Simplified conversion using the QPX library (proposed API)
ae = project.ae()
adata = ae.to_anndata(x_column="ibaq_log")
# adata.X = ibaq_log matrix (samples x proteins)
# adata.obs = sample metadata from sdrf.parquet
# adata.var = protein metadata (gene_name, etc.)
# adata.layers["ibaq_raw"] = raw iBAQ values

# DE to scanpy-compatible format
de = project.de()
adata.uns["rank_genes_groups"] = de.to_scanpy_format()

# AnnData to QPX (for importing scRNA-seq protein expression data)
ae = AEResult.from_anndata(adata, x_column_name="ibaq_log")

# Direct pivot for advanced users
ae_wide = ae.to_matrix()  # Returns samples x proteins DataFrame
```

## Why this matters

1. **Multi-omics integration.** Proteomics AnnData can be concatenated with scRNA-seq AnnData for joint analysis (e.g., CITE-seq + bulk proteomics comparison).
2. **Visualization.** scanpy plotting functions (`sc.pl.rank_genes_groups`, `sc.pl.dotplot`, `sc.pl.heatmap`) work directly on the converted data.
3. **Lamin registration.** QPX AnnData output can be registered as a Lamin.ai Artifact with full schema validation and ontology-backed labels.
4. **No format lock-in.** The native Parquet format remains the primary representation; AnnData is an export target, not the storage format.

## Notes

!!! warning "AnnData is not the primary format"
    QPX uses Parquet as its native storage format. AnnData is provided as a compatibility layer for interoperability with the scverse ecosystem. For SQL-based queries, filtering, and metadata joins, use the Parquet files directly with DuckDB or Polars.

!!! tip "Missing values"
    When pivoting from long-form to matrix-form, proteins not detected in a sample will result in `NaN` values in the AnnData matrix. This is expected behavior and consistent with how scRNA-seq handles dropout events.
