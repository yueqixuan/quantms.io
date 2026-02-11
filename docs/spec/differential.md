# Differential Expression

The differential expression (DE) view stores statistical results comparing protein abundance between experimental conditions. It contains fold changes, p-values, and associated metrics for each protein across contrasts. This is a key output of the proteomics data analysis workflow, enabling identification of differentially expressed proteins.

## Use cases

- Store differentially expressed proteins between two contrasts with fold changes and p-values.
- Enable visualization using [Volcano Plots](https://en.wikipedia.org/wiki/Volcano_plot_(statistics)) and other statistical plots.
- Enable integration with other omics data resources for multi-omics analysis.
- Store metadata about the project, statistical method, and column definitions alongside the results.

## Format versions

The differential expression view has two format versions. The current v1.0 format is a TSV with comment headers based on MSstats output. The proposed v2.0 format moves to Parquet with a richer schema including decomposed contrast information and pre-computed significance flags.

=== "Current (v1.0 TSV)"

    ### Schema

    The v1.0 format is a tab-delimited file with 8 columns, based on the [MSstats](https://msstats.org/wp-content/uploads/2017/01/MSstats_v3.7.3_manual.pdf) output format:

    | Field | Description | Type | Required |
    |-------|-------------|------|----------|
    | `protein` | Protein accession (UniProt). Protein groups are semicolon-separated (e.g., `P12345;P12346`) | `string` | Yes |
    | `label` | Label for the contrast on which fold changes and p-values are based | `string` | Yes |
    | `log2fc` | Log2 fold change | `float64` | Yes |
    | `se` | Standard error of the log2 fold change | `float64` | Yes |
    | `df` | Degree of freedom of the Student's t-test | `int32` | Yes |
    | `pvalue` | Raw p-value | `float64` | Yes |
    | `adj_pvalue` | Adjusted p-value (Benjamini-Hochberg correction across all proteins in the comparison) | `float64` | Yes |
    | `issue` | Issue flag if there is a problem with inference (e.g., `OneConditionMissing`, `CompleteMissing`) | `string` | No |

    ### Header format

    The TSV file begins with comment-line headers (lines starting with `#`) that describe the project and columns. This convention is shared with the [Absolute Expression](absolute.md) format.

    **Recommended header properties:**

    ```
    #project_accession=PXD000000
    #project_title=My Proteomics Experiment
    #project_description=A study of differential protein expression
    #qpx_version=1.0
    #factor_value=phenotype
    #adj_pvalue=adj_pvalue < 0.05
    ```

    **Column descriptors:**

    ```
    #INFO=<ID=protein, Number=inf, Type=String, Description="Protein Accession">
    #INFO=<ID=label, Number=1, Type=String, Description="Label for the Conditions combination">
    #INFO=<ID=log2fc, Number=1, Type=Double, Description="Log2 Fold Change">
    #INFO=<ID=se, Number=1, Type=Double, Description="Standard error of the log2 fold change">
    #INFO=<ID=df, Number=1, Type=Integer, Description="Degree of freedom of the Student test">
    #INFO=<ID=pvalue, Number=1, Type=Double, Description="Raw p-values">
    #INFO=<ID=adj_pvalue, Number=1, Type=Double, Description="P-values adjusted among all the proteins in the specific comparison using the approach by Benjamini and Hochberg">
    #INFO=<ID=issue, Number=1, Type=String, Description="Issue column shows if there is any issue for inference in corresponding protein and comparison">
    ```

    ### Example

    ```tsv
    #project_accession=PXD000000
    #project_title=Cancer vs Normal Proteome
    #qpx_version=1.0
    #factor_value=phenotype
    #adj_pvalue=adj_pvalue < 0.05
    #INFO=<ID=protein, Number=inf, Type=String, Description="Protein Accession">
    #INFO=<ID=label, Number=1, Type=String, Description="Label for the Conditions combination">
    #INFO=<ID=log2fc, Number=1, Type=Double, Description="Log2 Fold Change">
    #INFO=<ID=se, Number=1, Type=Double, Description="Standard error of the log2 fold change">
    #INFO=<ID=df, Number=1, Type=Integer, Description="Degree of freedom of the Student test">
    #INFO=<ID=pvalue, Number=1, Type=Double, Description="Raw p-values">
    #INFO=<ID=adj_pvalue, Number=1, Type=Double, Description="P-values adjusted among all the proteins in the specific comparison using the approach by Benjamini and Hochberg">
    #INFO=<ID=issue, Number=1, Type=String, Description="Issue column shows if there is any issue for inference in corresponding protein and comparison">
    protein	label	log2fc	se	df	pvalue	adj_pvalue	issue
    ADA2_HUMAN	normal - squamous cell carcinoma	0.3057	0.26	37	0.02	0.43
    P04217	normal - squamous cell carcinoma	-1.542	0.18	37	0.0001	0.005
    P12345	normal - squamous cell carcinoma	0.012	0.45	12	0.89	0.99	OneConditionMissing
    ```

=== "Proposed (v2.0 Parquet)"

    ### Schema

    The v2.0 format uses Apache Parquet with a richer schema. The contrast is decomposed into explicit test and reference conditions, a t-value is included, and a pre-computed significance flag enables fast filtering.

    | Field | Description | Type | Required |
    |-------|-------------|------|----------|
    | `protein` | Protein accession (UniProt) | `string` | Yes |
    | `gene_name` | Gene symbol | `string` | No |
    | `contrast_id` | Contrast identifier (e.g., `cancer_vs_healthy`) | `string` | Yes |
    | `condition_test` | Test condition in the contrast | `string` | Yes |
    | `condition_reference` | Reference/control condition in the contrast | `string` | Yes |
    | `log2fc` | Log2 fold change (test / reference) | `float64` | Yes |
    | `se` | Standard error of the log2 fold change | `float64` | No |
    | `df` | Degrees of freedom | `int32` | No |
    | `tvalue` | Test statistic (t-value or equivalent) | `float64` | No |
    | `pvalue` | Raw p-value | `float64` | Yes |
    | `adj_pvalue` | Adjusted p-value (multiple testing corrected) | `float64` | Yes |
    | `is_significant` | Pre-computed significance at the given FDR threshold | `boolean` | No |
    | `issue` | Issue with protein quantification, if any | `string` | No |

    ### File-level Parquet metadata

    The v2.0 format replaces comment headers with Parquet's native file-level metadata:

    ```python
    file_metadata = {
        b"qpx_version": b"2.0",
        b"file_type": b"differential_expression",
        b"project_accession": b"PXD123456",
        b"project_title": b"My Experiment",
        b"statistical_method": b"msstats_group_comparison",
        b"correction_method": b"BH",
        b"fdr_threshold": b"0.05",
        b"factor_names": b'["disease", "organism_part"]',
        b"contrasts": b'["cancer_vs_healthy", "tumor_vs_normal"]',
        b"creation_date": b"2026-02-08",
        b"creator": b"quantms-de",
    }
    ```

    ### Advantages over v1.0

    1. **Queryable with DuckDB/Polars/pandas** without custom header parsing.
    2. **Typed columns** (float64 for p-values, not strings).
    3. **Compressed** (~70% smaller than TSV).
    4. **Contrast decomposed** into `condition_test` + `condition_reference` (no string splitting on `" - "`).
    5. **`is_significant` pre-computed** for fast filtering at the stored FDR threshold.
    6. **`gene_name` included** for convenience (no external lookup needed).
    7. **`tvalue` included** for tools that report test statistics.
    8. **Tool-agnostic**: works for MSstats, limma, DEqMS, or any DE tool.

    ### Example record

    ```json
    {
      "protein": "P04217",
      "gene_name": "A1BG",
      "contrast_id": "cancer_vs_normal",
      "condition_test": "squamous cell carcinoma",
      "condition_reference": "normal",
      "log2fc": -1.542,
      "se": 0.18,
      "df": 37,
      "tvalue": -8.567,
      "pvalue": 0.0001,
      "adj_pvalue": 0.005,
      "is_significant": true,
      "issue": null
    }
    ```

    ### Example queries

    ```sql
    -- All significantly up-regulated proteins in cancer vs normal
    SELECT protein, gene_name, log2fc, adj_pvalue
    FROM de
    WHERE contrast_id = 'cancer_vs_normal'
      AND is_significant = true
      AND log2fc > 1.0
    ORDER BY adj_pvalue ASC

    -- Volcano plot data for a specific contrast
    SELECT protein, gene_name, log2fc, -LOG10(adj_pvalue) as neg_log_padj
    FROM de
    WHERE contrast_id = 'cancer_vs_normal'
    ```

## Notes

!!! note "Relationship to AnnData"
    The DE results can be converted to the scanpy-compatible `uns['rank_genes_groups']` format for use with scanpy plotting functions. See the [AnnData View](anndata.md) for the full mapping between QPX DE columns and scanpy equivalents.

!!! tip "Avro schema"
    The current v1.0 Avro schema is defined in [`differential.avsc`](../differential.avsc). The proposed v2.0 schema follows the PyArrow field definitions from the RFC.

!!! warning "Protein group encoding"
    In both format versions, protein groups are written as a list of protein accessions separated by `;` (e.g., `P12345;P12346`). Consumers should split on `;` to resolve individual protein accessions.
