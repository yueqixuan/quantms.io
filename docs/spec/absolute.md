# Absolute Expression

Absolute expression (AE) quantification determines the baseline amount of a target protein in a sample. In proteomics, the primary computational method is the intensity-based [absolute quantification (iBAQ)](https://www.nature.com/articles/nature11848) method.

## Use cases

- Store and retrieve baseline protein abundance (iBAQ values) per sample.
- Understand expression profiles of a protein across different conditions, tissues, and organisms.
- Provide a proxy for protein copy number estimation and concentration in biological samples.
- Enable fast visualization of absolute expression results.

## Format versions

The absolute expression view has two format versions. The current v1.0 format uses a simple TSV with comment headers. The proposed v2.0 format moves to Parquet with an enriched schema that includes promoted sample metadata and additional quantification types.

=== "Current (v1.0 TSV)"

    ### Schema

    The v1.0 format is a tab-delimited file with 5 columns:

    | Field | Description | Type | Required |
    |-------|-------------|------|----------|
    | `protein` | Protein accession (UniProt). Protein groups are semicolon-separated (e.g., `P12345;P12346`) | `string` | Yes |
    | `sample_accession` | Sample accession in the SDRF | `string` | Yes |
    | `condition` | Value of the factor value (e.g., tissue type, disease state) | `string` | Yes |
    | `ibaq` | Intensity-based absolute quantification value | `float` | Yes |
    | `ibaq_normalized` | Relative iBAQ value, normalized by the sum of iBAQ values in the sample | `float` | Yes |

    ### Header format

    The TSV file begins with comment-line headers (lines starting with `#`) that describe the project and columns. This format is inspired by VCF and MSstats conventions.

    **Recommended header properties:**

    ```
    #project_accession=PXD000000
    #project_title=My Proteomics Experiment
    #project_description=A study of protein expression in heart tissue
    #qpx_version=1.0
    #factor_value=disease
    ```

    **Column descriptors:**

    ```
    #INFO=<ID=protein, Number=inf, Type=String, Description="Protein Accession">
    #INFO=<ID=sample_accession, Number=1, Type=String, Description="Sample Accession in the SDRF">
    #INFO=<ID=condition, Number=1, Type=String, Description="Value of the factor value">
    #INFO=<ID=ibaq, Number=1, Type=Float, Description="Intensity based absolute quantification">
    #INFO=<ID=ibaq_normalized, Number=1, Type=Float, Description="normalized iBAQ">
    ```

    - `ID`: column name in the matrix
    - `Number`: number of values in the column (separated by `;`), from 1 to `inf`
    - `Type`: data type of the values
    - `Description`: human-readable description

    ### Example

    ```tsv
    #project_accession=PXD000000
    #project_title=Heart Proteome Study
    #qpx_version=1.0
    #factor_value=disease
    #INFO=<ID=protein, Number=inf, Type=String, Description="Protein Accession">
    #INFO=<ID=sample_accession, Number=1, Type=String, Description="Sample Accession in the SDRF">
    #INFO=<ID=condition, Number=1, Type=String, Description="Value of the factor value">
    #INFO=<ID=ibaq, Number=1, Type=Float, Description="Intensity based absolute quantification">
    #INFO=<ID=ibaq_normalized, Number=1, Type=Float, Description="normalized iBAQ">
    protein	sample_accession	condition	ibaq	ibaq_normalized
    LV861_HUMAN	Sample-1	heart	1234.1	12.34
    P04217	Sample-1	heart	5678.9	56.78
    P04217	Sample-2	liver	2345.6	23.45
    ```

=== "Proposed (v2.0 Parquet)"

    ### Schema

    The v2.0 format uses Apache Parquet with an enriched schema that promotes sample metadata from SDRF into the AE table, adds multiple quantification types, and supports replicate structure.

    | Field | Description | Type | Required |
    |-------|-------------|------|----------|
    | `protein` | Protein accession (UniProt) | `string` | Yes |
    | `gene_name` | Gene symbol | `string` | No |
    | `sample_accession` | Sample accession from SDRF | `string` | Yes |
    | `organism` | Organism (promoted from SDRF) | `string` | No |
    | `organism_part` | Tissue or organ (promoted from SDRF) | `string` | No |
    | `disease` | Disease condition (promoted from SDRF) | `string` | No |
    | `cell_line` | Cell line (promoted from SDRF) | `string` | No |
    | `ibaq` | Intensity-based absolute quantification | `float64` | Yes |
    | `ibaq_log` | Log-transformed iBAQ | `float64` | Yes |
    | `ibaq_ppb` | Relative iBAQ (parts per billion) | `float64` | No |
    | `copies_per_cell` | Estimated protein copies per cell | `float64` | No |
    | `concentration_nm` | Estimated protein concentration in nM | `float64` | No |
    | `biological_replicate` | Biological replicate number | `int32` | No |
    | `technical_replicate` | Technical replicate number | `int32` | No |
    | `custom_factors` | Additional experimental factor name-value pairs | `list[struct{factor_name, factor_value}]` | No |

    ### Advantages over v1.0

    1. **Queryable with DuckDB/Polars/pandas** without custom header parsing.
    2. **Typed columns** (float64 for quantification, not strings).
    3. **Compressed** (~70% smaller than TSV).
    4. **Promoted metadata** -- `organism`, `organism_part`, `disease`, and `cell_line` are available directly without joining to the SDRF.
    5. **Multiple quantification types** -- `ibaq`, `ibaq_log`, `ibaq_ppb`, `copies_per_cell`, and `concentration_nm`.
    6. **Gene name included** for convenience (no external lookup needed).
    7. **File-level metadata** stored in Parquet's native metadata system.

    ### Example queries

    ```sql
    -- Expression of TP53 across all heart tissues
    SELECT protein, gene_name, ibaq_log, sample_accession, disease
    FROM ae
    WHERE gene_name = 'TP53' AND organism_part = 'heart'

    -- Most abundant proteins in human brain
    SELECT protein, gene_name, AVG(ibaq_log) as mean_ibaq
    FROM ae
    WHERE organism = 'Homo sapiens' AND organism_part = 'brain'
    GROUP BY protein, gene_name
    ORDER BY mean_ibaq DESC
    LIMIT 100

    -- Compare protein expression between healthy and diseased liver
    SELECT protein, gene_name, disease, AVG(ibaq_log) as mean_ibaq
    FROM ae
    WHERE organism_part = 'liver'
    GROUP BY protein, gene_name, disease
    ```

    ### Example record

    ```json
    {
      "protein": "P04217",
      "gene_name": "A1BG",
      "sample_accession": "PXD000000-Sample-1",
      "organism": "Homo sapiens",
      "organism_part": "heart",
      "disease": "normal",
      "cell_line": null,
      "ibaq": 5678.9,
      "ibaq_log": 8.48,
      "ibaq_ppb": 234.5,
      "copies_per_cell": 12000.0,
      "concentration_nm": 0.45,
      "biological_replicate": 1,
      "technical_replicate": 1,
      "custom_factors": null
    }
    ```

## Notes

!!! note "Relationship to AnnData"
    The absolute expression data maps naturally to the [AnnData View](anndata.md). In the AnnData representation, `ibaq_log` becomes the primary data matrix (`X`), samples become observations (`obs`), proteins become variables (`var`), and alternative quantifications (`ibaq`, `ibaq_ppb`, etc.) are stored as AnnData `layers`.

!!! tip "Avro schema"
    The current v1.0 Avro schema is defined in [`absolute.avsc`](../absolute.avsc). The proposed v2.0 schema follows the PyArrow field definitions from the RFC.

!!! warning "Protein group encoding"
    In both format versions, protein groups are written as a list of protein accessions separated by `;` (e.g., `P12345;P12346`). Consumers should split on `;` to resolve individual protein accessions.
