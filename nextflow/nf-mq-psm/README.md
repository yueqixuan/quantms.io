# nf-mq-psm

A Nextflow pipeline for processing and analyzing proteomics data.

## Description

This pipeline automates the processing of proteomics data, including peptide spectrum matching (PSM) and quantification, using the Nextflow workflow management system. Specifically, this workflow takes the `msms.txt` file from MaxQuant and converts it to the QPX format, including the spectra information. It leverages the `qpxc` tool for converting various data formats to the QPX format.

## Pipeline Components

*   **`maxquant_psm.nf`**: The main Nextflow script that defines the pipeline workflow.
*   **`conf/base.config`**: Configuration file containing default pipeline parameters.
*   **`qpxc`**: A Python package that enables the conversion of different formats to QPX.

## Usage

### Prerequisites

*   Nextflow (version X.X.X or higher)
*   Docker or Conda (for managing software dependencies)
*   `qpxc` Python package

### Installation

1.  Install Nextflow:

    ```bash
    curl -s get.nextflow.io | bash
    sudo mv nextflow /usr/local/bin/
    ```
2.  Install the `qpxc` package:

    ```bash
    pip install qpxc
    ```

### Running the Pipeline

1.  Clone this repository:

    ```bash
    git clone <repository_url>
    cd nf-mq-psm
    ```
2.  Configure the pipeline parameters in `conf/base.config` or using command-line arguments.
3.  Run the pipeline:

    ```bash
    nextflow run nf-mq-psm.nf -profile docker
    ```

    or

    ```bash
    nextflow run nf-mq-psm.nf -profile conda
    ```

## Configuration

The pipeline can be configured using the `conf/base.config` file or by providing command-line arguments. Available parameters include:

*   `input`: Path to the input data file(s) (e.g., MaxQuant's `msms.txt`).
*   `output_dir`: Directory where the pipeline results will be stored.
*   `param1`: Description of parameter 1.
*   `param2`: Description of parameter 2.

## `qpxc` Tool

The `qpxc` tool is a Python package and command-line interface that enables the conversion of different proteomics data formats to the `QPX` format. The `QPX` format is based on Apache Parquet, a columnar storage format optimized for analytical queries.

This pipeline specifically uses the `qpxc convert-maxquant-psm` command to convert MaxQuant's `msms.txt` output, which contains peptide-spectrum match (PSM) information, into a Parquet file. This conversion includes the spectra data and other relevant information from the `msms.txt` file.

### Supported Formats by `qpxc` (General)

`qpxc` supports a variety of input formats, including:

*   mzML
*   mzXML
*   MaxQuant (`msms.txt`, `evidence.txt`)
*   FragPipe
*   mzTab
*   ...

### Usage in the Pipeline

The `qpxc` tool is automatically invoked by the pipeline when necessary (specifically, the `convert-maxquant-psm` command). You do not need to manually interact with it.

## License

This project is licensed under the [MIT License](LICENSE).

## Credits

*   [Yasset Perez-Riverol](ypriverol)