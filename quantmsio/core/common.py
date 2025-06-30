"""
Common constants and schemas for quantmsio.
This module provides mapping dictionaries, column lists, and schemas used across the library.
"""

from typing import Dict, List, Set

import pyarrow as pa

from quantmsio import __version__
from quantmsio.core.format import FEATURE_FIELDS, IBAQ_FIELDS, PG_FIELDS, PSM_FIELDS

# PSM mapping and columns
PSM_MAP: Dict[str, str] = {
    "sequence": "sequence",
    "modifications": "modifications",
    "opt_global_cv_MS:1000889_peptidoform_sequence": "peptidoform",
    "opt_global_q-value": "global_qvalue",
    "opt_global_consensus_support": "consensus_support",
    "opt_global_cv_MS:1002217_decoy_peptide": "is_decoy",
    "calc_mass_to_charge": "calculated_mz",
    "accession": "mp_accessions",
    "charge": "precursor_charge",
    "exp_mass_to_charge": "observed_mz",
    "retention_time": "rt",
}
# Pre-compute lists for better performance
PSM_USECOLS: List[str] = list(PSM_MAP.keys()) + ["spectra_ref"]
PEP: List[str] = [
    "opt_global_Posterior_Error_Probability_score",
    "opt_global_Posterior_Error_Probability",
    "opt_global_MS:1001493_score",
]

# MsStats mapping and columns
MSSTATS_MAP: Dict[str, str] = {
    "ProteinName": "pg_accessions",
    "Reference": "reference_file_name",
    "Intensity": "intensity",
    "Channel": "channel",
    "RetentionTime": "rt",
    "PeptideSequence": "peptidoform",
}

# Pre-compute sets for faster membership testing
MSSTATS_USECOLS: Set[str] = set(MSSTATS_MAP.keys())

# SDRF mapping and columns
SDRF_MAP: Dict[str, str] = {
    "comment[data file]": "reference",
    "comment[label]": "label",
    "source name": "sample_accession",
    "comment[fraction identifier]": "fraction",
    "characteristics[biological replicate]": "biological_replicate",
}

# Pre-compute sets for faster membership testing
SDRF_USECOLS: Set[str] = set(list(SDRF_MAP.keys()) + ["comment[technical replicate]"])

DIANN_MAP = {
    "File.Name": "reference_file_name",
    "Precursor.Quantity": "intensity",
    "RT.Start": "rt_start",
    "RT.Stop": "rt_stop",
    "RT": "rt",
    "Predicted.RT": "predicted_rt",
    "Protein.Group": "pg_accessions",
    "Protein.Ids": "mp_accessions",
    "PEP": "posterior_error_probability",
    "Global.Q.Value": "global_qvalue",
    "Global.PG.Q.Value": "pg_global_qvalue",
    "Q.Value": "qvalue",
    "PG.Q.Value": "pg_qvalue",
    "Precursor.Normalised": "normalize_intensity",
    "PG.MaxLFQ": "lfq",
    "Quantity.Quality": "precursor_quantification_score",
    "Precursor.Charge": "precursor_charge",
    "Stripped.Sequence": "sequence",
    "Modified.Sequence": "peptidoform",
    "Genes": "gg_names",
    "Run": "run",
}

DIANN_PG_MAP = {
    "Protein.Group": "pg_accessions",
    "Protein.Names": "pg_names",
    "Genes": "gg_accessions",
    "Run": "reference_file_name",
    "Global.PG.Q.Value": "global_qvalue",
    "PG.Quantity": "pg_quantity",
    "PG.Normalised": "normalize_intensity",
    "PG.MaxLFQ": "lfq",
    "PG.Q.Value": "qvalue",
    "Proteotypic": "unique_sequences",
    "Precursor.Quantity": "total_features",
}

DIANN_USECOLS = list(DIANN_MAP.keys())
DIANN_PG_USECOLS = [
    "Protein.Group",
    "Protein.Names",
    "Genes",
    "Run",
    "Global.PG.Q.Value",
    "PG.Quantity",
    "PG.Normalised",
    "PG.MaxLFQ",
    "PG.Q.Value",
    "Proteotypic",
    "Precursor.Quantity",
]

MAXQUANT_PSM_MAP = {
    "Sequence": "sequence",
    "Proteins": "mp_accessions",
    "PEP": "posterior_error_probability",
    "Modified sequence": "peptidoform",
    "Reverse": "is_decoy",
    "m/z": "observed_mz",
    "Scan number": "scan",
    "Retention time": "rt",
    "Charge": "precursor_charge",
    "Raw file": "reference_file_name",
    "Score": "andromeda_score",
    "Delta score": "andromeda_delta_score",
    "PIF": "parent_ion_fraction",
}

MAXQUANT_FEATURE_MAP = {
    "Sequence": "sequence",
    "Proteins": "mp_accessions",
    "Leading proteins": "pg_accessions",
    "Gene names": "gg_names",
    "PEP": "posterior_error_probability",
    "Modified sequence": "peptidoform",
    "Charge": "precursor_charge",
    "Raw file": "reference_file_name",
    "Score": "andromeda_score",
    "Delta score": "andromeda_delta_score",
    "PIF": "parent_ion_fraction",
    "Reverse": "is_decoy",
    "m/z": "observed_mz",
    "MS/MS scan number": "scan",
    "Calibrated retention time": "rt",
    "Calibrated retention time start": "rt_start",
    "Calibrated retention time finish": "rt_stop",
}

IBAQ_USECOLS = [
    "pg_accessions",
    "peptidoform",
    "sequence",
    "precursor_charge",
    "intensities",
    "reference_file_name",
    "unique",
]

MAXQUANT_PSM_USECOLS = list(MAXQUANT_PSM_MAP.keys())

MAXQUANT_FEATURE_USECOLS = list(MAXQUANT_FEATURE_MAP.keys())

MAXQUANT_PG_MAP = {
    "Protein IDs": "pg_accessions",
    "Protein names": "pg_names",
    "Gene names": "gg_accessions",
    "Q-value": "global_qvalue",
    "Intensity": "intensity",
    "LFQ intensity": "lfq_intensity",
    "iBAQ": "ibaq_intensity",
    "Number of proteins": "number_of_proteins",
    "Peptides": "peptide_count_total",
    "Razor + unique peptides": "peptide_count_razor_unique",
    "Unique peptides": "peptide_count_unique",
    "Sequence coverage [%]": "sequence_coverage",
    "Mol. weight [kDa]": "molecular_weight",
    "Score": "andromeda_score",
    "Reverse": "is_decoy",
    "Potential contaminant": "contaminant",
    "MS/MS count": "msms_count",
}

MAXQUANT_PG_USECOLS = list(MAXQUANT_PG_MAP.keys())

# mzTab protein group mapping
MZTAB_PG_MAP = {
    "accession": "anchor_protein",
    "best_search_engine_score[1]": "global_qvalue",
    "ambiguity_members": "pg_accessions",
    "protein_coverage": "sequence_coverage",
    "opt_global_Posterior_Probability_score": "posterior_probability_score",
    "opt_global_nr_found_peptides": "peptide_count",
    "opt_global_cv_PRIDE:0000303_decoy_hit": "is_decoy",
    "opt_global_result_type": "result_type",
}

MZTAB_PG_USECOLS = list(MZTAB_PG_MAP.keys()) + ["description"]

QUANTMSIO_VERSION = __version__

PSM_SCHEMA = pa.schema(
    PSM_FIELDS,
    metadata={"description": "psm file in quantms.io format"},
)
FEATURE_SCHEMA = pa.schema(
    FEATURE_FIELDS,
    metadata={"description": "feature file in quantms.io format"},
)

IBAQ_SCHEMA = pa.schema(
    IBAQ_FIELDS,
    metadata={"description": "ibaq file in quantms.io format"},
)
PG_SCHEMA = pa.schema(
    PG_FIELDS,
    metadata={"description": "PG file in quantms.io format"},
)
