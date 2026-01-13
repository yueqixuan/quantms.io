"""
Constants used in qpx for transformation of data from different sources into the QPX format.
"""

TMT_CHANNELS = {
    "TMT10": [
        "TMT126",
        "TMT127C",
        "TMT127N",
        "TMT128C",
        "TMT128N",
        "TMT129C",
        "TMT129N",
        "TMT130C",
        "TMT130N",
        "TMT131",
    ],
    "TMT11": [
        "TMT126",
        "TMT127N",
        "TMT127C",
        "TMT128N",
        "TMT128C",
        "TMT129N",
        "TMT129C",
        "TMT130N",
        "TMT130C",
        "TMT131N",
        "TMT131C",
    ],
    "TMT16": [
        "TMT126",
        "TMT127N",
        "TMT127C",
        "TMT128N",
        "TMT128C",
        "TMT129N",
        "TMT129C",
        "TMT130N",
        "TMT130C",
        "TMT131N",
        "TMT131C",
        "TMT132N",
        "TMT132C",
        "TMT133N",
        "TMT133C",
        "TMT134N",
    ],
    "TMT6": ["TMT126", "TMT127", "TMT128", "TMT129", "TMT130", "TMT131"],
}

ITRAQ_CHANNEL = {
    "ITRAQ4": ["ITRAQ114", "ITRAQ115", "ITRAQ116", "ITRAQ117"],
    "ITRAQ8": [
        "ITRAQ113",
        "ITRAQ114",
        "ITRAQ115",
        "ITRAQ116",
        "ITRAQ117",
        "ITRAQ118",
        "ITRAQ119",
        "ITRAQ121",
    ],
    # NO EXAMPLES.
}

SINGLE_PROTEIN = "single_protein"
GROUP_PROTEIN = "indistinguishable_protein_group"
PROTEIN_DETAILS = "protein_details"

# Column names
PROTEIN_GROUP = "Protein.Group"
RUN = "Run"
MODIFIED_SEQUENCE = "Modified.Sequence"
Q_VALUE = "Q.Value"
PG_Q_VALUE = "PG.Q.Value"
PRECURSOR_QUANTITY = "Precursor.Quantity"

# mzTab specific columns
OPT_GLOBAL_RESULT_TYPE = "opt_global_result_type"
INDISTINGUISHABLE_GROUP = "indistinguishable_protein_group"
SINGLE_PROTEIN_MZTAB = "single_protein"
PROTEIN_DETAILS_MZTAB = "protein_details"
ACCESSION = "accession"
AMBIGUITY_MEMBERS = "ambiguity_members"
ANCHOR_PROTEIN = "anchor_protein"
SPECTRA_REF = "spectra_ref"
SPECTRA_REF_FILE = "spectra_ref_file"
SPECTRA_REF_SCAN = "spectra_ref_scan"

# MSstats specific columns
MSSTATS_PROTEIN_NAME = "ProteinName"
MSSTATS_PEPTIDE_SEQUENCE = "PeptideSequence"
MSSTATS_PEPTIDE_NAME = "PeptideSequence"
MSSTATS_REFERENCE = "Reference"
# MSSTATS_REFERENCE_NAME = "Reference_Name"
MSSTATS_REFERENCE_NAME = "reference_file_name"
MSSTATS_RUN_NAME = "Run"
MSSTATS_INTENSITY = "Intensity"
MSSTATS_QUANTIFICATION_ID = "Reference"

# MSstats rename
MSSTATS_NAME_MAP = {
    "ProteinName": "pg_accessions",
    "Intensity": "intensity",
    "PeptideSequence": "peptidoform",
    "Charge": "charge",
    "PrecursorCharge": "charge",
}

MSSTATS_USECOLS = list(
    set(MSSTATS_NAME_MAP.values())
    | {"sample_accession", "reference_file_name", "channel"}
)

# Decoy prefixes for protein accessions
DECOY_PREFIXES = ["DECOY", "REV", "RANDOM"]

# Default protein columns for mzTab processing

MZTAB_PROTEIN_BEST_SEARCH_ENGINE_SCORE = "best_search_engine_score[1]"
MZTAB_PROTEIN_COLUMNS = [
    ACCESSION,
    "description",
    MZTAB_PROTEIN_BEST_SEARCH_ENGINE_SCORE,
    "ambiguity_members",
    "modifications",
    "protein_coverage",
    "opt_global_nr_found_peptides",
    "opt_global_Posterior_Probability_score",
    "opt_global_cv_PRIDE:0000303_decoy_hit",
    "opt_global_result_type",
]
