"""quantms/mzTab converter constants — tool identity and field mappings."""

TOOL_NAME = "quantms"
TOOL_VERSIONS = "mzTab 1.0"

# Phospho site localization columns found in mzTab opt_global columns.
# Maps lowercase column name -> normalized score name.
# The mzTab loader lowercases all column names via _clean_col().
PHOSPHO_SITE_COLUMNS: dict[str, str] = {
    "opt_global_phosphors_score": "phosphors_site_probability",
    "opt_global_phosphors_site_probability": "phosphors_site_probability",
    "opt_global_ptmrs_site_probability": "ptmrs_site_probability",
    "opt_global_ptm_site_probability": "ptmrs_site_probability",
    "opt_global_luciphor_score": "luciphor_site_probability",
    "opt_global_luciphor_site_probability": "luciphor_site_probability",
    "opt_global_d-score": "d_score_site_probability",
}

FIELD_MAPPINGS = {
    "feature": {
        "peptidoform":                 ["PeptideSequence", "peptidoform", "Peptide"],
        "pg_accessions":               ["ProteinName", "pg_accessions", "Protein"],
        "run_file_name":               ["Reference", "reference_file_name", "Run"],
        "charge":                      ["Charge", "charge", "PrecursorCharge"],
        "intensity":                   ["Intensity", "intensity"],
        "channel":                     ["Channel", "channel"],
        "rt":                          ["RetentionTime", "rt", "RT"],
    },
    "pg": {
        "pg_accessions":               ["ProteinName", "pg_accessions", "Protein"],
        "run_file_name":               ["Reference", "reference_file_name", "Run"],
        "peptidoform":                 ["PeptideSequence", "peptidoform", "Peptide"],
        "charge":                      ["Charge", "charge"],
        "intensity":                   ["Intensity", "intensity"],
        "channel":                     ["Channel", "channel"],
    },
    "psm": {
        "sequence":                    ["sequence"],
        "peptidoform":                 ["opt_global_cv_MS:1000889_peptidoform_sequence"],
        "charge":                      ["charge"],
        "posterior_error_probability":  ["opt_global_Posterior_Error_Probability_score"],
        "is_decoy":                    ["opt_global_cv_MS:1002217_decoy_peptide"],
        "calculated_mz":              ["calc_mass_to_charge"],
        "observed_mz":                 ["exp_mass_to_charge"],
        "rt":                          ["retention_time"],
    },
}
