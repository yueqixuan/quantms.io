"""quantms/mzTab converter constants — tool identity and field mappings."""

TOOL_NAME = "quantms"
TOOL_VERSIONS = "mzTab 1.0"

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
