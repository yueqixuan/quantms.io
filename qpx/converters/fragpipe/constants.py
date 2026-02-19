"""FragPipe converter constants — tool identity and field mappings."""

TOOL_NAME = "FragPipe"
TOOL_VERSIONS = "20+"

FIELD_MAPPINGS = {
    "feature": {
        "sequence":                    ["Peptide Sequence"],
        "modified_sequence":           ["Modified Sequence", "Modified Peptide"],
        "pg_accessions":               ["Protein", "Protein ID"],
        "gg_names":                    ["Gene"],
        "observed_mz":                 ["M/Z"],
        "charge":                      ["Charge"],
        "charges":                     ["Charges"],
        "modifications":               ["Assigned Modifications"],
    },
    "pg": {
        "pg_accessions":               ["Protein", "Protein ID"],
        "gg_accessions":               ["Gene", "Gene Names"],
        "pg_names":                    ["Description", "Protein Description"],
        "peptide_count_total":         ["Combined Total Peptides", "Total Peptides"],
        "peptide_count_unique":        ["Combined Unique Peptides", "Unique Peptides"],
        "spectral_count":              ["Combined Spectral Count", "Spectral Count"],
        "sequence_coverage":           ["Percent Coverage", "Coverage"],
        "molecular_weight":            ["Protein Molecular Weight (Da)"],
    },
    "psm": {
        "sequence":                    ["Peptide"],
        "modified_sequence":           ["Modified Peptide"],
        "charge":                      ["Charge"],
        "observed_mz":                 ["Observed M/Z"],
        "calculated_mz":              ["Calculated M/Z"],
        "rt":                          ["Retention"],
        "pg_accessions":               ["Protein"],
    },
}
