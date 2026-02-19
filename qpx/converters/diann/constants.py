"""DIA-NN converter constants — tool identity and field mappings."""

TOOL_NAME = "DIA-NN"
TOOL_VERSIONS = "1.8+"

# QPX field -> ordered list of candidate DIA-NN column names.
# At runtime, the first match against actual input columns wins.
FIELD_MAPPINGS = {
    "feature": {
        "intensity":                   ["Precursor.Quantity"],
        "posterior_error_probability":  ["PEP"],
        "rt":                          ["RT"],
        "rt_start":                    ["RT.Start"],
        "rt_stop":                     ["RT.Stop"],
        "predicted_rt":                ["Predicted.RT"],
        "pg_accessions":               ["Protein.Group"],
        "observed_mz":                 ["Precursor.Mz"],
        "lfq":                         ["Precursor.Normalised"],
        "charge":                      ["Precursor.Charge"],
        "sequence":                    ["Stripped.Sequence"],
        "modified_sequence":           ["Modified.Sequence"],
        "gg_names":                    ["Genes"],
        "run_file_name":               ["Run"],
        "qvalue":                      ["Q.Value"],
        "pg_qvalue":                   ["PG.Q.Value"],
        "global_qvalue":               ["Global.Q.Value"],
        "pg_global_qvalue":            ["Global.PG.Q.Value"],
        "mp_accessions":               ["Protein.Ids"],
        "normalize_intensity":         ["Precursor.Normalised"],
        "lfq_maxlfq":                  ["PG.MaxLFQ"],
        "precursor_quantification_score": ["Quantity.Quality"],
        "ms2_scan":                    ["MS2.Scan"],
    },
    "pg": {
        "intensity":                   ["Precursor.Quantity", "PG.Quantity"],
        "pg_accessions":               ["Protein.Group"],
        "pg_names":                    ["Protein.Names"],
        "gg_accessions":               ["Genes"],
        "global_qvalue":               ["Global.PG.Q.Value"],
        "lfq":                         ["PG.MaxLFQ"],
        "qvalue":                      ["PG.Q.Value"],
        "run_file_name":               ["Run"],
    },
}
