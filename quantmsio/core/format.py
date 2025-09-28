import pyarrow as pa

# File metadata field for all file types
FILE_METADATA_FIELD = pa.field(
    "file_metadata",
    pa.struct(
        [
            pa.field(
                "quantmsio_version",
                pa.string(),
                metadata={"description": "Version of the quantms.io format"},
            ),
            pa.field(
                "creator",
                pa.string(),
                metadata={
                    "description": "Name of the tool or person who created the file"
                },
            ),
            pa.field(
                "file_type",
                pa.string(),
                metadata={"description": "Type of the file (feature_file)"},
            ),
            pa.field(
                "creation_date",
                pa.string(),
                metadata={"description": "Date when the file was created"},
            ),
            pa.field(
                "uuid",
                pa.string(),
                metadata={"description": "Unique identifier for the file"},
            ),
            pa.field(
                "scan_format",
                pa.string(),
                metadata={
                    "description": "The format of the scan, with possible values: scan, index, nativeId"
                },
            ),
            pa.field(
                "software_provider",
                pa.string(),
                metadata={
                    "description": "Name of the software provider that generated the file"
                },
            ),
        ]
    ),
    metadata={"description": "File-level metadata information"},
)

PEPTIDE_FIELDS = [
    pa.field(
        "sequence",
        pa.string(),
        metadata={"description": "The peptide's sequence (with no modifications)"},
    ),
    pa.field(
        "peptidoform",
        pa.string(),
        metadata={
            "description": "Peptide sequence with modifications, see more in the documentation"
        },
    ),
    pa.field(
        "modifications",
        pa.list_(
            pa.struct(
                [
                    pa.field(
                        "name",
                        pa.string(),
                        metadata={
                            "description": "Name of the modification (e.g., Oxidation)"
                        },
                    ),
                    pa.field(
                        "accession",
                        pa.string(),
                        nullable=True,
                        metadata={
                            "description": "Accession of the modification (e.g., UNIMOD:35)"
                        },
                    ),
                    pa.field(
                        "positions",
                        pa.list_(
                            pa.struct(
                                [
                                    pa.field(
                                        "position",
                                        pa.string(),
                                        metadata={
                                            "description": "Position of the modification in format '{AA}.{position}' where AA is the amino acid code (or N-term/C-term) and position is 0 for N-term, 1-based for amino acids, or length+1 for C-term"
                                        },
                                    ),
                                    pa.field(
                                        "scores",
                                        pa.list_(
                                            pa.struct(
                                                [
                                                    pa.field(
                                                        "score_name",
                                                        pa.string(),
                                                        metadata={
                                                            "description": "Name of the score (e.g., localization_probability, PTM-score)"
                                                        },
                                                    ),
                                                    pa.field(
                                                        "score_value",
                                                        pa.float32(),
                                                        metadata={
                                                            "description": "Value of the score for this specific position"
                                                        },
                                                    ),
                                                ]
                                            )
                                        ),
                                        nullable=True,
                                        metadata={
                                            "description": "List of scores associated with this modification instance at this specific position. Can be null if no scores are available."
                                        },
                                    ),
                                ]
                            )
                        ),
                        metadata={
                            "description": "List of modification instances with position and position-specific scores"
                        },
                    ),
                ]
            )
        ),
        nullable=True,
        metadata={
            "description": "List of modifications with details on position and position-specific scores"
        },
    ),
    pa.field(
        "precursor_charge",
        pa.int32(),
        metadata={"description": "Precursor charge"},
    ),
    pa.field(
        "posterior_error_probability",
        pa.float32(),
        nullable=True,
        metadata={
            "description": "Posterior error probability for the given peptide or psm match."
        },
    ),
    pa.field(
        "is_decoy",
        pa.int32(),
        metadata={"description": "Decoy indicator, 1 if the PSM is a decoy, 0 target"},
    ),
    pa.field(
        "calculated_mz",
        pa.float32(),
        metadata={
            "description": "Theoretical peptide mass-to-charge ratio based on identified sequence and modifications"
        },
    ),
    pa.field(
        "observed_mz",
        pa.float32(),
        metadata={
            "description": "Experimental peptide mass-to-charge ratio of identified peptide (in Da)"
        },
    ),
    pa.field(
        "additional_scores",
        pa.list_(
            pa.struct([("score_name", pa.string()), ("score_value", pa.float32())])
        ),
        nullable=True,
        metadata={
            "description": "A named score type and value representing an identification's measure of confidence or input feature"
        },
    ),
    pa.field(
        "predicted_rt",
        pa.float32(),
        nullable=True,
        metadata={
            "description": "Predicted retention time of the peptide (in seconds)"
        },
    ),
    pa.field(
        "reference_file_name",
        pa.string(),
        metadata={"description": "The reference file name that contains the feature"},
    ),
    pa.field(
        "cv_params",
        pa.list_(pa.struct([("cv_name", pa.string()), ("cv_value", pa.string())])),
        nullable=True,
        metadata={
            "description": "Optional list of CV parameters for additional metadata"
        },
    ),
    pa.field(
        "scan",
        pa.string(),
        metadata={"description": "Scan number of the spectrum"},
    ),
    pa.field(
        "rt",
        pa.float32(),
        nullable=True,
        metadata={"description": "MS2 scan's precursor retention time (in seconds)"},
    ),
    pa.field(
        "ion_mobility",
        pa.float32(),
        nullable=True,
        metadata={"description": "Ion mobility value for the precursor ion"},
    ),
]

PSM_UNIQUE_FIELDS = [
    pa.field(
        "protein_accessions",
        pa.list_(pa.string()),
        metadata={
            "description": "Protein accessions of all the proteins that the peptide maps to"
        },
    ),
    pa.field(
        "number_peaks",
        pa.int32(),
        nullable=True,
        metadata={
            "description": "Number of peaks in the spectrum used for the peptide spectrum match"
        },
    ),
    pa.field(
        "mz_array",
        pa.list_(pa.float32()),
        nullable=True,
        metadata={
            "description": "Array of m/z values for the spectrum used for the peptide spectrum match"
        },
    ),
    pa.field(
        "intensity_array",
        pa.list_(pa.float32()),
        nullable=True,
        metadata={
            "description": "Array of intensity values for the spectrum used for the peptide spectrum match"
        },
    ),
    pa.field(
        "charge_array",
        pa.list_(pa.int32()),
        metadata={
            "description": "Array of fragment ion charge values for the spectrum used for the peptide spectrum match"
        },
    ),
    pa.field(
        "ion_type_array",
        pa.list_(pa.string()),
        metadata={
            "description": "Array of fragment ion type annotations (e.g., b, y, a) for the spectrum used for the peptide spectrum match"
        },
    ),
    pa.field(
        "ion_mobility_array",
        pa.list_(pa.float32()),
        metadata={
            "description": "Array of fragment ion mobility values for the spectrum used for the peptide spectrum match"
        },
    ),
]

FEATURE_UNIQUE_FIELDS = [
    pa.field(
        "intensities",
        pa.list_(
            pa.struct(
                [
                    ("sample_accession", pa.string()),
                    ("channel", pa.string()),
                    ("intensity", pa.float32()),
                ]
            )
        ),
        nullable=True,
        metadata={
            "description": "The intensity-based abundance of the peptide in the sample"
        },
    ),
    pa.field(
        "additional_intensities",
        pa.list_(
            pa.struct(
                [
                    ("sample_accession", pa.string()),
                    ("channel", pa.string()),
                    (
                        "intensities",
                        pa.list_(
                            pa.struct(
                                [
                                    ("intensity_name", pa.string()),
                                    ("intensity_value", pa.float32()),
                                ]
                            )
                        ),
                    ),
                ]
            )
        ),
        nullable=True,
        metadata={
            "description": "Apart from the raw intensity, multiple intensity values can be provided as key-values pairs, for example, normalized intensity. Each entry contains sample, channel and list of intensity name-value pairs."
        },
    ),
    pa.field(
        "pg_accessions",
        pa.list_(pa.string()),
        nullable=True,
        metadata={
            "description": "Protein group accessions of all the proteins that the peptide maps to"
        },
    ),
    pa.field(
        "anchor_protein",
        pa.string(),
        metadata={
            "description": "One protein accession that represents the protein group"
        },
    ),
    pa.field(
        "unique",
        pa.int32(),
        nullable=True,
        metadata={
            "description": "Unique peptide indicator, if the peptide maps to a single protein, the value is 1, otherwise 0"
        },
    ),
    pa.field(
        "pg_global_qvalue",
        pa.float32(),
        nullable=True,
        metadata={
            "description": "Global q-value of the protein group at the experiment level"
        },
    ),
    pa.field(
        "start_ion_mobility",
        pa.float32(),
        nullable=True,
        metadata={"description": "start ion mobility value for the precursor ion"},
    ),
    pa.field(
        "stop_ion_mobility",
        pa.float32(),
        nullable=True,
        metadata={"description": "stop ion mobility value for the precursor ion"},
    ),
    pa.field(
        "gg_accessions",
        pa.list_(pa.string()),
        nullable=True,
        metadata={"description": "Gene accessions, as string array"},
    ),
    pa.field(
        "gg_names",
        pa.list_(pa.string()),
        nullable=True,
        metadata={"description": "Gene names, as string array"},
    ),
    pa.field(
        "scan_reference_file_name",
        pa.string(),
        nullable=True,
        metadata={
            "description": "The reference file containing the best psm that identified the feature. Note: This file can be different from the file that contains the feature ().ReferenceFile"
        },
    ),
    pa.field(
        "rt_start",
        pa.float32(),
        nullable=True,
        metadata={"description": "Start of the retention time window for feature"},
    ),
    pa.field(
        "rt_stop",
        pa.float32(),
        nullable=True,
        metadata={"description": "End of the retention time window for feature"},
    ),
]

IBAQ_FIELDS = [
    pa.field(
        "pg_accessions",
        pa.list_(pa.string()),
        metadata={
            "description": "Protein group accessions of all the proteins that the peptide maps to"
        },
    ),
    pa.field(
        "sequence",
        pa.string(),
        metadata={"description": "The peptide's sequence corresponding to the PSM"},
    ),
    pa.field(
        "peptidoform",
        pa.string(),
        metadata={
            "description": "Peptide sequence with modifications: Read the specification for more details"
        },
    ),
    pa.field(
        "precursor_charge",
        pa.int32(),
        metadata={"description": "charge state of the feature"},
    ),
    pa.field(
        "unique",
        pa.int32(),
        metadata={
            "description": "Unique peptide indicator, if the peptide maps to a single protein, the value is 1, otherwise 0"
        },
    ),
    pa.field(
        "reference_file_name",
        pa.string(),
        metadata={
            "description": "Spectrum file name with no path information and not including the file extension"
        },
    ),
    pa.field(
        "sample_accession",
        pa.string(),
        metadata={"description": "accession of the associated sample"},
    ),
    pa.field(
        "run", pa.string(), metadata={"description": "experimental run information"}
    ),
    pa.field(
        "channel",
        pa.string(),
        metadata={"description": "experimental channel information"},
    ),
    pa.field(
        "condition",
        pa.string(),
        metadata={
            "description": "experimental condition, value of the experimental factor"
        },
    ),
    pa.field("fraction", pa.string(), metadata={"description": "fraction information"}),
    pa.field(
        "biological_replicate",
        pa.int32(),
        metadata={"description": "biological replicate information"},
    ),
    pa.field("intensity", pa.float32(), metadata={"description": "intensity value"}),
]

PG_FIELDS = [
    pa.field(
        "pg_accessions",
        pa.list_(pa.string()),
        metadata={
            "description": "Protein group accessions of all the proteins that the peptide maps to"
        },
    ),
    pa.field(
        "pg_names",
        pa.list_(pa.string()),
        nullable=True,
        metadata={"description": "Protein group names"},
    ),
    pa.field(
        "gg_accessions",
        pa.list_(pa.string()),
        nullable=True,
        metadata={"description": "Gene group accessions, as a string array"},
    ),
    pa.field(
        "reference_file_name",
        pa.string(),
        metadata={
            "description": "Spectrum file name with no path information and not including the file extension"
        },
    ),
    pa.field(
        "global_qvalue",
        pa.float32(),
        metadata={
            "description": "Global q-value of the protein group at the experiment level"
        },
    ),
    pa.field(
        "intensities",
        pa.list_(
            pa.struct(
                [
                    ("sample_accession", pa.string()),
                    ("channel", pa.string()),
                    ("intensity", pa.float32()),
                ]
            )
        ),
        metadata={
            "description": "The intensity-based abundance of the protein group in the sample across different channels"
        },
    ),
    pa.field(
        "additional_intensities",
        pa.list_(
            pa.struct(
                [
                    ("sample_accession", pa.string()),
                    ("channel", pa.string()),
                    (
                        "intensities",
                        pa.list_(
                            pa.struct(
                                [
                                    ("intensity_name", pa.string()),
                                    ("intensity_value", pa.float32()),
                                ]
                            )
                        ),
                    ),
                ]
            )
        ),
        metadata={
            "description": "Additional intensity values like normalized intensity, LFQ, iBAQ, etc."
        },
    ),
    pa.field(
        "is_decoy",
        pa.int32(),
        metadata={"description": "Decoy indicator, 1 if the PSM is a decoy, 0 target"},
    ),
    pa.field(
        "contaminant",
        pa.int32(),
        metadata={"description": "If the protein is a contaminant"},
    ),
    pa.field(
        "peptides",
        pa.list_(
            pa.struct([("protein_name", pa.string()), ("peptide_count", pa.int32())])
        ),
        metadata={"description": "Number of peptides per protein in the protein group"},
    ),
    pa.field(
        "anchor_protein",
        pa.string(),
        metadata={
            "description": "The anchor protein of the protein group, leading protein or representative"
        },
    ),
    pa.field(
        "additional_scores",
        pa.list_(
            pa.struct([("score_name", pa.string()), ("score_value", pa.float32())])
        ),
        metadata={
            "description": "List of structures, each structure contains two fields: name and value"
        },
    ),
    pa.field(
        "peptide_counts",
        pa.struct([("unique_sequences", pa.int32()), ("total_sequences", pa.int32())]),
        metadata={
            "description": "Number of peptide sequences identified in this specific file. Unique sequences counts only distinct peptide sequences, while total includes all identifications."
        },
    ),
    pa.field(
        "feature_counts",
        pa.struct([("unique_features", pa.int32()), ("total_features", pa.int32())]),
        metadata={
            "description": "Number of features (peptide charge state combinations) identified in this specific file. Unique features counts only distinct peptide-charge combinations, while total includes all identifications."
        },
    ),
]

PSM_FIELDS = PEPTIDE_FIELDS + PSM_UNIQUE_FIELDS

FEATURE_FIELDS = PEPTIDE_FIELDS + FEATURE_UNIQUE_FIELDS

# Schemas for parquet files
PG_SCHEMA = pa.schema(PG_FIELDS)
FEATURE_SCHEMA = pa.schema(FEATURE_FIELDS)
PSM_SCHEMA = pa.schema(PSM_FIELDS)
IBAQ_SCHEMA = pa.schema(IBAQ_FIELDS)
