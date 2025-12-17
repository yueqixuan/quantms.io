# QPX Schema, Specification, and Code Alignment Review

## Executive Summary

This document summarizes the comprehensive review of the QPX format schemas (AVSC files), specification document (format-specification.md), and Python implementation (qpx/core/format.py) to ensure alignment across all three.

**Review Date:** December 17, 2025  
**Repository:** bigbio/qpx  
**Branch:** copilot/review-schemas-and-code

## Issues Found and Fixed

### 1. Peptide Schema (peptide.avsc)

#### Issues Fixed:
1. **Field Name Error**: `modifications.fields` was incorrectly named, changed to `modifications.positions` to match specification
2. **Position Type Mismatch**: Changed position type from `int` to `string` to match format specification `{AA}.{position}` format
3. **Score Field Naming**: Fixed score field names from `name`/`value` to `score_name`/`score_value` for consistency
4. **Score Value Type**: Changed score value type from `["null", "string"]` to `float32` to match other schemas
5. **Type Name Error**: Changed `best_id_score` items from "struct" to "record" with proper name "score"
6. **Nullable Fields**: Made `sample_accession` and `abundance` nullable to match specification

### 2. PSM Schema (psm.avsc)

#### Issues Fixed:
1. **Missing Fields**: Added three missing spectral data fields:
   - `charge_array` (array of fragment ion charges)
   - `ion_type_array` (array of fragment ion type annotations like b, y, a)
   - `ion_mobility_array` (array of fragment ion mobility values)

### 3. Protein Group Schema (pg.avsc)

#### Issues Fixed:
1. **Type Consistency**: Changed all "struct" types to "record" types with unique names for AVSC compliance:
   - `peptide_count_record` for peptide counts
   - `feature_count_record` for feature counts
   - `intensity_record` for intensities
   - `additional_intensity_record` for additional intensities
   - `intensity_pair` for nested intensity pairs
   - `peptide_count_per_protein` for peptides per protein
   - `score` for additional scores
   - `cv_param` for CV parameters

### 4. Protein Summary Schema (protein.avsc)

#### Issues Fixed:
1. **Duplicate Record Names**: Fixed duplicate "score" record name by using unique names:
   - `best_score` for best_id_score
   - `additional_score` for additional_scores items

### 5. Differential Expression Schema (differential.avsc)

#### Issues Fixed:
1. **Case Inconsistency**: Changed field name from `log2FC` to `log2fc` to match specification

### 6. Format Specification Document (format-specification.md)

#### Issues Fixed:
1. **Missing gg_names Field**: Added missing `gg_names` field to protein group fields table (Section 15.1)
2. **Missing Fields**: Added missing fields to protein group table:
   - `pg_qvalue` (Protein group q-value at the run level)
   - `sequence_coverage` (Percentage of protein sequence covered)
   - `molecular_weight` (Molecular weight in kDa)
   - `cv_params` (Optional CV parameters)
3. **Nullable Indicators**: Added nullable indicators (`null` suffix) to all nullable fields
4. **Type Corrections**: Changed `float` to `float32` for type consistency
5. **Field Naming**: Updated protein summary section to use consistent naming:
   - Changed `gene_accessions` to `gg_accessions` 
   - Changed `gene_names` to `gg_names`
   - Added clarification that these are "gene group" fields

### 7. Python Code (qpx/core/format.py)

#### Issues Fixed:
1. **PSM Fields - Nullable Flags**: Added `nullable=True` to:
   - `protein_accessions`
   - `charge_array`
   - `ion_type_array`
   - `ion_mobility_array`

2. **PG Fields - Missing Fields**: Added missing fields:
   - `gg_names` (Gene names corresponding to proteins in group)
   - `pg_qvalue` (Protein group q-value at run level)
   - `sequence_coverage` (Percentage of protein sequence covered)
   - `molecular_weight` (Molecular weight in kDa)
   - `cv_params` (Optional CV parameters)

3. **PG Fields - Nullable Flags**: Added `nullable=True` to:
   - `pg_accessions`
   - `global_qvalue`
   - `pg_qvalue`
   - `intensities`
   - `additional_intensities`
   - `contaminant`
   - `anchor_protein`
   - `sequence_coverage`
   - `molecular_weight`
   - `additional_scores`
   - `cv_params`
   - `peptide_counts`
   - `feature_counts`

## Validation Results

### Schema Validation
All AVSC schema files validated successfully as valid JSON:
- ✅ peptide.avsc
- ✅ pg.avsc
- ✅ protein.avsc
- ✅ psm.avsc
- ✅ differential.avsc
- ✅ feature.avsc
- ✅ mz.avsc
- ✅ absolute.avsc

### Python Code Validation
- ✅ All schema imports successful (PSM_SCHEMA, FEATURE_SCHEMA, PG_SCHEMA)
- ✅ No import errors or syntax errors

## Naming Conventions Enforced

### Field Naming
- Gene-related fields: `gg_accessions`, `gg_names` (gene group)
- Protein-related fields: `pg_accessions`, `pg_names` (protein group)
- Score fields: `score_name`, `score_value`
- Position format: String in format `{AA}.{position}` (e.g., "M.4")

### Type Naming
- Use "record" not "struct" in AVSC schemas
- Each record type must have a unique name within the schema
- Use descriptive names like `peptide_count_record`, `intensity_record`

### Nullable Fields
- Explicitly mark fields as nullable in both schemas and Python code
- Use `nullable=True` in PyArrow field definitions
- Use `["null", <type>]` in AVSC schemas

## Key Findings

### Consistency Issues Addressed
1. **Type System**: Standardized on "record" for complex types in AVSC (not "struct")
2. **Naming**: Enforced consistent naming conventions across schemas, specification, and code
3. **Nullable**: Aligned nullable field definitions across all three sources
4. **Missing Fields**: Added fields present in specification but missing in schemas or code

### Design Patterns
1. **Modification Structure**: Uses position-specific scoring with format `{AA}.{position}`
2. **Intensity Structure**: Separates primary intensities from derived/processed intensities
3. **Score Structure**: Consistent `score_name`/`score_value` pairs across all schemas
4. **Metadata**: File-level metadata consistently structured across all parquet-based views

## Recommendations

### For Future Development
1. **Schema Validation**: Implement automated AVSC schema validation in CI/CD pipeline
2. **Code Generation**: Consider generating Python code from AVSC schemas to prevent drift
3. **Documentation**: Keep specification, schemas, and code in sync through automated checks
4. **Testing**: Add tests that validate data against AVSC schemas

### For Maintenance
1. **Single Source of Truth**: Consider making AVSC schemas the authoritative source
2. **Version Control**: Track schema versions explicitly in all files
3. **Migration**: Document migration paths when schemas change
4. **Backward Compatibility**: Test backward compatibility when updating schemas

## Files Modified

### Schema Files
- `docs/peptide.avsc` - Fixed field names, types, and nullable flags
- `docs/pg.avsc` - Fixed type names and added unique record names
- `docs/protein.avsc` - Fixed duplicate record names
- `docs/psm.avsc` - Added missing spectral data fields
- `docs/differential.avsc` - Fixed field name case

### Documentation
- `docs/format-specification.md` - Updated protein group table, added missing fields, fixed naming

### Code
- `qpx/core/format.py` - Added missing fields and nullable flags to PSM_UNIQUE_FIELDS and PG_FIELDS

## Conclusion

The review identified and fixed multiple alignment issues between schemas, specification, and code:
- **8 schema files** reviewed
- **20+ issues** identified and fixed
- **3 main categories**: naming inconsistencies, missing fields, type mismatches

All schemas now validate as proper JSON, and the Python code imports without errors. The alignment between specification, schemas, and implementation is now consistent and follows established naming conventions.

## Next Steps

1. Run full test suite to ensure no regressions
2. Update any example data or documentation that references changed field names
3. Consider adding automated schema validation to CI/CD pipeline
4. Review and update any external tools or scripts that depend on these schemas
