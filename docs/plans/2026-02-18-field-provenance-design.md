# Field Provenance: Self-Documenting Column Mappings

**Date**: 2026-02-18
**Status**: Approved

## Problem

QPX normalizes tool-specific column names into consensus fields (e.g., DIA-NN's `Precursor.Quantity` becomes `intensity`, MaxQuant's `PEP` becomes `posterior_error_probability`). Once converted, users cannot trace a QPX field back to its original tool column name or determine whether a PSI-MS CV term exists for it. The `tool-mappings.md` doc captures this knowledge statically, but it is not embedded in the data.

## Goal

Make QPX datasets self-documenting: a user reading `ontology.parquet` can see, for every non-trivial field, the original tool column name, the tool identity, and the PSI-MS CV term (if one exists). During conversion, warn when no CV term is found.

## Design Decisions

- **Scope**: Non-trivial fields only (intensities, scores, PEP, q-values, m/z, RT variants, protein accessions). Obvious mappings like `Sequence` -> `sequence` are excluded.
- **Storage**: Extend `ontology.parquet` with two new columns (`source_column_name`, `source_tool`).
- **Missing CV terms**: Warn during conversion, still write the ontology entry with `ontology_accession = null`.
- **Version handling**: Ordered candidate lists per field. Runtime detection picks the first match from the input data.
- **Single source of truth**: Adapters reference `constants.py` for column names instead of hardcoding strings.

## Schema Change: ontology.parquet

Two new nullable string fields:

| Field | Type | Description |
|-------|------|-------------|
| `source_column_name` | string, null | Original column name in tool output (e.g., `Precursor.Quantity`) |
| `source_tool` | string, null | Tool name (e.g., `DIA-NN`) |

Example row:

```
field_name:          intensity
source_column_name:  Precursor.Quantity
source_tool:         DIA-NN
ontology_name:       MS1 feature area
ontology_accession:  MS:1002498
ontology_source:     MS
ontology_version:    4.1.235
view:                feature
description:         Primary precursor intensity from DIA-NN Precursor.Quantity
```

When no CV term exists:

```
field_name:          lfq
source_column_name:  Precursor.Normalised
source_tool:         DIA-NN
ontology_name:       null
ontology_accession:  null
ontology_source:     null
view:                feature
description:         DIA-NN normalised precursor intensity (no CV term found)
```

## Converter constants.py

Each converter package gets a `constants.py` with tool identity and field mappings. Field mappings use ordered candidate lists to handle version differences.

```python
# qpx/converters/diann/constants.py

TOOL_NAME = "DIA-NN"
TOOL_VERSIONS = "1.8+"

# QPX field -> ordered list of candidate tool column names.
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
    },
    "pg": {
        "intensity":                   ["PG.Quantity"],
        "pg_accessions":               ["Protein.Group"],
        "global_qvalue":               ["Global.PG.Q.Value"],
        "lfq":                         ["PG.MaxLFQ"],
    },
}
```

Version variance example (hypothetical):

```python
FIELD_MAPPINGS = {
    "feature": {
        "rt": ["RT", "RetentionTime"],  # RT in v1.8+, RetentionTime in v1.7
    },
}
```

## Runtime Flow

```
constants.py          adapter.py                BaseConverter
    |                     |                          |
 FIELD_MAPPINGS  -->  resolve_columns()     -->  _build_field_ontology_entries()
 (candidates)        tries candidates            (CV lookup via _FIELD_CV_MAP + OBO)
                     against actual input         (warn if no CV term)
                     records which found                |
                          |                      write_ontology()
                     uses resolved names               |
                     in data extraction         ontology.parquet
                     (no hardcoded strings)    (scores + field provenance)
```

### resolve_columns()

A shared utility (in `BaseConverter` or a helper module) that takes `FIELD_MAPPINGS[view]` and the actual input column set, returns `dict[str, str]` mapping QPX field -> resolved tool column name.

```python
def resolve_columns(
    field_mappings: dict[str, list[str]],
    available_columns: set[str],
) -> dict[str, str]:
    """Resolve QPX fields to actual tool column names.

    Returns dict mapping QPX field name -> resolved tool column name.
    Skips fields where no candidate matches.
    """
    resolved = {}
    for qpx_field, candidates in field_mappings.items():
        for candidate in candidates:
            if candidate in available_columns:
                resolved[qpx_field] = candidate
                break
    return resolved
```

### _build_field_ontology_entries()

New method on `BaseConverter` that merges resolved mappings with CV term lookup:

```python
def _build_field_ontology_entries(
    self,
    view: str,
    resolved_mappings: dict[str, str],
    tool_name: str,
) -> list[dict]:
    """Build ontology entries for resolved field mappings."""
    entries = []
    for qpx_field, source_column in resolved_mappings.items():
        cv_info = _FIELD_CV_MAP.get(qpx_field) or _lookup_from_obo(qpx_field)
        if cv_info is None:
            self.logger.warning(
                "No CV term found for field '%s' (source: %s.%s)",
                qpx_field, tool_name, source_column,
            )
        entries.append({
            "field_name": qpx_field,
            "source_column_name": source_column,
            "source_tool": tool_name,
            "ontology_name": cv_info["ontology_name"] if cv_info else None,
            "ontology_accession": cv_info["ontology_accession"] if cv_info else None,
            "ontology_source": cv_info["ontology_source"] if cv_info else None,
            "ontology_version": ontology_version,
            "view": view,
            "description": cv_info["description"] if cv_info else f"{tool_name} {source_column} (no CV term)",
        })
    return entries
```

### write_ontology() (renamed from write_score_ontology)

Combines score entries + field provenance entries + any extra entries:

```python
def write_ontology(self, output_path, view, resolved_mappings, tool_name, extra_entries=None):
    entries = score_ontology_entries(self._discovered_scores, view=view)
    entries.extend(self._build_field_ontology_entries(view, resolved_mappings, tool_name))
    if extra_entries:
        entries.extend(extra_entries)
    # ... write to parquet
```

## Adapter Changes

Adapters stop hardcoding column names and instead use the resolved mapping:

```python
# Before (hardcoded):
pg_acc_raw = str(row.get("Protein.Group", ""))
intensity_val = safe_float(row.get("Precursor.Quantity"))

# After (from constants):
pg_acc_raw = str(row.get(resolved["pg_accessions"], ""))
intensity_val = safe_float(row.get(resolved["intensity"]))
```

## Files Changed

| File | Change |
|------|--------|
| `qpx/converters/diann/constants.py` | New: tool name, versions, field mappings |
| `qpx/converters/maxquant/constants.py` | New: tool name, versions, field mappings |
| `qpx/converters/fragpipe/constants.py` | New: tool name, versions, field mappings |
| `qpx/converters/quantms/constants.py` | New: tool name, versions, field mappings |
| `qpx/core/data/schemas/ontology.yaml` | Add `source_column_name` and `source_tool` fields |
| `qpx/core/scores.py` | Extend `_FIELD_CV_MAP`, add CV entries for intensity/lfq terms; extend `field_ontology_entries()` signature |
| `qpx/converters/base.py` | Add `resolve_columns()`, `_build_field_ontology_entries()`; rename `write_score_ontology()` -> `write_ontology()` |
| `qpx/converters/diann/feature_adapter.py` | Use `resolved` dict instead of hardcoded column names |
| `qpx/converters/diann/pg_adapter.py` | Use `resolved` dict instead of hardcoded column names |
| `qpx/converters/diann/converter.py` | Pass resolved mappings to `write_ontology()` |
| `qpx/converters/maxquant/feature_adapter.py` | Use `resolved` dict instead of hardcoded column names |
| `qpx/converters/maxquant/pg_adapter.py` | Use `resolved` dict instead of hardcoded column names |
| `qpx/converters/maxquant/converter.py` | Pass resolved mappings to `write_ontology()` |
| `qpx/converters/fragpipe/feature_adapter.py` | Use `resolved` dict instead of hardcoded column names |
| `qpx/converters/fragpipe/pg_adapter.py` | Use `resolved` dict instead of hardcoded column names |
| `qpx/converters/fragpipe/converter.py` | Pass resolved mappings to `write_ontology()` |
| `qpx/converters/quantms/feature_adapter.py` | Use `resolved` dict instead of hardcoded column names |
| `qpx/converters/quantms/pg_adapter.py` | Use `resolved` dict instead of hardcoded column names |
| `qpx/converters/quantms/converter.py` | Pass resolved mappings to `write_ontology()` |
| `docs/spec/ontology.md` | Document new fields |
| `docs/spec/tool-mappings.md` | Reference ontology.parquet as authoritative source |
| Tests | Unit tests for `resolve_columns()`, integration tests for ontology output |
