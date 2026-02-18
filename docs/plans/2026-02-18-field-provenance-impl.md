# Field Provenance Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Embed tool-specific column provenance into ontology.parquet so QPX datasets are self-documenting — every non-trivial field traces back to its original tool column name and PSI-MS CV term.

**Architecture:** Each converter package gets a `constants.py` declaring tool identity and ordered candidate column names per QPX field. A shared `resolve_columns()` utility resolves candidates against actual input data at runtime. `BaseConverter` builds ontology entries from the resolved mappings, merging with CV term lookup. Adapters reference constants instead of hardcoding column strings.

**Tech Stack:** Python, PyArrow, DuckDB, YAML schemas, pytest

**Design doc:** `docs/plans/2026-02-18-field-provenance-design.md`

---

### Task 1: Extend ontology.yaml schema

**Files:**
- Modify: `qpx/core/data/schemas/ontology.yaml`
- Test: `tests/unit/test_schema.py`

**Step 1: Write the failing test**

Add to `tests/unit/test_schema.py`:

```python
def test_ontology_schema_has_source_columns():
    """Verify ontology schema includes source_column_name and source_tool."""
    from qpx.core.data import OntologySchema
    schema = OntologySchema.get_arrow_schema()
    field_names = set(schema.names)
    assert "source_column_name" in field_names
    assert "source_tool" in field_names
```

**Step 2: Run test to verify it fails**

Run: `eval "$(/Users/yperez/miniconda3/bin/conda shell.bash hook 2>/dev/null)" && conda activate qpx && pytest tests/unit/test_schema.py::test_ontology_schema_has_source_columns -v`

Expected: FAIL — fields not yet in schema

**Step 3: Add new fields to ontology.yaml**

In `qpx/core/data/schemas/ontology.yaml`, add after the `description` field:

```yaml
  source_column_name:
    type: string
    doc: "Original column name in the tool output (e.g., Precursor.Quantity)"
  source_tool:
    type: string
    doc: "Tool name that produced this field (e.g., DIA-NN)"
```

**Step 4: Run test to verify it passes**

Run: `eval "$(/Users/yperez/miniconda3/bin/conda shell.bash hook 2>/dev/null)" && conda activate qpx && pytest tests/unit/test_schema.py::test_ontology_schema_has_source_columns -v`

Expected: PASS

**Step 5: Run all schema tests to check nothing breaks**

Run: `eval "$(/Users/yperez/miniconda3/bin/conda shell.bash hook 2>/dev/null)" && conda activate qpx && pytest tests/unit/test_schema.py tests/unit/test_writers.py -q --timeout=30`

Expected: All pass

**Step 6: Commit**

```bash
git add qpx/core/data/schemas/ontology.yaml tests/unit/test_schema.py
git commit -m "feat: add source_column_name and source_tool to ontology schema"
```

---

### Task 2: Add resolve_columns() utility and test

**Files:**
- Modify: `qpx/converters/base.py`
- Test: `tests/unit/test_resolve_columns.py` (new)

**Step 1: Write the failing test**

Create `tests/unit/test_resolve_columns.py`:

```python
"""Tests for resolve_columns() utility."""
import pytest


def test_resolve_columns_basic():
    from qpx.converters.base import resolve_columns
    mappings = {
        "intensity": ["Precursor.Quantity"],
        "rt": ["RT"],
    }
    available = {"Precursor.Quantity", "RT", "Run", "Sequence"}
    result = resolve_columns(mappings, available)
    assert result == {"intensity": "Precursor.Quantity", "rt": "RT"}


def test_resolve_columns_ordered_fallback():
    from qpx.converters.base import resolve_columns
    mappings = {
        "rt": ["RT", "RetentionTime"],
    }
    # Only the fallback name exists
    available = {"RetentionTime", "Sequence"}
    result = resolve_columns(mappings, available)
    assert result == {"rt": "RetentionTime"}


def test_resolve_columns_first_match_wins():
    from qpx.converters.base import resolve_columns
    mappings = {
        "rt": ["RT", "RetentionTime"],
    }
    # Both exist — first candidate wins
    available = {"RT", "RetentionTime"}
    result = resolve_columns(mappings, available)
    assert result == {"rt": "RT"}


def test_resolve_columns_missing_skipped():
    from qpx.converters.base import resolve_columns
    mappings = {
        "intensity": ["Precursor.Quantity"],
        "predicted_rt": ["Predicted.RT"],
    }
    available = {"Precursor.Quantity"}
    result = resolve_columns(mappings, available)
    assert result == {"intensity": "Precursor.Quantity"}
    assert "predicted_rt" not in result


def test_resolve_columns_empty():
    from qpx.converters.base import resolve_columns
    result = resolve_columns({}, {"col1", "col2"})
    assert result == {}
```

**Step 2: Run test to verify it fails**

Run: `eval "$(/Users/yperez/miniconda3/bin/conda shell.bash hook 2>/dev/null)" && conda activate qpx && pytest tests/unit/test_resolve_columns.py -v`

Expected: FAIL — `resolve_columns` not found

**Step 3: Implement resolve_columns()**

Add to `qpx/converters/base.py` before the `BaseConverter` class:

```python
def resolve_columns(
    field_mappings: dict[str, list[str]],
    available_columns: set[str],
) -> dict[str, str]:
    """Resolve QPX fields to actual tool column names.

    For each QPX field, tries candidates in order and returns the first match.
    Skips fields where no candidate matches the available columns.

    Args:
        field_mappings: QPX field name -> ordered list of candidate column names.
        available_columns: Set of column names present in the input data.

    Returns:
        Dict mapping QPX field name -> resolved tool column name.
    """
    resolved = {}
    for qpx_field, candidates in field_mappings.items():
        for candidate in candidates:
            if candidate in available_columns:
                resolved[qpx_field] = candidate
                break
    return resolved
```

**Step 4: Run test to verify it passes**

Run: `eval "$(/Users/yperez/miniconda3/bin/conda shell.bash hook 2>/dev/null)" && conda activate qpx && pytest tests/unit/test_resolve_columns.py -v`

Expected: All 5 PASS

**Step 5: Commit**

```bash
git add qpx/converters/base.py tests/unit/test_resolve_columns.py
git commit -m "feat: add resolve_columns() utility for field mapping resolution"
```

---

### Task 3: Extend _FIELD_CV_MAP and refactor field_ontology_entries()

**Files:**
- Modify: `qpx/core/scores.py` (lines 267-331)
- Test: `tests/unit/test_field_ontology.py` (new)

**Step 1: Write the failing test**

Create `tests/unit/test_field_ontology.py`:

```python
"""Tests for field_ontology_entries() with source provenance."""
import pytest


def test_field_ontology_entries_with_source_provenance():
    from qpx.core.scores import field_ontology_entries
    resolved = {"intensity": "Precursor.Quantity", "rt": "RT"}
    entries = field_ontology_entries(
        view="feature",
        resolved_mappings=resolved,
        tool_name="DIA-NN",
    )
    # Should have entries for resolved fields that have CV mappings
    rt_entries = [e for e in entries if e["field_name"] == "rt"]
    assert len(rt_entries) == 1
    assert rt_entries[0]["source_column_name"] == "RT"
    assert rt_entries[0]["source_tool"] == "DIA-NN"
    assert rt_entries[0]["ontology_accession"] == "MS:1000016"
    assert rt_entries[0]["view"] == "feature"


def test_field_ontology_entries_missing_cv_still_written():
    from qpx.core.scores import field_ontology_entries
    resolved = {"lfq": "Precursor.Normalised"}
    entries = field_ontology_entries(
        view="feature",
        resolved_mappings=resolved,
        tool_name="DIA-NN",
    )
    # lfq has no CV term but should still produce an entry
    assert len(entries) >= 1
    lfq_entries = [e for e in entries if e["field_name"] == "lfq"]
    assert len(lfq_entries) == 1
    assert lfq_entries[0]["source_column_name"] == "Precursor.Normalised"
    assert lfq_entries[0]["source_tool"] == "DIA-NN"
    # ontology_accession should be None (no CV term)
    assert lfq_entries[0]["ontology_accession"] is None


def test_field_ontology_entries_backward_compat():
    """Old call signature (no resolved_mappings) still works."""
    from qpx.core.scores import field_ontology_entries
    entries = field_ontology_entries(view="psm")
    # Should still produce entries for _FIELD_CV_MAP fields
    field_names = {e["field_name"] for e in entries}
    assert "posterior_error_probability" in field_names
    assert "rt" in field_names
    # Old entries should have source_column_name = None
    for e in entries:
        assert "source_column_name" in e
        assert "source_tool" in e
```

**Step 2: Run test to verify it fails**

Run: `eval "$(/Users/yperez/miniconda3/bin/conda shell.bash hook 2>/dev/null)" && conda activate qpx && pytest tests/unit/test_field_ontology.py -v`

Expected: FAIL — `field_ontology_entries()` doesn't accept `resolved_mappings` parameter

**Step 3: Extend _FIELD_CV_MAP and refactor field_ontology_entries()**

In `qpx/core/scores.py`:

1. Add `"intensity"` to `_FIELD_CV_MAP` (after line 297):

```python
    "intensity": {
        "ontology_name": "MS1 feature area",
        "ontology_accession": "MS:1002498",
        "ontology_source": "MS",
        "description": "Primary precursor intensity",
    },
```

2. Replace `field_ontology_entries()` (lines 306-331) with:

```python
def field_ontology_entries(
    view: str = "psm",
    resolved_mappings: dict[str, str] | None = None,
    tool_name: str | None = None,
) -> list[dict]:
    """Generate ontology.parquet rows for field-to-CV mappings.

    When *resolved_mappings* is provided, produces entries with
    ``source_column_name`` and ``source_tool`` populated from the
    actual converter resolution. Fields in the resolved mapping that
    have no CV term in ``_FIELD_CV_MAP`` still produce an entry
    (with ``ontology_accession=None``).

    For backward compatibility, calling without *resolved_mappings*
    emits entries for all ``_FIELD_CV_MAP`` fields with null source info.

    Args:
        view: The QPX view name (e.g. ``"psm"``, ``"feature"``).
        resolved_mappings: QPX field -> resolved tool column name.
        tool_name: Tool name string (e.g. ``"DIA-NN"``).

    Returns:
        List of dicts matching the ``OntologySchema``.
    """
    try:
        ontology_version = _get_ontology().version
    except Exception:
        ontology_version = None

    entries: list[dict] = []

    if resolved_mappings is not None:
        # New path: produce entries for each resolved field
        for field_name, source_column in sorted(resolved_mappings.items()):
            cv_info = _FIELD_CV_MAP.get(field_name)
            if cv_info is None:
                # Try OBO lookup
                cv_info = _lookup_from_obo(field_name)
            if cv_info is None:
                logger.warning(
                    "No CV term found for field '%s' (source: %s.%s)",
                    field_name,
                    tool_name or "unknown",
                    source_column,
                )
            entries.append({
                "field_name": field_name,
                "ontology_name": cv_info["ontology_name"] if cv_info else None,
                "ontology_accession": cv_info.get("ontology_accession") if cv_info else None,
                "ontology_source": cv_info.get("ontology_source") if cv_info else None,
                "ontology_version": ontology_version,
                "view": view,
                "description": cv_info["description"] if cv_info else f"{tool_name or 'unknown'} {source_column} (no CV term)",
                "source_column_name": source_column,
                "source_tool": tool_name,
            })
    else:
        # Backward-compatible path: emit _FIELD_CV_MAP entries
        for field_name, info in sorted(_FIELD_CV_MAP.items()):
            entries.append({
                "field_name": field_name,
                "ontology_name": info["ontology_name"],
                "ontology_accession": info["ontology_accession"],
                "ontology_source": info["ontology_source"],
                "ontology_version": ontology_version,
                "view": view,
                "description": info["description"],
                "source_column_name": None,
                "source_tool": None,
            })

    return entries
```

3. Also update `score_ontology_entries()` (lines 365-401) to include the two new keys with `None` values so all ontology entries have a consistent shape:

After `"description": info["description"],` add:
```python
                "source_column_name": None,
                "source_tool": None,
```

4. Similarly update `modification_ontology_entries()` (lines 334-362) — add the two `None` keys to the dict.

**Step 4: Run test to verify it passes**

Run: `eval "$(/Users/yperez/miniconda3/bin/conda shell.bash hook 2>/dev/null)" && conda activate qpx && pytest tests/unit/test_field_ontology.py -v`

Expected: All 3 PASS

**Step 5: Run existing tests to check nothing breaks**

Run: `eval "$(/Users/yperez/miniconda3/bin/conda shell.bash hook 2>/dev/null)" && conda activate qpx && pytest tests/unit/ -q --timeout=30 -k "not ae_view"`

Expected: All pass

**Step 6: Commit**

```bash
git add qpx/core/scores.py tests/unit/test_field_ontology.py
git commit -m "feat: extend field_ontology_entries() with source provenance support"
```

---

### Task 4: Create DIA-NN constants.py

**Files:**
- Create: `qpx/converters/diann/constants.py`
- Test: `tests/unit/test_diann_constants.py` (new)

**Step 1: Write the failing test**

Create `tests/unit/test_diann_constants.py`:

```python
"""Tests for DIA-NN converter constants."""


def test_diann_constants_structure():
    from qpx.converters.diann.constants import TOOL_NAME, TOOL_VERSIONS, FIELD_MAPPINGS
    assert TOOL_NAME == "DIA-NN"
    assert isinstance(TOOL_VERSIONS, str)
    assert "feature" in FIELD_MAPPINGS
    assert "pg" in FIELD_MAPPINGS


def test_diann_field_mappings_are_lists():
    from qpx.converters.diann.constants import FIELD_MAPPINGS
    for view, fields in FIELD_MAPPINGS.items():
        for qpx_field, candidates in fields.items():
            assert isinstance(candidates, list), (
                f"FIELD_MAPPINGS['{view}']['{qpx_field}'] must be a list"
            )
            assert len(candidates) > 0


def test_diann_key_fields_present():
    from qpx.converters.diann.constants import FIELD_MAPPINGS
    feature = FIELD_MAPPINGS["feature"]
    assert "intensity" in feature
    assert "posterior_error_probability" in feature
    assert "rt" in feature
    assert "pg_accessions" in feature
```

**Step 2: Run test to verify it fails**

Run: `eval "$(/Users/yperez/miniconda3/bin/conda shell.bash hook 2>/dev/null)" && conda activate qpx && pytest tests/unit/test_diann_constants.py -v`

Expected: FAIL — module not found

**Step 3: Create constants.py**

Create `qpx/converters/diann/constants.py`:

```python
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
```

**Step 4: Run test to verify it passes**

Run: `eval "$(/Users/yperez/miniconda3/bin/conda shell.bash hook 2>/dev/null)" && conda activate qpx && pytest tests/unit/test_diann_constants.py -v`

Expected: All 3 PASS

**Step 5: Commit**

```bash
git add qpx/converters/diann/constants.py tests/unit/test_diann_constants.py
git commit -m "feat: add DIA-NN converter constants with field mappings"
```

---

### Task 5: Create MaxQuant constants.py

**Files:**
- Create: `qpx/converters/maxquant/constants.py`

**Step 1: Create constants.py**

Create `qpx/converters/maxquant/constants.py`:

```python
"""MaxQuant converter constants — tool identity and field mappings."""

TOOL_NAME = "MaxQuant"
TOOL_VERSIONS = "2.x"

FIELD_MAPPINGS = {
    "feature": {
        "sequence":                    ["Sequence"],
        "modified_sequence":           ["Modified sequence"],
        "charge":                      ["Charge"],
        "run_file_name":               ["Raw file"],
        "is_decoy":                    ["Reverse"],
        "scan":                        ["MS/MS scan number"],
        "observed_mz":                 ["m/z"],
        "rt":                          ["Calibrated retention time"],
        "rt_start":                    ["Calibrated retention time start"],
        "rt_stop":                     ["Calibrated retention time finish"],
        "posterior_error_probability":  ["PEP"],
        "ion_mobility":                ["1/K0"],
        "pg_accessions":               ["Leading proteins"],
        "anchor_protein":              ["Leading razor protein"],
        "gg_names":                    ["Gene names"],
        "intensity":                   ["Intensity"],
        "andromeda_score":             ["Score"],
        "andromeda_delta_score":       ["Delta score"],
    },
    "psm": {
        "sequence":                    ["Sequence"],
        "modified_sequence":           ["Modified sequence"],
        "charge":                      ["Charge"],
        "run_file_name":               ["Raw file"],
        "is_decoy":                    ["Reverse"],
        "scan":                        ["Scan number", "MS/MS scan number"],
        "observed_mz":                 ["m/z"],
        "rt":                          ["Retention time"],
        "posterior_error_probability":  ["PEP"],
        "andromeda_score":             ["Score"],
        "andromeda_delta_score":       ["Delta score"],
    },
    "pg": {
        "pg_accessions":               ["Protein IDs"],
        "pg_names":                    ["Protein names"],
        "gg_accessions":               ["Gene names"],
        "anchor_protein":              ["Majority protein IDs"],
        "global_qvalue":               ["Q-value"],
        "is_decoy":                    ["Reverse"],
        "contaminant":                 ["Potential contaminant"],
        "sequence_coverage":           ["Sequence coverage [%]"],
        "molecular_weight":            ["Mol. weight [kDa]"],
        "peptide_count_total":         ["Peptides"],
        "peptide_count_unique":        ["Unique peptides"],
        "peptide_count_razor":         ["Razor + unique peptides"],
        "andromeda_score":             ["Score"],
        "intensity":                   ["Intensity"],
    },
}
```

**Step 2: Verify import works**

Run: `eval "$(/Users/yperez/miniconda3/bin/conda shell.bash hook 2>/dev/null)" && conda activate qpx && python -c "from qpx.converters.maxquant.constants import TOOL_NAME, FIELD_MAPPINGS; print(TOOL_NAME, len(FIELD_MAPPINGS))"`

Expected: `MaxQuant 3`

**Step 3: Commit**

```bash
git add qpx/converters/maxquant/constants.py
git commit -m "feat: add MaxQuant converter constants with field mappings"
```

---

### Task 6: Create FragPipe constants.py

**Files:**
- Create: `qpx/converters/fragpipe/constants.py`

**Step 1: Create constants.py**

```python
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
```

**Step 2: Verify import works**

Run: `eval "$(/Users/yperez/miniconda3/bin/conda shell.bash hook 2>/dev/null)" && conda activate qpx && python -c "from qpx.converters.fragpipe.constants import TOOL_NAME; print(TOOL_NAME)"`

Expected: `FragPipe`

**Step 3: Commit**

```bash
git add qpx/converters/fragpipe/constants.py
git commit -m "feat: add FragPipe converter constants with field mappings"
```

---

### Task 7: Create quantms constants.py

**Files:**
- Create: `qpx/converters/quantms/constants.py`

**Step 1: Create constants.py**

```python
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
```

**Step 2: Verify import works**

Run: `eval "$(/Users/yperez/miniconda3/bin/conda shell.bash hook 2>/dev/null)" && conda activate qpx && python -c "from qpx.converters.quantms.constants import TOOL_NAME; print(TOOL_NAME)"`

Expected: `quantms`

**Step 3: Commit**

```bash
git add qpx/converters/quantms/constants.py
git commit -m "feat: add quantms converter constants with field mappings"
```

---

### Task 8: Update BaseConverter with write_ontology()

**Files:**
- Modify: `qpx/converters/base.py` (lines 128-161)
- Test: `tests/unit/test_resolve_columns.py` (extend)

**Step 1: Write the failing test**

Add to `tests/unit/test_resolve_columns.py`:

```python
def test_base_converter_write_ontology_exists():
    """write_ontology method should exist on BaseConverter."""
    from qpx.converters.base import BaseConverter
    assert hasattr(BaseConverter, "write_ontology")


def test_base_converter_write_score_ontology_still_works():
    """write_score_ontology should still exist for backward compat."""
    from qpx.converters.base import BaseConverter
    assert hasattr(BaseConverter, "write_score_ontology")
```

**Step 2: Run test to verify it fails**

Run: `eval "$(/Users/yperez/miniconda3/bin/conda shell.bash hook 2>/dev/null)" && conda activate qpx && pytest tests/unit/test_resolve_columns.py::test_base_converter_write_ontology_exists -v`

Expected: FAIL — `write_ontology` not found

**Step 3: Add write_ontology() to BaseConverter**

In `qpx/converters/base.py`, add a new method after `write_score_ontology()`. Keep `write_score_ontology()` for backward compatibility:

```python
    def write_ontology(
        self,
        output_path: str | Path,
        view: str = "psm",
        resolved_mappings: dict[str, str] | None = None,
        tool_name: str | None = None,
        extra_entries: list[dict] | None = None,
    ) -> Path | None:
        """Write ontology.parquet with score + field provenance entries.

        Combines:
        1. Score ontology entries (from _discovered_scores)
        2. Field provenance entries (from resolved column mappings)
        3. Any extra entries passed by the converter

        Args:
            output_path: Path for the ontology Parquet file.
            view: QPX view name.
            resolved_mappings: QPX field -> resolved tool column name.
            tool_name: Tool name string (e.g. "DIA-NN").
            extra_entries: Additional ontology entries to include.

        Returns:
            Path to the written file, or None if no entries to write.
        """
        from qpx.core.scores import field_ontology_entries

        entries = score_ontology_entries(self._discovered_scores, view=view)
        if resolved_mappings:
            entries.extend(
                field_ontology_entries(
                    view=view,
                    resolved_mappings=resolved_mappings,
                    tool_name=tool_name,
                )
            )
        if extra_entries:
            entries.extend(extra_entries)
        if not entries:
            return None

        from qpx.writers.ontology import OntologyWriter

        output_path = Path(output_path)
        with OntologyWriter(output_path, creator="qpx") as writer:
            writer.write_batch(entries)
        self.logger.info(
            "Wrote %d ontology entries (%d scores, %d field mappings) to %s",
            len(entries),
            len(self._discovered_scores),
            len(resolved_mappings or {}),
            output_path,
        )
        return output_path
```

**Step 4: Run test to verify it passes**

Run: `eval "$(/Users/yperez/miniconda3/bin/conda shell.bash hook 2>/dev/null)" && conda activate qpx && pytest tests/unit/test_resolve_columns.py -v`

Expected: All PASS

**Step 5: Commit**

```bash
git add qpx/converters/base.py tests/unit/test_resolve_columns.py
git commit -m "feat: add write_ontology() to BaseConverter with field provenance"
```

---

### Task 9: Update DIA-NN converter to use constants and write_ontology()

This is the first converter migration. It serves as the reference pattern for the others.

**Files:**
- Modify: `qpx/converters/diann/converter.py`
- Modify: `qpx/converters/diann/feature_adapter.py`
- Modify: `qpx/converters/diann/pg_adapter.py`

**Step 1: Update converter.py to use write_ontology()**

In `qpx/converters/diann/converter.py`:

1. Add import: `from qpx.converters.diann.constants import TOOL_NAME, FIELD_MAPPINGS`
2. Add import: `from qpx.converters.base import resolve_columns`
3. In `convert_features()`: After loading the report table into DuckDB, resolve columns:
   ```python
   cols = {c[0] for c in self._conn.execute("SELECT column_name FROM information_schema.columns WHERE table_name='report'").fetchall()}
   self._resolved_feature = resolve_columns(FIELD_MAPPINGS.get("feature", {}), cols)
   ```
4. In `write_ontology()`: Replace `field_ontology_entries(view="feature")` call with passing `resolved_mappings=self._resolved_feature, tool_name=TOOL_NAME` to the base `write_ontology()` method or combine all resolved mappings.

**Step 2: Update feature_adapter.py SQL aliases to use constants**

In `qpx/converters/diann/feature_adapter.py`:

The DIA-NN adapter builds SQL queries that alias DIA-NN columns to QPX names. The SQL column references (e.g., `"Precursor.Quantity" AS intensity`) should be constructed from the constants. However, since DIA-NN uses SQL aliases, the `_build_feature_record()` reads the aliased names, not the raw DIA-NN names. The change here is:

1. Import: `from qpx.converters.diann.constants import FIELD_MAPPINGS`
2. Use `FIELD_MAPPINGS["feature"]` to build the `_DIANN_REPORT_COLS` list dynamically (extract all unique candidate column names) and to construct SQL aliases.

**Important**: The DIA-NN adapter's SQL aliasing pattern means adapters read QPX-aliased names after the SQL transform. The constants provide the original column names for the SQL aliases AND for ontology provenance. The `_build_feature_record()` method doesn't need to change its `row.get()` calls since it reads the aliased names.

**Step 3: Update pg_adapter.py similarly**

Same pattern — import constants, use for SQL construction and column resolution.

**Step 4: Run DIA-NN converter tests**

Run: `eval "$(/Users/yperez/miniconda3/bin/conda shell.bash hook 2>/dev/null)" && conda activate qpx && pytest tests/converters/test_diann_converter.py -v --timeout=60 -k "test_key_columns or test_schema_validation" 2>&1 | head -30`

Expected: Tests still pass (functional behavior unchanged)

**Step 5: Commit**

```bash
git add qpx/converters/diann/
git commit -m "refactor: DIA-NN converter uses constants.py for column names"
```

---

### Task 10: Update MaxQuant converter to use constants

**Files:**
- Modify: `qpx/converters/maxquant/converter.py`
- Modify: `qpx/converters/maxquant/feature_adapter.py`
- Modify: `qpx/converters/maxquant/pg_adapter.py`

Same pattern as Task 9. MaxQuant reads raw column names directly via `row.get()`, so the change is more straightforward:

1. Import constants
2. Resolve columns against the actual evidence.txt / proteinGroups.txt columns
3. Replace hardcoded strings in `row.get()` with `resolved["field_name"]`
4. Pass resolved mappings to `write_ontology()`

**Step 1: Update feature_adapter.py**

Replace hardcoded column names with constants references. Example:

```python
# Before:
pg_acc_raw = str(row.get("Leading proteins", ""))

# After:
from qpx.converters.maxquant.constants import FIELD_MAPPINGS
# (resolved dict passed to _transform_row)
pg_acc_raw = str(row.get(resolved.get("pg_accessions", "Leading proteins"), ""))
```

**Step 2: Run MaxQuant converter tests**

Run: `eval "$(/Users/yperez/miniconda3/bin/conda shell.bash hook 2>/dev/null)" && conda activate qpx && pytest tests/converters/test_maxquant_converter.py -q --timeout=30`

Expected: All pass

**Step 3: Commit**

```bash
git add qpx/converters/maxquant/
git commit -m "refactor: MaxQuant converter uses constants.py for column names"
```

---

### Task 11: Update FragPipe converter to use constants

Same pattern as Task 10.

**Files:**
- Modify: `qpx/converters/fragpipe/converter.py`
- Modify: `qpx/converters/fragpipe/feature_adapter.py`
- Modify: `qpx/converters/fragpipe/pg_adapter.py`

**Commit:**

```bash
git add qpx/converters/fragpipe/
git commit -m "refactor: FragPipe converter uses constants.py for column names"
```

---

### Task 12: Update quantms converter to use constants

**Files:**
- Modify: `qpx/converters/quantms/converter.py`
- Modify: `qpx/converters/quantms/feature_adapter.py`
- Modify: `qpx/converters/quantms/pg_adapter.py`

This converter already has `_detect_msstats_columns()` and `_detect_column()` patterns. The change is to migrate these to use `FIELD_MAPPINGS` from constants.py via `resolve_columns()`, replacing the inline candidate lists.

**Commit:**

```bash
git add qpx/converters/quantms/
git commit -m "refactor: quantms converter uses constants.py for column names"
```

---

### Task 13: Integration test — ontology.parquet contains source provenance

**Files:**
- Test: `tests/unit/test_ontology_provenance.py` (new)

**Step 1: Write integration test**

```python
"""Integration test: ontology.parquet entries include source provenance."""
import tempfile
from pathlib import Path

import pyarrow.parquet as pq
import pytest


def test_ontology_has_source_provenance(tmp_path):
    """Write ontology entries with provenance and verify Parquet output."""
    from qpx.core.scores import field_ontology_entries
    from qpx.writers.ontology import OntologyWriter

    resolved = {"intensity": "Precursor.Quantity", "rt": "RT"}
    entries = field_ontology_entries(
        view="feature",
        resolved_mappings=resolved,
        tool_name="DIA-NN",
    )
    assert len(entries) > 0

    out = tmp_path / "ontology.parquet"
    with OntologyWriter(out, creator="qpx") as writer:
        writer.write_batch(entries)

    table = pq.read_table(out)
    assert "source_column_name" in table.column_names
    assert "source_tool" in table.column_names

    df = table.to_pandas()
    rt_rows = df[df["field_name"] == "rt"]
    assert len(rt_rows) == 1
    assert rt_rows.iloc[0]["source_column_name"] == "RT"
    assert rt_rows.iloc[0]["source_tool"] == "DIA-NN"
    assert rt_rows.iloc[0]["ontology_accession"] == "MS:1000016"
```

**Step 2: Run test**

Run: `eval "$(/Users/yperez/miniconda3/bin/conda shell.bash hook 2>/dev/null)" && conda activate qpx && pytest tests/unit/test_ontology_provenance.py -v`

Expected: PASS (all infrastructure in place from Tasks 1-3)

**Step 3: Commit**

```bash
git add tests/unit/test_ontology_provenance.py
git commit -m "test: add integration test for ontology source provenance"
```

---

### Task 14: Update docs

**Files:**
- Modify: `docs/spec/ontology.md`
- Modify: `docs/spec/tool-mappings.md`

**Step 1: Update ontology.md**

Add the two new fields to the field table and add an example showing source provenance.

**Step 2: Update tool-mappings.md**

Add a note that `ontology.parquet` is the authoritative, machine-readable source for field mappings, while `tool-mappings.md` serves as a human-readable reference.

**Step 3: Commit**

```bash
git add docs/spec/ontology.md docs/spec/tool-mappings.md
git commit -m "docs: document source provenance fields in ontology spec"
```

---

### Task 15: Final verification

**Step 1: Run full test suite**

Run: `eval "$(/Users/yperez/miniconda3/bin/conda shell.bash hook 2>/dev/null)" && conda activate qpx && pytest tests/unit/ tests/integration/ -q --timeout=30 -k "not ae_view and not s3"`

Expected: All pass

**Step 2: Run converter tests**

Run: `eval "$(/Users/yperez/miniconda3/bin/conda shell.bash hook 2>/dev/null)" && conda activate qpx && pytest tests/converters/test_diann_converter.py tests/converters/test_maxquant_converter.py -q --timeout=60`

Expected: All pass (may have fixture timeouts for long OBO lookups — these are pre-existing)
