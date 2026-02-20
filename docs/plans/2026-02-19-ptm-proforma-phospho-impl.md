# ProForma/PTM Completion & Phospho Site Score Capture

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Ensure all converters properly build ProForma strings, populate the `modifications` field, and capture phospho site localization scores where available.

**Architecture:** Each converter already has `to_proforma()` in its constants. We extend `from_proforma()` with an optional `site_scores` parameter, wire modifications into all adapters that currently hardcode `None`, build ProForma for mzIdentML, and read phospho site columns from quantms/MaxQuant.

**Tech Stack:** Python, DuckDB, lxml (mzIdentML), existing ptm.py shared module

---

### Task 1: Extend `from_proforma()` to accept site scores

**Files:**
- Modify: `qpx/converters/ptm.py:130-206`
- Test: `tests/converters/test_ptm.py`

Add an optional `site_scores` parameter to `from_proforma()`. When provided as
`dict[int, list[dict]]` mapping position -> list of score dicts, merge them into
the `modification_position.scores` field instead of hardcoding `None`.

**Implementation:**
- Add parameter `site_scores: Optional[dict[int, list[dict]]] = None` to `from_proforma()`
- At line 196, replace `"scores": None` with `"scores": site_scores.get(position) if site_scores else None`
- Add tests: `test_from_proforma_with_site_scores` and `test_from_proforma_no_site_scores_unchanged`

---

### Task 2: Build ProForma in mzIdentML PSM adapter

**Files:**
- Modify: `qpx/converters/mzidentml/psm_adapter.py`
- Test: `tests/converters/test_mzidentml_psm.py` (or inline in existing test)

The mzIdentML adapter parses modifications from XML into structured dicts (line 134-174)
but never builds a ProForma string (line 387: `peptidoform = sequence`).

**Implementation:**
- Import `build_proforma` from `qpx.converters.ptm`
- Add a `_build_peptidoform()` helper that takes `sequence` and `modifications` list
  and builds ProForma using accession strings
- Normalize PSI-MOD accessions: if `accession` starts with `MOD:`, try mass lookup via `mass_to_unimod()`
- Wire into `_build_one_psm()` replacing the TODO at line 387

---

### Task 3: Populate `modifications` in MaxQuant PSM adapter

**Files:**
- Modify: `qpx/converters/maxquant/psm_adapter.py`

Currently `modifications` is always `None` (line 263). The adapter already calls
`to_proforma()` to build the peptidoform string.

**Implementation:**
- Import `from_proforma` from `qpx.converters.ptm`
- After building `peptidoform`, call `from_proforma(peptidoform, sequence)` to get structured mods
- Replace `"modifications": None` with the result

---

### Task 4: Populate `modifications` in quantms feature adapter

**Files:**
- Modify: `qpx/converters/quantms/feature_adapter.py`

Currently `modifications` is always `None` (line 383). The feature adapter has
access to `modifications_meta` from mzTab but never uses it for feature records.

**Implementation:**
- Import `from_proforma` from `qpx.converters.ptm`
- In `_build_feature_record()`, call `from_proforma(peptidoform, sequence, meta=modifications_meta)`
  when `peptidoform` != `sequence` (indicating modifications present)
- Store `modifications_meta` as an instance attribute so `_build_feature_record` can access it
- Replace `"modifications": None` with the result

---

### Task 5: Capture phospho site scores in quantms PSM adapter

**Files:**
- Modify: `qpx/converters/quantms/psm_adapter.py`
- Modify: `qpx/converters/quantms/constants.py`

The mzTab PSM table may contain `opt_global` columns with phospho site probabilities.
Common patterns: `opt_global_phosphors_score`, `opt_global_ptmrs_site_probability`,
`opt_global_luciphor_score`.

**Implementation:**
- Add `PHOSPHO_SITE_COLUMNS` mapping to quantms constants.py:
  ```python
  PHOSPHO_SITE_COLUMNS = {
      "opt_global_phosphors_score": "phosphors_site_probability",
      "opt_global_ptmrs_site_probability": "ptmrs_site_probability",
      "opt_global_luciphor_score": "luciphor_site_probability",
  }
  ```
- In `_transform_row()`, after building modifications, scan the PSM row for any
  opt_global phospho columns present in `actual_cols`
- Parse per-position probabilities and build `site_scores` dict
- Pass to `from_proforma(..., site_scores=site_scores)`
- Also add these as PSM-level additional_scores when a single value is present

---

### Task 6: Capture phospho site scores in MaxQuant PSM adapter

**Files:**
- Modify: `qpx/converters/maxquant/psm_adapter.py`
- Modify: `qpx/converters/maxquant/constants.py`

MaxQuant's msms.txt has `Phospho (STY) Probabilities` column with per-position
probability encoding like `_PEPTIDES(0.5)T(0.3)IDEK_`.

**Implementation:**
- Add `Phospho (STY) Probabilities` and `Phospho (STY) Score Diffs` to `_MQ_PSM_USECOLS`
- Add `parse_phospho_probabilities()` to maxquant constants.py:
  parses the MaxQuant probability string to extract `{position: probability}` mapping
- In `_transform_row()`, parse the probability column, build `site_scores` dict,
  and pass to `from_proforma(peptidoform, sequence, site_scores=site_scores)`
- Replace `"modifications": None` with the result

---

### Task 7: Add phospho score entries to `scores.py`

**Files:**
- Modify: `qpx/core/scores.py`

Add phospho-related scores to `_BUILTIN_SCORES` for proper ontology mapping.

**Implementation:**
- Add entries for: `phosphors_site_probability`, `ptmrs_site_probability`,
  `luciphor_site_probability`, `phospho_sty_probability`
- All are `higher_better: True` (higher probability = more confident localization)
- Add `andromeda_score` and `andromeda_delta_score` if not already present
