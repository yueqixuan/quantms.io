# mzIdentML Converter Improvements

## Date: 2026-02-15

## Overview

Four improvements to the mzIdentML conversion pipeline: PEP extraction, RT parsing,
MGF index fallback, and PRIDE metadata enrichment.

## Task 1: Extract posterior_error_probability (PEP)

**Files:** `qpx/core/cv_terms.py`, `qpx/converters/mzidentml/psm_adapter.py`

- Add `CV_PEP = "MS:1001493"` and `CV_PEP_GLOBAL = "MS:1002352"` to cv_terms.py
- In `_build_one_psm()`, scan SII cvParams for MS:1001493 or MS:1002352
- If found, populate the dedicated `posterior_error_probability` field (float)
- The score still also lands in `additional_scores` for backward compatibility
- Add PEP accessions to `SKIP_SCORE_ACCESSIONS` so the dedicated field is the canonical location
  (actually NO — keep it in additional_scores too, just also populate the dedicated field)

**Tests:** Add test in `tests/converters/test_mzidentml_converter.py` with synthetic mzid
containing MS:1001493 cvParam. Verify both `posterior_error_probability` field and
`additional_scores` are populated.

## Task 2: Extract retention time (RT)

**Files:** `qpx/core/cv_terms.py`, `qpx/converters/mzidentml/psm_adapter.py`,
`qpx/converters/mzidentml/mgf_parser.py`, `qpx/converters/mzidentml/converter.py`

**Two sources, priority order:**

1. **mzIdentML** — `MS:1000016` (scan start time) on `SpectrumIdentificationResult` cvParams.
   If `unitName="minute"`, multiply by 60. Store in parsed dict and assign to PSM `rt` field.

2. **MGF** — Parse `RTINSECONDS=` header in `mgf_parser.py`, store in spectrum dict.
   In `converter._attach_spectra()`, fill `rt` from MGF only if not already set from mzIdentML.

**Tests:** Synthetic mzid with RT cvParam. MGF with RTINSECONDS. Verify RT populates correctly,
and that mzIdentML RT takes priority over MGF RT.

## Task 3: MGF index fallback for scan matching

**Files:** `qpx/converters/mzidentml/mgf_parser.py`, `qpx/converters/mzidentml/converter.py`

- `MgfSpectraIndex` adds `_by_position: dict[int, dict]` tracking 0-based insertion order
- New method `get_spectrum_by_index(index: int) -> dict | None`
- `converter._attach_spectra()` tries scan-based lookup first, then index-based fallback
- Log how many were matched by scan vs. by index

**Tests:** Create MGF with non-standard scan numbering. Verify index fallback works.

## Task 4: PRIDE metadata enrichment

**Files:** `qpx/core/pride.py` (NEW), `qpx/dataset.py`

### `qpx/core/pride.py`

```python
def fetch_pride_metadata(accession: str) -> dict:
    """Fetch project metadata from PRIDE REST API v3.

    Returns dict with: project_title, project_description, pubmed_id,
    doi, organisms, instruments, keywords, submission_date, publication_date.
    """
```

Uses `urllib.request.urlopen()` to call
`https://www.ebi.ac.uk/pride/ws/archive/v3/projects/{accession}`.

PRIDE API v3 JSON fields mapped:
- `title` → `project_title`
- `projectDescription` → `project_description`
- `references[0].pubmedID` → `pubmed_id` (as string)
- `doi` → `doi`
- `organisms[].name` → `organisms` (list of strings)
- `instruments[].name` → `instruments` (list of strings)
- `keywords` → `keywords`
- `submissionDate` → `submission_date`
- `publicationDate` → `publication_date`

### `Dataset.enrich_from_pride()`

```python
def enrich_from_pride(self, project_accession: str | None = None) -> dict:
    """Fetch PRIDE metadata and update dataset.parquet."""
```

1. Read project_accession from existing dataset.parquet if not provided
2. Call `fetch_pride_metadata(accession)`
3. Update dataset.parquet with project_title, project_description, pubmed_id
4. Return the fetched metadata dict

### CLI integration

Add `--enrich-pride` flag to `qpxc convert mzidentml`. When set with `--project-accession`,
auto-enrich after conversion.

**Tests:** Mock the PRIDE API response, verify Dataset.enrich_from_pride() updates fields.
Integration test with real PXD054720 data.
