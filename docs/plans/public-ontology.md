# PublicOntology: Bionty-style CV Term Registry for QPX

## Context

QPX uses **12+ ontology sources** — scores, modifications (UNIMOD), instruments (PSI-MS), organisms (NCBITaxon), tissues (UBERON), diseases (MONDO), cell types (CL), enzymes, dissociation methods, cross-linkers (XLMOD), and more. Currently these terms are resolved ad-hoc through `CVTermRegistry` in `qpx/core/obo.py` with JSON caching.

**Problems identified by code audit:**
- `write_score_ontology()` exists in `BaseConverter` but is **never called** — ontology.parquet is never generated
- SDRF converter sets instrument/enzyme/dissociation `accession: ""` — empty accessions in run.parquet
- No offline fallback — network failure + expired cache = silent empty lookups
- No versioning — no way to know which PSI-MS version was used
- Hardcoded `_BUILTIN_SCORES` dict requires code changes + pip release for new tools

**Goal:** A bionty-style `PublicOntology` system that:
- Stores ontology data **outside the code** in a repo-level `/data` folder
- Ships a **bundled copy** with the pip package for offline/first-use
- **Auto-checks** for newer versions in the repo (no release needed)
- Uses **DuckDB** backend (consistent with rest of QPX)
- Compatible with **laminDB** artifact registration

---

## 1. Architecture: Three-Layer Resolution (Hybrid Bundle + Auto-Update)

```
Priority 1: Repo checkout    data/ontology/psi_ms.parquet       (dev working on repo)
Priority 2: Auto-updated     ~/.qpx/ontology/psi_ms.parquet     (pip users, auto-refreshed)
Priority 3: Bundled fallback  qpx/core/ontology/data/psi_ms.parquet  (offline/first-use)
```

### Layer 1 — Repo checkout (developers)

Walk up from the package's `__file__` looking for `data/ontology/` at the repo root. If found, read Parquet files directly. Developers always get the latest committed data without any network access.

### Layer 2 — Auto-updated cache (pip users)

For pip-installed users, Parquet files are cached to `~/.qpx/ontology/`. The library checks `versions.yaml` from GitHub every 24 hours (blocking, 5-second timeout). If a newer version is available, the new Parquet is downloaded before proceeding. On network failure, the existing cached copy is used with a warning.

### Layer 3 — Bundled fallback (offline / first use)

Small Parquet files (~450 KB total for Phase 1) are shipped inside the Python package at `qpx/core/ontology/data/`. These are always available, work offline immediately, and never require network access. This is the last-resort fallback.

### Resolution flow

```
PublicOntology("psi_ms")
│
├── 1. Repo checkout? Walk up from __file__ for data/ontology/psi_ms.parquet
│     └── Found → use it (done)
│
├── 2. Cache exists? ~/.qpx/ontology/psi_ms.parquet
│     ├── YES + version check due (>24h):
│     │     ├── Fetch versions.yaml (blocking, 5s timeout)
│     │     │    ├── New version → download new Parquet → use it
│     │     │    └── Same version → update timestamp → use cached
│     │     └── Network error → warn, use cached
│     ├── YES + not due → use cached (fast path)
│     └── NO → try to download from repo
│           ├── Success → cache + use
│           └── Fail → fall through to bundled
│
└── 3. Bundled: qpx/core/ontology/data/psi_ms.parquet (always available)
```

### Key design choices

- **Blocking version check** — always block until fresh. 5s timeout prevents hangs.
- **`auto_update=False`** — parameter to skip version checking entirely (CI, air-gapped, tests).
- **No new dependencies** — uses `urllib.request`, DuckDB, PyArrow (all already in QPX).
- **DuckDB backend** — consistent with rest of QPX. Scales to Phase 2 ontologies (UBERON 15K, MONDO 30K terms).

---

## 2. Schema Changes

### run.yaml — simplify to match sample pattern

Run metadata fields change from `ontology_term` struct to plain strings. The CV accession resolution moves to `ontology.parquet`.

```yaml
# BEFORE                              # AFTER
instrument:                           instrument:
  type: ontology_term         →         type: string
enzymes:                              enzymes:
  type: "list<ontology_term>" →         type: "list<string>"
dissociation_method:                  dissociation_method:
  type: ontology_term         →         type: string
modification_parameters:              modification_parameters:
  type: "list<modification_param>"      type: "list<modification_param>"  # UNCHANGED
additional_terms:                     # REMOVED (was never populated)
  type: "list<ontology_term>"
```

Modifications keep the struct because they have position, fixed/variable, and other fields that can't be reduced to a plain string.

### ontology.yaml — add ontology_version column

```yaml
fields:
  field_name:
    type: string
    required: true
    doc: "snake_case QPX field or score name, or term value for run/sample metadata"
  ontology_name:
    type: string
    doc: "Proper ontology term name"
  ontology_accession:
    type: string
    doc: "Ontology accession identifier (e.g., MS:1002252)"
  ontology_source:
    type: string
    doc: "Ontology prefix (MS, UBERON, UNIMOD, etc.)"
  ontology_version:
    type: string
    doc: "Version of the ontology used for resolution (e.g., 4.1.235)"
  view:
    type: string
    required: true
    doc: "QPX view this field belongs to (psm, feature, pg, run, sample)"
  description:
    type: string
    doc: "Human-readable description"
```

PK stays `[field_name, view]`.

### How ontology.parquet is populated

| field_name | ontology_name | ontology_accession | ontology_source | ontology_version | view | description |
|---|---|---|---|---|---|---|
| percolator_score | percolator:score | MS:1001492 | MS | 4.1.235 | psm | PSM-level score |
| Q Exactive HF | Q Exactive HF | MS:1002523 | MS | 4.1.235 | run | Mass spectrometer |
| Trypsin | Trypsin | MS:1001251 | MS | 4.1.235 | run | Proteolytic enzyme |
| HCD | beam-type collision... | MS:1000422 | MS | 4.1.235 | run | Fragmentation method |

Score field_names are normalized (snake_case). Instrument/enzyme/dissociation field_names use the original value as-is.

---

## 3. File Layout

### Repo-level data (NOT shipped in pip package)

```
data/ontology/
├── versions.yaml             # Current versions of all ontologies
├── psi_ms.parquet            # PSI-MS CV terms (~400 KB)
└── pride_cv.parquet          # PRIDE CV terms (~50 KB)
```

### Package code + bundled data (shipped in pip package)

```
qpx/core/ontology/
├── __init__.py               # from .registry import PublicOntology
├── sources.yaml              # Declares URLs, OBO sources, check interval
├── registry.py               # PublicOntology class
├── build.py                  # Maintainer build script (OBO → Parquet)
└── data/                     # Bundled Parquet files (fallback)
    ├── psi_ms.parquet        # ~400 KB
    └── pride_cv.parquet      # ~50 KB
```

### User cache (auto-managed)

```
~/.qpx/ontology/
├── versions.yaml             # Copy of repo versions.yaml (last fetched)
├── psi_ms.parquet            # Cached copy
├── pride_cv.parquet          # Cached copy
└── cache_meta.json           # When versions were last checked
```

---

## 4. sources.yaml — Package Configuration

```yaml
repo_base_url: "https://raw.githubusercontent.com/bigbio/quantms.io/main/data/ontology"
cache_dir: "~/.qpx/ontology"
check_interval_hours: 24
version_check_timeout_seconds: 5

ontologies:
  psi_ms:
    display_name: "PSI-MS Controlled Vocabulary"
    obo_url: "https://raw.githubusercontent.com/HUPO-PSI/psi-ms-CV/master/psi-ms.obo"
    source_prefix: "MS"

  pride_cv:
    display_name: "PRIDE Controlled Vocabulary"
    obo_url: "https://raw.githubusercontent.com/PRIDE-Archive/pride-ontology/master/pride_cv.obo"
    source_prefix: "PRIDE"

  # --- Phase 2 ---
  # unimod:
  #   display_name: "UNIMOD Modifications"
  #   obo_url: "https://www.unimod.org/obo/unimod.obo"
  #   source_prefix: "UNIMOD"
```

---

## 5. terms.parquet Schema (the bundled/cached ontology data)

| Column          | Arrow Type      | Description                               |
|-----------------|-----------------|-------------------------------------------|
| accession       | string (PK)     | CV accession (e.g., MS:1001492)           |
| name            | string          | Original CV term name                     |
| normalized_name | string          | snake_case normalized name                |
| definition      | string          | Term definition from OBO                  |
| source          | string          | "MS", "PRIDE", "UNIMOD", etc.             |
| is_score        | bool            | Whether this is a score term              |
| higher_better   | bool (nullable) | Score direction (null = unknown/non-score)|
| is_a            | list\<string\>  | Parent term accessions                    |
| synonyms        | list\<string\>  | Synonym strings from OBO                  |
| is_obsolete     | bool            | Whether term is obsolete                  |

Parquet metadata (footer):
```
qpx_ontology_source: "psi_ms"
qpx_ontology_version: "4.1.235"
qpx_built_at: "2025-03-15T10:30:00"
```

---

## 6. PublicOntology API

```python
class PublicOntology:
    """Bionty-style public ontology backed by Parquet + DuckDB.

    Three-layer resolution: repo checkout > auto-updated cache > bundled.
    Works offline immediately via bundled Parquet files.

    Usage:
        ms = PublicOntology("psi_ms")
        ms.search("percolator")
        ms.lookup("MS:1001492")
        ms.lookup_normalized("percolator_score")
        ms.validate(["percolator_score", "unknown"])
        ms.standardize("percolator:score")
        ms.scores()
        ms.df()
        ms.version   # "4.1.235"
    """
```

| Method | Returns | Description |
|--------|---------|-------------|
| `__init__(source_name, auto_update=True)` | — | Load from best-available source |
| `search(query, field="name", top_k=10)` | `pd.DataFrame` | DuckDB contains/LIKE search |
| `lookup(accession)` | `dict or None` | Exact accession lookup |
| `lookup_name(name)` | `dict or None` | Case-insensitive name lookup |
| `lookup_normalized(name)` | `dict or None` | Exact normalized name lookup |
| `validate(names, field="normalized_name")` | `pd.Series[bool]` | Check which names exist |
| `standardize(raw_name)` | `str or None` | Normalize + lookup → accession |
| `scores()` | `pd.DataFrame` | is_score=True terms |
| `df(include_obsolete=False)` | `pd.DataFrame` | Full table |
| `version` | `str` | Loaded version from Parquet metadata |
| `check_update()` | `str or None` | Check repo for newer version |
| `update()` | `PublicOntology` | Force download latest from repo |
| `from_obo(path, source)` | `PublicOntology` | Build from local OBO file |
| `from_parquet(path)` | `PublicOntology` | Load from any Parquet path |

---

## 7. Converter Integration

### Wire up write_score_ontology() in all converters

Currently `write_score_ontology()` is defined in `BaseConverter` but never called. Each converter orchestrator will call it at the end of `convert()`:

- `QuantMSConverter.compose()` → writes ontology.parquet with PSM/feature/PG scores
- `DiaNNConverter.convert_features()` / `convert_pg()` → writes ontology.parquet
- `MaxQuantConverter.convert()` → writes ontology.parquet
- `FragPipeConverter.convert()` → writes ontology.parquet

### SDRF converter: write run/sample ontology entries

The SDRF converter writes plain strings to run.parquet (instrument, enzymes, dissociation_method). After writing, it uses `PublicOntology("psi_ms")` to resolve these values to CV accessions and writes the mappings to ontology.parquet:

```python
# In SdrfConverter, after writing run.parquet:
ms = PublicOntology("psi_ms", auto_update=False)
entries = []
for instrument_name in unique_instruments:
    term = ms.lookup_name(instrument_name)
    if term:
        entries.append({
            "field_name": instrument_name,
            "ontology_name": term["name"],
            "ontology_accession": term["accession"],
            "ontology_source": term["source"],
            "ontology_version": ms.version,
            "view": "run",
            "description": term["definition"],
        })
# Write entries to ontology.parquet
```

---

## 8. Entity Types Coverage

| QPX Entity          | Ontology    | Phase | Uses In Codebase |
|---------------------|-------------|-------|------------------|
| Score               | PSI-MS      | 1     | scores.py, all PSM/feature/PG adapters |
| Instrument          | PSI-MS      | 1     | run.yaml, sdrf.py |
| Enzyme              | PSI-MS      | 1     | run.yaml, sdrf.py |
| DissociationMethod  | PSI-MS      | 1     | run.yaml, sdrf.py |
| FileFormat          | PSI-MS      | 1     | mzML/mzIdentML references |
| PRIDE Term          | PRIDE CV    | 1     | PRIDE-specific metadata |
| Modification        | UNIMOD      | 2     | psm_adapter.py, modifications.md |
| Tissue              | UBERON      | 2     | sample.yaml, sdrf.py |
| CellType            | CL          | 2     | sample.yaml, sdrf.py |
| Disease             | MONDO       | 2     | sample.yaml, sdrf.py |
| Organism            | NCBITaxon   | 2     | sample.yaml, sdrf.py |
| CrossLinker         | XLMOD       | 2     | mzidentml psm_adapter |

---

## 9. Key Design Decisions

1. **Hybrid Bundle + Auto-Update** — three layers: repo checkout > cache > bundled fallback.
2. **Offline-first** — bundled Parquet files (~450 KB) always work, no network needed.
3. **Auto version check** — every 24h, fetches `versions.yaml` (~200 bytes). Blocking with 5s timeout.
4. **No release needed** — maintainer pushes new Parquet to repo, users get it on next check.
5. **`auto_update=False`** — CI/offline use, skip all version checking.
6. **DuckDB backend** — consistent with rest of QPX. Scales to Phase 2 ontologies.
7. **No new dependencies** — uses DuckDB, PyArrow, urllib (all already in QPX).
8. **LaminDB compatible** — standard Parquet files, `ln.Artifact()` ready.
9. **run.yaml simplified** — instrument/enzyme/dissociation become plain strings (like sample.yaml). Modifications keep struct.
10. **ontology_version in ontology.parquet** — each CV mapping records which ontology version was used.
11. **CV resolution in ontology.parquet** — accessions live in ontology.parquet, not in run/sample/psm files. Clean separation of raw data vs. semantic grounding.
12. **Wire up converters** — all converters call `write_score_ontology()`. SDRF converter writes instrument/enzyme/dissociation mappings to ontology.parquet.
