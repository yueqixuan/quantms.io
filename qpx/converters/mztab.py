"""Shared mzTab loading logic for DuckDB-based converters.

This module provides helpers that parse an mzTab file (metadata, PSM,
protein, and MSstats sections) and load them into DuckDB tables.  The
quantms adapters compose these helpers rather than inheriting from the
legacy ``MzTabIndexer``.
"""

from __future__ import annotations

import gzip
import logging
import re
from pathlib import Path

import duckdb
import pandas as pd

logger = logging.getLogger(__name__)

# Prefixes that identify mzTab line types
_METADATA_PREFIX = "MTD"
_PROTEIN_HEADER_PREFIX = "PRH"
_PROTEIN_LINE_PREFIX = "PRT"
_PSM_HEADER_PREFIX = "PSH"
_PSM_LINE_PREFIX = "PSM"


def _open_mztab(path: str):
    """Open an mzTab file (plain or gzipped)."""
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return open(path, "r", encoding="utf-8")


# ------------------------------------------------------------------
# Public helpers
# ------------------------------------------------------------------


def load_mztab_sections(
    conn: duckdb.DuckDBPyConnection,
    mztab_path: str,
) -> None:
    """Parse an mzTab file and load metadata, proteins, and PSMs into DuckDB.

    After calling this function the connection will contain:
        * ``metadata``  -- two-column table (key TEXT, value TEXT)
        * ``proteins``  -- protein section with dynamic columns
        * ``psms``      -- PSM section with dynamic columns

    Args:
        conn: An open DuckDB connection (in-memory or persistent).
        mztab_path: Path to the mzTab file (plain-text or ``.gz``).
    """
    metadata_rows: list[tuple[str, str]] = []
    protein_header: list[str] | None = None
    protein_rows: list[list[str]] = []
    psm_header: list[str] | None = None
    psm_rows: list[list[str]] = []

    with _open_mztab(mztab_path) as fh:
        for line in fh:
            line = line.rstrip("\n\r")
            if not line:
                continue

            parts = line.split("\t")
            prefix = parts[0] if parts else ""

            if prefix == _METADATA_PREFIX and len(parts) >= 3:
                metadata_rows.append((parts[1], parts[2]))

            elif prefix == _PROTEIN_HEADER_PREFIX:
                protein_header = [_clean_col(c) for c in parts[1:]]

            elif prefix == _PROTEIN_LINE_PREFIX and protein_header is not None:
                protein_rows.append(parts[1:])

            elif prefix == _PSM_HEADER_PREFIX:
                psm_header = [_clean_col(c) for c in parts[1:]]

            elif prefix == _PSM_LINE_PREFIX and psm_header is not None:
                psm_rows.append(parts[1:])

    # Load metadata
    meta_df = pd.DataFrame(metadata_rows, columns=["key", "value"])
    conn.execute("DROP TABLE IF EXISTS metadata")
    conn.from_df(meta_df).create("metadata")

    # Load proteins
    if protein_header and protein_rows:
        prot_df = pd.DataFrame(protein_rows, columns=protein_header)
        conn.execute("DROP TABLE IF EXISTS proteins")
        conn.from_df(prot_df).create("proteins")
    else:
        conn.execute("CREATE TABLE IF NOT EXISTS proteins (accession VARCHAR)")

    # Load PSMs
    if psm_header and psm_rows:
        psm_df = pd.DataFrame(psm_rows, columns=psm_header)
        conn.execute("DROP TABLE IF EXISTS psms")
        conn.from_df(psm_df).create("psms")
    else:
        conn.execute("CREATE TABLE IF NOT EXISTS psms (sequence VARCHAR)")

    logger.info(f"mzTab loaded: {len(metadata_rows)} metadata rows, {len(protein_rows)} proteins, {len(psm_rows)} PSMs")


def load_msstats(
    conn: duckdb.DuckDBPyConnection,
    msstats_path: str,
) -> None:
    """Load an MSstats input CSV/TSV into DuckDB as the ``msstats`` table.

    Args:
        conn: An open DuckDB connection that already has ``metadata``.
        msstats_path: Path to the MSstats input file.
    """
    conn.execute(f"""
        CREATE OR REPLACE TABLE msstats AS
        SELECT * FROM read_csv_auto('{msstats_path}',
            header=true, auto_detect=true,
            quote='"', delim=',', null_padding=true,
            max_line_size=10000000)
        """)
    count = conn.execute("SELECT COUNT(*) FROM msstats").fetchone()[0]
    logger.info(f"Loaded {count:,} MSstats rows from {msstats_path}")


def extract_ms_runs(conn: duckdb.DuckDBPyConnection) -> dict[int, str]:
    """Return a mapping of MS run index to file stem from metadata.

    Example: ``{1: 'sample_01', 2: 'sample_02'}``
    """
    rows = conn.execute("""
        SELECT key, value FROM metadata
        WHERE key LIKE 'ms_run[%]-location'
        """).fetchall()

    ms_runs: dict[int, str] = {}
    for key, value in rows:
        m = re.search(r"\[(\d+)\]", key)
        if m:
            idx = int(m.group(1))
            ms_runs[idx] = Path(value).stem
    return ms_runs


def extract_modifications(conn: duckdb.DuckDBPyConnection) -> dict:
    """Parse modification metadata from the mzTab metadata table.

    Returns a dict keyed by modification accession, with values
    ``(name, list[site], list[position])``.
    """
    rows = conn.execute("""
        SELECT key, value FROM metadata
        WHERE key LIKE '%_mod[%'
           OR key LIKE '%variable_mod[%'
        ORDER BY key
        """).fetchall()

    temp: dict[str, dict] = {}
    for key, value in rows:
        base_key = key.split("-")[0]
        if base_key not in temp:
            temp[base_key] = {
                "name": None,
                "accession": None,
                "site": None,
                "position": None,
            }

        if "-site" in key:
            temp[base_key]["site"] = value.strip()
        elif "-position" in key:
            temp[base_key]["position"] = value.strip()
        else:
            # Main modification line -- parse [CV, accession, name, ...]
            vals = value.replace("[", "").replace("]", "").split(",")
            if len(vals) >= 3:
                temp[base_key]["accession"] = vals[1].strip()
                temp[base_key]["name"] = vals[2].strip()

    # Group by accession
    modifications: dict = {}
    for entry in temp.values():
        acc = entry["accession"]
        if acc and acc not in modifications:
            modifications[acc] = (
                entry["name"],
                [entry["site"]],
                [entry["position"]],
            )
        elif acc:
            modifications[acc][1].append(entry["site"])
            modifications[acc][2].append(entry["position"])

    return modifications


def extract_score_names(
    conn: duckdb.DuckDBPyConnection,
) -> dict[str, dict[int, str]]:
    """Parse search engine score metadata.

    Returns::

        {"psms": {1: "percolator:score", ...},
         "peptides": {...},
         "proteins": {...}}
    """
    rows = conn.execute("""
        SELECT key, value FROM metadata
        WHERE key LIKE '%search_engine_score%'
        """).fetchall()

    result: dict[str, dict[int, str]] = {
        "psms": {},
        "peptides": {},
        "proteins": {},
    }
    for key, value in rows:
        m = re.search(r"\[(\d+)\]", key)
        if not m:
            continue
        idx = int(m.group(1))
        # Extract the name term (third comma-separated value)
        parts = value.replace("[", "").replace("]", "").split(",")
        name = parts[2].strip() if len(parts) >= 3 else value.strip()

        if key.startswith("psm_"):
            result["psms"][idx] = name
        elif key.startswith("peptide_"):
            result["peptides"][idx] = name
        elif key.startswith("protein_"):
            result["proteins"][idx] = name

    return result


# ------------------------------------------------------------------
# Private helpers
# ------------------------------------------------------------------


def _clean_col(name: str) -> str:
    """Normalise a column name from the mzTab header."""
    return name.strip().lower().replace(" ", "_")
