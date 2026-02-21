"""FragPipe converter constants — tool identity, field mappings, and PTM parsing."""

from __future__ import annotations

import re
from typing import Optional

from qpx.converters.ptm import build_proforma, from_proforma, mass_to_unimod

TOOL_NAME = "FragPipe"
TOOL_VERSIONS = "20+"

FIELD_MAPPINGS = {
    "feature": {
        "sequence": ["Peptide Sequence"],
        "modified_sequence": ["Modified Sequence", "Modified Peptide"],
        "pg_accessions": ["Protein", "Protein ID"],
        "gg_names": ["Gene"],
        "observed_mz": ["M/Z"],
        "charge": ["Charge"],
        "charges": ["Charges"],
        "modifications": ["Assigned Modifications"],
    },
    "pg": {
        "pg_accessions": ["Protein", "Protein ID"],
        "gg_accessions": ["Gene", "Gene Names"],
        "pg_names": ["Description", "Protein Description"],
        "peptide_count_total": ["Combined Total Peptides", "Total Peptides"],
        "peptide_count_unique": ["Combined Unique Peptides", "Unique Peptides"],
        "spectral_count": ["Combined Spectral Count", "Spectral Count"],
        "sequence_coverage": ["Percent Coverage", "Coverage"],
        "molecular_weight": ["Protein Molecular Weight (Da)"],
    },
    "psm": {
        "sequence": ["Peptide"],
        "modified_sequence": ["Modified Peptide"],
        "charge": ["Charge"],
        "observed_mz": ["Observed M/Z"],
        "calculated_mz": ["Calculated M/Z"],
        "rt": ["Retention"],
        "pg_accessions": ["Protein"],
    },
}

# ---------------------------------------------------------------------------
# PTM parsing: FragPipe "Assigned Modifications" -> ProForma
# ---------------------------------------------------------------------------

_FP_MOD_TOKEN_RE = re.compile(r"(?:(\d+)([A-Z])|(N-term))\(([^)]+)\)")


def to_proforma(assigned_mods: str, sequence: str) -> str:
    """Build ProForma from a FragPipe sequence and Assigned Modifications string.

    Args:
        assigned_mods: Assigned Modifications string, e.g. ``5M(15.9949)`` or
            ``N-term(42.0106), 5M(15.9949)``.
        sequence: Stripped peptide sequence.

    Returns:
        ProForma string, e.g. ``PEPTM[UNIMOD:35]IDEK``.
    """
    if not sequence:
        return ""
    if not assigned_mods or not assigned_mods.strip():
        return sequence

    pos_mods: dict[int, str] = {}
    nterm_tag: Optional[str] = None

    for token in assigned_mods.split(","):
        token = token.strip()
        if not token:
            continue

        m = _FP_MOD_TOKEN_RE.match(token)
        if not m:
            continue

        pos_str, _aa, nterm, mass_str = m.groups()
        try:
            mass = float(mass_str)
        except ValueError:
            continue

        unimod_id = mass_to_unimod(mass)
        if unimod_id is not None:
            tag = f"UNIMOD:{unimod_id}"
        else:
            tag = f"+{mass}" if mass >= 0 else str(mass)

        if nterm:
            nterm_tag = tag
        elif pos_str is not None:
            pos_mods[int(pos_str)] = tag

    mods: list[tuple[int, str]] = []
    if nterm_tag:
        mods.append((0, nterm_tag))
    for pos, tag in pos_mods.items():
        mods.append((pos, tag))

    return build_proforma(sequence, mods)


def to_modifications(assigned_mods: str, sequence: str) -> list[dict] | None:
    """Parse modifications from FragPipe Assigned Modifications format.

    Args:
        assigned_mods: Assigned Modifications string, e.g. ``5M(15.9949)``.
        sequence: Stripped peptide sequence.

    Returns:
        List of modification dicts per QPX schema, or ``None`` if unmodified.
    """
    proforma = to_proforma(assigned_mods, sequence)
    return from_proforma(proforma, sequence, meta=None)
