"""MaxQuant converter constants — tool identity, field mappings, and ProForma."""

from __future__ import annotations

import re
from functools import lru_cache
from typing import Optional

from qpx.converters.ptm import build_proforma

TOOL_NAME = "MaxQuant"
TOOL_VERSIONS = "2.x"

FIELD_MAPPINGS = {
    "feature": {
        "sequence": ["Sequence"],
        "modified_sequence": ["Modified sequence"],
        "charge": ["Charge"],
        "run_file_name": ["Raw file"],
        "is_decoy": ["Reverse"],
        "scan": ["MS/MS scan number"],
        "observed_mz": ["m/z"],
        "rt": ["Calibrated retention time"],
        "rt_start": ["Calibrated retention time start"],
        "rt_stop": ["Calibrated retention time finish"],
        "posterior_error_probability": ["PEP"],
        "ion_mobility": ["1/K0"],
        "pg_accessions": ["Leading proteins"],
        "anchor_protein": ["Leading razor protein"],
        "gg_names": ["Gene names"],
        "intensity": ["Intensity"],
        "andromeda_score": ["Score"],
        "andromeda_delta_score": ["Delta score"],
        "mass": ["Mass"],
    },
    "psm": {
        "sequence": ["Sequence"],
        "modified_sequence": ["Modified sequence"],
        "charge": ["Charge"],
        "run_file_name": ["Raw file"],
        "is_decoy": ["Reverse"],
        "scan": ["Scan number", "MS/MS scan number"],
        "observed_mz": ["m/z"],
        "rt": ["Retention time"],
        "posterior_error_probability": ["PEP"],
        "andromeda_score": ["Score"],
        "andromeda_delta_score": ["Delta score"],
    },
    "pg": {
        "pg_accessions": ["Protein IDs"],
        "pg_names": ["Protein names"],
        "gg_accessions": ["Gene names"],
        "anchor_protein": ["Majority protein IDs"],
        "global_qvalue": ["Q-value"],
        "is_decoy": ["Reverse"],
        "contaminant": ["Potential contaminant"],
        "sequence_coverage": ["Sequence coverage [%]"],
        "molecular_weight": ["Mol. weight [kDa]"],
        "peptide_count_total": ["Peptides"],
        "peptide_count_unique": ["Unique peptides"],
        "peptide_count_razor": ["Razor + unique peptides"],
        "andromeda_score": ["Score"],
        "intensity": ["Intensity"],
    },
}

# ---------------------------------------------------------------------------
# MaxQuant modification name -> UNIMOD accession mapping
# ---------------------------------------------------------------------------

MQ_MOD_TO_UNIMOD: dict[str, int] = {
    "oxidation (m)": 35,
    "acetyl (protein n-term)": 1,
    "carbamidomethyl (c)": 4,
    "phospho (sty)": 21,
    "phospho (s)": 21,
    "phospho (t)": 21,
    "phospho (y)": 21,
    "deamidation (nq)": 7,
    "deamidation (n)": 7,
    "deamidation (q)": 7,
    "gln->pyro-glu": 28,
    "methyl (kr)": 34,
    "methyl (k)": 34,
    "methyl (r)": 34,
    "dimethyl (kr)": 36,
    "dimethyl (k)": 36,
    "dimethyl (r)": 36,
    "trimethyl (k)": 37,
    "gg (k)": 121,
    "formyl (protein n-term)": 122,
    "carbamyl (protein n-term)": 5,
    "tmt6plex (k)": 737,
    "tmt6plex (n-term)": 737,
    "tmtpro (k)": 730,
    "tmtpro (n-term)": 730,
    "itraq8plex (k)": 385,
    "itraq8plex (n-term)": 385,
    "pyro-glu from e": 214,
}

MQ_SHORT_TO_UNIMOD: dict[str, int] = {
    "ox": 35,
    "ac": 1,
    "ph": 21,
    "de": 7,
    "gl": 28,
    "me": 34,
}


def _resolve_mq_mod(tag: str) -> str:
    """Resolve a MaxQuant modification name to ``UNIMOD:N`` or return as-is.

    Handles full names (``Oxidation (M)``), short forms (``ac``),
    and already-resolved ``UNIMOD:`` prefixes.
    """
    m = re.match(r"(?i)unimod:(\d+)", tag)
    if m:
        return f"UNIMOD:{m.group(1)}"

    uid = MQ_MOD_TO_UNIMOD.get(tag.lower())
    if uid is None:
        uid = MQ_SHORT_TO_UNIMOD.get(tag.lower())
    if uid is not None:
        return f"UNIMOD:{uid}"

    return tag


@lru_cache(maxsize=8192)
def to_proforma(modified_sequence: Optional[str]) -> str:
    """Convert a MaxQuant Modified sequence to ProForma notation.

    MaxQuant wraps sequences in underscores and uses parenthetical
    modification names (possibly nested)::

        _PEPTM(Oxidation (M))IDEK_         -> PEPTM[UNIMOD:35]IDEK
        _(Acetyl (Protein N-term))PEPTIDEK_ -> [UNIMOD:1]-PEPTIDEK
        _(ac)PEPTIDEK_                      -> [UNIMOD:1]-PEPTIDEK
    """
    if not isinstance(modified_sequence, str) or not modified_sequence:
        return ""

    seq = modified_sequence.strip("_")
    if not seq:
        return ""

    mods: list[tuple[int, str]] = []
    plain_chars: list[str] = []
    i = 0
    n = len(seq)
    seq_pos = 0

    while i < n:
        if seq[i] == "(":
            # Walk forward tracking parenthesis depth to handle nested parens
            depth = 1
            j = i + 1
            while j < n and depth > 0:
                if seq[j] == "(":
                    depth += 1
                elif seq[j] == ")":
                    depth -= 1
                j += 1
            mod_name = seq[i + 1 : j - 1]
            tag = _resolve_mq_mod(mod_name)
            position = 0 if seq_pos == 0 else seq_pos
            mods.append((position, tag))
            i = j
        else:
            plain_chars.append(seq[i])
            seq_pos += 1
            i += 1

    plain_seq = "".join(plain_chars)
    return build_proforma(plain_seq, mods)


def parse_phospho_probabilities(
    prob_string: str,
) -> dict[int, list[dict]] | None:
    """Parse MaxQuant ``Phospho (STY) Probabilities`` into per-position scores.

    MaxQuant encodes per-position phosphorylation probabilities using the same
    underscore-wrapped format as Modified sequence, but with floating-point
    probabilities in parentheses instead of modification names::

        _PEPTIDES(0.95)T(0.03)Y(0.02)K_

    Only positions with non-zero probability are included in the output.

    Args:
        prob_string: The ``Phospho (STY) Probabilities`` column value.

    Returns:
        Dict mapping 1-indexed position -> list of score dicts, or ``None``
        if no probabilities found.
    """
    if not isinstance(prob_string, str) or not prob_string:
        return None

    seq = prob_string.strip("_")
    if not seq:
        return None

    site_scores: dict[int, list[dict]] = {}
    seq_pos = 0
    i = 0
    n = len(seq)

    while i < n:
        if seq[i] == "(":
            end = seq.index(")", i)
            prob_str = seq[i + 1 : end]
            try:
                prob = float(prob_str)
            except ValueError:
                i = end + 1
                continue
            if prob > 0:
                site_scores[seq_pos] = [
                    {
                        "score_name": "phospho_sty_probability",
                        "score_value": prob,
                        "higher_better": True,
                    }
                ]
            i = end + 1
        else:
            seq_pos += 1
            i += 1

    return site_scores if site_scores else None
