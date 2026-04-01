"""DIA-NN converter constants — tool identity, field mappings, and PTM parsing."""

from __future__ import annotations

import re
from functools import lru_cache

from qpx.converters.ptm import build_proforma, from_proforma

TOOL_NAME = "DIA-NN"
TOOL_VERSIONS = "1.8+"

# QPX field -> ordered list of candidate DIA-NN column names.
# At runtime, the first match against actual input columns wins.
FIELD_MAPPINGS = {
    "feature": {
        "intensity": ["Precursor.Quantity"],
        "posterior_error_probability": ["PEP"],
        "rt": ["RT"],
        "rt_start": ["RT.Start"],
        "rt_stop": ["RT.Stop"],
        "predicted_rt": ["Predicted.RT"],
        "pg_accessions": ["Protein.Group"],
        "observed_mz": ["Precursor.Mz"],
        "lfq": ["Precursor.Normalised"],
        "charge": ["Precursor.Charge"],
        "sequence": ["Stripped.Sequence"],
        "modified_sequence": ["Modified.Sequence"],
        "gg_names": ["Genes"],
        "run_file_name": ["Run"],
        "qvalue": ["Q.Value"],
        "pg_qvalue": ["PG.Q.Value"],
        "global_qvalue": ["Global.Q.Value"],
        "pg_global_qvalue": ["Global.PG.Q.Value"],
        "mp_accessions": ["Protein.Ids"],
        "normalize_intensity": ["Precursor.Normalised"],
        "ms1_area": ["Ms1.Area"],
        "lfq_maxlfq": ["PG.MaxLFQ"],
        "precursor_quantification_score": ["Quantity.Quality"],
        "ms2_scan": ["MS2.Scan"],
        "ion_mobility": ["IM"],
        "mass_error_ppm": ["Mass.Error (ppm)"],
        "ion_mobility_start": ["IM.Start"],
        "ion_mobility_stop": ["IM.Stop"],
        # Scores
        "cscore": ["CScore"],
        "evidence": ["Evidence"],
        "spectrum_similarity": ["Spectrum.Similarity"],
        "ms1_profile_corr": ["Ms1.Profile.Corr"],
        "mass_evidence": ["Mass.Evidence"],
        "averagine": ["Averagine"],
        "quantity_quality": ["Quantity.Quality"],
        "lib_qvalue": ["Lib.Q.Value"],
        "lib_pg_qvalue": ["Lib.PG.Q.Value"],
        "protein_qvalue": ["Protein.Q.Value"],
        "decoy_cscore": ["Decoy.CScore"],
        "decoy_evidence": ["Decoy.Evidence"],
        "translated_quality": ["Translated.Quality"],
        "translated_qvalue": ["Translated.Q.Value"],
        # CV params
        "proteotypic": ["Proteotypic"],
        "irt": ["iRT"],
        "predicted_irt": ["Predicted.iRT"],
        "predicted_im": ["Predicted.IM"],
        "predicted_iim": ["Predicted.iIM"],
        "iim": ["iIM"],
        # Decoy flag (explicit column; fallback to prefix heuristic if absent)
        "decoy": ["Decoy"],
    },
    "pg": {
        "intensity": ["Precursor.Quantity", "PG.Quantity"],
        "pg_accessions": ["Protein.Group"],
        "pg_names": ["Protein.Names"],
        "gg_accessions": ["Genes"],
        "global_qvalue": ["Global.PG.Q.Value"],
        "gg_qvalue": ["GG.Q.Value"],
        "lfq": ["PG.MaxLFQ"],
        "qvalue": ["PG.Q.Value"],
        "run_file_name": ["Run"],
    },
}

# ---------------------------------------------------------------------------
# PTM parsing: DIA-NN Modified.Sequence -> ProForma
# ---------------------------------------------------------------------------


@lru_cache(maxsize=8192)
def to_proforma(modified_sequence: str) -> str:
    """Convert a DIA-NN Modified.Sequence to ProForma notation.

    DIA-NN encodes modifications with parentheses::

        PEPTM(UniMod:35)IDE       -> PEPTM[UNIMOD:35]IDE
        (UniMod:1)PEPTIDE         -> [UNIMOD:1]-PEPTIDE

    Args:
        modified_sequence: The ``Modified.Sequence`` value from DIA-NN output.

    Returns:
        ProForma string, e.g. ``PEPTM[UNIMOD:35]IDEK``.
    """
    if not modified_sequence:
        return ""

    mods: list[tuple[int, str]] = []
    plain_chars: list[str] = []
    i = 0
    n = len(modified_sequence)
    seq_pos = 0

    while i < n:
        if modified_sequence[i] == "(":
            end = modified_sequence.index(")", i)
            inner = modified_sequence[i + 1 : end]
            normalised = re.sub(r"(?i)unimod:", "UNIMOD:", inner)
            position = 0 if seq_pos == 0 else seq_pos
            mods.append((position, normalised))
            i = end + 1
        else:
            plain_chars.append(modified_sequence[i])
            seq_pos += 1
            i += 1

    plain_seq = "".join(plain_chars)
    return build_proforma(plain_seq, mods)


def to_modifications(modified_sequence: str, sequence: str) -> list[dict] | None:
    """Parse modifications from a DIA-NN Modified.Sequence string.

    Converts to ProForma first, then delegates to the shared ``from_proforma``
    parser to produce the QPX modification list.

    Args:
        modified_sequence: The ``Modified.Sequence`` value from DIA-NN output.
        sequence: Stripped peptide sequence (no modification annotations).

    Returns:
        List of modification dicts per QPX schema, or ``None`` if unmodified.
    """
    proforma = to_proforma(modified_sequence)
    return from_proforma(proforma, sequence, meta=None)
