"""consensusXML -> QPX pg (protein group) records — identification only.

Protein groups come from the OpenMS protein-inference graph: indistinguishable
protein groups when present, else each protein hit as a singleton group. Peptide
counts come from the peptide→protein evidence links.

**No protein intensity is emitted** (interim): the consensusXML has no
protein-level abundance. Rows carry ``label``/``intensity`` null (identification-
only) until OpenMS ``-out_qpx`` provides the authoritative protein quant. One row
per ``(protein group, grouped_runs)`` — the quantification unit(s), from the SDRF
when given, otherwise the whole run set as a single unit.
"""

from __future__ import annotations

from collections import defaultdict
from typing import Optional

from qpx.converters.channel_labels import fraction_groups_from_sdrf
from qpx.converters.openms_consensus.feature_adapter import _run_stem


def _protein_to_peptides(cm) -> dict[str, set[str]]:
    """Map each protein accession to the set of peptide sequences supporting it."""
    acc_to_pep: dict[str, set[str]] = defaultdict(set)

    def _collect(pid):
        for hit in pid.getHits():
            seq = hit.getSequence().toUnmodifiedString()
            for ev in hit.getPeptideEvidences():
                acc = ev.getProteinAccession()
                if acc:
                    acc_to_pep[acc].add(seq)

    for pid in cm.getUnassignedPeptideIdentifications():
        _collect(pid)
    for cf in cm:
        for pid in cf.getPeptideIdentifications():
            _collect(pid)
    return acc_to_pep


def _is_decoy_accession(acc: str) -> bool:
    return str(acc).upper().startswith(("DECOY", "REV_", "RANDOM_"))


def consensus_protein_groups_to_records(
    consensusxml_path: str,
    sdrf_path: Optional[str] = None,
) -> list[dict]:
    """Return identification-only QPX pg record dicts (null intensity)."""
    import pyopenms as oms

    cm = oms.ConsensusMap()
    oms.ConsensusXMLFile().load(str(consensusxml_path), cm)

    all_runs = sorted({_run_stem(cm.getColumnHeaders()[i].filename) for i in cm.getColumnHeaders()})
    # grouped_runs unit per run, from the SDRF (fractions grouped); else one unit.
    run_to_grouped = fraction_groups_from_sdrf(sdrf_path)
    if run_to_grouped:
        units = {tuple(v) for v in run_to_grouped.values()}
    else:
        units = {tuple(all_runs)}

    acc_to_pep = _protein_to_peptides(cm)

    # Build groups: indistinguishable protein groups if present, else singletons.
    groups: list[list[str]] = []
    prot_ids = cm.getProteinIdentifications()
    if prot_ids and prot_ids[0].getIndistinguishableProteins():
        for grp in prot_ids[0].getIndistinguishableProteins():
            accs = [a.decode() if isinstance(a, bytes) else str(a) for a in grp.accessions]
            if accs:
                groups.append(accs)
    elif prot_ids:
        for hit in prot_ids[0].getHits():
            acc = hit.getAccession()
            groups.append([acc.decode() if isinstance(acc, bytes) else str(acc)])

    records: list[dict] = []
    for accs in groups:
        anchor = accs[0]
        peptide_seqs: set[str] = set()
        for acc in accs:
            peptide_seqs |= acc_to_pep.get(acc, set())
        n_pep = len(peptide_seqs)
        is_decoy = all(_is_decoy_accession(a) for a in accs)
        for unit in units:
            records.append(
                {
                    "pg_accessions": list(accs),
                    "anchor_protein": anchor,
                    "grouped_runs": list(unit),
                    "label": None,  # interim: identification-only, no protein quant
                    "intensity": None,
                    "is_decoy": is_decoy,
                    "peptides": [{"protein_name": a, "peptide_count": n_pep} for a in accs],
                }
            )
    return records
