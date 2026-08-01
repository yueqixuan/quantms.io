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


def _protein_maps(cm) -> tuple[dict[str, set[str]], dict[str, set[str]]]:
    """Return (accession -> peptide sequences, accession -> runs identified in).

    The per-protein run set lets pg rows be emitted only for the quantification
    units where the protein was actually identified, rather than every unit.
    """
    headers = cm.getColumnHeaders()
    map_run = {i: _run_stem(headers[i].filename) for i in headers}
    acc_to_pep: dict[str, set[str]] = defaultdict(set)
    acc_to_runs: dict[str, set[str]] = defaultdict(set)

    def _accession(ev) -> str | None:
        acc = ev.getProteinAccession()
        if not acc:
            return None
        return acc.decode() if isinstance(acc, bytes) else str(acc)

    def _collect(pid, runs):
        for hit in pid.getHits():
            seq = hit.getSequence().toUnmodifiedString()
            for ev in hit.getPeptideEvidences():
                acc = _accession(ev)
                if acc:
                    acc_to_pep[acc].add(seq)
                    acc_to_runs[acc] |= runs

    # Assigned IDs: the runs are the consensus feature's member maps.
    for cf in cm:
        cf_runs = {map_run.get(sub.getMapIndex()) for sub in cf.getFeatureList() if sub.getIntensity() > 0}
        cf_runs.discard(None)
        for pid in cf.getPeptideIdentifications():
            _collect(pid, cf_runs)
    # Unassigned IDs: the run comes from map_index / id_merge_index.
    for pid in cm.getUnassignedPeptideIdentifications():
        run = None
        for key in ("map_index", "id_merge_index"):
            if pid.metaValueExists(key):
                run = map_run.get(int(pid.getMetaValue(key)))
                break
        _collect(pid, {run} if run else set())
    return acc_to_pep, acc_to_runs


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

    acc_to_pep, acc_to_runs = _protein_maps(cm)

    # Build groups: indistinguishable protein groups if present, else singletons.
    # Also capture per-accession decoy flag (target_decoy meta) and q-value (the
    # protein hit score when the identification score IS a q-value).
    groups: list[list[str]] = []
    acc_decoy: dict[str, bool] = {}
    acc_qvalue: dict[str, float] = {}
    prot_ids = cm.getProteinIdentifications()
    if prot_ids:
        prot = prot_ids[0]
        score_is_qvalue = str(prot.getScoreType() or "").lower() in ("q-value", "qvalue", "fdr")
        for hit in prot.getHits():
            acc = hit.getAccession()
            acc = acc.decode() if isinstance(acc, bytes) else str(acc)
            if hit.metaValueExists("target_decoy"):
                acc_decoy[acc] = "decoy" in str(hit.getMetaValue("target_decoy")).lower()
            if score_is_qvalue and hit.getScore() is not None:
                acc_qvalue[acc] = float(hit.getScore())
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
        # Prefer the target_decoy meta; fall back to the accession prefix.
        is_decoy = all(acc_decoy.get(a, _is_decoy_accession(a)) for a in accs)
        qvals = [acc_qvalue[a] for a in accs if a in acc_qvalue]
        global_qvalue = min(qvals) if qvals else None
        # Only the quantification units where this group was actually identified
        # (its peptides appear in a run of that unit) — not every unit.
        group_runs: set[str] = set()
        for acc in accs:
            group_runs |= acc_to_runs.get(acc, set())
        group_units = [unit for unit in units if group_runs.intersection(unit)] or list(units)
        for unit in group_units:
            records.append(
                {
                    "pg_accessions": list(accs),
                    "anchor_protein": anchor,
                    "grouped_runs": list(unit),
                    "label": None,  # interim: identification-only, no protein quant
                    "intensity": None,
                    "global_qvalue": global_qvalue,
                    "is_decoy": is_decoy,
                    "peptides": [{"protein_name": a, "peptide_count": n_pep} for a in accs],
                }
            )
    return records
