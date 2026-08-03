"""consensusXML -> QPX pg (protein group) records.

Protein groups come from the OpenMS protein-inference graph: indistinguishable
protein groups when present, else each protein hit as a singleton group. Peptide
counts come from the peptide→protein evidence links.

**Interim protein intensity**: the consensusXML has no protein-level abundance
(that lived in the mzTab), so until OpenMS ``-out_qpx`` provides the authoritative
number we roll one up ourselves — the **unnormalized sum of the group's unique
peptides** (quantms ``unique_peptides`` policy, no normalization), per
``(protein group, grouped_runs unit, label)``. It is *not* the authoritative
quant, so every quantified row carries a ``quantification_method`` cv_param; set
``top=3`` to mirror the ProteomicsLFQ/IsobaricWorkflow default instead of summing
all peptides. Rows stay ``intensity``-null where a group has no unique-peptide
signal. One row per ``(protein group, grouped_runs, label)`` — the quantification
unit(s) from the SDRF when given, otherwise the whole run set as a single unit.
"""

from __future__ import annotations

import re
from collections import defaultdict
from typing import Optional

from qpx.converters.channel_labels import fraction_groups_from_sdrf
from qpx.converters.openms_consensus.feature_adapter import (
    _canonical_channel,
    _run_stem,
    load_consensus_map,
    to_proforma,
)

_GENE_RE = re.compile(r"GN=([^\s]+)")


def _protein_maps(cm) -> tuple[dict[str, set[str]], dict[str, set[str]], dict[str, set[tuple]]]:
    """Return (accession -> peptide sequences, -> runs identified in, -> features).

    The per-protein run set lets pg rows be emitted only for the quantification
    units where the protein was actually identified; the feature set (distinct
    peptidoform+charge) feeds ``feature_counts``.
    """
    headers = cm.getColumnHeaders()
    map_run = {i: _run_stem(headers[i].filename) for i in headers}
    acc_to_pep: dict[str, set[str]] = defaultdict(set)
    acc_to_runs: dict[str, set[str]] = defaultdict(set)
    acc_to_feat: dict[str, set[tuple]] = defaultdict(set)

    def _accession(ev) -> str | None:
        acc = ev.getProteinAccession()
        if not acc:
            return None
        return acc.decode() if isinstance(acc, bytes) else str(acc)

    def _collect(pid, runs):
        for hit in pid.getHits():
            seq_obj = hit.getSequence()
            seq = seq_obj.toUnmodifiedString()
            feat = (to_proforma(seq_obj), int(hit.getCharge() or 0))
            for ev in hit.getPeptideEvidences():
                acc = _accession(ev)
                if acc:
                    acc_to_pep[acc].add(seq)
                    acc_to_feat[acc].add(feat)
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
    return acc_to_pep, acc_to_runs, acc_to_feat


def _is_decoy_accession(acc: str) -> bool:
    return str(acc).upper().startswith(("DECOY", "REV_", "RANDOM_"))


def _is_contaminant(acc: str) -> bool:
    return "CONTAM" in str(acc).upper()


def _acc_str(acc) -> str:
    return acc.decode() if isinstance(acc, bytes) else str(acc)


def _protein_hit_meta(prot) -> tuple[dict[str, bool], dict[str, float], dict[str, str]]:
    """Per-accession decoy flag, q-value (when the protein score IS a q-value), gene."""
    score_is_qvalue = str(prot.getScoreType() or "").lower() in ("q-value", "qvalue", "fdr")
    acc_decoy: dict[str, bool] = {}
    acc_qvalue: dict[str, float] = {}
    acc_gene: dict[str, str] = {}
    for hit in prot.getHits():
        acc = _acc_str(hit.getAccession())
        if hit.metaValueExists("target_decoy"):
            acc_decoy[acc] = "decoy" in str(hit.getMetaValue("target_decoy")).lower()
        if score_is_qvalue and hit.getScore() is not None:
            acc_qvalue[acc] = float(hit.getScore())
        gene = _GENE_RE.search(str(hit.getDescription() or "")) if hasattr(hit, "getDescription") else None
        if gene:
            acc_gene[acc] = gene.group(1)
    return acc_decoy, acc_qvalue, acc_gene


def _build_groups(prot) -> list[list[str]]:
    """Protein groups: indistinguishable groups plus a singleton for every hit
    not covered by one.

    ``getIndistinguishableProteins()`` only returns the *grouped* proteins, so a
    hit that OpenMS left ungrouped would be lost if we returned those groups
    alone. Append singletons for the uncovered hits (this also reproduces the
    all-singletons behaviour when there are no indistinguishable groups).
    """
    groups: list[list[str]] = []
    covered: set[str] = set()
    for grp in prot.getIndistinguishableProteins():
        accs = [_acc_str(a) for a in grp.accessions]
        if accs:
            groups.append(accs)
            covered.update(accs)
    for hit in prot.getHits():
        acc = _acc_str(hit.getAccession())
        if acc not in covered:
            groups.append([acc])
            covered.add(acc)
    return groups


def _merge_protein_ids(cm) -> tuple[dict[str, bool], dict[str, float], dict[str, str], list[list[str]]]:
    """Merge every ProteinIdentification's metadata + groups from the inference graph.

    A merged multi-run consensusXML carries one ProteinIdentification per run, so
    process them all: merge per-accession decoy/gene, keep the best (min) q-value,
    and concatenate groups while deduplicating equivalent accession sets.
    """
    acc_decoy: dict[str, bool] = {}
    acc_qvalue: dict[str, float] = {}
    acc_gene: dict[str, str] = {}
    groups: list[list[str]] = []
    seen_groups: set[frozenset] = set()
    for prot in cm.getProteinIdentifications():
        decoy, qvalue, gene = _protein_hit_meta(prot)
        acc_decoy.update(decoy)
        acc_gene.update(gene)
        for acc, qv in qvalue.items():
            acc_qvalue[acc] = min(acc_qvalue[acc], qv) if acc in acc_qvalue else qv
        for grp in _build_groups(prot):
            key = frozenset(grp)
            if key not in seen_groups:
                seen_groups.add(key)
                groups.append(grp)
    return acc_decoy, acc_qvalue, acc_gene, groups


def _map_info(cm) -> dict[int, tuple[str, str]]:
    """Map index -> (run_file_name, canonical label), matching the feature adapter.

    ``label-free`` (pyopenms spelling, hyphen) collapses every map to ``LFQ``;
    labeled maps use the canonical isobaric channel label.
    """
    headers = cm.getColumnHeaders()
    is_labeled = cm.getExperimentType() != "label-free"
    return {i: (_run_stem(headers[i].filename), _canonical_channel(headers[i].label) if is_labeled else "LFQ") for i in headers}


def _peptide_intensities(cm, map_info: dict[int, tuple[str, str]]):
    """Return ``((peptide, run, label) -> summed intensity, peptide -> {accessions})``.

    Peptides are keyed by unmodified sequence (matching :func:`_protein_maps`);
    the accession set drives the unique-to-group test. Feature intensities are
    summed per (peptide, run, label) — charge states / peptidoforms of one
    sequence roll up together.
    """
    pep_intensity: dict[tuple[str, str, str], float] = defaultdict(float)
    pep_accs: dict[str, set[str]] = defaultdict(set)
    for cf in cm:
        pids = cf.getPeptideIdentifications()
        if not pids or not pids[0].getHits():
            continue
        hit = pids[0].getHits()[0]
        seq = hit.getSequence().toUnmodifiedString()
        for ev in hit.getPeptideEvidences():
            acc = ev.getProteinAccession()
            if acc:
                pep_accs[seq].add(acc.decode() if isinstance(acc, bytes) else str(acc))
        for sub in cf.getFeatureList():
            inten = float(sub.getIntensity())
            if inten <= 0:
                continue
            run, label = map_info.get(sub.getMapIndex(), (None, None))
            if run is not None:
                pep_intensity[(seq, run, label)] += inten
    return pep_intensity, pep_accs


def _protein_intensity(group_peps, group_accs, unit, label, pep_intensity, pep_accs, top) -> Optional[float]:
    """Interim protein-group intensity for one ``(unit, label)``.

    Sum of the group's **unique** peptides — those mapping only to proteins in the
    group (the quantms ``unique_peptides`` policy) — with **no normalization**.
    ``top > 0`` keeps only the N most-abundant peptides (quantms/ProteomicsLFQ
    default is 3); ``top = 0`` sums all. Returns ``None`` when the group has no
    unique-peptide signal in this unit+label (stays identification-only).
    """
    abundances = []
    for pep in group_peps:
        if not pep_accs.get(pep, set()).issubset(group_accs):
            continue  # shared outside the group -> excluded by unique_peptides
        ab = sum(pep_intensity.get((pep, run, label), 0.0) for run in unit)
        if ab > 0:
            abundances.append(ab)
    if not abundances:
        return None
    abundances.sort(reverse=True)
    if top and top > 0:
        abundances = abundances[:top]
    return sum(abundances)


def consensus_protein_groups_to_records(
    consensusxml_path: str | None = None,
    sdrf_path: Optional[str] = None,
    cm=None,
    top: int = 0,
) -> list[dict]:
    """Return QPX pg record dicts: one per (protein group, unit, label).

    ``intensity`` is an **interim, unnormalized total** — the sum of the group's
    unique-peptide feature intensities for that unit+label (quantms
    ``unique_peptides`` policy, no normalization), until OpenMS ``-out_qpx``
    provides the authoritative protein quant. ``top`` bounds the peptides used
    (``0`` = all; set to 3 to mirror the ProteomicsLFQ/IsobaricWorkflow default).
    Each quantified row is stamped with a ``quantification_method`` cv_param so it
    is never mistaken for the authoritative number. Pass either
    ``consensusxml_path`` (loaded here) or an already-loaded ``cm``.
    """
    cm = cm if cm is not None else load_consensus_map(consensusxml_path)

    map_info = _map_info(cm)
    all_runs = sorted({run for run, _ in map_info.values()})
    # grouped_runs unit per run, from the SDRF (fractions grouped); else one unit.
    run_to_grouped = fraction_groups_from_sdrf(sdrf_path)
    units = {tuple(v) for v in run_to_grouped.values()} if run_to_grouped else {tuple(all_runs)}
    # One pg row per (protein group, unit, label): the isobaric channels, or "LFQ".
    labels = sorted({label for _, label in map_info.values()})

    acc_to_pep, acc_to_runs, acc_to_feat = _protein_maps(cm)
    acc_decoy, acc_qvalue, acc_gene, groups = _merge_protein_ids(cm)
    pep_intensity, pep_accs = _peptide_intensities(cm, map_info)
    quant_method = "unnormalized_unique_peptide_sum" if not top else f"unnormalized_unique_peptide_top{top}_sum"

    records: list[dict] = []
    for accs in groups:
        anchor = accs[0]
        group_accs = set(accs)
        peptide_seqs: set[str] = set()
        feats: set[tuple] = set()
        for acc in accs:
            peptide_seqs |= acc_to_pep.get(acc, set())
            feats |= acc_to_feat.get(acc, set())
        n_pep = len(peptide_seqs)
        n_feat = len(feats)
        # Prefer the target_decoy meta; fall back to the accession prefix.
        is_decoy = all(acc_decoy.get(a, _is_decoy_accession(a)) for a in accs)
        qvals = [acc_qvalue[a] for a in accs if a in acc_qvalue]
        global_qvalue = min(qvals) if qvals else None
        genes = [acc_gene[a] for a in accs if a in acc_gene] or None
        # Only the quantification units where this group was actually identified
        # (its peptides appear in a run of that unit) — not every unit.
        group_runs: set[str] = set()
        for acc in accs:
            group_runs |= acc_to_runs.get(acc, set())
        group_units = [unit for unit in units if group_runs.intersection(unit)] or list(units)
        for unit in group_units:
            for label in labels:
                intensity = _protein_intensity(peptide_seqs, group_accs, unit, label, pep_intensity, pep_accs, top)
                cv_params = [{"cv_name": "quantification_method", "cv_value": quant_method}] if intensity is not None else None
                records.append(
                    {
                        "pg_accessions": list(accs),
                        "anchor_protein": anchor,
                        "grouped_runs": list(unit),
                        "label": label,
                        # interim unnormalized total; null when the group has no
                        # unique-peptide signal in this unit+label (see cv_params).
                        "intensity": intensity,
                        "global_qvalue": global_qvalue,
                        "is_decoy": is_decoy,
                        "contaminant": any(_is_contaminant(a) for a in accs),
                        "gg_accessions": genes,
                        "gg_names": genes,
                        "peptide_counts": {"unique_sequences": n_pep, "total_sequences": n_pep},
                        "feature_counts": {"unique_features": n_feat, "total_features": n_feat},
                        "peptides": [{"protein_name": a, "peptide_count": len(acc_to_pep.get(a, set()))} for a in accs],
                        "cv_params": cv_params,
                    }
                )
    return records
