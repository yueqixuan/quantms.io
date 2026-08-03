"""Orchestrate consensusXML (+ SDRF) -> a QPX dataset (feature/psm/pg [+ run/sample]).

Interim quantms path while OpenMS ``-out_qpx`` is pre-1.1. pg carries an interim
unnormalized unique-peptide-sum intensity (until ``-out_qpx`` provides the real
quant). run/sample come from the SDRF (reusing :class:`SdrfConverter`) when an
SDRF is provided; when it is, the consensusXML channels are also checked against
the SDRF ``comment[label]`` and mismatches are logged as warnings.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional

from qpx.converters.openms_consensus.feature_adapter import (
    check_channels_vs_sdrf,
    consensus_features_to_records,
    load_consensus_map,
)
from qpx.converters.openms_consensus.pg_adapter import consensus_protein_groups_to_records
from qpx.converters.openms_consensus.psm_adapter import consensus_psms_to_records
from qpx.writers.feature import FeatureWriter
from qpx.writers.pg import PgWriter
from qpx.writers.psm import PsmWriter

_log = logging.getLogger(__name__)

_STRUCTURE_ALL = ("feature", "psm", "pg", "run", "sample")


def _validate_structures(structures: tuple[str, ...], sdrf_path: Optional[str]) -> None:
    """Reject unknown structure names and run/sample requested without an SDRF."""
    unknown = [s for s in structures if s not in _STRUCTURE_ALL]
    if unknown:
        raise ValueError(f"Unknown structure(s) {unknown}; valid values are {list(_STRUCTURE_ALL)}")
    if ("run" in structures or "sample" in structures) and not sdrf_path:
        raise ValueError("The 'run'/'sample' structures require an SDRF (pass sdrf_path)")


class OpenMSConsensusConverter:  # pylint: disable=too-few-public-methods
    """consensusXML + SDRF -> QPX views.

    A single-entry orchestrator (``convert``) — the interim counterpart to the
    other converter classes, kept as a class for call-site symmetry with them.
    """

    def convert(
        self,
        consensusxml_path: str,
        output_folder: str,
        output_prefix: str = "openms",
        sdrf_path: Optional[str] = None,
        structures: tuple[str, ...] = _STRUCTURE_ALL,
        creator: str = "openms-consensus",
        pg_top: int = 0,
    ) -> dict[str, Path]:
        """Write the requested QPX views and return ``{structure: parquet path}``.

        ``structures`` selects which of feature/psm/pg/run/sample to emit. pg
        carries an interim unnormalized unique-peptide-sum intensity; ``pg_top``
        bounds the peptides used (0 = all; 3 mirrors ProteomicsLFQ/IsobaricWorkflow).
        """
        _validate_structures(structures, sdrf_path)

        out = Path(output_folder)
        out.mkdir(parents=True, exist_ok=True)
        written: dict[str, Path] = {}

        # Parse the consensusXML once and share the in-memory map across the
        # feature/psm/pg adapters (they only read it) — the load dominates the
        # convert cost, so re-parsing per adapter would triple it.
        cm = None
        if {"feature", "psm", "pg"}.intersection(structures):
            cm = load_consensus_map(consensusxml_path)
            # Ground-truth check: the channels read from the consensusXML maps must
            # match the SDRF comment[label] set — warn on any inconsistency (e.g. a
            # mis-declared plex or the wrong SDRF) rather than silently trusting one.
            if sdrf_path:
                for msg in check_channels_vs_sdrf(cm, sdrf_path):
                    _log.warning("consensusXML/SDRF channel mismatch: %s", msg)

        if "feature" in structures:
            recs = consensus_features_to_records(cm=cm)
            path = out / f"{output_prefix}.feature.parquet"
            with FeatureWriter(str(path), creator=creator) as w:
                if recs:
                    w.write_batch(recs)
            written["feature"] = path

        if "psm" in structures:
            recs = consensus_psms_to_records(cm=cm)
            path = out / f"{output_prefix}.psm.parquet"
            with PsmWriter(str(path), creator=creator) as w:
                if recs:
                    w.write_batch(recs)
            written["psm"] = path

        if "pg" in structures:
            recs = consensus_protein_groups_to_records(sdrf_path=sdrf_path, cm=cm, top=pg_top)
            path = out / f"{output_prefix}.pg.parquet"
            with PgWriter(str(path), creator=creator) as w:
                if recs:
                    w.write_batch(recs)
            written["pg"] = path

        if "run" in structures or "sample" in structures:
            from qpx.converters.sdrf import SdrfConverter

            with SdrfConverter() as sdrf_conv:
                sdrf_conv.convert(
                    sdrf_path=sdrf_path,
                    sample_output=str(out / f"{output_prefix}.sample.parquet"),
                    run_output=str(out / f"{output_prefix}.run.parquet"),
                )
            # Record only the structures the caller actually requested.
            if "run" in structures:
                written["run"] = out / f"{output_prefix}.run.parquet"
            if "sample" in structures:
                written["sample"] = out / f"{output_prefix}.sample.parquet"

        return written
