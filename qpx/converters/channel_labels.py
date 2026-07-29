"""Shared quantification channel-label resolution for QuantMS/OpenMS QPX.

Both QPX paths for QuantMS output — the mzTab+MSstats converter
(:mod:`qpx.converters.quantms`) and the enrichment of OpenMS ``-out_qpx``
(:mod:`qpx.converters.openms`) — must put the *same* canonical reporter labels
into ``intensities[].label`` so TMT/iTRAQ/LFQ results are consistent and join
to ``sample``/``run`` the same way DIA-NN (quantmsdiann) output does.

The canonical vocabulary is the shared sdrf-pipelines ``channel_map`` (single
source of truth for label spellings, incl. the TMT10 ch10 = ``TMT131`` vs
TMT11 ch10 = ``TMT131N`` distinction). The plex is resolved, most-solid first,
from the SDRF ``comment[label]`` set (ground truth) and only then from the
channel indices present in the data.
"""

from __future__ import annotations

from typing import Iterable, Optional

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
from sdrf_pipelines.converters.channel_map import CHANNEL_LABELS, normalize_label
from sdrf_pipelines.converters.openms.utils import infer_itraqplex, infer_tmtplex

__all__ = [
    "read_sdrf_labels",
    "experiment_type_from_labels",
    "resolve_channel_labels",
    "relabel_intensities_parquet",
    "normalize_label",
]

_FAMILY_PREFIX = {"TMT": "tmt", "ITRAQ": "itraq", "iTRAQ": "itraq"}


def read_sdrf_labels(sdrf_path: Optional[str]) -> Optional[set[str]]:
    """Return the distinct ``comment[label]`` values declared in the SDRF.

    This is the experiment's ground-truth channel set. Returns ``None`` when no
    SDRF is available or it has no label column.
    """
    if not sdrf_path:
        return None
    try:
        df = pd.read_csv(sdrf_path, sep="\t", dtype=str)
    except Exception:
        return None
    col = next((c for c in df.columns if c.strip().lower() == "comment[label]"), None)
    if col is None:
        return None
    labels = {str(v).strip() for v in df[col].dropna() if str(v).strip()}
    return labels or None


def experiment_type_from_labels(sdrf_labels: Optional[set[str]]) -> str:
    """Classify the labeling type (``"TMT"`` / ``"iTRAQ"`` / ``"LFQ"``) from SDRF labels."""
    if not sdrf_labels:
        return "LFQ"
    upper = {label.upper() for label in sdrf_labels}
    if any("TMT" in label for label in upper):
        return "TMT"
    if any("ITRAQ" in label for label in upper):
        return "iTRAQ"
    return "LFQ"


def resolve_channel_labels(
    experiment_type: str,
    sdrf_labels: Optional[set[str]] = None,
    channel_indices: Optional[Iterable] = None,
) -> dict[int, str]:
    """Resolve the ``{1-based channel index -> canonical label}`` map for a run.

    Resolution order (most-solid first):

    1. From the SDRF ``comment[label]`` set via sdrf-pipelines
       ``infer_tmtplex`` / ``infer_itraqplex`` — robust even when a channel is
       empty in the data.
    2. Fallback: the smallest plex of the reagent family whose channel count
       covers the highest index present in ``channel_indices``.

    Returns ``{}`` for non-isobaric or unresolvable input (caller falls back to
    the raw value / ``"LFQ"``).
    """
    prefix = _FAMILY_PREFIX.get(experiment_type)
    if prefix is None:
        return {}

    # 1. SDRF-declared plex (ground truth).
    if sdrf_labels:
        canonical = {normalize_label(label) for label in sdrf_labels if label}
        infer = infer_tmtplex if prefix == "tmt" else infer_itraqplex
        try:
            return CHANNEL_LABELS[infer(canonical)]
        except (ValueError, KeyError):
            pass

    # 2. Fallback: infer plex from the channel indices present in the data.
    indices: set[int] = set()
    for value in channel_indices or ():
        try:
            indices.add(int(float(value)))
        except (TypeError, ValueError):
            continue
    if not indices:
        return {}
    max_index = max(indices)

    family = {plex: labels for plex, labels in CHANNEL_LABELS.items() if plex.startswith(prefix)}
    for plex in sorted(family, key=lambda p: len(family[p])):
        if len(family[plex]) >= max_index:
            return family[plex]
    return {}


def _relabel_entries(rows, channel_labels: dict[int, str], is_lfq: bool):
    """Relabel the ``label`` field of a list<struct> intensities column.

    OpenMS ``-out_qpx`` writes the *filename* into feature ``intensities[].label``
    and a bare 1-based index (``"1"``, ``"2"``, ...) into pg ``intensities[].label``.
    Recover the channel index from a numeric label, else from the array position
    (channels are emitted in index order), then map to the canonical label.
    LFQ collapses to the single ``"LFQ"`` label.
    """
    out = []
    for row in rows:
        if row is None:
            out.append(None)
            continue
        new_row = []
        for position, entry in enumerate(row):
            label = entry.get("label")
            if is_lfq:
                new_label = "LFQ"
            else:
                text = str(label) if label is not None else ""
                index = int(text) if text.isdigit() else position + 1
                new_label = channel_labels.get(index, normalize_label(text) if text else text)
            updated = dict(entry)
            updated["label"] = new_label
            new_row.append(updated)
        out.append(new_row)
    return out


def relabel_intensities_parquet(
    src_path: str,
    dst_path: str,
    channel_labels: dict[int, str],
    is_lfq: bool,
    columns: tuple[str, ...] = ("intensities", "additional_intensities"),
) -> None:
    """Rewrite ``src_path`` to ``dst_path`` with canonical channel labels.

    Streams by row group so memory stays bounded on large feature tables.
    Non-intensity columns pass through untouched.
    """
    parquet = pq.ParquetFile(src_path)
    writer: Optional[pq.ParquetWriter] = None
    try:
        for group in range(parquet.num_row_groups):
            table = parquet.read_row_group(group)
            for column in columns:
                if column not in table.column_names:
                    continue
                field_index = table.schema.get_field_index(column)
                original = table.column(column)
                relabeled = pa.array(
                    _relabel_entries(original.to_pylist(), channel_labels, is_lfq),
                    type=original.type,
                )
                table = table.set_column(field_index, column, relabeled)
            if writer is None:
                writer = pq.ParquetWriter(dst_path, table.schema)
            writer.write_table(table)
    finally:
        if writer is not None:
            writer.close()
