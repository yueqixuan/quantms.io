"""Update QPX dataset metadata from a revised SDRF file.

Re-generates sample.parquet (and optionally run.parquet metadata fields)
from a new SDRF, while protecting fields that affect identification or
quantification.

Protected fields (changes blocked with warning):
    - instrument, enzymes, modification_parameters, fraction
    - label / sample-run channel mapping
    - comment[data file] (run identity)

Allowed fields (sample metadata updates):
    - disease, organism_part, cell_type, cell_line, sex, age,
      developmental_stage, ancestry, individual, sample_description
    - Any additional characteristics[*] columns
"""

from __future__ import annotations

import logging
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq

logger = logging.getLogger(__name__)

# Fields in run.parquet that affect identification/quantification.
# If these differ between old and new SDRF, we warn and skip.
PROTECTED_RUN_FIELDS = {
    "instrument",
    "enzymes",
    "modification_parameters",
    "fraction",
    "dissociation_method",
}

# Fields checked at the sample-run mapping level.
# If the mapping of samples to runs changes, we block.
PROTECTED_MAPPING_FIELDS = {
    "comment[data file]",
    "comment[label]",
    "source name",
}


def _check_protected_changes(
    old_sdrf_path: str | Path | None,
    new_sdrf_path: str | Path,
) -> list[str]:
    """Compare old vs new SDRF for protected field changes.

    Returns list of warning messages. Empty list means safe to proceed.
    If old_sdrf_path is None, skip comparison (first-time update).
    """
    if old_sdrf_path is None:
        return []

    import pandas as pd

    old = pd.read_csv(old_sdrf_path, sep="\t")
    new = pd.read_csv(new_sdrf_path, sep="\t")

    old.columns = [c.lower().strip() for c in old.columns]
    new.columns = [c.lower().strip() for c in new.columns]

    warnings = []

    # Check that data file references haven't changed
    old_files = set(old.get("comment[data file]", []))
    new_files = set(new.get("comment[data file]", []))
    if old_files and old_files != new_files:
        added = new_files - old_files
        removed = old_files - new_files
        if removed:
            warnings.append(
                f"BLOCKED: {len(removed)} data files removed from SDRF — "
                f"this would break run references. Removed: {sorted(removed)[:5]}"
            )
        if added:
            warnings.append(
                f"WARNING: {len(added)} new data files in SDRF — "
                f"these won't have corresponding feature/PG data."
            )

    # Check sample-run mapping (source name ↔ data file ↔ label)
    for col in ["comment[label]"]:
        if col in old.columns and col in new.columns:
            old_map = set(
                zip(
                    old.get("source name", []),
                    old.get("comment[data file]", []),
                    old[col],
                )
            )
            new_map = set(
                zip(
                    new.get("source name", []),
                    new.get("comment[data file]", []),
                    new[col],
                )
            )
            if old_map != new_map:
                warnings.append(
                    f"BLOCKED: sample-run-label mapping changed — "
                    f"this would invalidate intensity assignments."
                )

    # Check protected run-level fields
    for col_pattern in [
        "comment[instrument]",
        "comment[cleavage agent details]",
        "comment[fraction identifier]",
        "comment[dissociation method]",
    ]:
        if col_pattern in old.columns and col_pattern in new.columns:
            old_vals = set(old[col_pattern].dropna().unique())
            new_vals = set(new[col_pattern].dropna().unique())
            if old_vals != new_vals:
                warnings.append(
                    f"BLOCKED: '{col_pattern}' changed from {old_vals} to {new_vals} — "
                    f"this affects identification/quantification."
                )

    # Check modification parameters
    mod_cols_old = [c for c in old.columns if "modification parameters" in c]
    mod_cols_new = [c for c in new.columns if "modification parameters" in c]
    if mod_cols_old or mod_cols_new:
        old_mods = set()
        for c in mod_cols_old:
            old_mods.update(old[c].dropna().unique())
        new_mods = set()
        for c in mod_cols_new:
            new_mods.update(new[c].dropna().unique())
        if old_mods != new_mods:
            warnings.append(
                "BLOCKED: modification parameters changed — "
                "this would invalidate PSM/feature identifications."
            )

    return warnings


def update_metadata(
    dataset_path: Path,
    new_sdrf_path: Path,
    old_sdrf_path: Path | None = None,
    force: bool = False,
) -> dict:
    """Update sample.parquet and run.parquet metadata from a new SDRF.

    Args:
        dataset_path: QPX dataset directory.
        new_sdrf_path: Path to the updated SDRF TSV file.
        old_sdrf_path: Path to the original SDRF (for safety checks).
            If None, skips protected-field validation.
        force: If True, apply changes even with protected-field warnings.

    Returns:
        Dict with stats: {"samples_updated", "runs_updated", "warnings"}.
    """
    from qpx.converters.sdrf import SdrfConverter

    # Step 1: Check for protected field changes
    check_warnings = _check_protected_changes(old_sdrf_path, new_sdrf_path)
    blocked = [w for w in check_warnings if w.startswith("BLOCKED")]

    if blocked and not force:
        for w in check_warnings:
            logger.warning(w)
        raise ValueError(
            f"{len(blocked)} protected field(s) changed in the new SDRF. "
            f"Use --force to override (data integrity not guaranteed)."
        )

    for w in check_warnings:
        if w.startswith("WARNING"):
            logger.warning(w)

    # Step 2: Re-generate sample.parquet and run.parquet from new SDRF
    sample_path = dataset_path / "quantms.sample.parquet"
    run_path = dataset_path / "quantms.run.parquet"

    # Back up originals
    if sample_path.exists():
        backup = dataset_path / "quantms.sample.parquet.bak"
        if not backup.exists():
            import shutil

            shutil.copy2(sample_path, backup)
            logger.info("Backed up original sample.parquet → %s", backup.name)

    if run_path.exists():
        backup = dataset_path / "quantms.run.parquet.bak"
        if not backup.exists():
            import shutil

            shutil.copy2(run_path, backup)
            logger.info("Backed up original run.parquet → %s", backup.name)

    # Step 3: Convert new SDRF
    converter = SdrfConverter()
    converter.convert(
        sdrf_path=str(new_sdrf_path),
        sample_output=str(sample_path),
        run_output=str(run_path),
    )

    # Count results
    n_samples = len(pq.read_table(str(sample_path)))
    n_runs = len(pq.read_table(str(run_path)))

    logger.info(
        "Updated metadata: %d samples, %d runs from %s",
        n_samples, n_runs, new_sdrf_path.name,
    )

    return {
        "samples_updated": n_samples,
        "runs_updated": n_runs,
        "warnings": check_warnings,
    }
