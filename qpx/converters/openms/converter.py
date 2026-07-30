"""OpenMS QPX converter — enriches ProteomicsLFQ ``-out_qpx`` output.

OpenMS ``ProteomicsLFQ -out_qpx`` already produces QPX-compliant
``psm.parquet``, ``feature.parquet``, and ``pg.parquet``.  This converter
validates those files, copies them to the output folder, and generates the
missing metadata tables (``run``, ``sample``, ``ontology``, ``provenance``,
``dataset``) from an SDRF file.
"""

from __future__ import annotations

import logging
import shutil
from pathlib import Path
from typing import Optional

import pyarrow as pa
import pyarrow.parquet as pq

from qpx._version import __version__
from qpx.converters.channel_labels import (
    channel_labels_from_consensusxml,
    experiment_type_from_labels,
    fraction_groups_from_consensusxml,
    read_sdrf_labels,
    relabel_intensities_parquet,
    resolve_channel_labels,
)
from qpx.converters.orchestrator import BaseOrchestrator
from qpx.converters.sdrf import SdrfConverter
from qpx.core.constants import FEATURE, ONTOLOGY, PG, PSM, RUN, SAMPLE
from qpx.core.data.loader import load_schema
from qpx.core.scores import score_ontology_entries
from qpx.writers.base import parquet_write_options

logger = logging.getLogger(__name__)

# Mapping from QPX view name to the schema loader name
_VIEW_SCHEMAS = {
    PSM: "psm",
    FEATURE: "feature",
    PG: "pg",
}

# There is no PSI-MS CV term for OpenMS's experimental-design ``fraction_group``
# grouping key, so it is captured under a local quantms CV accession placeholder
# until a PSI-MS term is minted (bigbio/qpx#221, OpenMS#9817). The qpx cv_param
# struct carries only (cv_name, cv_value); the local accession is documented here
# and recorded via ``cv_name`` so the value round-trips unambiguously.
FRACTION_GROUP_CV_ACCESSION = "quantms:fraction_group"
FRACTION_GROUP_CV_NAME = "fraction_group"

# Run-file extensions stripped when matching consensusXML ``<map name>`` values
# against pg ``grouped_runs`` / feature ``run_file_name`` (which -out_qpx may
# write with or without an extension).
_RUN_EXTENSIONS = (".mzml", ".raw", ".mzxml", ".d", ".wiff", ".dia")


def _run_key_candidates(name: str | None) -> list[str]:
    """Return match candidates for a run/map filename (raw, basename, stem)."""
    if not name:
        return []
    text = str(name).strip()
    if not text:
        return []
    candidates = [text]
    base = text.rsplit("/", 1)[-1].rsplit("\\", 1)[-1]
    if base and base not in candidates:
        candidates.append(base)
    lower = base.lower()
    for ext in _RUN_EXTENSIONS:
        if lower.endswith(ext):
            stem = base[: -len(ext)]
            if stem and stem not in candidates:
                candidates.append(stem)
            break
    return candidates


def _fraction_group_lookup(fraction_groups: dict[str, dict[str, str]]) -> dict[str, str]:
    """Flatten ``{map_name -> design}`` into ``{run_key -> fraction_group}``.

    Each map name contributes several keys (raw, basename, extension-stripped)
    so a ``grouped_runs`` / ``run_file_name`` value matches regardless of whether
    -out_qpx kept the run-file extension.
    """
    lookup: dict[str, str] = {}
    for name, design in fraction_groups.items():
        fraction_group = design.get("fraction_group")
        if not fraction_group:
            continue
        for key in _run_key_candidates(name):
            lookup.setdefault(key, fraction_group)
    return lookup


def _fraction_group_for(run_names, lookup: dict[str, str]) -> str | None:
    """First fraction_group matching any of ``run_names`` (str or iterable)."""
    if run_names is None:
        return None
    names = [run_names] if isinstance(run_names, str) else list(run_names)
    for name in names:
        for key in _run_key_candidates(name):
            if key in lookup:
                return lookup[key]
    return None


def _add_fraction_group_cv_params(
    path: Path,
    run_column: str,
    fraction_groups: dict[str, dict[str, str]],
    compression: str,
) -> int:
    """Append a ``fraction_group`` cv_param to each row of a pg/feature parquet.

    ``run_column`` names the column identifying the row's run(s): ``grouped_runs``
    (a list, pg) or ``run_file_name`` (a scalar, feature). The row's OpenMS
    ``fraction_group`` is looked up from the consensusXML design and appended to
    the existing ``cv_params`` list (existing cv_params are preserved, duplicates
    of ``fraction_group`` are not re-added). Best-effort: rows whose run does not
    resolve to a fraction_group are left untouched. Streams by row group so
    memory stays bounded. Returns the number of rows annotated.
    """
    lookup = _fraction_group_lookup(fraction_groups)
    if not lookup:
        return 0

    parquet = pq.ParquetFile(str(path))
    schema = parquet.schema_arrow
    if run_column not in schema.names or "cv_params" not in schema.names:
        return 0

    tmp = path.with_suffix(".fgroup.tmp")
    writer = pq.ParquetWriter(str(tmp), schema, **parquet_write_options(schema, compression))
    annotated = 0
    try:
        for group in range(parquet.num_row_groups):
            table = parquet.read_row_group(group)
            run_values = table.column(run_column).to_pylist()
            cv_values = table.column("cv_params").to_pylist()
            new_cv = []
            for run_names, cv_params in zip(run_values, cv_values):
                fraction_group = _fraction_group_for(run_names, lookup)
                params = list(cv_params) if cv_params else []
                if fraction_group is not None and not any((p or {}).get("cv_name") == FRACTION_GROUP_CV_NAME for p in params):
                    params.append({"cv_name": FRACTION_GROUP_CV_NAME, "cv_value": str(fraction_group)})
                    annotated += 1
                new_cv.append(params or None)
            field_index = table.schema.get_field_index("cv_params")
            cv_array = pa.array(new_cv, type=table.column("cv_params").type)
            table = table.set_column(field_index, "cv_params", cv_array)
            writer.write_table(table)
    finally:
        writer.close()
    tmp.replace(path)
    return annotated


def _discover_parquet(qpx_dir: Path) -> dict[str, Path]:
    """Discover QPX parquet files in a directory.

    Looks for files matching ``*.psm.parquet``, ``*.feature.parquet``,
    and ``*.pg.parquet``.  Returns a dict mapping view name to path.
    """
    found: dict[str, Path] = {}
    for view in (PSM, FEATURE, PG):
        matches = sorted(qpx_dir.glob(f"*.{view}.parquet"))
        if matches:
            found[view] = matches[0]
            if len(matches) > 1:
                logger.warning(
                    "Multiple %s parquet files found in %s; using %s",
                    view,
                    qpx_dir,
                    matches[0].name,
                )
    return found


def _collect_score_names(table_path: Path) -> set[str]:
    """Read ``additional_scores`` from a parquet file and collect score names."""
    table = pq.read_table(str(table_path), columns=["additional_scores"])
    names: set[str] = set()
    col = table.column("additional_scores")
    for row in col.to_pylist():
        if row:
            for score in row:
                if score and "name" in score and score["name"]:
                    names.add(score["name"])
    return names


def _validate_core(discovered: dict[str, Path]) -> None:
    """Validate each discovered parquet file against its QPX schema."""
    for view, path in discovered.items():
        schema = load_schema(_VIEW_SCHEMAS[view])
        table = pq.read_table(str(path))
        result = schema.validate_full(table)
        if not result.is_valid:
            errors = "; ".join(i.message for i in result.errors)
            raise ValueError(f"Validation failed for {path.name}: {errors}")
        logger.info(
            "%s: %s (%d rows, %d errors, %d warnings)",
            view,
            path.name,
            table.num_rows,
            len(result.errors),
            len(result.warnings),
        )


def _copy_core(
    discovered: dict[str, Path],
    output_folder: Path,
    output_prefix: str,
    channel_labels: Optional[dict[int, str]] = None,
    is_lfq: bool | None = None,
    compression: str = "zstd",
) -> dict[str, Path]:
    """
    Copy core parquet files to the output directory.

    ``feature`` and ``pg`` carry ``intensities[].label`` — OpenMS ``-out_qpx``
    writes the run filename (feature) or a bare channel index (pg) there, so
    those two are rewritten with canonical channel labels while copying when
    the experiment type is known. ``is_lfq=None`` preserves their source labels
    because no SDRF evidence is available. The rest pass through untouched.
    """
    output_paths: dict[str, Path] = {}
    for view, src_path in discovered.items():
        dst = output_folder / f"{output_prefix}.{view}.parquet"
        if view in (FEATURE, PG) and is_lfq is not None:
            if src_path.resolve() == dst.resolve():
                tmp = dst.with_suffix(".relabel.tmp")
                relabel_intensities_parquet(
                    str(src_path),
                    str(tmp),
                    channel_labels or {},
                    is_lfq,
                    compression=compression,
                )
                tmp.replace(dst)
            else:
                relabel_intensities_parquet(
                    str(src_path),
                    str(dst),
                    channel_labels or {},
                    is_lfq,
                    compression=compression,
                )
            logger.info("Relabeled channel labels in %s -> %s", src_path.name, dst.name)
        elif src_path.resolve() != dst.resolve():
            shutil.copy2(str(src_path), str(dst))
            logger.info("Copied %s -> %s", src_path.name, dst.name)
        else:
            logger.info("Skipping copy for %s (same location)", view)
        output_paths[view] = dst
    return output_paths


class OpenMSConverter(BaseOrchestrator):
    """Orchestrate OpenMS QPX enrichment to a full QPX dataset.

    Usage::

        converter = OpenMSConverter(
            qpx_dir="path/to/openms_qpx_output",
            sdrf_path="metadata.sdrf.tsv",
        )
        converter.convert(
            output_folder="./qpx_full",
            output_prefix="openms",
        )
    """

    def __init__(
        self,
        qpx_dir: str | Path,
        sdrf_path: str | Path | None = None,
        consensusxml_path: str | Path | None = None,
        compression: str = "zstd",
    ):
        """Initialize the OpenMS QPX converter.

        Parameters
        ----------
        qpx_dir : str or Path
            Directory containing OpenMS ``-out_qpx`` parquet files.
        sdrf_path : str, Path or None
            Optional SDRF metadata file for sample/run generation.
        consensusxml_path : str, Path or None
            Optional OpenMS ``.consensusXML`` (the ``-out_cxml`` companion of
            ``-out_qpx``). When given, its ColumnHeaders provide the
            authoritative channel count/order for relabeling isobaric channels;
            otherwise the plex is resolved from the SDRF + data indices.
        compression : str
            Parquet compression codec (default ``zstd``).

        """
        self.qpx_dir = Path(qpx_dir)
        self.sdrf_path = str(sdrf_path) if sdrf_path else None
        self.consensusxml_path = str(consensusxml_path) if consensusxml_path else None
        self._compression = compression

    def convert(
        self,
        output_folder: str | Path,
        output_prefix: str = "openms",
        project_accession: str | None = None,
    ) -> None:
        """Run the full enrichment pipeline.

        1. Discover and validate core QPX parquet files
        2. Copy core files to output folder (if needed)
        3. Convert SDRF to run.parquet + sample.parquet
        4. Collect score names for ontology
        5. Write ontology.parquet, provenance.parquet, dataset.parquet
        """
        output_folder = Path(output_folder)
        output_folder.mkdir(parents=True, exist_ok=True)

        discovered = self.discover_and_validate()

        # Resolve canonical channel labels so the OpenMS -out_qpx run-filename /
        # bare-index labels become TMT126.. / LFQ, consistent with the mzTab and
        # DIA-NN (quantmsdiann) QPX paths. Prefer the consensusXML ColumnHeaders
        # (authoritative channel count/order) when available; otherwise resolve
        # from the SDRF-declared plex.
        sdrf_labels = read_sdrf_labels(self.sdrf_path)
        experiment_type = experiment_type_from_labels(sdrf_labels) if sdrf_labels else None
        channel_labels = {}
        if self.consensusxml_path and experiment_type:
            channel_labels = channel_labels_from_consensusxml(self.consensusxml_path, experiment_type, sdrf_labels)
            if channel_labels:
                logger.info("Resolved %d channels from consensusXML", len(channel_labels))
        if not channel_labels and experiment_type:
            channel_labels = resolve_channel_labels(experiment_type, sdrf_labels)
        is_lfq = experiment_type == "LFQ" if experiment_type else None
        if experiment_type is None:
            logger.info("No SDRF channel labels available; preserving OpenMS intensity labels")
        output_paths = _copy_core(
            discovered,
            output_folder,
            output_prefix,
            channel_labels,
            is_lfq,
            self._compression,
        )

        # Capture OpenMS's experimental-design ``fraction_group`` (the
        # replicate/fraction grouping key) from the consensusXML as a cv_param on
        # pg + feature rows. Best-effort: no-op for isobaric / pre-design
        # consensusXML files (no fraction_group UserParams) and when the file is
        # absent. See bigbio/qpx#221, OpenMS#9817.
        self._capture_fraction_groups(output_paths)

        ontology_entries = self._convert_sdrf(output_folder, output_prefix)
        ontology_entries.extend(self._collect_ontology(output_paths))

        self._write_ontology(output_folder, output_prefix, ontology_entries)
        structures = list(discovered.keys())
        provenance_records = self._build_provenance(structures, experiment_type)
        self._write_provenance(output_folder, output_prefix, provenance_records)
        self._write_dataset(
            output_folder,
            output_prefix,
            project_accession,
            software_name=provenance_records[0]["tool_name"],
            software_version=None,
            provenance_records=provenance_records,
        )
        self._write_mudata(output_folder, output_prefix)
        logger.info("OpenMS QPX enrichment complete -> %s", output_folder)

    def discover_and_validate(self) -> dict[str, Path]:
        """Discover and validate core QPX parquet files."""
        discovered = _discover_parquet(self.qpx_dir)
        if not discovered:
            raise FileNotFoundError(
                f"No QPX parquet files (*.psm.parquet, *.feature.parquet, *.pg.parquet) found in {self.qpx_dir}"
            )
        logger.info(
            "Discovered %d core QPX file(s): %s",
            len(discovered),
            ", ".join(f"{v}={p.name}" for v, p in discovered.items()),
        )
        _validate_core(discovered)
        return discovered

    def _capture_fraction_groups(self, output_paths: dict[str, Path]) -> None:
        """Annotate pg + feature rows with the consensusXML ``fraction_group``.

        Best-effort and side-effect free when no consensusXML is supplied or the
        file declares no ``fraction_group`` UserParams (isobaric / pre-design).
        """
        if not self.consensusxml_path:
            return
        fraction_groups = fraction_groups_from_consensusxml(self.consensusxml_path)
        if not fraction_groups:
            return
        logger.info("Resolved fraction_group for %d run(s) from consensusXML", len(fraction_groups))
        for view, run_column in ((PG, "grouped_runs"), (FEATURE, "run_file_name")):
            path = output_paths.get(view)
            if path is None or not path.exists():
                continue
            annotated = _add_fraction_group_cv_params(path, run_column, fraction_groups, self._compression)
            if annotated:
                logger.info("Added fraction_group cv_param to %d %s row(s)", annotated, view)

    def _convert_sdrf(
        self,
        output_folder: Path,
        output_prefix: str,
    ) -> list[dict]:
        """Convert SDRF to run + sample, returning ontology entries."""
        if not self.sdrf_path:
            logger.warning("No SDRF path provided — skipping sample/run conversion")
            return []
        with SdrfConverter(compression=self._compression) as sdrf_conv:
            sdrf_conv.convert(
                sdrf_path=self.sdrf_path,
                sample_output=str(output_folder / f"{output_prefix}.sample.parquet"),
                run_output=str(output_folder / f"{output_prefix}.run.parquet"),
            )
            entries = list(sdrf_conv.run_ontology_entries())
        logger.info("SDRF conversion complete (sample + run)")
        return entries

    @staticmethod
    def _collect_ontology(output_paths: dict[str, Path]) -> list[dict]:
        """Collect score-based ontology entries from core tables."""
        entries: list[dict] = []
        all_scores: set[str] = set()
        for view, path in output_paths.items():
            try:
                scores = _collect_score_names(path)
            except (KeyError, pa.ArrowInvalid):
                logger.debug("No additional_scores in %s", path.name)
                continue
            if scores:
                entries.extend(score_ontology_entries(scores, view=view))
                all_scores |= scores
        if all_scores:
            logger.info("Discovered scores: %s", sorted(all_scores))
        return entries

    @staticmethod
    def _build_provenance(structures: list[str], experiment_type: str | None = None) -> list[dict]:
        """Build provenance records for OpenMS + QPX conversion."""
        is_isobaric = experiment_type in {"TMT", "iTRAQ"}
        step_name = "isobaric_quantification" if is_isobaric else "label_free_quantification"
        tool_name = "OpenMS/IsobaricWorkflow" if is_isobaric else "OpenMS/ProteomicsLFQ"
        return [
            {
                "step_order": 1,
                "step_category": "quantification",
                "step_name": step_name,
                "tool_name": tool_name,
                "tool_version": None,
                "tool_uri": None,
                "parameters": None,
                "config": None,
                "output_views": structures,
            },
            {
                "step_order": 2,
                "step_category": "format_conversion",
                "step_name": "openms_qpx_enrichment",
                "tool_name": "qpx",
                "tool_version": __version__,
                "tool_uri": None,
                "parameters": None,
                "config": None,
                "output_views": [SAMPLE, RUN, ONTOLOGY],
            },
        ]
