"""MaxQuant orchestrator — composes PSM, Feature, and PG adapters."""

import logging
from datetime import datetime
from pathlib import Path

from qpx._version import __version__
from qpx.core.scores import score_ontology_entries, field_ontology_entries
from qpx.converters.base import resolve_columns
from qpx.converters.maxquant.constants import TOOL_NAME, FIELD_MAPPINGS
from qpx.converters.maxquant.psm_adapter import MaxQuantPsmAdapter
from qpx.converters.maxquant.feature_adapter import MaxQuantFeatureAdapter
from qpx.converters.maxquant.pg_adapter import MaxQuantPgAdapter

logger = logging.getLogger(__name__)


class MaxQuantConverter:
    """Orchestrate full MaxQuant conversion to QPX format."""

    def __init__(self, memory_limit_gb=None):
        self._memory = f"{int(memory_limit_gb)}GB" if memory_limit_gb else "16GB"
        self._resolved_mappings: dict[str, str] = {}

    def convert(
        self,
        output_folder,
        msms_file=None,
        evidence_file=None,
        protein_groups_file=None,
        sdrf_file=None,
        output_prefix=None,
        structures=None,
        protein_file=None,
        batch_size=100_000,
        n_workers=None,
        spectral_data=False,
        standardized_intensities=False,
        project_accession=None,
    ):
        output_folder = Path(output_folder)
        output_folder.mkdir(parents=True, exist_ok=True)
        prefix = output_prefix or "maxquant"

        if structures is None:
            structures = []
            if msms_file:
                structures.append("psm")
            if evidence_file:
                structures.append("feature")
            if protein_groups_file:
                structures.append("pg")

        ontology_entries: list[dict] = []

        # Convert SDRF for sample/run metadata (if SDRF has required columns)
        if sdrf_file:
            try:
                from qpx.converters.sdrf import SdrfConverter

                with SdrfConverter() as sdrf_conv:
                    sdrf_conv.convert(
                        sdrf_path=str(sdrf_file),
                        sample_output=str(output_folder / f"{prefix}.sample.parquet"),
                        run_output=str(output_folder / f"{prefix}.run.parquet"),
                    )
                    ontology_entries.extend(sdrf_conv.run_ontology_entries())
                logger.info("SDRF conversion complete")
            except Exception as exc:
                logger.warning(
                    "SDRF conversion skipped (incomplete SDRF?): %s", exc
                )
                # Clean up any corrupt partial files left by the failed write
                for suffix in (".sample.parquet", ".run.parquet"):
                    corrupt = output_folder / f"{prefix}{suffix}"
                    if corrupt.exists():
                        corrupt.unlink()
                        logger.debug("Removed corrupt %s", corrupt)

        if "psm" in structures and msms_file:
            with MaxQuantPsmAdapter(duckdb_memory=self._memory) as adapter:
                adapter.convert(
                    msms_path=str(msms_file),
                    output_path=str(output_folder / f"{prefix}.psm.parquet"),
                    chunksize=batch_size,
                    spectral_data=spectral_data,
                )
                ontology_entries.extend(
                    score_ontology_entries(adapter._discovered_scores, view="psm")
                )
                # Accumulate resolved column mappings for provenance
                if hasattr(adapter, "_resolved"):
                    self._resolved_mappings.update(adapter._resolved)
            logger.info("MaxQuant PSM conversion complete")

        if "feature" in structures and evidence_file:
            with MaxQuantFeatureAdapter(duckdb_memory=self._memory) as adapter:
                adapter.convert(
                    evidence_path=str(evidence_file),
                    output_path=str(output_folder / f"{prefix}.feature.parquet"),
                    sdrf_path=str(sdrf_file) if sdrf_file else None,
                    chunksize=batch_size,
                )
                ontology_entries.extend(
                    score_ontology_entries(adapter._discovered_scores, view="feature")
                )
                if hasattr(adapter, "_resolved"):
                    self._resolved_mappings.update(adapter._resolved)
            logger.info("MaxQuant feature conversion complete")

        if "pg" in structures and protein_groups_file:
            with MaxQuantPgAdapter(duckdb_memory=self._memory) as adapter:
                adapter.convert(
                    protein_groups_path=str(protein_groups_file),
                    output_path=str(output_folder / f"{prefix}.pg.parquet"),
                    sdrf_path=str(sdrf_file) if sdrf_file else None,
                    chunksize=batch_size,
                )
                ontology_entries.extend(
                    score_ontology_entries(adapter._discovered_scores, view="pg")
                )
                if hasattr(adapter, "_resolved"):
                    self._resolved_mappings.update(adapter._resolved)
            logger.info("MaxQuant PG conversion complete")

        # Add field-level CV term entries with source provenance
        ontology_entries.extend(
            field_ontology_entries(
                view="psm",
                resolved_mappings=self._resolved_mappings,
                tool_name=TOOL_NAME,
            )
        )

        # Write combined ontology.parquet
        if ontology_entries:
            from qpx.writers.ontology import OntologyWriter

            onto_path = output_folder / f"{prefix}.ontology.parquet"
            with OntologyWriter(onto_path, creator="qpx") as writer:
                writer.write_batch(ontology_entries)
            logger.info(
                "Wrote %d ontology entries to %s", len(ontology_entries), onto_path
            )

        # Write provenance.parquet
        self._write_provenance(output_folder, prefix, structures)

        # Write dataset.parquet
        self._write_dataset(output_folder, prefix, project_accession)

    # ------------------------------------------------------------------
    # Provenance
    # ------------------------------------------------------------------

    @staticmethod
    def _write_provenance(
        output_folder: Path, prefix: str, structures: list[str]
    ) -> None:
        """Write provenance.parquet with MaxQuant + QPX conversion steps."""
        from qpx.writers.provenance import ProvenanceWriter

        records = [
            {
                "step_order": 1,
                "step_category": "database_search",
                "step_name": "spectrum_identification",
                "tool_name": "MaxQuant",
                "tool_version": None,
                "tool_uri": None,
                "parameters": None,
                "config": None,
                "output_views": structures,
            },
            {
                "step_order": 2,
                "step_category": "format_conversion",
                "step_name": "maxquant_to_qpx",
                "tool_name": "qpx",
                "tool_version": __version__,
                "tool_uri": None,
                "parameters": None,
                "config": None,
                "output_views": structures + ["sample", "run", "ontology"],
            },
        ]
        prov_path = output_folder / f"{prefix}.provenance.parquet"
        with ProvenanceWriter(prov_path, creator="qpx") as writer:
            writer.write_batch(records)
        logger.info("Wrote %d provenance steps to %s", len(records), prov_path)

    # ------------------------------------------------------------------
    # Dataset metadata
    # ------------------------------------------------------------------

    @staticmethod
    def _write_dataset(
        output_folder: Path,
        prefix: str,
        project_accession: str | None,
    ) -> None:
        """Write dataset.parquet with project-level metadata."""
        from qpx.writers.dataset import DatasetWriter

        record = {
            "project_accession": project_accession or "unknown",
            "project_title": None,
            "project_description": None,
            "pubmed_id": None,
            "software_name": "MaxQuant",
            "software_version": None,
            "creation_date": datetime.now().isoformat(),
            "file_checksums": None,
            "file_row_counts": None,
            "file_sizes_bytes": None,
            "total_structures": None,
            "packaged_at": None,
        }
        ds_path = output_folder / f"{prefix}.dataset.parquet"
        with DatasetWriter(ds_path, creator="qpx") as writer:
            writer.write_batch([record])
        logger.info("Wrote dataset metadata to %s", ds_path)
