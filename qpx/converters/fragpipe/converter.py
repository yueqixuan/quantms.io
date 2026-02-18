"""FragPipe orchestrator — composes PSM, Feature, and PG adapters."""

import logging
from datetime import datetime
from pathlib import Path

from qpx._version import __version__
from qpx.core.scores import score_ontology_entries, field_ontology_entries
from qpx.converters.fragpipe.psm_adapter import FragPipePsmAdapter
from qpx.converters.fragpipe.pg_adapter import FragPipePgAdapter

logger = logging.getLogger(__name__)


class FragPipeConverter:
    """Orchestrate full FragPipe conversion to QPX format."""

    def __init__(self, output_directory=None):
        self._output_dir = Path(output_directory) if output_directory else None

    def convert(
        self,
        psm_file=None,
        ion_file=None,
        peptide_file=None,
        pg_file=None,
        sdrf_file=None,
        output_prefix=None,
        batch_size=500_000,
        output_folder=None,
        project_accession=None,
    ):
        out = Path(output_folder) if output_folder else self._output_dir
        if out is None:
            raise ValueError("output_folder or output_directory required")
        out.mkdir(parents=True, exist_ok=True)
        prefix = output_prefix or "fragpipe"

        ontology_entries: list[dict] = []
        produced_structures: list[str] = []

        if psm_file:
            with FragPipePsmAdapter() as adapter:
                adapter.convert(
                    psm_path=str(psm_file),
                    output_path=str(out / f"{prefix}.psm.parquet"),
                    chunksize=batch_size,
                )
                ontology_entries.extend(
                    score_ontology_entries(adapter._discovered_scores, view="psm")
                )
            produced_structures.append("psm")
            logger.info("FragPipe PSM conversion complete")

        if ion_file or peptide_file:
            from qpx.converters.fragpipe.feature_adapter import FragPipeFeatureAdapter

            feature_path = str(ion_file or peptide_file)
            with FragPipeFeatureAdapter() as adapter:
                adapter.convert(
                    feature_path=feature_path,
                    output_path=str(out / f"{prefix}.feature.parquet"),
                    sdrf_path=str(sdrf_file) if sdrf_file else None,
                    chunksize=batch_size,
                )
                ontology_entries.extend(
                    score_ontology_entries(adapter._discovered_scores, view="feature")
                )
            produced_structures.append("feature")
            logger.info("FragPipe feature conversion complete")

        if pg_file:
            with FragPipePgAdapter() as adapter:
                adapter.convert(
                    protein_path=str(pg_file),
                    output_path=str(out / f"{prefix}.pg.parquet"),
                    chunksize=batch_size,
                )
                ontology_entries.extend(
                    score_ontology_entries(adapter._discovered_scores, view="pg")
                )
            produced_structures.append("pg")
            logger.info("FragPipe PG conversion complete")

        if sdrf_file:
            from qpx.converters.sdrf import SdrfConverter

            with SdrfConverter() as conv:
                conv.convert(
                    sdrf_path=str(sdrf_file),
                    sample_output=str(out / f"{prefix}.sample.parquet"),
                    run_output=str(out / f"{prefix}.run.parquet"),
                )
                ontology_entries.extend(conv.run_ontology_entries())
            logger.info("SDRF conversion complete")

        # Add field-level CV term entries
        ontology_entries.extend(field_ontology_entries(view="psm"))

        # Write combined ontology.parquet
        if ontology_entries:
            from qpx.writers.ontology import OntologyWriter

            onto_path = out / f"{prefix}.ontology.parquet"
            with OntologyWriter(onto_path, creator="qpx") as writer:
                writer.write_batch(ontology_entries)
            logger.info(
                "Wrote %d ontology entries to %s", len(ontology_entries), onto_path
            )

        # Write provenance.parquet
        self._write_provenance(out, prefix, produced_structures)

        # Write dataset.parquet
        self._write_dataset(out, prefix, project_accession)

    # ------------------------------------------------------------------
    # Provenance
    # ------------------------------------------------------------------

    @staticmethod
    def _write_provenance(
        output_folder: Path, prefix: str, structures: list[str]
    ) -> None:
        """Write provenance.parquet with FragPipe/MSFragger + QPX conversion steps."""
        from qpx.writers.provenance import ProvenanceWriter

        records = [
            {
                "step_order": 1,
                "step_category": "database_search",
                "step_name": "spectrum_identification",
                "tool_name": "FragPipe/MSFragger",
                "tool_version": None,
                "tool_uri": None,
                "parameters": None,
                "config": None,
                "output_views": structures,
            },
            {
                "step_order": 2,
                "step_category": "format_conversion",
                "step_name": "fragpipe_to_qpx",
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
            "software_name": "FragPipe/MSFragger",
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
