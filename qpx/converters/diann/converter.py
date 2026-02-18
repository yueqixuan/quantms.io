"""DIA-NN orchestrator — composes Feature and PG adapters."""

import logging
from datetime import datetime
from pathlib import Path

from qpx._version import __version__
from qpx.core.scores import score_ontology_entries, field_ontology_entries
from qpx.converters.diann.feature_adapter import DiannFeatureAdapter
from qpx.converters.diann.pg_adapter import DiannPgAdapter

logger = logging.getLogger(__name__)


class DiaNNConverter:
    """Orchestrate full DIA-NN conversion to QPX format."""

    def __init__(
        self,
        report_path,
        sdrf_path,
        duckdb_max_memory=None,
        duckdb_threads=None,
    ):
        self.report_path = str(report_path)
        self.sdrf_path = str(sdrf_path)
        self._memory = duckdb_max_memory or "16GB"
        self._threads = duckdb_threads or 4
        self._ontology_entries: list[dict] = []

    def convert_features(
        self,
        mzml_info_folder,
        qvalue_threshold=0.05,
        output_folder=".",
        output_prefix=None,
        protein_file=None,
        partitions=None,
        batch_size=100,
    ):
        output_folder = Path(output_folder)
        prefix = output_prefix or "diann"
        with DiannFeatureAdapter(
            duckdb_memory=self._memory, duckdb_threads=self._threads
        ) as adapter:
            adapter.convert(
                diann_report=self.report_path,
                sdrf_path=self.sdrf_path,
                mzml_info_folder=str(mzml_info_folder),
                output_path=str(output_folder / f"{prefix}.feature.parquet"),
                qvalue_threshold=qvalue_threshold,
            )
            self._ontology_entries.extend(
                score_ontology_entries(adapter._discovered_scores, view="feature")
            )
        logger.info("DIA-NN feature conversion complete")

    def convert_pg(
        self,
        pg_matrix_path,
        output_folder=".",
        output_prefix=None,
        batch_size=100,
        standardized_intensities=False,
    ):
        output_folder = Path(output_folder)
        prefix = output_prefix or "diann"
        with DiannPgAdapter(
            duckdb_memory=self._memory, duckdb_threads=self._threads
        ) as adapter:
            adapter.convert(
                diann_report=self.report_path,
                pg_matrix_path=str(pg_matrix_path),
                sdrf_path=self.sdrf_path,
                output_path=str(output_folder / f"{prefix}.pg.parquet"),
            )
            self._ontology_entries.extend(
                score_ontology_entries(adapter._discovered_scores, view="pg")
            )
        logger.info("DIA-NN PG conversion complete")

    def write_ontology(self, output_folder: str | Path, prefix: str = "diann") -> None:
        """Write combined ontology.parquet with all accumulated entries."""
        # Add field-level CV term entries
        self._ontology_entries.extend(field_ontology_entries(view="feature"))

        if not self._ontology_entries:
            return
        from qpx.writers.ontology import OntologyWriter

        output_folder = Path(output_folder)
        onto_path = output_folder / f"{prefix}.ontology.parquet"
        with OntologyWriter(onto_path, creator="qpx") as writer:
            writer.write_batch(self._ontology_entries)
        logger.info(
            "Wrote %d ontology entries to %s",
            len(self._ontology_entries),
            onto_path,
        )

    def write_provenance(
        self, output_folder: str | Path, prefix: str = "diann"
    ) -> None:
        """Write provenance.parquet with DIA-NN + QPX conversion steps."""
        from qpx.writers.provenance import ProvenanceWriter

        output_folder = Path(output_folder)
        records = [
            {
                "step_order": 1,
                "step_category": "quantification",
                "step_name": "precursor_quantification",
                "tool_name": "DIA-NN",
                "tool_version": None,
                "tool_uri": None,
                "parameters": None,
                "config": None,
                "output_views": ["feature", "pg"],
            },
            {
                "step_order": 2,
                "step_category": "format_conversion",
                "step_name": "diann_to_qpx",
                "tool_name": "qpx",
                "tool_version": __version__,
                "tool_uri": None,
                "parameters": [
                    {"key": "report_path", "value": Path(self.report_path).name},
                ],
                "config": None,
                "output_views": ["feature", "pg", "sample", "run", "ontology"],
            },
        ]
        prov_path = output_folder / f"{prefix}.provenance.parquet"
        with ProvenanceWriter(prov_path, creator="qpx") as writer:
            writer.write_batch(records)
        logger.info("Wrote %d provenance steps to %s", len(records), prov_path)

    def write_dataset(
        self,
        output_folder: str | Path,
        prefix: str = "diann",
        project_accession: str | None = None,
    ) -> None:
        """Write dataset.parquet with project-level metadata."""
        from qpx.writers.dataset import DatasetWriter

        output_folder = Path(output_folder)
        record = {
            "project_accession": project_accession or "unknown",
            "project_title": None,
            "project_description": None,
            "pubmed_id": None,
            "software_name": "DIA-NN",
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
