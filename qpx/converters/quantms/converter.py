"""QuantMS orchestrator — composes PSM, Feature, and PG adapters."""

import logging
from datetime import datetime
from pathlib import Path

import duckdb

from qpx._version import __version__
from qpx.core.scores import (
    score_ontology_entries,
    field_ontology_entries,
    modification_ontology_entries,
)
from qpx.converters.sdrf import SdrfConverter
from qpx.converters.mztab import (
    load_mztab_sections,
    load_msstats,
    extract_modifications,
)
from qpx.converters.quantms.psm_adapter import QuantmsPsmAdapter
from qpx.converters.quantms.feature_adapter import QuantmsFeatureAdapter
from qpx.converters.quantms.pg_adapter import QuantmsPgAdapter

logger = logging.getLogger(__name__)


class QuantMSConverter:
    """Orchestrate full QuantMS mzTab conversion to QPX format.

    Delegates to individual adapters based on the requested structures.
    """

    def __init__(
        self,
        mztab_path,
        sdrf_file,
        msstats_file=None,
        database_path=None,
    ):
        self.mztab_path = str(mztab_path)
        self.sdrf_file = str(sdrf_file)
        self.msstats_file = str(msstats_file) if msstats_file else None
        self.database_path = str(database_path) if database_path else None

    def convert(
        self,
        output_folder,
        output_prefix="quantms",
        structures=None,
        compute_ibaq=True,
        project_accession=None,
    ):
        output_folder = Path(output_folder)
        output_folder.mkdir(parents=True, exist_ok=True)

        if structures is None:
            structures = ["psm", "feature", "pg"]

        ontology_entries: list[dict] = []

        # Always convert SDRF for sample/run metadata
        with SdrfConverter() as sdrf_conv:
            sdrf_conv.convert(
                sdrf_path=self.sdrf_file,
                sample_output=str(output_folder / f"{output_prefix}.sample.parquet"),
                run_output=str(output_folder / f"{output_prefix}.run.parquet"),
            )
            ontology_entries.extend(sdrf_conv.run_ontology_entries())
            logger.info("SDRF conversion complete")

        # Shared DuckDB connection — parse mzTab and MSstats once
        needs_mztab = any(s in structures for s in ("psm", "feature", "pg"))
        if not needs_mztab:
            self._write_dataset(
                output_folder, output_prefix, project_accession, []
            )
            return

        shared_conn = duckdb.connect(":memory:")
        shared_conn.execute("SET memory_limit='16GB'")
        shared_conn.execute("SET threads=4")
        provenance_records: list[dict] = []
        try:
            load_mztab_sections(shared_conn, self.mztab_path)
            if self.msstats_file:
                load_msstats(shared_conn, self.msstats_file)

            if "psm" in structures:
                with QuantmsPsmAdapter(conn=shared_conn) as adapter:
                    adapter.convert(
                        mztab_path=self.mztab_path,
                        output_path=str(output_folder / f"{output_prefix}.psm.parquet"),
                    )
                    ontology_entries.extend(
                        score_ontology_entries(adapter._discovered_scores, view="psm")
                    )
                    logger.info("PSM conversion complete")

            if "feature" in structures and self.msstats_file:
                with QuantmsFeatureAdapter(conn=shared_conn) as adapter:
                    adapter.convert(
                        mztab_path=self.mztab_path,
                        msstats_path=self.msstats_file,
                        output_path=str(output_folder / f"{output_prefix}.feature.parquet"),
                    )
                    ontology_entries.extend(
                        score_ontology_entries(adapter._discovered_scores, view="feature")
                    )
                    logger.info("Feature conversion complete")

            if "pg" in structures and self.msstats_file:
                with QuantmsPgAdapter(conn=shared_conn) as adapter:
                    adapter.convert(
                        mztab_path=self.mztab_path,
                        msstats_path=self.msstats_file,
                        output_path=str(output_folder / f"{output_prefix}.pg.parquet"),
                        compute_ibaq=compute_ibaq,
                    )
                    ontology_entries.extend(
                        score_ontology_entries(adapter._discovered_scores, view="pg")
                    )
                    logger.info("PG conversion complete")

            # PTM ontology entries from mzTab modification metadata
            modifications_meta = extract_modifications(shared_conn)
            ontology_entries.extend(
                modification_ontology_entries(modifications_meta, view="psm")
            )

            # Field-level CV term entries
            ontology_entries.extend(field_ontology_entries(view="psm"))

            # Build provenance from mzTab metadata
            provenance_records = self._build_provenance(shared_conn, self.mztab_path)
        finally:
            shared_conn.close()

        # Write combined ontology.parquet
        if ontology_entries:
            from qpx.writers.ontology import OntologyWriter

            onto_path = output_folder / f"{output_prefix}.ontology.parquet"
            with OntologyWriter(onto_path, creator="qpx") as writer:
                writer.write_batch(ontology_entries)
            logger.info(
                "Wrote %d ontology entries to %s", len(ontology_entries), onto_path
            )

        # Write provenance.parquet
        if provenance_records:
            from qpx.writers.provenance import ProvenanceWriter

            prov_path = output_folder / f"{output_prefix}.provenance.parquet"
            with ProvenanceWriter(prov_path, creator="qpx") as writer:
                writer.write_batch(provenance_records)
            logger.info(
                "Wrote %d provenance steps to %s",
                len(provenance_records),
                prov_path,
            )

        # Write dataset.parquet
        self._write_dataset(
            output_folder, output_prefix, project_accession, provenance_records
        )

    # ------------------------------------------------------------------
    # Provenance extraction from mzTab metadata
    # ------------------------------------------------------------------

    @staticmethod
    def _build_provenance(
        conn: duckdb.DuckDBPyConnection,
        mztab_path: str | None = None,
    ) -> list[dict]:
        """Extract processing provenance from mzTab metadata.

        Reads ``software[N]`` entries from the metadata table to identify
        which search engines and tools were used.
        """
        import re

        provenance_records: list[dict] = []
        step_order = 0

        try:
            rows = conn.execute(
                """
                SELECT key, value FROM metadata
                WHERE key LIKE 'software[%'
                  AND key NOT LIKE '%-%'
                ORDER BY key
                """
            ).fetchall()
        except Exception:
            rows = []

        for key, value in rows:
            step_order += 1
            # Parse CV term: [CV, accession, name, value]
            parts = value.replace("[", "").replace("]", "").split(",")
            tool_name = parts[2].strip() if len(parts) >= 3 else value.strip()

            # Try to extract version from the same metadata entry
            version = None
            if len(parts) >= 4 and parts[3].strip():
                version = parts[3].strip()

            # Determine step category from software index
            m = re.search(r"\[(\d+)\]", key)
            idx = int(m.group(1)) if m else step_order

            provenance_records.append({
                "step_order": idx,
                "step_category": "database_search",
                "step_name": "spectrum_identification",
                "tool_name": tool_name,
                "tool_version": version,
                "tool_uri": None,
                "parameters": None,
                "config": None,
                "output_views": ["psm", "feature", "pg"],
            })

        # Add QPX conversion step
        step_order = max((r["step_order"] for r in provenance_records), default=0) + 1
        params = None
        if mztab_path:
            params = [{"key": "mztab_path", "value": Path(mztab_path).name}]
        provenance_records.append({
            "step_order": step_order,
            "step_category": "format_conversion",
            "step_name": "mztab_to_qpx",
            "tool_name": "qpx",
            "tool_version": __version__,
            "tool_uri": None,
            "parameters": params,
            "config": None,
            "output_views": ["psm", "feature", "pg", "sample", "run", "ontology"],
        })

        return provenance_records

    # ------------------------------------------------------------------
    # Dataset metadata
    # ------------------------------------------------------------------

    @staticmethod
    def _write_dataset(
        output_folder: Path,
        output_prefix: str,
        project_accession: str | None,
        provenance_records: list[dict],
    ) -> None:
        """Write dataset.parquet with project-level metadata."""
        from qpx.writers.dataset import DatasetWriter

        sw_name = None
        sw_version = None
        if provenance_records:
            sw_name = provenance_records[0].get("tool_name")
            sw_version = provenance_records[0].get("tool_version")

        dataset_record = {
            "project_accession": project_accession or "unknown",
            "project_title": None,
            "project_description": None,
            "pubmed_id": None,
            "software_name": sw_name or "OpenMS/quantms",
            "software_version": sw_version or "quantms",
            "creation_date": datetime.now().isoformat(),
            "file_checksums": None,
            "file_row_counts": None,
            "file_sizes_bytes": None,
            "total_structures": None,
            "packaged_at": None,
        }
        ds_path = output_folder / f"{output_prefix}.dataset.parquet"
        with DatasetWriter(ds_path, creator="qpx") as writer:
            writer.write_batch([dataset_record])
        logger.info("Wrote dataset metadata to %s", ds_path)
