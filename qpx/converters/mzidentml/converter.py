"""MzIdentML full-dataset converter.

Orchestrates production of a complete QPX dataset directory from an
mzIdentML file and optional MGF spectra file.  Composes the existing
:class:`MzIdentMLPsmAdapter` for PSM extraction and adds pepmap,
provenance, ontology, and dataset metadata structures.
"""

from __future__ import annotations

import gzip
import logging
from datetime import datetime
from pathlib import Path

try:
    from lxml import etree
except ImportError:
    etree = None  # type: ignore[assignment]

from qpx._version import __version__
from qpx.converters.mzidentml.psm_adapter import MzIdentMLPsmAdapter
from qpx.core.scores import score_ontology_entries, field_ontology_entries
from qpx.writers.dataset import DatasetWriter
from qpx.writers.ontology import OntologyWriter
from qpx.writers.pepmap import PepMapWriter
from qpx.writers.provenance import ProvenanceWriter
from qpx.writers.psm import PsmWriter

logger = logging.getLogger(__name__)


class MzIdentMLConverter:
    """Convert mzIdentML (+ optional MGF) to a complete QPX dataset."""

    def __init__(self, compression: str = "zstd"):
        self._compression = compression

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def convert(
        self,
        mzid_path: str | Path,
        output_folder: str | Path,
        output_prefix: str = "mzidentml",
        mgf_path: str | Path | None = None,
        include_spectra: bool = False,
        project_accession: str | None = None,
    ) -> Path:
        """Convert mzIdentML to QPX dataset.

        Args:
            mzid_path: Path to ``.mzid`` or ``.mzid.gz`` file.
            output_folder: Directory for output files.
            output_prefix: File prefix (e.g. ``"PXD054720"``).
            mgf_path: Optional MGF file for spectra attachment.
            include_spectra: If ``True`` and *mgf_path* given, attach
                spectra arrays to PSM records.
            project_accession: PRIDE / ProteomeXchange accession.

        Returns:
            Path to *output_folder*.
        """
        if etree is None:
            raise ImportError(
                "lxml is required for mzIdentML conversion. "
                "Install it with: pip install qpx[mzidentml]"
            )
        mzid_path = Path(mzid_path)
        output_folder = Path(output_folder)
        output_folder.mkdir(parents=True, exist_ok=True)

        # 1. Parse mzIdentML and build PSM records via the existing adapter
        adapter = MzIdentMLPsmAdapter()
        try:
            parsed = adapter._parse_mzidentml(str(mzid_path))
            records = adapter._build_psm_records(parsed)
        finally:
            adapter.close()

        if not records:
            logger.warning("No PSM records produced from %s", mzid_path)
            return output_folder

        # 2. Optionally enrich PSMs with spectra from MGF
        if include_spectra and mgf_path is not None:
            records = self._attach_spectra(records, mgf_path)

        # 3. Track scores for ontology generation
        adapter._track_scores(records)
        discovered_scores = adapter.get_discovered_scores()

        # 4. Write PSM parquet
        psm_path = output_folder / f"{output_prefix}.psm.parquet"
        with PsmWriter(psm_path, creator="mzidentml", compression=self._compression) as writer:
            writer.write_batch(records)
        logger.info("Wrote %d PSMs to %s", len(records), psm_path)

        # 5. Write pepmap parquet
        pepmap_path = output_folder / f"{output_prefix}.pepmap.parquet"
        pepmap_records = self._build_pepmap(records)
        if pepmap_records:
            with PepMapWriter(pepmap_path, creator="mzidentml", compression=self._compression) as writer:
                writer.write_batch(pepmap_records)
            logger.info(
                "Wrote %d pepmap entries to %s", len(pepmap_records), pepmap_path
            )

        # 6. Write provenance parquet
        provenance_path = output_folder / f"{output_prefix}.provenance.parquet"
        provenance_records = self._build_provenance(mzid_path)
        if provenance_records:
            with ProvenanceWriter(provenance_path, creator="mzidentml", compression=self._compression) as writer:
                writer.write_batch(provenance_records)
            logger.info(
                "Wrote %d provenance steps to %s",
                len(provenance_records),
                provenance_path,
            )

        # 7. Write ontology parquet
        ontology_path = output_folder / f"{output_prefix}.ontology.parquet"
        ontology_entries = score_ontology_entries(discovered_scores, view="psm")
        ontology_entries.extend(field_ontology_entries(view="psm"))
        if ontology_entries:
            with OntologyWriter(ontology_path, creator="mzidentml", compression=self._compression) as writer:
                writer.write_batch(ontology_entries)
            logger.info(
                "Wrote %d ontology entries to %s", len(ontology_entries), ontology_path
            )

        # 8. Write dataset metadata parquet
        dataset_path = output_folder / f"{output_prefix}.dataset.parquet"
        # Extract software info from provenance for metadata
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
            "software_name": sw_name or "qpx",
            "software_version": sw_version or __version__,
            "creation_date": datetime.now().isoformat(),
            "file_checksums": None,
            "file_row_counts": None,
            "file_sizes_bytes": None,
            "total_structures": None,
            "packaged_at": None,
        }
        with DatasetWriter(dataset_path, creator="mzidentml", compression=self._compression) as writer:
            writer.write_batch([dataset_record])
        logger.info("Wrote dataset metadata to %s", dataset_path)

        return output_folder

    # ------------------------------------------------------------------
    # Spectra attachment
    # ------------------------------------------------------------------

    @staticmethod
    def _attach_spectra(records: list[dict], mgf_path: str | Path) -> list[dict]:
        """Attach mz_array and intensity_array from MGF to matching PSM records."""
        from qpx.converters.mzidentml.mgf_parser import MgfSpectraIndex

        mgf_index = MgfSpectraIndex(mgf_path)
        matched_by_scan = 0
        matched_by_index = 0
        for rec in records:
            scans = rec.get("scan") or []
            if not scans:
                continue
            scan = scans[0]
            spectrum = mgf_index.get_spectrum(scan)
            if spectrum is None:
                # Fallback: try 0-based index lookup
                spectrum = mgf_index.get_spectrum_by_index(scan)
                if spectrum is not None:
                    matched_by_index += 1
            else:
                matched_by_scan += 1

            if spectrum is not None:
                rec["mz_array"] = spectrum["mz_array"]
                rec["intensity_array"] = spectrum["intensity_array"]
                # Fill RT from MGF only if not already set from mzIdentML
                if rec.get("rt") is None and spectrum.get("rt") is not None:
                    rec["rt"] = spectrum["rt"]

        total_matched = matched_by_scan + matched_by_index
        logger.info(
            "Attached spectra to %d / %d PSMs from MGF index (%d spectra, %d by scan, %d by index)",
            total_matched,
            len(records),
            len(mgf_index),
            matched_by_scan,
            matched_by_index,
        )
        return records

    # ------------------------------------------------------------------
    # Pepmap construction
    # ------------------------------------------------------------------

    @staticmethod
    def _build_pepmap(records: list[dict]) -> list[dict]:
        """Derive deduplicated peptide-protein map from PSM records."""
        seen: set[tuple[str, str]] = set()
        pepmap_records: list[dict] = []
        seq_proteins: dict[str, set[str]] = {}

        for psm in records:
            seq = psm.get("sequence")
            if not seq:
                continue
            proteins = psm.get("protein_accessions") or []
            for prot in proteins:
                key = (seq, prot)
                if key not in seen:
                    seen.add(key)
                    pepmap_records.append(
                        {
                            "sequence": seq,
                            "protein_accession": prot,
                            "protein_name": None,
                            "gene_name": None,
                            "start_position": None,
                            "end_position": None,
                            "is_unique": None,  # computed below
                            "is_proteotypic": None,
                        }
                    )
                # Track all proteins per sequence for is_unique
                seq_proteins.setdefault(seq, set()).add(prot)

        # Compute is_unique: True if the sequence maps to exactly 1 protein
        for rec in pepmap_records:
            n_prots = len(seq_proteins.get(rec["sequence"], set()))
            rec["is_unique"] = n_prots == 1

        return pepmap_records

    # ------------------------------------------------------------------
    # Provenance extraction
    # ------------------------------------------------------------------

    @staticmethod
    def _build_provenance(mzid_path: str | Path) -> list[dict]:
        """Extract processing provenance from mzIdentML AnalysisSoftware elements."""
        mzid_path = Path(mzid_path)

        if str(mzid_path).endswith(".gz"):
            with gzip.open(mzid_path, "rb") as fh:
                tree = etree.parse(fh)
        else:
            tree = etree.parse(str(mzid_path))

        root = tree.getroot()
        ns = root.nsmap.get(None, "")

        provenance_records: list[dict] = []
        step_order = 0

        for sw in root.iter(f"{{{ns}}}AnalysisSoftware"):
            step_order += 1

            name_attr = sw.get("name", "")
            version = sw.get("version")

            # Try to get the proper software name from SoftwareName/cvParam
            tool_name = name_attr
            software_name_el = sw.find(f"{{{ns}}}SoftwareName")
            if software_name_el is not None:
                cv = software_name_el.find(f"{{{ns}}}cvParam")
                if cv is not None:
                    tool_name = cv.get("name", name_attr)

            provenance_records.append(
                {
                    "step_order": step_order,
                    "step_category": "database_search",
                    "step_name": "spectrum_identification",
                    "tool_name": tool_name,
                    "tool_version": version,
                    "tool_uri": sw.get("uri"),
                    "parameters": None,
                    "config": None,
                    "output_views": ["psm"],
                }
            )

        return provenance_records
