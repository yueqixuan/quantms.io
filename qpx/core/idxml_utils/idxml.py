import logging
from typing import Optional, List, Dict, Any, Iterator
from pathlib import Path

import pandas as pd
import pyarrow as pa
import pyopenms as oms

from qpx.core.common import PSM_SCHEMA
from qpx.core.openms import OpenMSHandler
from qpx.utils.file_utils import ParquetBatchWriter
from qpx.utils.pride_utils import (
    generate_scan_number,
    standardize_protein_string_accession,
)
from qpx.utils.mztab_utils import parse_peptidoform_openms


class IdXmlPsm:
    """IdXML PSM (Peptide-Spectrum Match) processor.

    This class processes PSM data from idXML files using OpenMS functionality,
    providing PSM-specific functionality for converting idXML to parquet format.
    """

    def __init__(self, idxml_path: str, spectral_data: bool = False):
        """Initialize IdXML PSM processor.

        Args:
            idxml_path: Path to idXML file
        """
        self.idxml_path = Path(idxml_path)
        self.logger = logging.getLogger(f"{__name__}.{self.__class__.__name__}")

        self._spectral_data = spectral_data

        prot_ids = []
        pep_ids = []
        oms.IdXMLFile().load(idxml_path, prot_ids, pep_ids)
        self.prot_ids = prot_ids
        self.pep_ids = pep_ids
        self.openms_handler = OpenMSHandler()
        self.logger.info(
            f"Loaded {len(self.prot_ids)} proteins and {len(self.pep_ids)} peptides from {idxml_path}"
        )

        self.search_metadata = self._extract_search_metadata()

        self.psm_records = []
        for batch in self.iter_idxml_psms():
            self.psm_records.extend(batch)
        self.logger.info(f"Loaded {len(self.psm_records)} PSM records")

    def _extract_search_metadata(self) -> Dict[str, Any]:
        """Extract search metadata from protein identifications."""
        metadata = {}

        if not self.prot_ids:
            return metadata

        prot_id = self.prot_ids[0]

        # Extract basic search engine information
        metadata.update(self._extract_search_engine_info(prot_id))

        # Extract search parameters
        search_params = self._get_search_parameters(prot_id)
        if search_params:
            metadata.update(self._extract_search_parameters(search_params))

        # Ensure search_engine is always present
        if "search_engine" not in metadata:
            metadata["search_engine"] = "unknown"

        self.logger.info(f"Extracted search metadata: {metadata}")
        return metadata

    def _extract_search_engine_info(self, prot_id) -> Dict[str, Any]:
        """Extract search engine name and version."""
        metadata = {}

        # Extract search engine name
        try:
            search_engine = prot_id.getSearchEngine()
            if search_engine:
                metadata["search_engine"] = search_engine
        except (AttributeError, RuntimeError) as e:
            self.logger.debug(f"Could not extract search engine: {e}")
        except Exception as e:
            self.logger.warning(f"Unexpected error extracting search engine: {e}")

        # Extract search engine version
        try:
            version = prot_id.getSearchEngineVersion()
            if version:
                metadata["search_engine_version"] = version
        except (AttributeError, RuntimeError) as e:
            self.logger.debug(f"Could not extract search engine version: {e}")
        except Exception as e:
            self.logger.warning(
                f"Unexpected error extracting search engine version: {e}"
            )

        return metadata

    def _get_search_parameters(self, prot_id):
        """Safely get search parameters object."""
        try:
            return prot_id.getSearchParameters()
        except (AttributeError, RuntimeError) as e:
            self.logger.debug(f"Could not access search parameters: {e}")
        except Exception as e:
            self.logger.warning(f"Unexpected error accessing search parameters: {e}")
        return None

    def _extract_search_parameters(self, search_params) -> Dict[str, Any]:
        """Extract detailed search parameters."""
        metadata = {}

        # Extract enzyme information
        try:
            enzyme = search_params.getEnzyme()
            if enzyme:
                metadata["enzyme"] = enzyme.getName()
        except (AttributeError, RuntimeError) as e:
            self.logger.debug(f"Could not extract enzyme: {e}")
        except Exception as e:
            self.logger.warning(f"Unexpected error extracting enzyme: {e}")

        # Extract missed cleavages
        try:
            missed_cleavages = search_params.getMissedCleavages()
            metadata["missed_cleavages"] = str(missed_cleavages)
        except (AttributeError, RuntimeError) as e:
            self.logger.debug(f"Could not extract missed cleavages: {e}")
        except Exception as e:
            self.logger.warning(f"Unexpected error extracting missed cleavages: {e}")

        # Extract precursor tolerance
        try:
            precursor_tol = search_params.getPrecursorMassTolerance()
            metadata["precursor_tolerance"] = f"{precursor_tol:.6f}"
        except (AttributeError, RuntimeError) as e:
            self.logger.debug(f"Could not extract precursor tolerance: {e}")
        except Exception as e:
            self.logger.warning(f"Unexpected error extracting precursor tolerance: {e}")

        # Extract fragment tolerance
        try:
            fragment_tol = search_params.getFragmentMassTolerance()
            metadata["peak_mass_tolerance"] = f"{fragment_tol:.6f}"
        except (AttributeError, RuntimeError) as e:
            self.logger.debug(f"Could not extract fragment tolerance: {e}")
        except Exception as e:
            self.logger.warning(f"Unexpected error extracting fragment tolerance: {e}")

        return metadata

    def iter_idxml_psms(
        self, batch_size: int = 1000, spectral_data: bool = False
    ) -> Iterator[List[Dict[str, Any]]]:
        """
        Memory-efficient iterator for reading PSMs from idXML file.

        Args:
            idxml_path: Path to idXML file
            batch_size: Number of PSMs to yield in each batch

        Yields:
            List of PSM dictionaries, with each batch containing up to batch_size PSMs

        Raises:
            RuntimeError: If file cannot be loaded or contains no data
        """

        current_batch = []

        # Process peptide identifications in batches
        for peptide_id in self.pep_ids:
            # Get spectrum information
            rt = peptide_id.getRT() if peptide_id.hasRT() else 0.0
            mz = peptide_id.getMZ() if peptide_id.hasMZ() else 0.0

            # Get spectrum reference from metadata if available
            spectrum_ref = (
                peptide_id.getMetaValue("spectrum_reference")
                if peptide_id.metaValueExists("spectrum_reference")
                else None
            )

            # Get run information
            run = self.openms_handler._get_run(self.prot_ids, peptide_id)

            # Get hits for this peptide identification
            hits = peptide_id.getHits()

            # Get peptide identification metadata
            for hit in hits:
                # Get sequence and peptidoform
                unmodified_sequence = hit.getSequence().toUnmodifiedString()
                peptidoform = parse_peptidoform_openms(hit.getSequence().toString())

                # Calculate theoretical m/z
                try:
                    calc_mz = (
                        hit.getSequence().getMonoWeight(
                            oms.Residue.ResidueType.Full, hit.getCharge()
                        )
                        / hit.getCharge()
                    )
                except (AttributeError, RuntimeError, ZeroDivisionError) as e:
                    self.logger.debug(
                        f"Could not calculate theoretical m/z for sequence {unmodified_sequence}: {e}"
                    )
                    calc_mz = 0.0
                except Exception as e:
                    self.logger.warning(
                        f"Unexpected error calculating m/z for sequence {unmodified_sequence}: {e}"
                    )
                    calc_mz = 0.0

                # Get ion mobility if available
                ion_mobility = (
                    float(im)
                    if (im := peptide_id.getMetaValue("IM")) is not None
                    else None
                )

                # Extract protein accessions
                protein_accessions = [
                    acc.decode() for acc in hit.extractProteinAccessionsSet()
                ]

                # Check if PSM is decoy
                is_decoy = self.openms_handler._is_decoy(hit)

                # Get additional scores and metadata
                additional_scores = self.openms_handler._extract_additional_scores(hit)
                modifications = (
                    self.openms_handler._extract_modifications_from_sequence(
                        hit.getSequence()
                    )
                )

                psm_record = {
                    "sequence": unmodified_sequence,
                    "peptidoform": peptidoform,
                    "charge": hit.getCharge(),
                    "rt": rt,
                    "mz": mz,
                    "spectrum_ref": spectrum_ref
                    or f"scan={self.pep_ids.index(peptide_id)}",
                    "scan": run,
                    "modifications": modifications,
                    "score": hit.getScore(),
                    "rank": hit.getRank() + 1,  # Convert to 1-based
                    "calculated_mz": calc_mz,
                    "protein_accessions": protein_accessions,
                    "is_decoy": is_decoy,
                    "additional_scores": additional_scores,
                    "cv_params": [],
                    "ion_mobility": ion_mobility,
                }

                current_batch.append(psm_record)

                # Yield batch when it reaches batch_size
                if len(current_batch) >= batch_size:
                    yield current_batch
                    current_batch = []

        # Yield any remaining PSMs
        if current_batch:
            yield current_batch

    def iter_psm_table(self, chunksize: int = 1000000) -> Iterator[pa.Table]:
        """Iterate over PSM table in chunks.

        Args:
            chunksize: Number of rows to process in each chunk

        Yields:
            PyArrow Table containing PSM data with proper struct types
        """
        for psm_batch in self.iter_idxml_psms(chunksize):
            # Transform each record to PyArrow-compatible format
            transformed_records = []
            for record in psm_batch:
                transformed_record = self._transform_psm_record_to_arrow(record)
                if transformed_record:
                    transformed_records.append(transformed_record)

            # Convert to PyArrow Table with proper schema
            if transformed_records:
                table = self._create_arrow_table(transformed_records)
                yield table
            else:
                # Return empty table with correct schema
                yield pa.table([], schema=PSM_SCHEMA)

    def _transform_psm_record_to_arrow(
        self, record: Dict[str, Any]
    ) -> Optional[Dict[str, Any]]:
        """Transform a single PSM record to PyArrow-compatible format.

        Args:
            record: Dictionary containing PSM data

        Returns:
            Transformed record dictionary with PyArrow structs or None if invalid
        """
        try:
            # Create a new record with PSM schema columns
            transformed_record = {}

            # Basic fields
            transformed_record["sequence"] = record.get("sequence", "")
            transformed_record["peptidoform"] = record.get("peptidoform", "")
            transformed_record["precursor_charge"] = int(record.get("charge", 0))
            transformed_record["rt"] = (
                float(record.get("rt", 0.0)) if record.get("rt") else None
            )
            transformed_record["observed_mz"] = (
                float(record.get("mz", 0.0)) if record.get("mz") else None
            )
            transformed_record["calculated_mz"] = (
                float(record.get("calculated_mz", 0.0))
                if record.get("calculated_mz")
                else None
            )

            # Protein accessions
            protein_accessions = record.get("protein_accessions", [])
            if protein_accessions:
                # Standardize protein accessions
                standardized_accessions = []
                for acc in protein_accessions:
                    standardized_acc = standardize_protein_string_accession(
                        acc, is_sorted=True
                    )
                    standardized_accessions.append(standardized_acc)
                transformed_record["protein_accessions"] = standardized_accessions
            else:
                transformed_record["protein_accessions"] = []

            # Modifications
            modifications = record.get("modifications", [])
            transformed_record["modifications"] = modifications if modifications else []

            # Scores
            additional_scores = record.get("additional_scores", [])

            # Add search engine score if available
            if record.get("score") is not None:
                search_engine = self.search_metadata.get("search_engine", "unknown")
                additional_scores.append(
                    {
                        "score_name": f"{search_engine}_score",
                        "score_value": float(record["score"]),
                    }
                )

            transformed_record["additional_scores"] = (
                additional_scores if additional_scores else []
            )

            # Decoy flag
            transformed_record["is_decoy"] = int(record.get("is_decoy", 0))

            # Generate scan number from spectrum reference
            spectrum_ref = record.get("spectrum_ref", "")
            if spectrum_ref:
                transformed_record["scan"] = generate_scan_number(spectrum_ref)
                transformed_record["reference_file_name"] = (
                    self._extract_reference_file_name(spectrum_ref)
                )
            else:
                transformed_record["scan"] = "unknown"
                transformed_record["reference_file_name"] = "unknown"

            # CV parameters
            cv_params = record.get("cv_params", [])

            # Add search metadata as CV parameters
            if self.search_metadata:
                # Add basic search parameters
                for key in [
                    "enzyme",
                    "missed_cleavages",
                    "precursor_tolerance",
                    "peak_mass_tolerance",
                ]:
                    if key in self.search_metadata:
                        cv_params.append(
                            {
                                "cv_name": f"search_{key}",
                                "cv_value": str(self.search_metadata[key]),
                            }
                        )

                # Add search engine info
                if "search_engine" in self.search_metadata:
                    cv_params.append(
                        {
                            "cv_name": "search_engine",
                            "cv_value": self.search_metadata["search_engine"],
                        }
                    )
                    if "search_engine_version" in self.search_metadata:
                        cv_params.append(
                            {
                                "cv_name": "search_engine_version",
                                "cv_value": self.search_metadata[
                                    "search_engine_version"
                                ],
                            }
                        )

            transformed_record["cv_params"] = cv_params

            # Set defaults for optional fields
            transformed_record["posterior_error_probability"] = None
            transformed_record["predicted_rt"] = None

            if self._spectral_data:
                transformed_record["ion_mobility"] = record.get("ion_mobility")
            else:
                transformed_record["ion_mobility"] = None

            transformed_record["number_peaks"] = None
            transformed_record["mz_array"] = None
            transformed_record["intensity_array"] = None
            transformed_record["charge_array"] = None
            transformed_record["ion_type_array"] = None
            transformed_record["ion_mobility_array"] = None

            return transformed_record

        except Exception as e:
            self.logger.error(f"=== DEBUG: Error transforming PSM record: {e} ===")
            self.logger.error(f"=== DEBUG: Record keys: {list(record.keys())} ===")
            self.logger.error(
                f"=== DEBUG: Record sample: {dict(list(record.items())[:5])} ==="
            )
            import traceback

            self.logger.error(
                f"=== DEBUG: Transformation traceback: {traceback.format_exc()} ==="
            )
            return None

    def _extract_reference_file_name(self, spectrum_ref: str) -> str:
        """Extract reference file name from spectrum reference.

        Args:
            spectrum_ref: Spectrum reference string

        Returns:
            Reference file name
        """
        if not spectrum_ref:
            return "unknown"

        # Try to extract file name from spectrum reference
        # This is a simplified approach - may need adjustment based on actual idXML format
        if ":" in spectrum_ref:
            parts = spectrum_ref.split(":")
            if len(parts) > 1:
                return parts[0]

        return spectrum_ref

    def _create_arrow_table(self, records: List[Dict[str, Any]]) -> pa.Table:
        """Create PyArrow Table from transformed records.

        Args:
            records: List of transformed record dictionaries

        Returns:
            PyArrow Table with proper schema
        """
        # Convert to DataFrame first for easier column handling
        df = pd.DataFrame(records)

        # Convert to PyArrow Table with schema
        table = pa.Table.from_pandas(df, schema=PSM_SCHEMA)
        return table

    def _create_file_metadata(self):
        """Create file metadata structure according to PSM avsc schema"""
        import uuid
        from datetime import datetime
        from qpx import __version__

        metadata = {
            "qpx_version": __version__,
            "creator": "QPX",
            "file_type": "psm_file",
            "creation_date": datetime.now().isoformat(),
            "uuid": str(uuid.uuid4()),
            "scan_format": "scan",  # Default scan format
            "software_provider": "QPX",
        }

        # Add search engine info if available
        if self.search_metadata:
            if "search_engine" in self.search_metadata:
                metadata["creator"] = self.search_metadata["search_engine"]
            if "search_engine_version" in self.search_metadata:
                metadata["software_provider"] = (
                    f"{self.search_metadata['search_engine']} {self.search_metadata['search_engine_version']}"
                )

        return metadata

    def convert_to_parquet(
        self,
        output_path: str,
        chunksize: int = 1000000,
        protein_file: Optional[str] = None,
    ) -> None:
        """Write PSM data to Parquet file using efficient batch writing.

        Args:
            output_path: Path to output Parquet file
            chunksize: Number of rows to process in each chunk
            protein_file: Optional protein filter file (not implemented for idXML)
        """
        logger = logging.getLogger("qpx.core.idxml.idxml")

        # Log input and output paths
        logger.info(f"=== DEBUG: convert_to_parquet started ===")
        logger.info(f"Output path: {output_path}")
        logger.info(f"Converting from idXML file: {self.idxml_path}")
        logger.info(f"Chunk size: {chunksize}")

        if protein_file:
            logger.warning("Protein filtering not implemented for idXML files")

        # Create file metadata for parquet file
        logger.info("=== DEBUG: Creating file metadata ===")
        file_metadata = self._create_file_metadata()
        logger.info(f"=== DEBUG: File metadata: {file_metadata} ===")

        # Initialize batch writer with file metadata
        logger.info("=== DEBUG: Initializing ParquetBatchWriter ===")
        batch_writer = ParquetBatchWriter(
            output_path, PSM_SCHEMA, file_metadata=file_metadata
        )
        logger.info("=== DEBUG: ParquetBatchWriter initialized ===")

        try:
            # Process data in chunks and write batches
            logger.info("=== DEBUG: Starting chunk processing ===")
            chunk_count = 0
            total_records_written = 0

            for table in self.iter_psm_table(chunksize=chunksize):
                chunk_count += 1
                transformed_records = table.to_pylist()
                logger.info(
                    f"=== DEBUG: Processing chunk {chunk_count}, records: {len(transformed_records)} ==="
                )

                # Write batch if we have records
                if transformed_records:
                    batch_writer.write_batch(transformed_records)
                    total_records_written += len(transformed_records)
                    logger.info(
                        f"=== DEBUG: Wrote {len(transformed_records)} records in chunk {chunk_count} ==="
                    )
                else:
                    logger.warning(f"=== DEBUG: No records in chunk {chunk_count} ===")

            logger.info(
                f"=== DEBUG: Processed {chunk_count} chunks, total records: {total_records_written} ==="
            )

        except Exception as e:
            logger.error(f"=== DEBUG: Error during parquet conversion: {e} ===")
            logger.error(f"=== DEBUG: Exception type: {type(e)} ===")
            import traceback

            logger.error(f"=== DEBUG: Traceback: {traceback.format_exc()} ===")
            raise
        finally:
            # Ensure final batch is written and writer is closed
            logger.info("=== DEBUG: Closing ParquetBatchWriter ===")
            batch_writer.close()
            logger.info("=== DEBUG: Parquet file closed successfully ===")
            logger.info("Parquet file closed successfully")

    def get_psm_count(self) -> int:
        """Get total number of PSMs.

        Returns:
            Total count of PSM entries
        """
        count = 0
        for batch in self.iter_idxml_psms():
            count += len(batch)
        return count
