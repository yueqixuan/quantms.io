import logging
from typing import Optional, List, Dict, Any, Tuple, Union
from pathlib import Path

import pandas as pd
import pyarrow as pa

from quantmsio.core.common import (
    OPENMS_IS_DECOY,
    OPENMS_PEPTIDOFORM_COLUMN,
    OPENMS_POSTERIOR_ERRORPROBABILITY,
    PSM_SCHEMA,
)
from quantmsio.core.openms import get_openms_score_name
from quantmsio.core.quantms.mztab import MzTabIndexer
from quantmsio.utils.file_utils import extract_protein_list, ParquetBatchWriter
from quantmsio.utils.pride_utils import (
    generate_scan_number,
    get_petidoform_msstats_notation,
    standardize_protein_string_accession,
)
from quantmsio.utils.mztab_utils import (
    extract_ms_runs_from_metadata,
    parse_pepidoform_with_modifications,
)


class Psm:
    """PSM (Peptide-Spectrum Match) processor using composition pattern.

    This class processes PSM data from an MzTabIndexer instance, providing
    PSM-specific functionality without inheriting from the indexer.
    """

    def __init__(self, mztab_indexer: MzTabIndexer, spectral_data: bool = False):
        """Initialize PSM processor with an MzTabIndexer instance.

        Args:
            mztab_indexer: An initialized MzTabIndexer instance
        """
        self._indexer = mztab_indexer
        self.logger = logging.getLogger(f"{__name__}.{self.__class__.__name__}")
        self._spectral_data = spectral_data

        # Initialize PSM-specific data from the indexer
        self._ms_runs = self._extract_ms_runs()

        # This dictionary could be big, then we need to optimize it
        self._protein_global_qvalue_map = self._get_protein_qvalue()
        self._modifications = self._get_metadata_modifications()
        self._score_names = self._get_score_names()

    def _extract_ms_runs(self) -> dict:
        """
        Extract MS runs from the mzTab metadata.

        This method uses the enhanced mztab_utils function to extract MS run
        information from the metadata DataFrame. It provides better error
        handling and more robust parsing than the previous implementation.

        Returns:
            Dictionary mapping MS run IDs to file paths
        """
        try:
            metadata_df = self._indexer.get_metadata()
            ms_runs = extract_ms_runs_from_metadata(metadata_df)

            self.logger.debug(f"Extracted {len(ms_runs)} MS runs from metadata")
            return ms_runs

        except Exception as e:
            self.logger.warning(f"Could not extract MS runs: {e}")
            return {}

    def _get_protein_qvalue(self) -> dict:
        """
        Get protein q-value mapping optimized for performance.

        Returns:
            Dictionary mapping protein accessions to q-values
        """
        try:
            proteins_df = self._indexer.get_protein_best_searchengine_score()

            # Early return if required columns don't exist
            if not (
                "accession" in proteins_df.columns
                and "opt_global_qvalue" in proteins_df.columns
            ):
                return {}

            # Check q-value status once and cache the result
            is_qvalue = self._indexer._is_protein_score_qvalue()

            # Transform the accession to standardize the protein string accession, sort True
            proteins_df["accession"] = proteins_df["accession"].apply(
                standardize_protein_string_accession, is_sorted=True
            )

            if is_qvalue:
                # Filter out rows with NaN values and create mapping efficiently
                valid_rows = proteins_df.dropna(
                    subset=["accession", "opt_global_qvalue"]
                )
                self.logger.info(
                    f"Valid rows, number of rows, removed rows: {len(proteins_df)}, {len(valid_rows)}, {len(proteins_df) - len(valid_rows)}"
                )
                return dict(
                    zip(valid_rows["accession"], valid_rows["opt_global_qvalue"])
                )
            else:
                # If not q-value, return empty dict (all values would be None anyway)
                return {}

        except Exception as e:
            self.logger.warning(f"Could not get protein q-value map: {e}")
            return {}

    def _get_score_names(self) -> dict:
        """
        Get score names from the indexer metadata using vectorized operations.

        This method extracts PSM, peptide, and protein search engine score names
        from the mzTab metadata. It uses pandas' string operations for efficient
        and readable parsing.

        Returns:
            A dictionary containing mappings of score indexes to canonical names
            for 'psms', 'peptides', and 'proteins'.
        """
        score_names = {"psms": {}, "peptides": {}, "proteins": {}}
        try:
            metadata_df = self._indexer.get_metadata()
            if (
                metadata_df.empty
                or "key" not in metadata_df.columns
                or "value" not in metadata_df.columns
            ):
                return score_names

            # Filter for rows containing search engine score information
            score_pattern = "search_engine_score"
            score_df = metadata_df[
                metadata_df["key"].str.contains(score_pattern, na=False)
            ].copy()

            if score_df.empty:
                return score_names

            # Extract score index and name term using vectorized operations
            score_df.loc[:, "mztab_score_name"] = (
                score_df["key"].str.extract(r"\[(\d+)\]").astype(int)
            )
            score_df.loc[:, "name_term"] = (
                score_df["value"].str.split(",").str[2].str.strip()
            )

            # Drop rows where extraction failed
            score_df.dropna(subset=["mztab_score_name", "name_term"], inplace=True)

            # Apply get_openms_score_name to get canonical names
            score_df.loc[:, "canonical_name"] = score_df["name_term"].apply(
                get_openms_score_name
            )

            # Populate the score_names dictionary
            for score_type in ["psm", "peptide", "protein"]:
                type_df = score_df[score_df["key"].str.startswith(f"{score_type}_")]
                if not type_df.empty:
                    score_names[f"{score_type}s"] = dict(
                        zip(type_df["mztab_score_name"], type_df["canonical_name"])
                    )

        except Exception as e:
            self.logger.warning(f"Could not extract score names: {e}")

        return score_names

    def _is_modification_row(self, row_key):
        """Check if a metadata row contains modification information."""
        return (
            "fixed_mod[" in row_key
            or "var_mod[" in row_key
            or "variable_mod[" in row_key
        )

    def _extract_modification_values(self, row_value):
        """Extract accession and name from row value."""
        values = row_value.replace("[", "").replace("]", "").split(",")
        if len(values) >= 3:
            return values[1].strip(), values[2].strip()
        return None, None

    def _parse_modification_row(self, row):
        """Parse a single metadata row to extract modification data."""
        if not self._is_modification_row(row["key"]):
            return None, None, None, None, None

        current_mod_index = row["key"].split("-")[0]
        site, position, accession, name = None, None, None, None

        if "site" in row["key"]:
            site = row["value"].strip()
        elif "position" in row["key"]:
            position = row["value"].strip()
        else:
            accession, name = self._extract_modification_values(row["value"])

        return current_mod_index, name, accession, site, position

    def _create_modification_entry(self, name, accession, site, position):
        """Create a new modification entry tuple."""
        return (name, accession, site, position)

    def _update_modification_entry(
        self, existing_entry, name, accession, site, position
    ):
        """Update an existing modification entry with new data."""
        updated_name = name if name is not None else existing_entry[0]
        updated_accession = accession if accession is not None else existing_entry[1]
        updated_site = site if site is not None else existing_entry[2]
        updated_position = position if position is not None else existing_entry[3]

        return (updated_name, updated_accession, updated_site, updated_position)

    def _group_modifications_by_accession(self, temp_modifications):
        """Convert temporary modifications to final format grouped by accession."""
        modifications = {}
        for name, accession, sites, positions in temp_modifications.values():
            if accession not in modifications:
                modifications[accession] = (name, [sites], [positions])
            else:
                modifications[accession][1].append(sites)
                modifications[accession][2].append(positions)
        return modifications

    def _process_metadata_rows(self, metadata_df):
        """Process metadata rows to extract modification information."""
        temp_modifications = {}

        for _, row in metadata_df.iterrows():
            current_mod_index, name, accession, site, position = (
                self._parse_modification_row(row)
            )

            if current_mod_index is None:
                continue

            if current_mod_index not in temp_modifications:
                temp_modifications[current_mod_index] = self._create_modification_entry(
                    name, accession, site, position
                )
            else:
                temp_modifications[current_mod_index] = self._update_modification_entry(
                    temp_modifications[current_mod_index],
                    name,
                    accession,
                    site,
                    position,
                )

        return temp_modifications

    def _get_metadata_modifications(self) -> dict:
        """Get modifications from metadata."""
        try:
            metadata_df = self._indexer.get_metadata()

            # Process metadata rows to extract modification information
            temp_modifications = self._process_metadata_rows(metadata_df)

            # Group modifications by accession
            modifications = self._group_modifications_by_accession(temp_modifications)

            return modifications
        except Exception as e:
            self.logger.warning(f"Could not get modifications: {e}")
            return {}

    def iter_psm_table(
        self, chunksize: int = 1000000, protein_str: Optional[str] = None
    ):
        """Iterate over PSM table in chunks.

        Args:
            chunksize: Number of rows to process in each chunk
            protein_str: Optional protein accession filter

        Yields:
            PyArrow Table containing PSM data with proper struct types
        """
        for df in self._indexer.stream_section("PSM", chunk_size=chunksize):
            if protein_str:
                df = df[df["accession"].str.contains(f"{protein_str}", na=False)]

            # Convert to list of dictionaries for safer transformations
            records = df.to_dict("records")

            # Transform each record individually to PyArrow-compatible format
            transformed_records = []
            for record in records:
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

    @staticmethod
    def search_best_protein_global_qvalue(
        accessions: Union[str, List[str]], qvalue_map: dict
    ) -> float | None:
        """
        Return the qvalue of the protein group that contains all the accessions. If not found, return None.
        """
        if not accessions or not qvalue_map:
            return None

        # Convert list to string if needed
        if isinstance(accessions, list):
            accessions = ";".join(sorted(accessions))

        # First, try the exact combination
        if accessions in qvalue_map:
            return qvalue_map[accessions]

        return None

    def _create_file_metadata(self):
        """Create file metadata structure according to PSM avsc schema"""
        import uuid
        from datetime import datetime
        from quantmsio import __version__

        return {
            "quantmsio_version": __version__,
            "creator": "quantms.io",
            "file_type": "psm_file",
            "creation_date": datetime.now().isoformat(),
            "uuid": str(uuid.uuid4()),
            "scan_format": "scan",  # Default scan format
            "software_provider": "quantms.io",
        }

    def _transform_psm_record_to_arrow(self, record: dict) -> Optional[dict]:
        """Transform a single PSM record to PyArrow-compatible format.

        Args:
            record: Dictionary containing PSM data

        Returns:
            Transformed record dictionary with PyArrow structs or None if invalid
        """
        try:
            # Create a new record with only PSM_MAP columns
            transformed_record = {}

            peptidoform, modification_details = self._parse_modifications_for_arrow(
                openms_peptidoform=record[OPENMS_PEPTIDOFORM_COLUMN],
                openms_modifications=record["modifications"],
                reference_modifications=self._modifications,
            )

            transformed_record["peptidoform"] = peptidoform
            transformed_record["modifications"] = modification_details
            transformed_record["sequence"] = record["sequence"]
            transformed_record["precursor_charge"] = int(record["charge"])
            transformed_record["rt"] = float(record["retention_time"])
            transformed_record["protein_accessions"] = (
                standardize_protein_string_accession(
                    record["accession"], is_sorted=True
                ).split(";")
            )

            # Generate additional scores as PyArrow structs
            additional_scores = self._generate_additional_scores_for_arrow(record)
            transformed_record["additional_scores"] = additional_scores

            ## Add the protein score to the additional scores if protein_accessions, exist in self._protein_global_qvalue_map
            best_protein_global_qvalue = self.search_best_protein_global_qvalue(
                transformed_record["protein_accessions"],
                self._protein_global_qvalue_map,
            )
            if best_protein_global_qvalue is not None:
                transformed_record["additional_scores"].append(
                    {
                        "score_name": "protein_global_qvalue",
                        "score_value": float(best_protein_global_qvalue),
                    }
                )

            # Generate CV parameters as PyArrow structs
            consensus_support = record.get("consensus_support")
            if consensus_support is None:
                consensus_support = record.get("opt_global_consensus_support")
            cv_params = self._generate_cv_params_for_arrow(consensus_support)
            transformed_record["cv_params"] = cv_params

            # Handle spectra_ref processing
            spectra_ref = record.get("spectra_ref")
            transformed_record["reference_file_name"] = self._get_reference_file_name(
                spectra_ref
            )

            transformed_record["scan"] = generate_scan_number(record["spectra_ref"])

            # Add posterior error probability
            if OPENMS_POSTERIOR_ERRORPROBABILITY in record:
                transformed_record["posterior_error_probability"] = float(
                    record[OPENMS_POSTERIOR_ERRORPROBABILITY]
                )
            else:
                transformed_record["posterior_error_probability"] = None

            # Add is_decoy
            if OPENMS_IS_DECOY in record:
                transformed_record["is_decoy"] = int(record[OPENMS_IS_DECOY])
            else:
                transformed_record["is_decoy"] = None

            # Add observed_mz
            if "exp_mass_to_charge" in record:
                transformed_record["observed_mz"] = float(record["exp_mass_to_charge"])
            else:
                transformed_record["observed_mz"] = None

            # Add calc_mass_to_charge
            if "calc_mass_to_charge" in record:
                transformed_record["calculated_mz"] = float(
                    record["calc_mass_to_charge"]
                )
            else:
                transformed_record["calculated_mz"] = None

            # Add predicted RT
            if "opt_global_predicted_rt" in record:
                transformed_record["predicted_rt"] = float(
                    record["opt_global_predicted_rt"]
                )
            else:
                transformed_record["predicted_rt"] = None

            # Non supported columns
            transformed_record["number_peaks"] = None
            transformed_record["ion_mobility"] = None
            transformed_record["mz_array"] = None
            transformed_record["intensity_array"] = None
            transformed_record["charge_array"] = None
            transformed_record["ion_type_array"] = None
            transformed_record["ion_mobility_array"] = None

            return transformed_record

        except Exception as e:
            self.logger.warning(f"Error transforming PSM record: {e}")
            return None

    def _parse_modifications_for_arrow(
        self,
        openms_peptidoform: str,
        openms_modifications: str,
        reference_modifications: dict,
    ) -> Tuple[str, List[Dict[str, Any]]]:
        """
        The following fucntion will take a peptiform in mzTab format from OpenMS like:
        PEPTIDE(Oxidation)R
        and return a list of modifications in the following format:
        [{"modification_name": "Oxidation", "fields": [{"position": 1, "localization_probability": 1.0}]}]
        In addition, it will return the peptidoform in ProForma format:
        PEPTIDE[Oxidation]R
        """
        peptidoform, modification_details = parse_pepidoform_with_modifications(
            openms_peptidoform=openms_peptidoform,
            openms_modifications=openms_modifications,
            reference_modifications=reference_modifications,
        )
        return peptidoform, modification_details

    def _generate_additional_scores_for_arrow(
        self, record: dict
    ) -> List[Dict[str, Any]]:
        """Generate additional scores as list of dicts for PyArrow structs.

        Args:
            record: Dictionary containing PSM data

        Returns:
            List of score dictionaries compatible with PyArrow structs
        """
        struct_list = []

        # Process score names
        for score_index, score_name in self._score_names["psms"].items():
            value = record[f"search_engine_score[{score_index}]"]
            if value is not None and value != "null":
                try:
                    struct = {
                        "score_name": score_name,
                        "score_value": float(value),
                    }
                    struct_list.append(struct)
                except (ValueError, TypeError):
                    self.logger.warning(
                        f"Could not convert score value '{value}' to float for score '{score_name}'"
                    )

        # Handle global_qvalue in the PSM table
        if (
            "opt_global_q-value" in record and record["opt_global_q-value"] is not None
        ) or ("global_qvalue" in record and record["global_qvalue"] is not None):
            try:
                global_qvalue = float(
                    record["opt_global_q-value"]
                    if "opt_global_q-value" in record
                    else record["global_qvalue"]
                )
                struct = {
                    "score_name": "global_qvalue",
                    "score_value": global_qvalue,
                }
                struct_list.append(struct)
            except (ValueError, TypeError):
                self.logger.warning(f"Could not convert global_qvalue to float")

        return struct_list

    def _generate_cv_params_for_arrow(
        self, consensus_support
    ) -> Optional[List[Dict[str, str]]]:
        """Generate CV parameters as list of dicts for PyArrow structs.

        Args:
            consensus_support: Consensus support value

        Returns:
            List of CV parameter dictionaries or None
        """
        cv_list = []
        if consensus_support and consensus_support != "null":
            struct = {
                "cv_name": "consesus_support",
                "cv_value": str(consensus_support),
            }
            cv_list.append(struct)

        return cv_list

    def _create_arrow_table(self, records: List[dict]) -> pa.Table:
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

    def _get_reference_file_name(self, spectra_ref) -> Optional[str]:
        """Get reference file name from spectra_ref.

        Args:
            spectra_ref: Spectra reference string

        Returns:
            Reference file name or None
        """
        if not spectra_ref or spectra_ref == "null":
            return None
        try:
            colon_idx = spectra_ref.index(":")
            ms_run_id_str = spectra_ref[:colon_idx]
            # Convert ms_run[X] format to numeric index
            ms_run_id = int(ms_run_id_str.split("[")[1].split("]")[0])
            return self._ms_runs.get(ms_run_id)
        except (ValueError, AttributeError, IndexError):
            return None

    def _extract_pep_columns(self) -> None:
        """Extract peptide columns from DuckDB."""
        try:
            # Get column names from DuckDB
            source = self._indexer._get_table_source(
                self._indexer._MZTAB_INDEXER_TABLE_PSMS
            )
            columns = self._indexer._duckdb.execute(
                f"PRAGMA table_info({source})"
            ).fetchall()
            self._pep_columns = [col[1] for col in columns]
        except Exception as e:
            self.logger.error(f"Could not extract peptide columns: {e}")
            self._pep_columns = []

    def extract_from_pep(self, chunksize: int = 2000000) -> dict:
        """Extract peptide information from DuckDB.

        Args:
            chunksize: Number of rows to process in each chunk

        Returns:
            Dictionary mapping peptide sequences to their best scores
        """
        self._extract_pep_columns()
        pep_map = {}

        try:
            # Query required columns
            source = self._indexer._get_table_source(
                self._indexer._MZTAB_INDEXER_TABLE_PSMS
            )
            query = f"""
            SELECT 
                sequence,
                modifications,
                charge,
                best_search_engine_score[1] as score,
                spectra_ref,
                scan
            FROM {source}
            """

            # Process in chunks
            offset = 0
            while True:
                chunk_query = f"{query} LIMIT {chunksize} OFFSET {offset}"
                df = self._indexer._duckdb.execute(chunk_query).df()

                if df.empty:
                    break

                # Process each row
                for _, row in df.iterrows():
                    peptidoform = get_petidoform_msstats_notation(
                        row["sequence"], row["modifications"], self._modifications
                    )
                    key = (peptidoform, row["charge"])

                    # Get scan number
                    scan_number = None
                    if pd.notna(row["spectra_ref"]):
                        scan_number = generate_scan_number(row["spectra_ref"])

                    # Get reference file name
                    ref_file = None
                    if pd.notna(row["spectra_ref"]):
                        try:
                            ms_run_id = row["spectra_ref"].split(":")[0]
                            ref_file = self._ms_runs.get(ms_run_id)
                        except (ValueError, AttributeError):
                            pass

                    # Update map with best score
                    if key not in pep_map or row["score"] < pep_map[key][0]:
                        pep_map[key] = [row["score"], ref_file, scan_number]

                offset += chunksize

        except Exception as e:
            self.logger.error(f"Error extracting peptide information: {e}")
            return {}

        return pep_map

    @staticmethod
    def slice(df: pd.DataFrame, partitions: list):
        cols = df.columns
        if not isinstance(partitions, list):
            raise Exception(f"{partitions} is not a list")
        if len(partitions) == 0:
            raise Exception(f"{partitions} is empty")
        for partion in partitions:
            if partion not in cols:
                raise Exception(f"{partion} does not exist")
        for key, df in df.groupby(partitions):
            yield key, df

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
            protein_file: Optional protein filter file
        """
        logger = logging.getLogger("quantmsio.core.psm")

        # Log input and output paths
        logger.info(f"Output path: {output_path}")
        logger.info(f"Converting Indexer from path: {self._indexer._database_path}")
        num_psms = self._indexer._get_num_psms()
        logger.info(f"Number of PSMs: {num_psms}")

        if self._spectral_data:
            logger.info(
                "Loading spectra information into quantms.io, but the required data is not included in quantms mzTab, so it will be empty."
            )
        else:
            logger.info("Spectra information will not be loaded into quantms.io")

        if protein_file:
            logger.info(f"Protein filter file: {protein_file}")

        protein_list: list = (
            extract_protein_list(protein_file) if protein_file else None
        )
        protein_str: Optional[str] = "|".join(protein_list) if protein_list else None

        # Create file metadata for parquet file
        file_metadata = self._create_file_metadata()

        # Initialize batch writer with file metadata
        batch_writer = ParquetBatchWriter(
            output_path, PSM_SCHEMA, file_metadata=file_metadata
        )

        try:
            # Process data in chunks and write batches
            for table in self.iter_psm_table(
                chunksize=chunksize, protein_str=protein_str
            ):
                transformed_records = table.to_pylist()
                # Write batch if we have records
                if transformed_records:
                    batch_writer.write_batch(transformed_records)

        finally:
            # Ensure final batch is written and writer is closed
            batch_writer.close()

            if Path(output_path).exists():
                self.logger.info(f"[Writer] Successfully wrote PSM to: {output_path}")

            # Clean up the temporary MzTabIndexer
            self._indexer.cleanup_duckdb()
