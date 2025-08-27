"""QuantMS-IO protein groups processing module."""

import logging
import os
import time
from pathlib import Path
from typing import Optional, Union

import pandas as pd
import pyarrow as pa

from quantmsio.core.format import PG_SCHEMA
from quantmsio.core.quantms.mztab import MzTab


class MzTabProteinGroups(MzTab):
    """Handle protein groups in mzTab format with optimized quantification."""

    def __init__(self, mztab_path: Union[Path, str]):
        # Initialize tracking lists first
        self._temp_files = []  # Track any temporary files created
        self._file_handles = []  # Track any open file handles

        super().__init__(mztab_path)
        self._protein_columns = self._extract_protein_columns()

    def __enter__(self):
        """Context manager entry."""
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Context manager exit with cleanup."""
        self.cleanup()

    def cleanup(self):
        """Clean up any temporary files and resources."""

        logger = logging.getLogger("quantmsio.core.mztab")

        # Close any open file handles
        for file_handle in self._file_handles:
            try:
                if hasattr(file_handle, "close") and not file_handle.closed:
                    file_handle.close()
                    logger.debug("Closed file handle")
            except Exception as e:
                logger.warning(f"Error closing file handle: {e}")
        self._file_handles.clear()

        # Clean up any temporary files
        for temp_file in self._temp_files:
            try:
                if os.path.exists(temp_file):
                    os.unlink(temp_file)
                    logger.debug(f"Deleted temporary file: {temp_file}")
            except Exception as e:
                logger.warning(f"Error deleting temporary file {temp_file}: {e}")
        self._temp_files.clear()

    def __del__(self):
        """Destructor to ensure cleanup."""
        try:
            self.cleanup()
        except Exception as e:
            logging.getLogger("quantmsio.core.quantms.pg").warning(
                f"Exception during __del__ cleanup: {e}"
            )

    def _extract_protein_columns(self):
        """Extract protein columns from mzTab header or first PRT line."""
        protein_columns = []
        try:
            # Use safe file opening with automatic cleanup
            with self._safe_file_open(self.mztab_path, "r") as file:
                for line in file:
                    if line.startswith("PRH"):
                        protein_columns = line.strip().split("\t")[1:]
                        break
                    elif line.startswith("PRT\t") and not protein_columns:
                        # Fallback: use first PRT line to determine column count
                        prt_parts = line.strip().split("\t")
                        # Generate default column names based on actual data
                        protein_columns = [f"col_{i}" for i in range(len(prt_parts))]
                        break
        except Exception as e:
            logging.getLogger("quantmsio.core.mztab").warning(
                f"Could not extract protein columns: {e}"
            )
        finally:
            # File handle is automatically tracked and will be cleaned up
            pass
        return protein_columns

    def iter_protein_groups_batch(
        self, chunksize: int = 1000000, protein_str: Optional[str] = None
    ):
        """Iterate over protein groups in chunks for memory efficiency.

        Args:
            chunksize: Number of rows per chunk
            protein_str: Optional protein filter string

        Yields:
            pd.DataFrame: Chunk of protein groups data
        """
        logger = logging.getLogger("quantmsio.core.mztab")

        try:
            protein_lines = []
            # Use safe file opening - file handle will be tracked for cleanup
            with self._safe_file_open(self.mztab_path, "r") as file:
                for line in file:
                    if line.startswith("PRT\t"):
                        parts = line.strip().split("\t")
                        # Remove the "PRT" identifier (first column) to get actual data
                        data_parts = parts[1:]

                        # Adjust column names to match actual data length
                        if len(data_parts) != len(self._protein_columns):
                            logger.warning(
                                f"Column mismatch: header has {len(self._protein_columns)} columns, data has {len(data_parts)} columns. Adjusting..."
                            )
                            if len(data_parts) > len(self._protein_columns):
                                # Add missing column names
                                while len(self._protein_columns) < len(data_parts):
                                    self._protein_columns.append(
                                        f"extra_col_{len(self._protein_columns)}"
                                    )
                            else:
                                # Truncate excess column names
                                self._protein_columns = self._protein_columns[
                                    : len(data_parts)
                                ]

                        protein_lines.append(data_parts)

                        # Yield chunk when we reach chunksize
                        if len(protein_lines) >= chunksize:
                            chunk = pd.DataFrame(
                                protein_lines, columns=self._protein_columns
                            )

                            # Filter by protein if specified
                            if protein_str and not chunk.empty:
                                mask = chunk.get("accession", pd.Series()).str.contains(
                                    protein_str, case=False, na=False
                                )
                                chunk = chunk[mask]

                            if not chunk.empty:
                                yield chunk
                            protein_lines = []

            # Yield remaining lines
            if protein_lines:
                chunk = pd.DataFrame(protein_lines, columns=self._protein_columns)

                # Filter by protein if specified
                if protein_str and not chunk.empty:
                    mask = chunk.get("accession", pd.Series()).str.contains(
                        protein_str, case=False, na=False
                    )
                    chunk = chunk[mask]

                if not chunk.empty:
                    yield chunk

        except Exception as e:
            logger.error(f"Error reading protein groups: {e}")
            # Return empty DataFrame as fallback
            yield pd.DataFrame(columns=self._protein_columns)

    def _is_protein_line(self, line_num: int) -> bool:
        """Check if a line is a protein data line."""
        try:
            # Use safe file opening with automatic cleanup
            with self._safe_file_open(self.mztab_path, "r") as file:
                for i, line in enumerate(file):
                    if i == line_num:
                        return line.startswith("PRT")
            return False
        except Exception:
            return False

    def _extract_protein_name(self, accession: str) -> str:
        """Extract protein name from accession."""
        if not accession or pd.isna(accession):
            return ""

        # Handle different accession formats
        if "|" in accession:
            # UniProt format: sp|P12345|PROT_HUMAN
            parts = accession.split("|")
            if len(parts) >= 3:
                return parts[2]  # PROT_HUMAN
            elif len(parts) >= 2:
                return parts[1]  # P12345

        return accession

    def _extract_gene_names(self, description: str) -> list:
        """Extract gene names from protein description."""
        import re

        if not description or pd.isna(description):
            return []

        gene_names = []

        # Look for GN= pattern (common in UniProt)
        gn_match = re.search(r"GN=([^\s]+)", description)
        if gn_match:
            gene_names.append(gn_match.group(1))

        # Look for Gene_Symbol= pattern
        gs_match = re.search(r"Gene_Symbol=([^\s]+)", description)
        if gs_match:
            gene_names.append(gs_match.group(1))

        return list(set(gene_names))  # Remove duplicates

    def _safe_int_conversion(self, value, default=0):
        """Safely convert value to int with fallback."""
        try:
            if pd.isna(value) or value in ["null", "NULL", "", None]:
                return default
            return int(float(value))
        except (ValueError, TypeError):
            return default

    def _safe_file_open(self, file_path, mode="r"):
        """Safely open files with proper resource tracking."""
        import gzip
        from contextlib import contextmanager

        @contextmanager
        def managed_file():
            file_handle = None
            try:
                # Handle both compressed and uncompressed files
                if str(file_path).endswith(".gz"):
                    # For gzipped files, use gzip.open with text mode
                    file_handle = gzip.open(
                        file_path, "rt" if "t" not in mode else mode, encoding="utf-8"
                    )
                else:
                    # For regular files
                    file_handle = open(file_path, mode, encoding="utf-8")

                # Track the file handle for cleanup (only if context manager doesn't close it)
                self._file_handles.append(file_handle)

                yield file_handle

            except Exception as e:
                logging.getLogger("quantmsio.core.mztab").error(
                    f"Error opening file {file_path}: {e}"
                )
                raise
            finally:
                if file_handle:
                    # Remove from tracking list since we're closing it here
                    try:
                        self._file_handles.remove(file_handle)
                    except ValueError as e:
                        logging.getLogger("quantmsio.core.quantms.pg").warning(
                            f"File handle already removed: {e}"
                        )

                    # Close the file
                    try:
                        if not file_handle.closed:
                            file_handle.close()
                    except Exception as e:
                        logging.getLogger("quantmsio.core.quantms.pg").warning(
                            f"Exception closing file handle: {e}"
                        )

        return managed_file()

    def _convert_to_parquet_format(self, df: pd.DataFrame) -> pa.Table:
        """Convert DataFrame to parquet format using PG_SCHEMA."""
        if df.empty:
            # For empty DataFrames, create an empty table with the correct schema
            return pa.Table.from_arrays(
                [pa.array([], type=field.type) for field in PG_SCHEMA], schema=PG_SCHEMA
            )
        else:
            return pa.Table.from_pandas(df, schema=PG_SCHEMA, preserve_index=False)

    def quantify_from_msstats_optimized(
        self,
        msstats_path: str,
        sdrf_path: str,
        compute_topn: bool = True,
        topn: int = 3,
        compute_ibaq: bool = True,
        file_num: int = 10,
        duckdb_max_memory: str = "16GB",
        duckdb_threads: int = 4,
    ) -> pd.DataFrame:
        """Optimized protein quantification using DuckDB SQL aggregation."""
        logger = logging.getLogger("quantmsio.core.mztab")
        logger.info(
            "[OPTIMIZED] Starting protein group quantification using DuckDB SQL"
        )

        total_start = time.time()

        # Step 1: Create protein groups table from mzTab (only once)
        logger.info("[SETUP] Creating protein groups table from mzTab...")
        pg_start = time.time()
        protein_groups_info = self._create_protein_groups_table_optimized()
        pg_time = time.time() - pg_start
        logger.info(
            f"[SETUP] Created protein groups table with {len(protein_groups_info)} entries in {pg_time:.2f}s"
        )

        # Step 2: Initialize MsstatsIN with DuckDB
        logger.info("[DATA] Loading msstats data with DuckDB...")
        msstats_start = time.time()
        from quantmsio.core.quantms.msstats_in import MsstatsIN

        with MsstatsIN(
            msstats_path, sdrf_path, duckdb_max_memory, duckdb_threads
        ) as msstats_in:
            # Set up optimized processing views in MsstatsIN
            msstats_in._setup_optimized_processing()

            # Get experiment type
            experiment_type = msstats_in.experiment_type
            logger.info(f"[INFO] Detected experiment type: {experiment_type}")

            # Step 3: Create joined view between msstats and protein groups
            self._create_msstats_protein_join_optimized(msstats_in, protein_groups_info)
            msstats_time = time.time() - msstats_start
            logger.info(f"[DATA] MsstatsIN setup completed in {msstats_time:.2f}s")

            # Step 4: Process in batches using SQL aggregation
            logger.info("[PROCESS] Processing protein quantification in batches...")
            process_start = time.time()

            expanded_rows = []
            processed_files = 0

            for file_batch in self._get_file_batches_optimized(msstats_in, file_num):
                batch_start = time.time()

                # Generate experiment-specific SQL
                sql = self._get_protein_aggregation_sql(experiment_type, file_batch)

                # Execute SQL aggregation
                try:
                    batch_results = msstats_in._duckdb.execute(sql).df()

                    logger.info(
                        f"[SQL] Returned {len(batch_results)} rows for batch {file_batch}"
                    )

                    if len(batch_results) > 0:
                        # Transform results to final schema
                        for _, row in batch_results.iterrows():
                            protein_row = self._create_optimized_protein_row(
                                row,
                                protein_groups_info,
                                experiment_type,
                                compute_topn,
                                topn,
                                compute_ibaq,
                            )
                            expanded_rows.append(protein_row)

                        logger.info(
                            f"[SUCCESS] Converted {len(batch_results)} SQL rows to {len(batch_results)} protein rows"
                        )
                    else:
                        logger.warning(
                            f"[WARNING] No data returned from SQL for files: {file_batch}"
                        )

                    processed_files += len(file_batch)
                    batch_time = time.time() - batch_start
                    logger.info(
                        f"[BATCH] Processed {len(file_batch)} files in {batch_time:.2f}s ({processed_files} total)"
                    )

                except Exception as e:
                    logger.warning(f"[ERROR] SQL batch failed: {e}, skipping batch")
                    continue

            process_time = time.time() - process_start
            logger.info(
                f"[PROCESS] Completed quantification processing in {process_time:.2f}s"
            )

            # Step 5: Convert to DataFrame
            logger.info("[CONVERT] Converting results to DataFrame...")
            df_start = time.time()

            if expanded_rows:
                result_df = pd.DataFrame(expanded_rows)
                logger.info(
                    f"[CONVERT] Created DataFrame with {len(result_df)} rows and {len(result_df.columns)} columns"
                )
            else:
                logger.warning(
                    "[WARNING] No data to convert - creating empty DataFrame"
                )
                result_df = pd.DataFrame(
                    columns=[
                        "pg_accessions",
                        "anchor_protein",
                        "pg_names",
                        "gg_accessions",
                        "reference_file_name",
                        "intensities",
                        "additional_intensities",
                        "is_decoy",
                        "contaminant",
                        "peptides",
                        "additional_scores",
                        "global_qvalue",
                        "molecular_weight",
                    ]
                )

            df_time = time.time() - df_start
            total_time = time.time() - total_start

            logger.info(f"[CONVERT] DataFrame conversion completed in {df_time:.2f}s")
            logger.info(
                f"[SUCCESS] OPTIMIZED quantification completed in {total_time:.2f}s total"
            )

            return result_df

        # Context manager automatically cleans up DuckDB resources
        # Cleanup any temporary files created during processing
        self.cleanup()

    def _create_protein_groups_table_optimized(self) -> dict:
        """Create protein groups lookup from mzTab protein section."""
        logger = logging.getLogger("quantmsio.core.mztab")
        protein_groups = {}

        try:
            # Parse protein groups section efficiently
            total_rows = 0
            for chunk in self.iter_protein_groups_batch():
                logger.debug(f"Processing chunk with {len(chunk)} rows")
                for _, row in chunk.iterrows():
                    total_rows += 1
                    result_type = row.get("opt_global_result_type", "single_protein")
                    accession = row.get("accession")

                    if pd.isna(accession) or not accession:
                        continue

                    # Determine pg_accessions based on result type
                    if result_type in ["single_protein", "unknown", "protein_details"]:
                        pg_accessions = [accession.strip()]
                        anchor_protein = accession.strip()
                    else:
                        # Skip other unknown types
                        continue

                    # Store protein group information
                    protein_groups[anchor_protein] = {
                        "result_type": result_type,
                        "pg_accessions": pg_accessions,
                        "anchor_protein": anchor_protein,
                        "description": row.get("description", ""),
                        "sequence_coverage": row.get("protein_coverage"),
                        "global_qvalue": row.get("best_search_engine_score[1]"),
                        "peptide_count": row.get("opt_global_nr_found_peptides"),
                        "is_decoy": self._safe_int_conversion(
                            row.get("opt_global_cv_PRIDE:0000303_decoy_hit", 0)
                        ),
                        "sequence_length": row.get("sequence_length", 0),
                        "pg_names": [
                            self._extract_protein_name(acc) for acc in pg_accessions
                        ],
                        "gg_accessions": self._extract_gene_names(
                            row.get("description", "")
                        ),
                    }

        except Exception as e:
            logger.error(f"Error creating protein groups table: {e}")

        logger.info(
            f"Created {len(protein_groups)} protein groups from {total_rows} total rows"
        )
        return protein_groups

    def _create_msstats_protein_join_optimized(
        self, msstats_in, protein_groups_info: dict
    ):
        """Create optimized join view between msstats and protein groups in DuckDB."""
        logger = logging.getLogger("quantmsio.core.mztab")

        try:
            # Create protein groups lookup table in DuckDB
            protein_data = []
            for anchor, info in protein_groups_info.items():
                # Create entries for all proteins in the group for lookup
                for protein in info["pg_accessions"]:
                    protein_data.append(
                        {
                            "protein_name": protein,
                            "anchor_protein": anchor,
                            "result_type": info["result_type"],
                            "pg_accessions": ";".join(info["pg_accessions"]),
                            "description": info["description"],
                            "sequence_coverage": info["sequence_coverage"],
                            "global_qvalue": info["global_qvalue"],
                            "peptide_count": info["peptide_count"],
                            "is_decoy": info["is_decoy"],
                            "sequence_length": info["sequence_length"],
                        }
                    )

            # Convert to DataFrame and load into DuckDB
            if protein_data:
                protein_df = pd.DataFrame(protein_data)
                msstats_in._duckdb.execute("DROP TABLE IF EXISTS protein_groups")
                msstats_in._duckdb.execute(
                    "CREATE TABLE protein_groups AS SELECT * FROM protein_df"
                )

                # Create index for better join performance
                msstats_in._duckdb.execute(
                    "CREATE INDEX IF NOT EXISTS idx_protein_name ON protein_groups(protein_name)"
                )

                # Simplified join - use exact match first, then fallback for unmatched
                msstats_in._duckdb.execute(
                    """
                    DROP VIEW IF EXISTS processed_msstats_with_pg;
                    CREATE VIEW processed_msstats_with_pg AS
                    SELECT
                        m.*,
                        COALESCE(pg.anchor_protein, m.pg_accessions_raw) as anchor_protein,
                        pg.result_type,
                        pg.pg_accessions,
                        pg.description,
                        pg.sequence_coverage,
                        pg.global_qvalue,
                        pg.peptide_count,
                        pg.is_decoy,
                        pg.sequence_length
                    FROM processed_msstats m
                    LEFT JOIN protein_groups pg ON m.pg_accessions_raw = pg.protein_name
                """
                )

                # Log statistics
                count = msstats_in._duckdb.execute(
                    "SELECT COUNT(*) FROM processed_msstats"
                ).fetchone()[0]
                matched_count = msstats_in._duckdb.execute(
                    "SELECT COUNT(*) FROM processed_msstats_with_pg WHERE anchor_protein IS NOT NULL"
                ).fetchone()[0]
                logger.info(
                    f"[DATA] Created processed_msstats view with {count} rows, {matched_count} matched to protein groups"
                )

            else:
                logger.warning(
                    "[WARNING] No protein groups data to create lookup table"
                )

        except Exception as e:
            logger.error(f"Error creating msstats protein join: {e}")
            raise

    def _get_file_batches_optimized(self, msstats_in, batch_size: int):
        """Get file batches for processing."""
        try:
            # Get unique files from the view
            files_df = msstats_in._duckdb.execute(
                "SELECT DISTINCT reference_file_name FROM processed_msstats_with_pg ORDER BY reference_file_name"
            ).df()

            unique_files = files_df["reference_file_name"].tolist()

            # Yield batches
            for i in range(0, len(unique_files), batch_size):
                yield unique_files[i : i + batch_size]

        except Exception as e:
            logging.getLogger("quantmsio.core.mztab").error(
                f"Error getting file batches: {e}"
            )
            return

    def _get_protein_aggregation_sql(
        self, experiment_type: str, file_batch: list
    ) -> str:
        """Generate SQL for protein aggregation based on experiment type."""

        # Simplified SQL without expensive ARRAY_AGG operations
        base_sql = f"""
        SELECT 
            anchor_protein,
            pg_accessions,
            result_type,
            description,
            sequence_coverage,
            global_qvalue,
            peptide_count,
            is_decoy,
            sequence_length,
            reference_file_name,
            channel,
            SUM(intensity) as total_intensity,
            COUNT(*) as feature_count,
            COUNT(DISTINCT peptidoform) as unique_peptide_count,
            MAX(intensity) as max_intensity,
            AVG(intensity) as avg_intensity
        FROM processed_msstats_with_pg 
        WHERE reference_file_name = ANY({file_batch})
        AND anchor_protein IS NOT NULL
        AND intensity > 0
        GROUP BY anchor_protein, pg_accessions, result_type, description, 
                 sequence_coverage, global_qvalue, peptide_count, is_decoy, 
                 sequence_length, reference_file_name, channel
        ORDER BY anchor_protein, reference_file_name, channel
        """

        return base_sql

    def _create_optimized_protein_row(
        self,
        row: pd.Series,
        protein_groups_info: dict,
        experiment_type: str,
        compute_topn: bool = True,
        topn: int = 3,
        compute_ibaq: bool = True,
    ) -> dict:
        """Create protein row from SQL aggregation results."""

        anchor_protein = row["anchor_protein"]
        file_name = row["reference_file_name"]
        channel = row["channel"]
        total_intensity = float(row["total_intensity"])
        max_intensity = float(row.get("max_intensity", total_intensity))
        avg_intensity = float(row.get("avg_intensity", total_intensity))
        unique_peptide_count = int(row.get("unique_peptide_count", 1))

        # Get protein group info
        if anchor_protein in protein_groups_info:
            pg_info = protein_groups_info[anchor_protein]
        else:
            # Fallback info
            pg_accessions = (
                row["pg_accessions"].split(";")
                if pd.notna(row["pg_accessions"])
                else [anchor_protein]
            )
            pg_info = {
                "pg_accessions": pg_accessions,
                "anchor_protein": anchor_protein,
                "pg_names": [self._extract_protein_name(acc) for acc in pg_accessions],
                "gg_accessions": self._extract_gene_names(row.get("description", "")),
                "is_decoy": self._safe_int_conversion(row.get("is_decoy", 0)),
                "sequence_length": self._safe_int_conversion(
                    row.get("sequence_length", 0)
                ),
                "global_qvalue": (
                    float(row.get("global_qvalue", 1.0))
                    if pd.notna(row.get("global_qvalue"))
                    else 1.0
                ),
                "peptide_count": self._safe_int_conversion(
                    row.get("peptide_count", unique_peptide_count)
                ),
                "sequence_coverage": (
                    row.get("sequence_coverage")
                    if pd.notna(row.get("sequence_coverage"))
                    else None
                ),
            }

        # Main intensity
        intensities = [
            {
                "sample_accession": file_name,
                "channel": channel,
                "intensity": total_intensity,
            }
        ]

        # Additional intensities - simplified calculations without peptide details
        additional_intensities = []

        if compute_topn:
            # TopN: Use max intensity as approximation (since we don't have individual peptide data)
            topn_intensity = max_intensity

            additional_intensities.append(
                {
                    "sample_accession": file_name,
                    "channel": channel,
                    "intensity": topn_intensity,
                    "intensity_type": "TopN",
                }
            )

        if compute_ibaq and pg_info["sequence_length"] > 0:
            # iBAQ: total intensity / sequence length
            ibaq_intensity = total_intensity / max(pg_info["sequence_length"], 1)

            additional_intensities.append(
                {
                    "sample_accession": file_name,
                    "channel": channel,
                    "intensity": ibaq_intensity,
                    "intensity_type": "iBAQ",
                }
            )

        # Create peptides structure
        peptides = []
        for protein in pg_info["pg_accessions"]:
            peptides.append(
                {
                    "protein_name": protein,
                    "peptide_count": max(
                        pg_info["peptide_count"], unique_peptide_count
                    ),
                }
            )

        # Create additional scores
        additional_scores = []
        if pg_info.get("sequence_coverage") is not None:
            additional_scores.append(
                {
                    "score_name": "sequence_coverage_percent",
                    "score_value": float(pg_info["sequence_coverage"]),
                }
            )
        additional_scores.append(
            {
                "score_name": "peptide_count",
                "score_value": float(
                    max(pg_info["peptide_count"], unique_peptide_count)
                ),
            }
        )
        additional_scores.append(
            {
                "score_name": "max_intensity",
                "score_value": max_intensity,
            }
        )
        additional_scores.append(
            {
                "score_name": "avg_intensity",
                "score_value": avg_intensity,
            }
        )

        return {
            "pg_accessions": pg_info["pg_accessions"],
            "anchor_protein": pg_info["anchor_protein"],
            "pg_names": pg_info["pg_names"],
            "gg_accessions": pg_info["gg_accessions"],
            "reference_file_name": file_name,
            "intensities": intensities,
            "additional_intensities": additional_intensities,
            "is_decoy": pg_info["is_decoy"],
            "contaminant": 0,  # mzTab doesn't have contaminant info
            "peptides": peptides,
            "additional_scores": additional_scores,
            "global_qvalue": pg_info["global_qvalue"],
            "peptide_counts": {
                "unique_sequences": max(pg_info["peptide_count"], unique_peptide_count),
                "total_sequences": max(pg_info["peptide_count"], unique_peptide_count),
            },
            "feature_counts": {
                "unique_features": int(row.get("feature_count", 0)),
                "total_features": int(row.get("feature_count", 0)),
            },
        }
