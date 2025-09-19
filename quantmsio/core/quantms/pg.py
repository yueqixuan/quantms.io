"""QuantMS-IO protein groups processing module."""

import logging
import os
import time
from typing import Optional

import pandas as pd
import pyarrow as pa
import re

from quantmsio.core.format import PG_SCHEMA
from quantmsio.core.quantms.mztab import MzTabIndexer
from quantmsio.utils.constants import MZTAB_PROTEIN_BEST_SEARCH_ENGINE_SCORE


class MzTabProteinGroups:
    """Handle protein groups in mzTab format with optimized quantification using composition pattern."""

    def __init__(self, mztab_indexer):
        """Initialize MzTabProteinGroups with an MzTabIndexer instance.

        Args:
            mztab_path: Path to mzTab file
        """
        self._indexer: MzTabIndexer = mztab_indexer

        # Initialize tracking lists first
        self._temp_files = []  # Track any temporary files created
        self._file_handles = []  # Track any open file handles
        # protein_details_df = self._indexer.get_protein_details()
        self.logger = logging.getLogger(f"{__name__}.{self.__class__.__name__}")

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

    def _is_protein_line(self, line_num: int) -> bool:
        """Check if a line is a protein data line."""
        try:
            # Use safe file opening with automatic cleanup
            with self._safe_file_open(self._mztab_path, "r") as file:
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

        accessions = accession.split(";")
        acc_list = list()
        for acc in accessions:
            # Handle different accession formats
            if "|" in acc:
                # UniProt format: sp|P12345|PROT_HUMAN
                parts = acc.split("|")
                if len(parts) >= 3:
                    acc_list.append(parts[2])  # PROT_HUMAN
                elif len(parts) >= 2:
                    acc_list.append(parts[1])  # P12345
        if acc_list:
            return ";".join(acc_list)

        return accession

    def _extract_gene_names(self, description: str) -> list:
        """Extract gene names from protein description."""

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
        compute_topn: bool = True,
        topn: int = 3,
        compute_ibaq: bool = True,
        file_num: int = 1,
    ) -> pd.DataFrame:
        """Optimized protein quantification using DuckDB SQL aggregation."""
        logger = logging.getLogger("quantmsio.core.quantms.pg")
        logger.info(
            "[OPTIMIZED] Starting protein group quantification using DuckDB SQL"
        )

        total_start = time.time()

        # Step 1: Loading protein groups table from MzTabIndexer
        logger.info("[SETUP] Loading protein groups table from MzTabIndexer...")
        pg_start = time.time()
        protein_groups_info = self._load_protein_groups_table_optimized()
        pg_time = time.time() - pg_start
        logger.info(
            f"[SETUP] Created protein groups table with {len(protein_groups_info)} entries in {pg_time:.2f}s"
        )

        # Step 2: Loading msstats data from MzTabIndexer for enhanced analysis
        logger.info("[DATA] Loading msstats data from MzTabIndexer...")
        msstats_start = time.time()

        # Step 3: Create joined view between msstats and protein groups
        self._create_msstats_protein_join_optimized(protein_groups_info)
        msstats_time = time.time() - msstats_start
        logger.info(f"[DATA] MzTabIndexer setup completed in {msstats_time:.2f}s")

        # Step 4: Process in batches using SQL aggregation
        logger.info("[PROCESS] Processing protein quantification in batches...")
        process_start = time.time()

        expanded_rows = []
        processed_files = 0

        for batch in self.get_sql_batchs():
            batch_start = time.time()

            try:
                batch_data = self.get_sql_batch_data(batch)

                logger.info(f"[SQL] Returned {len(batch_data)} rows for batch {batch}")

                if len(batch_data) > 0:
                    protein_row = self._create_optimized_protein_row(
                        batch_data,
                        compute_topn,
                        topn,
                        compute_ibaq,
                    )
                    expanded_rows.extend(protein_row)

                    logger.info(
                        f"[SUCCESS] Converted {len(batch_data)} SQL rows to protein rows"
                    )
                else:
                    logger.warning(f"[WARNING] No data returned from SQL")
                processed_files += len(batch_data)
                batch_time = time.time() - batch_start
                logger.info(
                    f"[BATCH] Processed {len(batch_data)} rows in {batch_time:.2f}s ({processed_files} total)"
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
            logger.warning("[WARNING] No data to convert - creating empty DataFrame")
            result_df = pd.DataFrame(
                columns=[
                    "pg_accessions",
                    "pg_names",
                    "gg_accessions",
                    "reference_file_name",
                    "peptide_counts",
                    "feature_counts",
                    "global_qvalue",
                    "intensities",
                    "additional_intensities",
                    "peptides",
                    "anchor_protein",
                    "is_decoy",
                    "contaminant",
                    "additional_scores",
                ]
            )

        df_time = time.time() - df_start
        total_time = time.time() - total_start

        logger.info(f"[CONVERT] DataFrame conversion completed in {df_time:.2f}s")
        logger.info(
            f"[SUCCESS] OPTIMIZED quantification completed in {total_time:.2f}s total"
        )

        # Clean up the temporary MzTabIndexer
        self._indexer.cleanup_duckdb()

        return result_df

    def _load_protein_groups_table_optimized(self):
        """Load protein groups lookup from mzTab protein section."""
        logger = logging.getLogger("quantmsio.core.quantms.pg")

        protein_data = list()

        try:
            protein_df = self._indexer.get_proteins()
            total_rows = 0

            for chunk in self.iter_in_chunks(protein_df):
                for _, row in chunk.iterrows():
                    total_rows += 1

                    result_type = row.get("opt_global_result_type", "single_protein")
                    accession = row.get("accession")

                    if pd.isna(accession) or not accession or accession == "null":
                        continue

                    # Skip other types
                    if result_type not in [
                        "single_protein",
                        "indistinguishable_protein_group",
                    ]:
                        continue

                    pg_accessions = accession.strip()
                    anchor_protein = pg_accessions
                    description = self.check_clean_null(row["description"], None)

                    protein_data.append(
                        {
                            "result_type": result_type,
                            "pg_accessions": pg_accessions,
                            "anchor_protein": anchor_protein,
                            "description": description,
                            "sequence_coverage": self.check_clean_null(
                                row["protein_coverage"], 0
                            ),
                            "global_qvalue": self.check_clean_null(
                                row[MZTAB_PROTEIN_BEST_SEARCH_ENGINE_SCORE], 0
                            ),
                            "peptide_count": self.check_clean_null(
                                row["opt_global_nr_found_peptides"], 0
                            ),
                            "is_decoy": self._safe_int_conversion(
                                self.check_clean_null(
                                    row.get("opt_global_cv_PRIDE:0000303_decoy_hit"), 0
                                )
                            ),
                            "sequence_length": self.check_clean_null(
                                row.get("sequence_length"), 0
                            ),
                            "pg_names": self._extract_protein_name(pg_accessions),
                            "gg_accessions": self._extract_gene_names(description),
                        }
                    )

        except Exception as e:
            logger.error(f"Error loading protein groups table: {e}")

        logger.info(
            f"Loaded {len(protein_data)} protein groups from {total_rows} total rows"
        )
        return protein_data

    def check_clean_null(self, value, fill_value):
        return value if pd.notna(value) and value != "null" else fill_value

    def iter_in_chunks(self, df, chunksize=100000):
        for start in range(0, len(df), chunksize):
            yield df.iloc[start : start + chunksize]

    def _create_msstats_protein_join_optimized(self, protein_groups_info):
        """Create optimized join view between msstats and protein groups in DuckDB."""
        logger = logging.getLogger("quantmsio.core.quantms.pg")

        try:
            if protein_groups_info:
                protein_df = pd.DataFrame(protein_groups_info)

                self._indexer._duckdb.execute("DROP TABLE IF EXISTS protein_groups")
                self._indexer._duckdb.from_df(protein_df).create("protein_groups")

                # Create index for better join performance
                self._indexer._duckdb.execute(
                    "CREATE INDEX IF NOT EXISTS idx_protein_name ON protein_groups(pg_accessions)"
                )

                # Get 'unique' from PSM table
                unique_peptide_df = self._indexer.get_unique_from_psm_table()

                self._indexer._duckdb.execute("DROP TABLE IF EXISTS unique_peptide")
                self._indexer._duckdb.from_df(unique_peptide_df).create(
                    "unique_peptide"
                )

                # Simplified join - use exact match first, then fallback for unmatched
                self._indexer._duckdb.execute(
                    """
                    DROP VIEW IF EXISTS processed_msstats_with_pg;
                    CREATE VIEW processed_msstats_with_pg AS
                    SELECT
                        m.*,
                        COALESCE(pg.anchor_protein, m.pg_accessions) AS anchor_protein,
                        pg.result_type,
                        pg.description,
                        pg.sequence_coverage,
                        pg.global_qvalue,
                        pg.peptide_count,
                        pg.is_decoy,
                        pg.sequence_length,
                        pg.pg_names,
                        pg.gg_accessions,
                        unpep."unique"
                    FROM msstats m
                    LEFT JOIN protein_groups pg ON m.pg_accessions = pg.pg_accessions
                    LEFT JOIN unique_peptide unpep ON m.pg_accessions = unpep.pg_accessions AND m.peptidoform = unpep.peptidoform
                """
                )

                # Log statistics
                count = self._indexer._duckdb.execute(
                    "SELECT COUNT(*) FROM msstats"
                ).fetchone()[0]
                matched_count = self._indexer._duckdb.execute(
                    "SELECT COUNT(*) FROM processed_msstats_with_pg WHERE anchor_protein IS NOT NULL"
                ).fetchone()[0]
                logger.info(
                    f"[DATA] Created msstats view with {count} rows, {matched_count} matched to protein groups"
                )

            else:
                logger.warning(
                    "[WARNING] No protein groups data to create lookup table"
                )

        except Exception as e:
            logger.error(f"Error creating msstats protein join: {e}")
            raise

    def _create_optimized_protein_row(
        self,
        batch_data,
        compute_topn: bool = True,
        topn: int = 3,
        compute_ibaq: bool = True,
    ):

        result = []
        for (anchor_protein, reference_file_name), group in batch_data.groupby(
            ["anchor_protein", "reference_file_name"]
        ):

            unique_group = group[group["unique"].isin([1, "1"])]
            peptide_count = {
                "unique_sequences": unique_group["peptidoform"].nunique(),
                "total_sequences": group["peptidoform"].nunique(),
            }
            feature_count = {
                "unique_features": len(
                    set(zip(unique_group["peptidoform"], unique_group["charge"]))
                ),
                "total_features": len(set(zip(group["peptidoform"], group["charge"]))),
            }

            if float(group["peptide_count"].iloc[0]) != 0:
                unique_peptide_count = float(group["peptide_count"].iloc[0])
            else:
                unique_peptide_count = group["peptidoform"].nunique()

            max_intensity = group["intensity"].max()
            avg_intensity = group["intensity"].mean()

            # TODO peptides
            peptides = [
                {
                    "protein_name": anchor_protein,
                    "peptide_count": unique_peptide_count,
                }
            ]

            intensities = []
            additional_intensities = []

            for channel, channel_group in group.groupby("channel"):

                # sample_accession
                sample_accessions = channel_group["sample_accession"].unique()
                if len(sample_accessions) != 1:
                    raise ValueError("number of sample_accession is more than 1")
                sample_accession = sample_accessions[0]

                # intensity sum
                sum_intensity = float(channel_group["intensity"].sum())

                intensities.append(
                    {
                        "sample_accession": sample_accession,
                        "channel": channel,
                        "intensity": sum_intensity,
                    }
                )

                # additional_intensities
                extra_intensities = []
                if compute_topn:
                    topn_intensity = float(channel_group["intensity"].max())
                    extra_intensities.append(
                        {
                            "intensity_name": "TopN",
                            "intensity_value": topn_intensity,
                        }
                    )
                if compute_ibaq:
                    seq_len = channel_group["sequence_length"].iloc[0]
                    if pd.notna(seq_len):
                        seq_len = float(seq_len)
                        ibaq_intensity = sum_intensity / max(seq_len, 1)
                        extra_intensities.append(
                            {
                                "intensity_name": "ibaq",
                                "intensity_value": ibaq_intensity,
                            }
                        )
                additional_intensities.append(
                    {
                        "sample_accession": sample_accession,
                        "channel": channel,
                        "intensities": extra_intensities,
                    }
                )

                # additional_scores
                extra_scores = []
                sequence_coverage = float(group["sequence_coverage"].iloc[0])
                if sequence_coverage != 0:
                    extra_scores.append(
                        {
                            "score_name": "sequence_coverage_percent",
                            "score_value": sequence_coverage,
                        }
                    )
                extra_scores.append(
                    {
                        "score_name": "peptide_count",
                        "score_value": float(
                            max(
                                float(group["peptide_count"].iloc[0]),
                                unique_peptide_count,
                            )
                        ),
                    }
                )
                extra_scores.append(
                    {
                        "score_name": "max_intensity",
                        "score_value": max_intensity,
                    }
                )
                extra_scores.append(
                    {
                        "score_name": "avg_intensity",
                        "score_value": avg_intensity,
                    }
                )

            result.append(
                {
                    "pg_accessions": group["pg_accessions"].iloc[0].split(";"),
                    "pg_names": group["pg_names"].iloc[0].split(";"),
                    "gg_accessions": group["gg_accessions"].iloc[0],
                    "reference_file_name": reference_file_name,
                    "peptide_counts": peptide_count,
                    "feature_counts": feature_count,
                    "global_qvalue": float(group["global_qvalue"].iloc[0]),
                    "intensities": intensities,
                    "additional_intensities": additional_intensities,
                    "peptides": peptides,
                    "anchor_protein": anchor_protein.split(";")[0],
                    "is_decoy": int(group["is_decoy"].iloc[0]),
                    "contaminant": 0,  # mzTab doesn't have contaminant info
                    "additional_scores": extra_scores,
                }
            )

        return result

    def get_sql_batchs(self, file_num: int = 2):

        logger = logging.getLogger("quantmsio.core.quantms.pg")

        try:
            query = """
                SELECT DISTINCT reference_file_name
                FROM processed_msstats_with_pg
                ORDER BY reference_file_name
            """
            files_df = self._indexer._duckdb.execute(query).df()
            unique_files = files_df["reference_file_name"].tolist()

            for i in range(0, len(unique_files), file_num):
                yield unique_files[i : i + file_num]

        except Exception as e:
            logger.error(f"Error getting file batches: {e}")
            return

    def get_sql_batch_data(self, file_batch):

        batch_query = f"""
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
                pg_names,
                gg_accessions,
                peptidoform,
                intensity,
                sample_accession,
                charge,
                "unique"
            FROM processed_msstats_with_pg
            WHERE reference_file_name = ANY({file_batch})
            AND anchor_protein IS NOT NULL
            AND intensity > 0
        """
        batch_data = self._indexer._duckdb.execute(batch_query).df()

        return batch_data
