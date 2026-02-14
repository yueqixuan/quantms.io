"""Dataset class — the central entry point for opening QPX datasets."""

from pathlib import Path

from qpx.core.engine import DuckDBEngine
from qpx.core.convert import QueryResult
from qpx.core.models.base import ValidationResult, ValidationIssue
from qpx.core.data import (
    Feature, PSM, PG, MzSpectra,
    Sample, Run, DatasetMeta, Ontology, Provenance,
    PeptideProteinMap,
    BaseStructure,
)


class Dataset:
    """
    A QPX dataset — a directory of related Parquet/H5AD files.

    Usage:
        ds = qpx.Dataset("PXD014414/")
        ds.feature.filter("charge > 2").to_df()
        ds.pg.join(ds.run).to_df()

    Partial access:
        ds = qpx.Dataset("path/", structures=["feature", "sample", "run"])
    """

    # Data structure registry: name → (Class, file_suffix)
    _STRUCTURE_REGISTRY = {
        "psm":        (PSM,          ".psm.parquet"),
        "feature":    (Feature,      ".feature.parquet"),
        "pg":         (PG,           ".pg.parquet"),
        "mz":         (MzSpectra,    ".mz.parquet"),
        "sample":     (Sample,       ".sample.parquet"),
        "run":        (Run,          ".run.parquet"),
        "dataset":    (DatasetMeta,  ".dataset.parquet"),
        "ontology":   (Ontology,     ".ontology.parquet"),
        "provenance": (Provenance,   ".provenance.parquet"),
        "peptide_protein_map": (PeptideProteinMap, ".peptide_protein_map.parquet"),
    }

    def __init__(
        self,
        path: str | Path,
        structures: list[str] | None = None,
        duckdb_memory: str = "16GB",
        duckdb_threads: int = 4,
        s3_config: dict | None = None,
    ):
        self._is_s3 = isinstance(path, str) and path.startswith("s3://")
        self.path = path if self._is_s3 else Path(path)
        self._engine = DuckDBEngine(
            memory_limit=duckdb_memory,
            threads=duckdb_threads,
            s3_config=s3_config if self._is_s3 else None,
        )

        self._structures: dict[str, BaseStructure] = {}
        self._discover_and_register(structures)

    def _discover_and_register(self, requested: list[str] | None):
        """Scan directory, find QPX files, register as DuckDB tables.

        Checks for single Parquet files first, then falls back to
        Hive-partitioned directories (e.g., feature/ with run_file_name= subdirs).
        For S3 paths, attempts to register each structure via S3 glob.
        """
        if self._is_s3:
            self._discover_s3(requested)
        else:
            self._discover_local(requested)

    def _discover_s3(self, requested: list[str] | None):
        """Register structures from S3 path."""
        for name, (cls, suffix) in self._STRUCTURE_REGISTRY.items():
            if requested and name not in requested:
                continue
            s3_glob = f"{self.path.rstrip('/')}/*{suffix}"
            try:
                self._engine.register_s3_parquet(name, s3_glob)
                self._structures[name] = cls(
                    engine=self._engine,
                    table_name=name,
                    file_path=f"{self.path}/{name}",
                )
            except Exception:
                pass  # Structure not present in S3

    def _discover_local(self, requested: list[str] | None):
        """Register structures from local filesystem."""
        for name, (cls, suffix) in self._STRUCTURE_REGISTRY.items():
            if requested and name not in requested:
                continue
            # Check for single file first
            matches = sorted(self.path.glob(f"*{suffix}"))
            if matches:
                file_path = matches[0]  # Take first match
                self._engine.register_parquet(name, file_path)
                self._structures[name] = cls(
                    engine=self._engine,
                    table_name=name,
                    file_path=file_path,
                )
            else:
                # Check for Hive-partitioned directory
                part_dir = self.path / name
                if part_dir.is_dir() and list(part_dir.glob("**/*.parquet")):
                    self._engine.register_partitioned_parquet(name, part_dir)
                    self._structures[name] = cls(
                        engine=self._engine,
                        table_name=name,
                        file_path=part_dir,
                    )

    # --- Data structure accessors (lazy, return None if not present) ---
    @property
    def psm(self) -> PSM | None:
        return self._structures.get("psm")

    @property
    def feature(self) -> Feature | None:
        return self._structures.get("feature")

    @property
    def pg(self) -> PG | None:
        return self._structures.get("pg")

    @property
    def mz(self) -> MzSpectra | None:
        return self._structures.get("mz")

    @property
    def sample(self) -> Sample | None:
        return self._structures.get("sample")

    @property
    def run(self) -> Run | None:
        return self._structures.get("run")

    @property
    def dataset_meta(self) -> DatasetMeta | None:
        return self._structures.get("dataset")

    @property
    def ontology(self) -> Ontology | None:
        return self._structures.get("ontology")

    @property
    def provenance(self) -> Provenance | None:
        return self._structures.get("provenance")

    @property
    def peptide_protein_map(self) -> PeptideProteinMap | None:
        return self._structures.get("peptide_protein_map")

    # --- View accessors (cached) ---
    @property
    def protein_view(self):
        if not hasattr(self, "_protein_view"):
            from qpx.views.api import ProteinView
            self._protein_view = ProteinView(self)
        return self._protein_view

    @property
    def peptide_view(self):
        if not hasattr(self, "_peptide_view"):
            from qpx.views.api import PeptideView
            self._peptide_view = PeptideView(self)
        return self._peptide_view

    @property
    def identification_summary(self):
        if not hasattr(self, "_identification_summary"):
            from qpx.views.api import IdentificationSummaryView
            self._identification_summary = IdentificationSummaryView(self)
        return self._identification_summary

    @property
    def run_summary(self):
        if not hasattr(self, "_run_summary"):
            from qpx.views.api import RunSummaryView
            self._run_summary = RunSummaryView(self)
        return self._run_summary

    @property
    def modification_view(self):
        if not hasattr(self, "_modification_view"):
            from qpx.views.api import ModificationView
            self._modification_view = ModificationView(self)
        return self._modification_view

    @property
    def qc_view(self):
        if not hasattr(self, "_qc_view"):
            from qpx.views.api import QualityControlView
            self._qc_view = QualityControlView(self)
        return self._qc_view

    # --- Cross-structure queries ---
    def sql(self, query: str) -> QueryResult:
        """Execute arbitrary SQL across registered structures."""
        return QueryResult(self._engine.execute(query))

    @property
    def available_structures(self) -> list[str]:
        return list(self._structures.keys())

    # --- Analysis helpers ---

    def _abundance_sql(self, level: str) -> str:
        """Build the long-form abundance SQL for a given level.

        Returns SQL that produces columns: sample_accession, feature_id, intensity.
        """
        if level == "protein":
            if self.pg is None or self.run is None:
                raise ValueError(
                    "level='protein' requires pg and run structures."
                )
            return """
            SELECT rs.sample_accession,
                   pg.anchor_protein AS feature_id,
                   i.intensity
            FROM pg,
                 run r,
                 UNNEST(r.samples) AS _t1(rs),
                 UNNEST(pg.intensities) AS _t2(i)
            WHERE pg.run_file_name = r.run_file_name
              AND i.label = rs.label
              AND pg.is_decoy = false
            """
        elif level == "peptide":
            if self.feature is None or self.run is None:
                raise ValueError(
                    "level='peptide' requires feature and run structures."
                )
            return """
            SELECT rs.sample_accession,
                   f.sequence AS feature_id,
                   SUM(i.intensity) AS intensity
            FROM feature f,
                 run r,
                 UNNEST(r.samples) AS _t1(rs),
                 UNNEST(f.intensities) AS _t2(i)
            WHERE f.run_file_name = r.run_file_name
              AND i.label = rs.label
              AND f.is_decoy = false
            GROUP BY rs.sample_accession, f.sequence
            """
        else:
            raise ValueError(f"level must be 'protein' or 'peptide', got '{level}'")

    def abundance(self, level: str = "protein") -> QueryResult:
        """Lazy long-form abundance query — scalable for large datasets.

        Returns a QueryResult that stays lazy until you materialize it.
        DuckDB handles memory management, so this works on datasets that
        don't fit in memory.

        Parameters
        ----------
        level : str
            "protein" (uses pg + run) or "peptide" (uses feature + run).

        Returns
        -------
        QueryResult
            Long-form table with columns: sample_accession, feature_id, intensity.
            Call .to_df(), .to_arrow(), .to_polars(), or iterate row-by-row.

        Examples
        --------
        Lazy iteration (constant memory):
            for row in ds.abundance("protein"):
                sample, protein, intensity = row

        Materialize when data fits in memory:
            df = ds.abundance("protein").to_df()

        Write directly to Parquet (out-of-core):
            ds.abundance("protein").to_arrow()  # Arrow is more compact than pandas
        """
        sql = self._abundance_sql(level)
        return QueryResult(self._engine.execute(sql))

    def design_matrix(
        self,
        level: str = "protein",
        value_col: str = "intensity",
        fillna: float | None = 0.0,
        output_path: str | Path | None = None,
    ) -> "pd.DataFrame | Path":
        """Pivot abundance data into a samples-by-features matrix.

        For small-to-medium datasets, returns a pandas DataFrame.
        For large datasets, pass *output_path* to write the pivot directly
        to Parquet via DuckDB (out-of-core, no full pandas materialization).

        For streaming access to large data without pivoting, use
        :meth:`abundance` instead.

        Parameters
        ----------
        level : str
            "protein" (uses pg + run) or "peptide" (uses feature + run).
        value_col : str
            Which intensity column to pivot (default: "intensity").
        fillna : float | None
            Value to fill missing cells. Use *None* to keep NaN.
            Ignored when *output_path* is set.
        output_path : str | Path | None
            If set, writes the pivoted matrix directly to this Parquet path
            using DuckDB PIVOT (out-of-core). Returns the Path instead of a
            DataFrame.

        Returns
        -------
        pd.DataFrame | Path
            DataFrame (in-memory) or Path (when output_path is set).
        """
        import pandas as pd

        sql = self._abundance_sql(level)

        if output_path is not None:
            # Out-of-core: DuckDB PIVOT → Parquet, never touches pandas
            output_path = Path(output_path)
            pivot_sql = f"""
            COPY (
                PIVOT ({sql}) ON feature_id USING SUM(intensity)
            ) TO '{output_path}' (FORMAT PARQUET)
            """
            self._engine.execute(pivot_sql)
            return output_path

        # In-memory path for small/medium datasets
        df = self._engine.execute(sql).fetchdf()

        if df.empty:
            return pd.DataFrame()

        matrix = df.pivot_table(
            index="sample_accession",
            columns="feature_id",
            values="intensity",
            aggfunc="sum",
        )
        matrix.columns.name = None

        if fillna is not None:
            matrix = matrix.fillna(fillna)

        return matrix

    # --- Validation ---
    def validate(
        self, structures: list[str] | None = None
    ) -> dict[str, ValidationResult]:
        """Validate the dataset or specific structures against their schemas.

        Parameters
        ----------
        structures : list[str] | None
            Structure names to validate.  If *None*, validates all available.

        Returns
        -------
        dict[str, ValidationResult]
            Mapping of structure name to its validation result.
        """
        targets = structures or self.available_structures
        results: dict[str, ValidationResult] = {}
        for name in targets:
            struct = self._structures.get(name)
            if struct is None:
                result = ValidationResult(structure=name)
                result.issues.append(ValidationIssue(
                    structure=name,
                    check="missing_structure",
                    severity="error",
                    column=None,
                    message=f"Structure '{name}' not found in dataset at {self.path}",
                ))
                results[name] = result
            else:
                results[name] = struct.validate()
        return results

    # --- Integrity ---
    def compute_integrity(self) -> dict:
        """Compute checksums, row counts, and file sizes for all Parquet files.

        Returns a dict of integrity fields suitable for writing to dataset.parquet.
        """
        import hashlib
        from datetime import datetime, timezone
        import pyarrow.parquet as pq

        checksums, row_counts, sizes = {}, {}, {}
        for f in sorted(self.path.glob("*.parquet")):
            name = f.name
            sizes[name] = f.stat().st_size
            sha = hashlib.sha256(f.read_bytes()).hexdigest()
            checksums[name] = sha
            try:
                row_counts[name] = pq.read_metadata(f).num_rows
            except Exception:
                row_counts[name] = -1

        return {
            "file_checksums": checksums,
            "file_row_counts": row_counts,
            "file_sizes_bytes": sizes,
            "total_structures": len(checksums),
            "packaged_at": datetime.now(timezone.utc).isoformat(),
        }

    def verify_integrity(self) -> dict[str, list[str]]:
        """Verify dataset files against stored integrity data.

        Returns dict with 'errors' and 'warnings' lists.
        """
        import hashlib

        errors, warnings = [], []
        if self.dataset_meta is None:
            errors.append("No dataset.parquet found — cannot verify integrity")
            return {"errors": errors, "warnings": warnings}

        meta_df = self.dataset_meta.to_df()
        if meta_df.empty or "file_checksums" not in meta_df.columns:
            warnings.append("No integrity data stored in dataset.parquet")
            return {"errors": errors, "warnings": warnings}

        stored_checksums = meta_df["file_checksums"].iloc[0]
        if not isinstance(stored_checksums, dict):
            warnings.append("file_checksums is null")
            return {"errors": errors, "warnings": warnings}

        # Skip dataset.parquet itself — writing integrity changes its own checksum
        dataset_suffix = self._STRUCTURE_REGISTRY["dataset"][1]
        for name, expected_sha in stored_checksums.items():
            if name.endswith(dataset_suffix):
                continue
            fpath = self.path / name
            if not fpath.exists():
                errors.append(f"Missing file: {name}")
                continue
            actual_sha = hashlib.sha256(fpath.read_bytes()).hexdigest()
            if actual_sha != expected_sha:
                errors.append(f"Checksum mismatch: {name}")

        return {"errors": errors, "warnings": warnings}

    # --- Write-back ---
    # Writer registry: name → (WriterClassName, file_suffix)
    _WRITER_REGISTRY = {
        "psm":        ("PsmWriter",       ".psm.parquet"),
        "feature":    ("FeatureWriter",    ".feature.parquet"),
        "pg":         ("PgWriter",         ".pg.parquet"),
        "mz":         ("MzWriter",         ".mz.parquet"),
        "sample":     ("SampleWriter",     ".sample.parquet"),
        "run":        ("RunWriter",        ".run.parquet"),
        "dataset":    ("DatasetWriter",    ".dataset.parquet"),
        "ontology":   ("OntologyWriter",   ".ontology.parquet"),
        "provenance": ("ProvenanceWriter", ".provenance.parquet"),
        "peptide_protein_map": ("PeptideProteinMapWriter", ".peptide_protein_map.parquet"),
    }

    def save_structure(
        self,
        data,
        structure: str,
        prefix: str | None = None,
    ) -> Path:
        """Write validated Parquet back into the dataset directory.

        Parameters
        ----------
        data : list[dict] | pd.DataFrame | pa.Table
            Records to write. Schema-validated by the QPX writer.
        structure : str
            Structure name (e.g., "feature", "pg", "sample").
        prefix : str | None
            File name prefix. Defaults to the dataset directory name.

        Returns
        -------
        Path
            Path to the written Parquet file.
        """
        import pyarrow as pa

        if structure not in self._WRITER_REGISTRY:
            raise ValueError(
                f"Unknown structure '{structure}'. "
                f"Valid: {list(self._WRITER_REGISTRY.keys())}"
            )

        writer_name, suffix = self._WRITER_REGISTRY[structure]

        import qpx.writers as writers_mod
        writer_cls = getattr(writers_mod, writer_name)

        file_prefix = prefix or self.path.name
        output_path = self.path / f"{file_prefix}{suffix}"

        with writer_cls(output_path) as writer:
            if isinstance(data, list):
                writer.write_batch(data)
            elif isinstance(data, pa.Table):
                # Cast to the writer's schema to handle nullability/metadata differences
                casted = data.cast(writer.arrow_schema)
                writer.write_table(casted)
            else:
                writer.write_dataframe(data)

        return output_path

    def save_anndata(self, adata, name: str) -> Path:
        """Write an AnnData object to the dataset directory as .h5ad.

        Parameters
        ----------
        adata : anndata.AnnData
            The AnnData object to save.
        name : str
            Base name for the file (without extension), e.g., "de_results".

        Returns
        -------
        Path
            Path to the written .h5ad file.
        """
        output_path = self.path / f"{name}.h5ad"
        adata.write_h5ad(output_path)
        return output_path

    # --- External registration (for mokume integration) ---
    def register_external(self, name: str, file_path: str | Path) -> None:
        """Register an external Parquet file (e.g., mokume output) in the DuckDB engine."""
        self._engine.register_parquet(name, file_path)

    # --- Lifecycle ---
    def refresh(self) -> None:
        """Re-scan the dataset directory for new or updated files.

        Clears all cached view instances so next access gets fresh data.
        Call this after writing new structures or AnnData files.
        """
        for attr in [
            "_protein_view", "_peptide_view", "_identification_summary",
            "_run_summary", "_modification_view", "_qc_view",
        ]:
            if hasattr(self, attr):
                delattr(self, attr)

        self._structures.clear()
        self._discover_and_register(None)

    def close(self):
        self._engine.close()

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()

    def __repr__(self):
        structs = ", ".join(self.available_structures)
        return f"Dataset('{self.path}', structures=[{structs}])"
