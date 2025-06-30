#!/usr/bin/env python3
"""
Comprehensive Test Suite for QuantMS Conversions
Tests PSM, Feature, and Protein Group conversions for LFQ and TMT datasets
"""

import os
import subprocess
import tempfile
import time
from pathlib import Path

import pyarrow.parquet as pq
import pytest


def get_workspace_root():
    """Get the workspace root directory"""
    return Path(__file__).parent.parent


def get_test_files():
    """Get paths to test files"""
    workspace_root = get_workspace_root()

    # Dataset paths
    lfq_dataset = workspace_root / "tests/examples/quantms/dda-lfq-full"
    tmt_dataset = workspace_root / "tests/examples/quantms/dda-plex-full"

    # File paths
    lfq_files = {
        "mztab": lfq_dataset / "PXD007683-LFQ.sdrf_openms_design_openms.mzTab.gz",
        "sdrf": lfq_dataset / "PXD007683-LFQ.sdrf.tsv",
        "msstats": lfq_dataset / "PXD007683-LFQ.sdrf_openms_design_msstats_in.csv.gz",
    }

    tmt_files = {
        "mztab": tmt_dataset / "PXD007683TMT.sdrf_openms_design_openms.mzTab.gz",
        "sdrf": tmt_dataset / "PXD007683-TMT.sdrf.tsv",
        "msstats": tmt_dataset / "PXD007683TMT.sdrf_openms_design_msstats_in.csv.gz",
    }

    return lfq_files, tmt_files


def run_command(cmd, description, workspace_root):
    """Run a command and return success/failure with timing"""
    print(f"\n[RUN] Running: {description}")
    print(f"Command: {' '.join(cmd)}")

    start_time = time.time()
    try:
        # Ensure UTF-8 encoding for subprocess
        env = os.environ.copy()
        env["PYTHONIOENCODING"] = "utf-8"

        result = subprocess.run(
            cmd,
            cwd=workspace_root,
            capture_output=True,
            text=True,
            env=env,
            encoding="utf-8",
            errors="replace",
            timeout=600,  # 10 minute timeout (TMT protein groups needs ~6 minutes)
        )
        end_time = time.time()
        duration = end_time - start_time

        if result.returncode == 0:
            print(f"[OK] Success in {duration:.2f}s")
            return True, duration, result.stdout, result.stderr
        else:
            print(f"[FAIL] Failed in {duration:.2f}s")
            print(f"Error: {result.stderr}")
            return False, duration, result.stdout, result.stderr

    except subprocess.TimeoutExpired:
        print("[TIMEOUT] Timeout after 10 minutes")
        return False, 600, "", "Timeout"
    except Exception as e:
        print(f"[ERROR] Exception: {str(e)}")
        return False, 0, "", str(e)


def analyze_output_file(file_path):
    """Analyze a parquet output file"""
    if not file_path.exists():
        return None

    try:
        df = pq.read_table(file_path).to_pandas()
        size_mb = file_path.stat().st_size / (1024 * 1024)

        return {
            "rows": len(df),
            "columns": len(df.columns),
            "size_mb": size_mb,
            "column_names": list(df.columns)[:10],  # First 10 columns
        }
    except Exception as e:
        return {"error": str(e)}


@pytest.mark.integration
def test_lfq_psm_conversion():
    """Test LFQ PSM conversion"""
    workspace_root = get_workspace_root()
    lfq_files, _ = get_test_files()

    # Skip if files don't exist
    if not lfq_files["mztab"].exists():
        pytest.skip(f"Test file not found: {lfq_files['mztab']}")

    with tempfile.TemporaryDirectory() as temp_dir:
        output_folder = Path(temp_dir)
        output_prefix = "LFQ_psm_test"

        cmd = [
            "python",
            "-m",
            "quantmsio.quantmsioc",
            "convert",
            "quantms-psm",
            "--input-file",
            str(lfq_files["mztab"]),
            "--output-folder",
            str(output_folder),
            "--output-prefix",
            output_prefix,
        ]

        success, duration, stdout, stderr = run_command(
            cmd, "LFQ PSM Conversion", workspace_root
        )

        # Find output file (with UUID)
        output_files = list(output_folder.glob(f"{output_prefix}-*.psm.parquet"))
        assert (
            len(output_files) == 1
        ), f"Expected 1 output file, found {len(output_files)}"

        output_file = output_files[0]
        analysis = analyze_output_file(output_file)

        assert success, f"Command failed: {stderr}"
        assert analysis is not None, "Failed to analyze output file"
        assert "error" not in analysis, f"Analysis error: {analysis.get('error')}"
        assert analysis["rows"] > 0, "No data rows in output"
        assert (
            analysis["columns"] >= 15
        ), f"Expected at least 15 columns, got {analysis['columns']}"


@pytest.mark.integration
def test_lfq_feature_conversion():
    """Test LFQ Feature conversion"""
    workspace_root = get_workspace_root()
    lfq_files, _ = get_test_files()

    # Skip if files don't exist
    required_files = ["mztab", "sdrf", "msstats"]
    for file_type in required_files:
        if not lfq_files[file_type].exists():
            pytest.skip(f"Test file not found: {lfq_files[file_type]}")

    with tempfile.TemporaryDirectory() as temp_dir:
        output_folder = Path(temp_dir)
        output_prefix = "LFQ_feature_test"

        cmd = [
            "python",
            "-m",
            "quantmsio.quantmsioc",
            "convert",
            "quantms-feature",
            "--input-file",
            str(lfq_files["mztab"]),
            "--sdrf-file",
            str(lfq_files["sdrf"]),
            "--msstats-file",
            str(lfq_files["msstats"]),
            "--output-folder",
            str(output_folder),
            "--output-prefix",
            output_prefix,
        ]

        success, duration, stdout, stderr = run_command(
            cmd, "LFQ Feature Conversion", workspace_root
        )

        # Find output file (with UUID)
        output_files = list(output_folder.glob(f"{output_prefix}-*.feature.parquet"))
        assert (
            len(output_files) == 1
        ), f"Expected 1 output file, found {len(output_files)}"

        output_file = output_files[0]
        analysis = analyze_output_file(output_file)

        assert success, f"Command failed: {stderr}"
        assert analysis is not None, "Failed to analyze output file"
        assert "error" not in analysis, f"Analysis error: {analysis.get('error')}"
        assert analysis["rows"] > 0, "No data rows in output"
        assert (
            analysis["columns"] >= 20
        ), f"Expected at least 20 columns, got {analysis['columns']}"


@pytest.mark.integration
def test_lfq_protein_groups_conversion():
    """Test LFQ Protein Groups conversion"""
    workspace_root = get_workspace_root()
    lfq_files, _ = get_test_files()

    # Skip if files don't exist
    required_files = ["mztab", "sdrf", "msstats"]
    for file_type in required_files:
        if not lfq_files[file_type].exists():
            pytest.skip(f"Test file not found: {lfq_files[file_type]}")

    with tempfile.TemporaryDirectory() as temp_dir:
        output_folder = Path(temp_dir)
        output_prefix = "LFQ_pg_test"

        cmd = [
            "python",
            "-m",
            "quantmsio.quantmsioc",
            "convert",
            "quantms-pg",
            "--input-file",
            str(lfq_files["mztab"]),
            "--msstats-file",
            str(lfq_files["msstats"]),
            "--sdrf-file",
            str(lfq_files["sdrf"]),
            "--output-folder",
            str(output_folder),
            "--output-prefix",
            output_prefix,
        ]

        success, duration, stdout, stderr = run_command(
            cmd, "LFQ Protein Groups Conversion", workspace_root
        )

        # Find output file (with UUID)
        output_files = list(output_folder.glob(f"{output_prefix}-*.pg.parquet"))
        assert (
            len(output_files) == 1
        ), f"Expected 1 output file, found {len(output_files)}"

        output_file = output_files[0]
        analysis = analyze_output_file(output_file)

        assert success, f"Command failed: {stderr}"
        assert analysis is not None, "Failed to analyze output file"
        assert "error" not in analysis, f"Analysis error: {analysis.get('error')}"
        assert analysis["rows"] > 0, "No data rows in output"
        assert (
            analysis["columns"] >= 10
        ), f"Expected at least 10 columns, got {analysis['columns']}"


@pytest.mark.integration
def test_tmt_psm_conversion():
    """Test TMT PSM conversion"""
    workspace_root = get_workspace_root()
    _, tmt_files = get_test_files()

    # Skip if files don't exist
    if not tmt_files["mztab"].exists():
        pytest.skip(f"Test file not found: {tmt_files['mztab']}")

    with tempfile.TemporaryDirectory() as temp_dir:
        output_folder = Path(temp_dir)
        output_prefix = "TMT_psm_test"

        cmd = [
            "python",
            "-m",
            "quantmsio.quantmsioc",
            "convert",
            "quantms-psm",
            "--input-file",
            str(tmt_files["mztab"]),
            "--output-folder",
            str(output_folder),
            "--output-prefix",
            output_prefix,
        ]

        success, duration, stdout, stderr = run_command(
            cmd, "TMT PSM Conversion", workspace_root
        )

        # Find output file (with UUID)
        output_files = list(output_folder.glob(f"{output_prefix}-*.psm.parquet"))
        assert (
            len(output_files) == 1
        ), f"Expected 1 output file, found {len(output_files)}"

        output_file = output_files[0]
        analysis = analyze_output_file(output_file)

        assert success, f"Command failed: {stderr}"
        assert analysis is not None, "Failed to analyze output file"
        assert "error" not in analysis, f"Analysis error: {analysis.get('error')}"
        assert analysis["rows"] > 0, "No data rows in output"
        assert (
            analysis["columns"] >= 15
        ), f"Expected at least 15 columns, got {analysis['columns']}"


@pytest.mark.integration
def test_tmt_feature_conversion():
    """Test TMT Feature conversion"""
    workspace_root = get_workspace_root()
    _, tmt_files = get_test_files()

    # Skip if files don't exist
    required_files = ["mztab", "sdrf", "msstats"]
    for file_type in required_files:
        if not tmt_files[file_type].exists():
            pytest.skip(f"Test file not found: {tmt_files[file_type]}")

    with tempfile.TemporaryDirectory() as temp_dir:
        output_folder = Path(temp_dir)
        output_prefix = "TMT_feature_test"

        cmd = [
            "python",
            "-m",
            "quantmsio.quantmsioc",
            "convert",
            "quantms-feature",
            "--input-file",
            str(tmt_files["mztab"]),
            "--sdrf-file",
            str(tmt_files["sdrf"]),
            "--msstats-file",
            str(tmt_files["msstats"]),
            "--output-folder",
            str(output_folder),
            "--output-prefix",
            output_prefix,
        ]

        success, duration, stdout, stderr = run_command(
            cmd, "TMT Feature Conversion", workspace_root
        )

        # Find output file (with UUID)
        output_files = list(output_folder.glob(f"{output_prefix}-*.feature.parquet"))
        assert (
            len(output_files) == 1
        ), f"Expected 1 output file, found {len(output_files)}"

        output_file = output_files[0]
        analysis = analyze_output_file(output_file)

        assert success, f"Command failed: {stderr}"
        assert analysis is not None, "Failed to analyze output file"
        assert "error" not in analysis, f"Analysis error: {analysis.get('error')}"
        assert analysis["rows"] > 0, "No data rows in output"
        assert (
            analysis["columns"] >= 20
        ), f"Expected at least 20 columns, got {analysis['columns']}"


@pytest.mark.integration
def test_tmt_protein_groups_conversion():
    """Test TMT Protein Groups conversion"""
    workspace_root = get_workspace_root()
    _, tmt_files = get_test_files()

    # Skip if files don't exist
    required_files = ["mztab", "sdrf", "msstats"]
    for file_type in required_files:
        if not tmt_files[file_type].exists():
            pytest.skip(f"Test file not found: {tmt_files[file_type]}")

    with tempfile.TemporaryDirectory() as temp_dir:
        output_folder = Path(temp_dir)
        output_prefix = "TMT_pg_test"

        cmd = [
            "python",
            "-m",
            "quantmsio.quantmsioc",
            "convert",
            "quantms-pg",
            "--input-file",
            str(tmt_files["mztab"]),
            "--msstats-file",
            str(tmt_files["msstats"]),
            "--sdrf-file",
            str(tmt_files["sdrf"]),
            "--output-folder",
            str(output_folder),
            "--output-prefix",
            output_prefix,
        ]

        success, duration, stdout, stderr = run_command(
            cmd, "TMT Protein Groups Conversion", workspace_root
        )

        # Find output file (with UUID)
        output_files = list(output_folder.glob(f"{output_prefix}-*.pg.parquet"))
        assert (
            len(output_files) == 1
        ), f"Expected 1 output file, found {len(output_files)}"

        output_file = output_files[0]
        analysis = analyze_output_file(output_file)

        assert success, f"Command failed: {stderr}"
        assert analysis is not None, "Failed to analyze output file"
        assert "error" not in analysis, f"Analysis error: {analysis.get('error')}"
        assert analysis["rows"] > 0, "No data rows in output"
        assert (
            analysis["columns"] >= 10
        ), f"Expected at least 10 columns, got {analysis['columns']}"
