"""Tests for mzIdentML parser module."""

import tempfile
import shutil
import unittest
from pathlib import Path

import pandas as pd

from qpx.core.mzidentml import MzIdentML


class TestMzIdentML(unittest.TestCase):
    """Test cases for MzIdentML parser."""

    def setUp(self):
        """Set up test fixtures."""
        self.test_data_dir = Path(__file__).parents[2] / "examples" / "mzidentml"
        self.test_mzid_file = self.test_data_dir / "test_sample.mzid"
        self.temp_dir = Path(tempfile.mkdtemp())

        if self.test_mzid_file.exists():
            self.parser = MzIdentML(self.test_mzid_file)
        else:
            self.parser = None

    def tearDown(self):
        """Clean up test fixtures."""
        if self.temp_dir.exists():
            shutil.rmtree(self.temp_dir)

    def test_initialization(self):
        """Test MzIdentML initialization and basic properties."""
        if self.parser is None:
            self.skipTest("Test mzIdentML file not found")

        self.assertIsInstance(self.parser, MzIdentML)
        self.assertEqual(self.parser.mzid_path, self.test_mzid_file)
        self.assertGreater(self.parser.get_psm_count(), 0)

    def test_psm_count(self):
        """Test PSM count matches expected value."""
        if self.parser is None:
            self.skipTest("Test mzIdentML file not found")

        # Test file has 3 PSMs
        self.assertEqual(self.parser.get_psm_count(), 3)

    def test_dataframe_conversion(self):
        """Test DataFrame conversion and required columns."""
        if self.parser is None:
            self.skipTest("Test mzIdentML file not found")

        df = self.parser.to_dataframe()
        self.assertIsInstance(df, pd.DataFrame)
        self.assertFalse(df.empty)

        # Check required columns exist
        required_columns = [
            "sequence",
            "peptidoform",
            "modifications",
            "precursor_charge",
            "posterior_error_probability",
            "is_decoy",
            "calculated_mz",
            "observed_mz",
            "additional_scores",
            "protein_accessions",
            "reference_file_name",
            "scan",
            "q_value",
        ]

        for col in required_columns:
            self.assertIn(col, df.columns, f"Missing column: {col}")

    def test_sequence_parsing(self):
        """Test peptide sequence parsing."""
        if self.parser is None:
            self.skipTest("Test mzIdentML file not found")

        df = self.parser.to_dataframe()

        # Check sequences from test file
        sequences = df["sequence"].tolist()
        self.assertIn("PEPTIDER", sequences)
        self.assertIn("TESTPEPTIDE", sequences)
        self.assertIn("OXIDIZEDM", sequences)

    def test_charge_state_parsing(self):
        """Test charge state parsing."""
        if self.parser is None:
            self.skipTest("Test mzIdentML file not found")

        df = self.parser.to_dataframe()

        # Test file has charges 2 and 3
        charges = df["precursor_charge"].tolist()
        self.assertIn(2, charges)
        self.assertIn(3, charges)

    def test_modification_parsing(self):
        """Test modification parsing."""
        if self.parser is None:
            self.skipTest("Test mzIdentML file not found")

        df = self.parser.to_dataframe()

        # Find row with Phospho modification
        phospho_row = df[df["sequence"] == "TESTPEPTIDE"]
        self.assertFalse(phospho_row.empty)

        mods = phospho_row.iloc[0]["modifications"]
        self.assertIsNotNone(mods)
        self.assertIsInstance(mods, list)
        self.assertEqual(len(mods), 1)
        self.assertEqual(mods[0]["name"], "Phospho")

        # Find row with Oxidation modification
        ox_row = df[df["sequence"] == "OXIDIZEDM"]
        self.assertFalse(ox_row.empty)

        mods = ox_row.iloc[0]["modifications"]
        self.assertIsNotNone(mods)
        self.assertEqual(mods[0]["name"], "Oxidation")

    def test_protein_accession_parsing(self):
        """Test protein accession parsing."""
        if self.parser is None:
            self.skipTest("Test mzIdentML file not found")

        df = self.parser.to_dataframe()

        # Check protein accessions
        for _, row in df.iterrows():
            accessions = row["protein_accessions"]
            self.assertIsNotNone(accessions)
            self.assertIsInstance(accessions, list)
            self.assertGreater(len(accessions), 0)

    def test_scan_number_extraction(self):
        """Test scan number extraction from spectrum ID."""
        if self.parser is None:
            self.skipTest("Test mzIdentML file not found")

        df = self.parser.to_dataframe()

        scans = df["scan"].tolist()
        self.assertIn("100", scans)
        self.assertIn("200", scans)
        self.assertIn("300", scans)

    def test_score_extraction(self):
        """Test score extraction from PSMs."""
        if self.parser is None:
            self.skipTest("Test mzIdentML file not found")

        df = self.parser.to_dataframe()

        # Check q-values are parsed
        for _, row in df.iterrows():
            q_value = row["q_value"]
            self.assertIsNotNone(q_value)
            self.assertGreater(q_value, 0)
            self.assertLess(q_value, 1)

    def test_parquet_conversion(self):
        """Test parquet file conversion."""
        if self.parser is None:
            self.skipTest("Test mzIdentML file not found")

        output_path = self.temp_dir / "test_output.parquet"
        self.parser.to_parquet(output_path)

        self.assertTrue(output_path.exists())
        self.assertGreater(output_path.stat().st_size, 0)

        # Verify parquet can be read back
        df = pd.read_parquet(output_path)
        self.assertEqual(len(df), 3)

    def test_decoy_detection(self):
        """Test decoy PSM detection."""
        if self.parser is None:
            self.skipTest("Test mzIdentML file not found")

        df = self.parser.to_dataframe()

        # Test file has no decoys
        for _, row in df.iterrows():
            self.assertEqual(row["is_decoy"], 0)

    def test_mz_values(self):
        """Test m/z value parsing."""
        if self.parser is None:
            self.skipTest("Test mzIdentML file not found")

        df = self.parser.to_dataframe()

        for _, row in df.iterrows():
            # Check experimental m/z
            self.assertIsNotNone(row["observed_mz"])
            self.assertGreater(row["observed_mz"], 0)

            # Check calculated m/z
            self.assertIsNotNone(row["calculated_mz"])
            self.assertGreater(row["calculated_mz"], 0)


class TestMzIdentMLEdgeCases(unittest.TestCase):
    """Test edge cases for MzIdentML parser."""

    def test_nonexistent_file(self):
        """Test handling of non-existent file."""
        with self.assertRaises(Exception):
            MzIdentML("/nonexistent/path/file.mzid")

    def test_empty_dataframe_handling(self):
        """Test handling when no PSMs are found."""
        # This test would require a specially crafted empty mzid file
        pass


class TestScanNumberExtraction(unittest.TestCase):
    """Test scan number extraction from various native ID formats."""

    def setUp(self):
        """Set up test fixtures."""
        self.test_data_dir = Path(__file__).parents[2] / "examples" / "mzidentml"
        self.test_mzid_file = self.test_data_dir / "test_sample.mzid"

        if self.test_mzid_file.exists():
            self.parser = MzIdentML(self.test_mzid_file)
        else:
            self.parser = None

    def test_thermo_scan_format(self):
        """Test Thermo scan=XXX format."""
        if self.parser is None:
            self.skipTest("Test mzIdentML file not found")

        result = self.parser._extract_scan_number(
            "controllerType=0 controllerNumber=1 scan=12345"
        )
        self.assertEqual(result, "12345")

    def test_waters_cycle_format(self):
        """Test Waters/Agilent cycle=XXX format."""
        if self.parser is None:
            self.skipTest("Test mzIdentML file not found")

        result = self.parser._extract_scan_number(
            "sample=1 period=1 cycle=1055 experiment=4"
        )
        self.assertEqual(result, "1055")

    def test_index_format(self):
        """Test index=XXX format."""
        if self.parser is None:
            self.skipTest("Test mzIdentML file not found")

        result = self.parser._extract_scan_number("index=500")
        self.assertEqual(result, "500")

    def test_spectrum_format(self):
        """Test spectrum=XXX format."""
        if self.parser is None:
            self.skipTest("Test mzIdentML file not found")

        result = self.parser._extract_scan_number("spectrum=999")
        self.assertEqual(result, "999")

    def test_trailing_number_format(self):
        """Test format with trailing number."""
        if self.parser is None:
            self.skipTest("Test mzIdentML file not found")

        result = self.parser._extract_scan_number("some_prefix_123")
        self.assertEqual(result, "123")

    def test_unknown_format_returns_original(self):
        """Test unknown format returns original string."""
        if self.parser is None:
            self.skipTest("Test mzIdentML file not found")

        result = self.parser._extract_scan_number("unknown_format")
        self.assertEqual(result, "unknown_format")


class TestMultiMzMLSupport(unittest.TestCase):
    """Test multi-mzML folder support."""

    def test_mzml_folder_parameter(self):
        """Test MzIdentML accepts mzml_folder parameter."""
        test_data_dir = Path(__file__).parents[2] / "examples" / "mzidentml"
        test_mzid_file = test_data_dir / "test_sample.mzid"

        if not test_mzid_file.exists():
            self.skipTest("Test mzIdentML file not found")

        # Should not raise error with mzml_folder parameter
        parser = MzIdentML(
            test_mzid_file,
            mzml_folder=test_data_dir,  # Use test dir (no mzML files expected)
            spectral_data=True,
        )
        self.assertIsNotNone(parser)
        self.assertEqual(parser._mzml_folder, test_data_dir)

    def test_mzml_file_and_folder_mutually_exclusive(self):
        """Test that mzml_path and mzml_folder can both be set (CLI validates)."""
        test_data_dir = Path(__file__).parents[2] / "examples" / "mzidentml"
        test_mzid_file = test_data_dir / "test_sample.mzid"

        if not test_mzid_file.exists():
            self.skipTest("Test mzIdentML file not found")

        # Both can be set at parser level (CLI validates mutual exclusion)
        parser = MzIdentML(
            test_mzid_file,
            mzml_path=test_mzid_file,  # Just a path, won't be used
            mzml_folder=test_data_dir,
            spectral_data=False,  # Disable to avoid errors
        )
        self.assertIsNotNone(parser)


if __name__ == "__main__":
    unittest.main()
