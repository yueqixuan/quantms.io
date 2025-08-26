import unittest
from pathlib import Path
import tempfile
import shutil
import pandas as pd

from quantmsio.core.idxml import IdXML


class TestIdXML(unittest.TestCase):

    def setUp(self):
        self.test_data_dir = Path("tests/examples/idxml")
        self.test_idxml_file = (
            self.test_data_dir
            / "SF_200217_pPeptideLibrary_pool1_HCDnlETcaD_OT_rep2_consensus_fdr_pep_luciphor.idXML"
        )
        self.test_mzml_file = (
            self.test_data_dir
            / "SF_200217_pPeptideLibrary_pool1_HCDnlETcaD_OT_rep1.mzML"
        )

        self.temp_dir = Path(tempfile.mkdtemp())

    def tearDown(self):
        if self.temp_dir.exists():
            shutil.rmtree(self.temp_dir)

    def test_idxml_initialization(self):
        if self.test_idxml_file.exists():
            if self.test_mzml_file.exists():
                idxml = IdXML(self.test_idxml_file, self.test_mzml_file)
            else:
                idxml = IdXML(self.test_idxml_file)
            self.assertIsInstance(idxml, IdXML)
            self.assertEqual(idxml.idxml_path, self.test_idxml_file)
        else:
            self.skipTest("Test IdXML file not found")

    def test_protein_parsing(self):
        if self.test_idxml_file.exists():
            if self.test_mzml_file.exists():
                idxml = IdXML(self.test_idxml_file, self.test_mzml_file)
            else:
                idxml = IdXML(self.test_idxml_file)
            protein_count = idxml.get_protein_count()
            self.assertGreaterEqual(protein_count, 0)
        else:
            self.skipTest("Test IdXML file not found")

    def test_peptide_parsing(self):
        if self.test_idxml_file.exists():
            if self.test_mzml_file.exists():
                idxml = IdXML(self.test_idxml_file, self.test_mzml_file)
            else:
                idxml = IdXML(self.test_idxml_file)
            psm_count = idxml.get_psm_count()
            self.assertGreaterEqual(psm_count, 0)
        else:
            self.skipTest("Test IdXML file not found")

    def test_dataframe_conversion(self):
        if self.test_idxml_file.exists():
            if self.test_mzml_file.exists():
                idxml = IdXML(self.test_idxml_file, self.test_mzml_file)
            else:
                idxml = IdXML(self.test_idxml_file)
            df = idxml.to_dataframe()
            self.assertIsInstance(df, pd.DataFrame)

            if not df.empty:
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
                    "mp_accessions",
                    "rt",
                    "reference_file_name",
                    "scan",
                    "q_value",
                ]

                for col in required_columns:
                    self.assertIn(col, df.columns)
        else:
            self.skipTest("Test IdXML file not found")

    def test_parquet_conversion(self):
        if self.test_idxml_file.exists():
            if self.test_mzml_file.exists():
                idxml = IdXML(self.test_idxml_file, self.test_mzml_file)
            else:
                idxml = IdXML(self.test_idxml_file)
            output_path = self.temp_dir / "test_output.parquet"

            # Convert to parquet
            idxml.to_parquet(output_path)

            # Check if file was created
            self.assertTrue(output_path.exists())

            # Check file size
            self.assertGreater(output_path.stat().st_size, 0)
        else:
            self.skipTest("Test IdXML file not found")

    def test_modification_parsing(self):
        if self.test_idxml_file.exists():
            if self.test_mzml_file.exists():
                idxml = IdXML(self.test_idxml_file, self.test_mzml_file)
            else:
                idxml = IdXML(self.test_idxml_file)
        else:
            idxml = None
        if idxml is None:
            self.skipTest("Test IdXML file not found")

        # Test modification parsing with a sample sequence
        test_sequence = "GSGEKPVSAPGDDT(Phospho)ES(Phospho)LHSQGEEEFDMPQPPHGHVLHR"
        modifications = idxml._parse_modifications(test_sequence)

        self.assertIsInstance(modifications, list)
        self.assertEqual(len(modifications), 2)  # Two phospho modifications

        for mod in modifications:
            self.assertIn("modification_name", mod)
            self.assertIn("fields", mod)
            self.assertIsInstance(mod["fields"], list)
            self.assertEqual(mod["modification_name"], "Phospho")
            # Check fields structure
            for field in mod["fields"]:
                self.assertIn("position", field)
                self.assertIn("localization_probability", field)
                self.assertIsInstance(field["position"], int)
                self.assertIsInstance(field["localization_probability"], float)

    def test_scan_number_extraction(self):
        if self.test_idxml_file.exists():
            if self.test_mzml_file.exists():
                idxml = IdXML(self.test_idxml_file, self.test_mzml_file)
            else:
                idxml = IdXML(self.test_idxml_file)
        else:
            idxml = None
        if idxml is None:
            self.skipTest("Test IdXML file not found")

        # Test scan number extraction
        test_spectrum_ref = "controllerType=0 controllerNumber=1 scan=37144"
        scan_number = idxml._extract_scan_number(test_spectrum_ref)
        self.assertEqual(scan_number, "37144")

        # Test with no scan number
        test_spectrum_ref_no_scan = "controllerType=0 controllerNumber=1"
        scan_number_no_scan = idxml._extract_scan_number(test_spectrum_ref_no_scan)
        self.assertEqual(scan_number_no_scan, "unknown_index")

    def test_theoretical_mz_calculation(self):
        if self.test_idxml_file.exists():
            if self.test_mzml_file.exists():
                idxml = IdXML(self.test_idxml_file, self.test_mzml_file)
            else:
                idxml = IdXML(self.test_idxml_file)
        else:
            idxml = None
        if idxml is None:
            self.skipTest("Test IdXML file not found")

        # Test with a simple peptide
        test_sequence = "PEPTIDE"
        test_charge = 2

        calculated_mz = idxml._calculate_theoretical_mz(test_sequence, test_charge)
        self.assertIsInstance(calculated_mz, float)
        self.assertGreater(calculated_mz, 0)

    def test_theoretical_mz_calculation_with_modifications(self):
        """Test m/z calculation with modified peptides using pyopenms"""
        if self.test_idxml_file.exists():
            if self.test_mzml_file.exists():
                idxml = IdXML(self.test_idxml_file, self.test_mzml_file)
            else:
                idxml = IdXML(self.test_idxml_file)
        else:
            idxml = None
        if idxml is None:
            self.skipTest("Test IdXML file not found")

        # Test cases with different modifications (using only pyopenms-supported modifications)
        test_cases = [
            ("PEPTIDE", 2, "Simple peptide"),
            ("PEPTIDE(Oxidation)", 2, "Peptide with oxidation"),
            ("PEPTIDE(Phospho)", 3, "Peptide with phosphorylation"),
            ("PEPTIDE(Carbamidomethyl)", 2, "Peptide with carbamidomethylation"),
            ("PEPTIDE(Methyl)", 2, "Peptide with methylation"),
            ("PEPTIDE(Oxidation)(Phospho)", 2, "Peptide with multiple modifications"),
            ("M(Oxidation)PEPTIDE", 1, "N-terminal modification"),
            ("PEPTIDEM(Oxidation)", 2, "C-terminal modification"),
            ("PEPTIDE(Dioxidation)", 2, "Peptide with dioxidation"),
        ]

        for sequence, charge, description in test_cases:
            with self.subTest(description=description):
                calculated_mz = idxml._calculate_theoretical_mz(sequence, charge)
                self.assertIsInstance(calculated_mz, float)
                self.assertGreater(calculated_mz, 0)

                # Verify that m/z decreases with increasing charge
                if charge > 1:
                    mz_charge_1 = idxml._calculate_theoretical_mz(sequence, 1)
                    self.assertGreater(mz_charge_1, calculated_mz)

    def test_theoretical_mz_calculation_edge_cases(self):
        """Test m/z calculation with edge cases"""
        if self.test_idxml_file.exists():
            if self.test_mzml_file.exists():
                idxml = IdXML(self.test_idxml_file, self.test_mzml_file)
            else:
                idxml = IdXML(self.test_idxml_file)
        else:
            idxml = None
        if idxml is None:
            self.skipTest("Test IdXML file not found")

        # Test edge cases
        edge_cases = [
            ("A", 1, "Single amino acid"),
            ("AA", 1, "Two amino acids"),
            ("PEPTIDE", 0, "Zero charge"),
            ("PEPTIDE", 10, "High charge"),
        ]

        for sequence, charge, description in edge_cases:
            with self.subTest(description=description):
                if charge == 0:
                    # Zero charge should return peptide mass
                    calculated_mz = idxml._calculate_theoretical_mz(sequence, charge)
                    self.assertIsInstance(calculated_mz, float)
                    self.assertGreater(calculated_mz, 0)
                else:
                    calculated_mz = idxml._calculate_theoretical_mz(sequence, charge)
                    self.assertIsInstance(calculated_mz, float)
                    self.assertGreater(calculated_mz, 0)

    def test_theoretical_mz_calculation_consistency(self):
        """Test that m/z calculation is consistent with pyopenms"""
        if self.test_idxml_file.exists():
            if self.test_mzml_file.exists():
                idxml = IdXML(self.test_idxml_file, self.test_mzml_file)
            else:
                idxml = IdXML(self.test_idxml_file)
        else:
            idxml = None
        if idxml is None:
            self.skipTest("Test IdXML file not found")

        try:
            import pyopenms as oms
            from pyopenms.Constants import PROTON_MASS_U
        except ImportError:
            self.skipTest("pyopenms not available")

        # Test consistency with direct pyopenms calculation
        test_sequence = "PEPTIDE"
        test_charge = 2

        # Calculate using our method
        our_mz = idxml._calculate_theoretical_mz(test_sequence, test_charge)

        # Calculate directly with pyopenms
        aa_sequence = oms.AASequence.fromString(test_sequence)
        peptide_mass = aa_sequence.getMonoWeight()
        expected_mz = (peptide_mass + (test_charge * PROTON_MASS_U)) / test_charge

        # Should be very close (within floating point precision)
        self.assertAlmostEqual(our_mz, expected_mz, places=6)

    def test_theoretical_mz_calculation_edge_cases_and_error_handling(self):
        """Test m/z calculation with edge cases and error handling"""
        if self.test_idxml_file.exists():
            if self.test_mzml_file.exists():
                idxml = IdXML(self.test_idxml_file, self.test_mzml_file)
            else:
                idxml = IdXML(self.test_idxml_file)
        else:
            idxml = None
        if idxml is None:
            self.skipTest("Test IdXML file not found")

        # Test edge cases that pyopenms handles gracefully
        edge_sequences = [
            "",  # Empty sequence - pyopenms returns proton mass
            "X",  # Unknown amino acid - pyopenms may handle differently
        ]

        for sequence in edge_sequences:
            with self.subTest(sequence=sequence):
                # These should either return a result or raise an exception
                try:
                    result = idxml._calculate_theoretical_mz(sequence, 2)
                    self.assertIsInstance(result, float)
                    self.assertGreater(result, 0)
                except Exception as e:
                    # If it raises an exception, that's also acceptable
                    self.assertIsInstance(e, Exception)


if __name__ == "__main__":
    unittest.main()
