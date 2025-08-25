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


if __name__ == "__main__":
    unittest.main()
