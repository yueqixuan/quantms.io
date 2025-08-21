"""
IdXML parser for quantmsio package.
This module provides functionality to parse OpenMS IdXML files and convert them to quantms.io PSM format.
"""

import xml.etree.ElementTree as ET
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
from pathlib import Path
from typing import Union, Optional, Dict, List, Tuple
import logging
import re

from quantmsio.core.common import PSM_SCHEMA
from quantmsio.core.format import PSM_FIELDS
from quantmsio.utils.pride_utils import generate_scan_number
from quantmsio.core.openms import OpenMSHandler


class IdXML:
    """
    Parser for OpenMS IdXML files.

    This class provides functionality to parse IdXML files and convert them to quantms.io PSM format.
    IdXML is an OpenMS format for storing peptide and protein identifications.
    """

    def __init__(
        self, idxml_path: Union[Path, str], mzml_path: Optional[Union[Path, str]] = None
    ):
        """
        Initialize the IdXML parser.

        :param idxml_path: Path to the IdXML file
        :param mzml_path: Optional path to the mzML file for attaching spectra
        """
        self.idxml_path = Path(idxml_path)
        self._mzml_path: Optional[Path] = Path(mzml_path) if mzml_path else None
        self._protein_map = {}
        self._peptide_identifications = []
        self._parse_protein_identifications()
        self._parse_peptide_identifications()
        if self._mzml_path is not None:
            self._attach_mzml_spectra()

    def _parse_protein_identifications(self) -> None:
        """
        Parse protein identifications from the IdXML file.
        Extracts protein accessions and their target/decoy status.
        """
        try:
            tree = ET.parse(self.idxml_path)
            root = tree.getroot()

            for protein_id in root.findall(".//ProteinIdentification"):
                for protein_hit in protein_id.findall(".//ProteinHit"):
                    accession = protein_hit.get("accession", "")
                    if accession:
                        is_decoy = 0
                        target_decoy_param = protein_hit.find(
                            './/UserParam[@name="target_decoy"]'
                        )
                        if target_decoy_param is not None:
                            is_decoy = (
                                1 if target_decoy_param.get("value") == "decoy" else 0
                            )
                        elif accession.startswith("DECOY_"):
                            is_decoy = 1

                        self._protein_map[accession] = {
                            "is_decoy": is_decoy,
                            "accession": accession,
                        }

        except ET.ParseError as e:
            logging.error(f"Error parsing IdXML file: {e}")
            raise
        except Exception as e:
            logging.error(
                f"Unexpected error while parsing protein identifications: {e}"
            )
            raise

    def _parse_peptide_identifications(self) -> None:
        """
        Parse peptide identifications from the IdXML file.
        Extracts peptide hits with their associated information.
        """
        try:
            tree = ET.parse(self.idxml_path)
            root = tree.getroot()

            for peptide_id in root.findall(".//PeptideIdentification"):
                mz = float(peptide_id.get("MZ", 0))
                rt = float(peptide_id.get("RT", 0))
                spectrum_ref = peptide_id.get("spectrum_reference", "")

                scan = self._extract_scan_number(spectrum_ref)

                for peptide_hit in peptide_id.findall(".//PeptideHit"):
                    peptide_data = self._parse_peptide_hit(
                        peptide_hit, mz, rt, scan, spectrum_ref
                    )
                    if peptide_data:
                        self._peptide_identifications.append(peptide_data)

        except ET.ParseError as e:
            logging.error(f"Error parsing IdXML file: {e}")
            raise
        except Exception as e:
            logging.error(
                f"Unexpected error while parsing peptide identifications: {e}"
            )
            raise

    def _parse_peptide_hit(
        self,
        peptide_hit: ET.Element,
        mz: float,
        rt: float,
        scan: str,
        spectrum_ref: str,
    ) -> Optional[Dict]:
        """
        Parse individual peptide hit information.

        :param peptide_hit: XML element containing peptide hit data
        :param mz: Precursor m/z value
        :param rt: Retention time value
        :param scan: Scan number from spectrum reference
        :return: Dictionary containing parsed peptide hit data
        """
        try:
            sequence = peptide_hit.get("sequence", "")
            if not sequence:
                return None

            charge = int(peptide_hit.get("charge", 1))
            score = float(peptide_hit.get("score", 0))

            protein_refs = peptide_hit.get("protein_refs", "")
            modifications = self._parse_modifications(sequence)

            additional_scores = []
            q_value = None
            posterior_error_probability = None

            for user_param in peptide_hit.findall(".//UserParam"):
                param_name = user_param.get("name", "")
                param_value = user_param.get("value", "")

                if param_name == "q-value":
                    try:
                        q_value = float(param_value)
                    except ValueError:
                        pass
                elif param_name == "Posterior Error Probability_score":
                    try:
                        posterior_error_probability = float(param_value)
                    except ValueError:
                        pass
                else:
                    try:
                        score_value = float(param_value)
                        additional_scores.append(
                            {"score_name": param_name, "score_value": score_value}
                        )
                    except ValueError:
                        pass

            is_decoy = 0
            target_decoy_param = peptide_hit.find('.//UserParam[@name="target_decoy"]')
            if target_decoy_param is not None:
                is_decoy = 1 if target_decoy_param.get("value") == "decoy" else 0

            calculated_mz = self._calculate_theoretical_mz(sequence, charge)

            protein_accessions = []
            if protein_refs:
                protein_accessions = [ref.strip() for ref in protein_refs.split(",")]

            clean_sequence = sequence
            if "(" in sequence:
                clean_sequence = re.sub(r"[\(\[].*?[\)\]]", "", sequence)

            return {
                "sequence": clean_sequence,
                "peptidoform": sequence,
                "modifications": modifications,
                "precursor_charge": charge,
                "posterior_error_probability": posterior_error_probability,
                "is_decoy": is_decoy,
                "calculated_mz": calculated_mz,
                "observed_mz": mz,
                "additional_scores": additional_scores,
                "mp_accessions": protein_accessions,
                "rt": rt,
                "reference_file_name": spectrum_ref,
                "scan": scan,
                "q_value": q_value,
            }

        except Exception as e:
            logging.warning(f"Error parsing peptide hit: {e}")
            return None

    def _parse_modifications(self, sequence: str) -> List[Dict]:
        """
        Parse modifications from peptide sequence.

        :param sequence: Peptide sequence with modification annotations
        :return: List of modification dictionaries
        """
        modifications = []
        aa_positions = []
        aa_count = 0

        i = 0
        while i < len(sequence):
            if sequence[i] == "(":
                j = i + 1
                while j < len(sequence) and sequence[j] != ")":
                    j += 1
                i = j + 1 if j < len(sequence) else len(sequence)
            elif sequence[i].isalpha():
                aa_count += 1
                aa_positions.append((i, aa_count))
                i += 1
            else:
                i += 1

        mod_pattern = r"\(([^)]+)\)"
        matches = list(re.finditer(mod_pattern, sequence))

        for match in matches:
            mod_name = match.group(1)
            position = match.start()

            aa_pos = 0
            for orig_index, aa_position in aa_positions:
                if orig_index < position:
                    aa_pos = aa_position
                else:
                    break

            modifications.append(
                {
                    "modification_name": mod_name,
                    "position": aa_pos,
                    "localization_probability": 1.0,
                }
            )

        return modifications

    def _extract_scan_number(self, spectrum_ref: str) -> str:
        """
        Extract scan number from spectrum reference string.

        :param spectrum_ref: Spectrum reference string
        :return: Extracted scan number
        """
        scan_match = re.search(r"scan=(\d+)", spectrum_ref)
        if scan_match:
            return scan_match.group(1)
        return "unknown_index"

    def _calculate_theoretical_mz(self, sequence: str, charge: int) -> float:
        """
        Calculate theoretical m/z for a peptide sequence.
        This is a simplified calculation.

        :param sequence: Peptide sequence
        :param charge: Charge state
        :return: Theoretical m/z value
        """
        aa_masses = {
            "A": 71.03711,
            "R": 156.10111,
            "N": 114.04293,
            "D": 115.02694,
            "C": 103.00919,
            "E": 129.04259,
            "Q": 128.05858,
            "G": 57.02146,
            "H": 137.05891,
            "I": 113.08406,
            "L": 113.08406,
            "K": 128.09496,
            "M": 131.04049,
            "F": 147.06841,
            "P": 97.05276,
            "S": 87.03203,
            "T": 101.04768,
            "W": 186.07931,
            "Y": 163.06333,
            "V": 99.06841,
        }

        peptide_mass = 0
        for aa in sequence:
            if aa in aa_masses:
                peptide_mass += aa_masses[aa]

        peptide_mass += 1.007825 + 17.00274

        if charge > 0:
            return (peptide_mass + (charge - 1) * 1.007825) / charge
        else:
            return peptide_mass

    def to_dataframe(self) -> pd.DataFrame:
        """
        Convert parsed peptide identifications to a pandas DataFrame.

        :return: DataFrame containing PSM data
        """
        if not self._peptide_identifications:
            return pd.DataFrame()

        df = pd.DataFrame(self._peptide_identifications)

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
            "consensus_support",
        ]

        # Spectra fields are optional but add them for schema consistency
        spectra_optional = ["number_peaks", "mz_array", "intensity_array"]
        required_columns.extend(
            [c for c in spectra_optional if c not in required_columns]
        )

        for col in required_columns:
            if col not in df.columns:
                df[col] = None

        return df

    def _attach_mzml_spectra(self) -> None:
        """
        Attach number_peaks, mz_array, intensity_array using the provided mzML file.
        Uses scan number as identifier.
        """
        try:
            handler = OpenMSHandler()
            for psm in self._peptide_identifications:
                scan_value = psm.get("scan")
                if not scan_value or scan_value == "unknown_index":
                    continue
                try:
                    scan_int = int(str(scan_value))
                except ValueError:
                    continue
                try:
                    num_peaks, mzs, intens = handler.get_spectrum_from_scan(
                        str(self._mzml_path), scan_int
                    )
                except Exception as e:
                    logging.warning(f"Error fetching spectrum (scan={scan_int}): {e}")
                    continue
                try:
                    psm["number_peaks"] = int(num_peaks)
                    psm["mz_array"] = [float(x) for x in mzs]
                    psm["intensity_array"] = [float(x) for x in intens]
                except Exception:
                    psm["number_peaks"] = None
                    psm["mz_array"] = None
                    psm["intensity_array"] = None
        except Exception as e:
            logging.warning(f"Failed to attach mzML spectra: {e}")

    def to_parquet(self, output_path: Union[Path, str]) -> None:
        """
        Convert IdXML data to parquet format and save to file.

        :param output_path: Output file path for parquet file
        """
        df = self.to_dataframe()
        if df.empty:
            logging.warning("No peptide identifications found to convert")
            return

        table = pa.Table.from_pandas(df)

        pq.write_table(table, output_path)
        logging.info(f"Successfully converted IdXML to parquet: {output_path}")

    def get_psm_count(self) -> int:
        """
        Get the total number of PSMs found in the IdXML file.

        :return: Number of PSMs
        """
        return len(self._peptide_identifications)

    def get_protein_count(self) -> int:
        """
        Get the total number of proteins found in the IdXML file.

        :return: Number of proteins
        """
        return len(self._protein_map)
