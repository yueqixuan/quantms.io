"""
IdXML parser for quantmsio package.
This module provides functionality to parse OpenMS IdXML files and convert them to quantms.io PSM format.
"""

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
from pathlib import Path
from typing import Union, Optional, Dict, List
import logging
import re

from quantmsio.core.openms import OpenMSHandler

try:
    import pyopenms as oms

    PYOPENMS_AVAILABLE = True
except ImportError:
    PYOPENMS_AVAILABLE = False
    logging.error(
        "pyopenms is required but not available. Please install pyopenms to use this parser."
    )


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
        if not PYOPENMS_AVAILABLE:
            raise ImportError(
                "pyopenms is required but not available. Please install pyopenms to use this parser."
            )
        self.idxml_path = Path(idxml_path)
        self._mzml_path: Optional[Path] = Path(mzml_path) if mzml_path else None
        self._protein_map = {}
        self._peptide_identifications = []
        self._parse_with_pyopenms()
        if self._mzml_path is not None:
            self._attach_mzml_spectra()

    def _parse_with_pyopenms(self) -> None:
        """Parse IdXML file using pyopenms"""
        try:
            protein_identifications, peptide_identifications = self._load_idxml_file()
            self._parse_proteins(protein_identifications)
            self._parse_peptides(peptide_identifications)
        except Exception as e:
            logging.error(f"Error parsing IdXML file with pyopenms: {e}")
            raise

    def _load_idxml_file(self) -> tuple[list, list]:
        """Load IdXML file and return protein and peptide identifications"""
        protein_identifications = []
        peptide_identifications = []
        oms.IdXMLFile().load(
            str(self.idxml_path), protein_identifications, peptide_identifications
        )
        return protein_identifications, peptide_identifications

    def _parse_proteins(self, protein_identifications: list) -> None:
        """Parse protein identifications"""
        for protein_id in protein_identifications:
            for protein_hit in protein_id.getHits():
                self._process_protein_hit(protein_hit)

    def _process_protein_hit(self, protein_hit) -> None:
        """Process single protein hit"""
        accession = protein_hit.getAccession()
        if accession:
            is_decoy = self._determine_protein_is_decoy(protein_hit, accession)
            self._protein_map[accession] = {
                "is_decoy": is_decoy,
                "accession": accession,
            }

    def _determine_protein_is_decoy(self, protein_hit, accession: str) -> int:
        """Determine if protein is decoy based on metadata or accession prefix"""
        is_decoy = 0
        try:
            target_decoy_value = protein_hit.getMetaValue("target_decoy")
            if target_decoy_value is not None:
                is_decoy = 1 if str(target_decoy_value).lower() == "decoy" else 0
        except (AttributeError, ValueError, TypeError):
            pass
        if is_decoy == 0 and accession.startswith("DECOY_"):
            is_decoy = 1

        return is_decoy

    def _parse_peptides(self, peptide_identifications: list) -> None:
        """Parse peptide identifications"""
        for peptide_id in peptide_identifications:
            self._process_peptide_identification(peptide_id)

    def _process_peptide_identification(self, peptide_id) -> None:
        """Process single peptide identification"""
        mz = peptide_id.getMZ()
        rt = peptide_id.getRT()
        spectrum_ref = self._extract_spectrum_reference(peptide_id)
        scan = self._extract_scan_number(spectrum_ref)
        reference_file_name = self._extract_reference_file_name(spectrum_ref)

        for peptide_hit in peptide_id.getHits():
            peptide_data = self._parse_peptide_hit_pyopenms(
                peptide_hit, mz, rt, scan, reference_file_name
            )
            if peptide_data:
                self._peptide_identifications.append(peptide_data)

    def _extract_spectrum_reference(self, peptide_id) -> str:
        """Extract spectrum reference from peptide identification"""
        try:
            return peptide_id.getSpectrumReference()
        except (AttributeError, ValueError, TypeError):
            try:
                return peptide_id.getMetaValue("spectrum_reference")
            except (AttributeError, ValueError, TypeError):
                return ""

    def _parse_peptide_hit_pyopenms(
        self,
        peptide_hit,
        mz: float,
        rt: float,
        scan: str,
        reference_file_name: str,
    ) -> Optional[Dict]:
        """Parse single peptide hit using pyopenms"""
        try:
            sequence = peptide_hit.getSequence().toString()
            if not sequence:
                return None

            charge = peptide_hit.getCharge()
            protein_accessions = self._extract_protein_accessions(peptide_hit)
            is_decoy = self._determine_is_decoy(peptide_hit, sequence)
            modifications = self._parse_modifications(sequence)
            additional_scores, q_value, posterior_error_probability = (
                self._extract_scores(peptide_hit)
            )
            calculated_mz = self._calculate_theoretical_mz(sequence, charge)
            clean_sequence = self._clean_sequence(sequence)

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
                "reference_file_name": reference_file_name,
                "scan": scan,
                "q_value": q_value,
            }

        except Exception as e:
            logging.warning(f"Error parsing peptide hit with pyopenms: {e}")
            return None

    def _extract_protein_accessions(self, peptide_hit) -> List[str]:
        """Extract protein accessions from peptide hit"""
        protein_accessions = []
        for protein_ref in peptide_hit.getPeptideEvidences():
            protein_accessions.append(protein_ref.getProteinAccession())
        return protein_accessions

    def _determine_is_decoy(self, peptide_hit, sequence: str) -> int:
        """Determine if peptide is decoy based on metadata or sequence prefix"""
        is_decoy = 0
        try:
            target_decoy_value = peptide_hit.getMetaValue("target_decoy")
            if target_decoy_value is not None:
                is_decoy = 1 if str(target_decoy_value).lower() == "decoy" else 0
        except (AttributeError, ValueError, TypeError):
            pass

        if is_decoy == 0 and sequence.startswith("DECOY_"):
            is_decoy = 1

        return is_decoy

    def _extract_scores(
        self, peptide_hit
    ) -> tuple[List[Dict], Optional[float], Optional[float]]:
        """Extract additional scores, q-value, and posterior error probability"""
        additional_scores = []
        q_value = None
        posterior_error_probability = None

        try:
            q_value = self._extract_meta_value(peptide_hit, "q-value")
            posterior_error_probability = self._extract_meta_value(
                peptide_hit, "Posterior Error Probability_score"
            )
            additional_scores = self._extract_additional_scores(peptide_hit)
        except (AttributeError, ValueError, TypeError):
            pass

        return additional_scores, q_value, posterior_error_probability

    def _extract_meta_value(self, peptide_hit, key: str) -> Optional[float]:
        """Extract and convert meta value to float"""
        try:
            value = peptide_hit.getMetaValue(key)
            if value is not None:
                return float(value)
        except (AttributeError, ValueError, TypeError):
            pass
        return None

    def _extract_additional_scores(self, peptide_hit) -> List[Dict]:
        """Extract additional scores from peptide hit"""
        additional_scores = []

        try:
            meta_keys = peptide_hit.getMetaValueKeys()
            important_scores = [
                "Luciphor_pep_score",
                "Luciphor_global_flr",
                "Luciphor_local_flr",
                "consensus_support",
                "search_engine_sequence",
                "target_decoy",
            ]

            # Extract important scores first
            for score_name in important_scores:
                score_value = self._extract_meta_value(peptide_hit, score_name)
                if score_value is not None:
                    additional_scores.append(
                        {"score_name": score_name, "score_value": score_value}
                    )

            # Extract other scores
            excluded_scores = [
                "q-value",
                "Posterior Error Probability_score",
            ] + important_scores
            for key in meta_keys:
                if key not in excluded_scores:
                    score_value = self._extract_meta_value(peptide_hit, key)
                    if score_value is not None:
                        additional_scores.append(
                            {"score_name": key, "score_value": score_value}
                        )
        except (AttributeError, ValueError, TypeError):
            # Fallback to known scores only
            known_scores = [
                "Luciphor_pep_score",
                "Luciphor_global_flr",
                "Luciphor_local_flr",
                "consensus_support",
                "search_engine_sequence",
                "target_decoy",
            ]

            for score_name in known_scores:
                score_value = self._extract_meta_value(peptide_hit, score_name)
                if score_value is not None:
                    additional_scores.append(
                        {"score_name": score_name, "score_value": score_value}
                    )

        return additional_scores

    def _clean_sequence(self, sequence: str) -> str:
        """Clean sequence by removing modification annotations"""
        if "(" in sequence:
            return re.sub(r"[\(\[].*?[\)\]]", "", sequence)
        return sequence

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
                    "fields": [
                        {
                            "position": aa_pos,
                            "localization_probability": 1.0,
                        }
                    ],
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

    def _extract_reference_file_name(self, spectrum_ref: str) -> str:
        """
        Extract reference file name from spectrum reference string.
        If no file name is found in spectrum_ref, use the mzML file name if available.

        :param spectrum_ref: Spectrum reference string
        :return: Reference file name
        """
        file_match = re.search(r"file=([^,\s]+)", spectrum_ref)
        if file_match:
            return file_match.group(1)

        if self._mzml_path:
            return self._mzml_path.stem

        return ""

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
