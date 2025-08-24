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
    logging.warning("pyopenms不可用，将回退到XML解析")


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
        
        if PYOPENMS_AVAILABLE:
            self._parse_with_pyopenms()
        else:
            self._parse_with_xml_fallback()
            
        if self._mzml_path is not None:
            self._attach_mzml_spectra()

    def _parse_with_pyopenms(self) -> None:
        """使用pyopenms解析IdXML文件"""
        try:
            # 使用OpenMS的IdXMLFile解析器
            protein_identifications = []
            peptide_identifications = []
            
            oms.IdXMLFile().load(str(self.idxml_path), protein_identifications, peptide_identifications)
            
            # 解析蛋白质信息
            for protein_id in protein_identifications:
                for protein_hit in protein_id.getHits():
                    accession = protein_hit.getAccession()
                    if accession:
                        is_decoy = 0
                        # 检查是否为decoy
                        if accession.startswith("DECOY_"):
                            is_decoy = 1
                        
                        self._protein_map[accession] = {
                            "is_decoy": is_decoy,
                            "accession": accession,
                        }
            
            # 解析肽段信息
            for peptide_id in peptide_identifications:
                mz = peptide_id.getMZ()
                rt = peptide_id.getRT()
                
                # 尝试获取光谱引用
                spectrum_ref = ""
                try:
                    # OpenMS API中获取光谱引用的方法
                    if hasattr(peptide_id, 'getSpectrumReference'):
                        spectrum_ref = peptide_id.getSpectrumReference()
                    elif hasattr(peptide_id, 'getMetaValue'):
                        spectrum_ref = peptide_id.getMetaValue("spectrum_reference")
                except:
                    pass
                
                scan = self._extract_scan_number(spectrum_ref)
                reference_file_name = self._extract_reference_file_name(spectrum_ref)
                
                for peptide_hit in peptide_id.getHits():
                    peptide_data = self._parse_peptide_hit_pyopenms(
                        peptide_hit, mz, rt, scan, spectrum_ref, reference_file_name
                    )
                    if peptide_data:
                        self._peptide_identifications.append(peptide_data)
                        
        except Exception as e:
            logging.error(f"使用pyopenms解析IdXML文件时出错: {e}")
            logging.info("回退到XML解析方法")
            self._parse_with_xml_fallback()

    def _parse_with_xml_fallback(self) -> None:
        """回退到XML解析方法"""
        import xml.etree.ElementTree as ET
        
        try:
            tree = ET.parse(self.idxml_path)
            root = tree.getroot()

            # 解析蛋白质信息
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

            # 解析肽段信息
            for peptide_id in root.findall(".//PeptideIdentification"):
                mz = float(peptide_id.get("MZ", 0))
                rt = float(peptide_id.get("RT", 0))
                spectrum_ref = peptide_id.get("spectrum_reference", "")

                scan = self._extract_scan_number(spectrum_ref)
                reference_file_name = self._extract_reference_file_name(spectrum_ref)

                for peptide_hit in peptide_id.findall(".//PeptideHit"):
                    peptide_data = self._parse_peptide_hit_xml(
                        peptide_hit, mz, rt, scan, spectrum_ref, reference_file_name
                    )
                    if peptide_data:
                        self._peptide_identifications.append(peptide_data)

        except ET.ParseError as e:
            logging.error(f"Error parsing IdXML file: {e}")
            raise
        except Exception as e:
            logging.error(
                f"Unexpected error while parsing protein identifications: {e}"
            )
            raise

    def _parse_peptide_hit_pyopenms(
        self,
        peptide_hit,
        mz: float,
        rt: float,
        scan: str,
        spectrum_ref: str,
        reference_file_name: str,
    ) -> Optional[Dict]:
        """使用pyopenms解析单个肽段命中"""
        try:
            sequence = peptide_hit.getSequence().toString()
            if not sequence:
                return None

            charge = peptide_hit.getCharge()
            score = peptide_hit.getScore()

            # 获取蛋白质引用
            protein_accessions = []
            for protein_ref in peptide_hit.getPeptideEvidences():
                protein_accessions.append(protein_ref.getProteinAccession())

            # 检查是否为decoy
            is_decoy = 0
            # 这里需要根据具体的OpenMS API来获取decoy信息
            # 暂时使用序列检查
            if sequence.startswith("DECOY_"):
                is_decoy = 1

            modifications = self._parse_modifications(sequence)

            additional_scores = []
            q_value = None
            posterior_error_probability = None

            # 尝试获取额外的分数信息
            try:
                # 尝试从OpenMS API获取q-value
                if hasattr(peptide_hit, 'getMetaValue'):
                    try:
                        q_value = peptide_hit.getMetaValue("q-value")
                        if q_value is not None:
                            q_value = float(q_value)
                    except:
                        pass
                    
                    try:
                        posterior_error_probability = peptide_hit.getMetaValue("Posterior Error Probability_score")
                        if posterior_error_probability is not None:
                            posterior_error_probability = float(posterior_error_probability)
                    except:
                        pass
                    
                    # 获取其他分数 - 改进版本
                    try:
                        # 尝试获取所有元数据键
                        meta_keys = peptide_hit.getMetaValueKeys()
                        
                        # 定义已知的重要分数名称，确保它们被包含
                        important_scores = [
                            "Luciphor_pep_score", "Luciphor_global_flr", "Luciphor_local_flr",
                            "consensus_support", "search_engine_sequence", "target_decoy"
                        ]
                        
                        # 首先添加已知的重要分数
                        for score_name in important_scores:
                            try:
                                score_value = peptide_hit.getMetaValue(score_name)
                                if score_value is not None:
                                    additional_scores.append(
                                        {"score_name": score_name, "score_value": float(score_value)}
                                    )
                            except:
                                pass
                        
                        # 然后添加其他所有分数
                        for key in meta_keys:
                            if key not in ["q-value", "Posterior Error Probability_score"] + important_scores:
                                try:
                                    score_value = peptide_hit.getMetaValue(key)
                                    if score_value is not None:
                                        additional_scores.append(
                                            {"score_name": key, "score_value": float(score_value)}
                                        )
                                except:
                                    pass
                    except:
                        # 如果getMetaValueKeys()失败，尝试直接获取已知的分数
                        known_scores = [
                            "Luciphor_pep_score", "Luciphor_global_flr", "Luciphor_local_flr",
                            "consensus_support", "search_engine_sequence", "target_decoy"
                        ]
                        
                        for score_name in known_scores:
                            try:
                                score_value = peptide_hit.getMetaValue(score_name)
                                if score_value is not None:
                                    additional_scores.append(
                                        {"score_name": score_name, "score_value": float(score_value)}
                                    )
                            except:
                                pass
            except:
                pass

            calculated_mz = self._calculate_theoretical_mz(sequence, charge)

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
                "reference_file_name": reference_file_name,
                "scan": scan,
                "q_value": q_value,
            }

        except Exception as e:
            logging.warning(f"Error parsing peptide hit with pyopenms: {e}")
            return None

    def _parse_peptide_hit_xml(
        self,
        peptide_hit,
        mz: float,
        rt: float,
        scan: str,
        spectrum_ref: str,
        reference_file_name: str,
    ) -> Optional[Dict]:
        """使用XML解析单个肽段命中（回退方法）"""
        try:
            sequence = peptide_hit.get("sequence", "")
            if not sequence:
                return None

            charge = int(peptide_hit.get("charge", 1))

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
                "reference_file_name": reference_file_name,
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

    def _extract_reference_file_name(self, spectrum_ref: str) -> str:
        """
        Extract reference file name from spectrum reference string.
        If no file name is found in spectrum_ref, use the mzML file name if available.

        :param spectrum_ref: Spectrum reference string
        :return: Reference file name
        """
        # Try to extract file name from spectrum reference
        file_match = re.search(r"file=([^,\s]+)", spectrum_ref)
        if file_match:
            return file_match.group(1)
        
        # If no file name in spectrum_ref and mzML path is available, use mzML file name
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
