import warnings
import logging
from typing import List, Dict, Any, Optional, Iterator, Union

import numpy as np
import pyopenms as oms
from pyopenms import SpectrumLookup
from pathlib import Path

from qpx.core.common import OPENMS_NAMES_MAP


class OpenMSHandler:
    def __init__(self) -> None:
        self._mzml_exp = None
        self._consensus_xml_path = None
        self._spec_lookup = None
        self.logger = logging.getLogger(f"{__name__}.{self.__class__.__name__}")

    def get_spectrum_from_scan(self, mzml_path: str, scan_number: int):
        """
        Get a spectrum from a mzML file using the scan number
        :param mzml_path: path to the mzML file
        :param scan_number: scan number
        :return: spectrum
        """
        if self._mzml_exp is None:
            self._mzml_exp = oms.MSExperiment()
            oms.MzMLFile().load(mzml_path, self._mzml_exp)
            self._spec_lookup = SpectrumLookup()
            self._spec_lookup.readSpectra(self._mzml_exp, "scan=(?<SCAN>\\d+)")
        try:
            index = self._spec_lookup.findByScanNumber(scan_number)
        except IndexError:
            message = (
                "scan_number" + str(scan_number) + "not found in file: " + mzml_path
            )
            warnings.warn(message, category=None, stacklevel=1, source=None)
            return [], []
        spectrum = self._mzml_exp.getSpectrum(index)
        spectrum_mz, spectrum_intensities = spectrum.get_peaks()
        return len(spectrum_mz), spectrum_mz, spectrum_intensities

    def get_intensity_map(
        self, consensusxml_path: str, experiment_type: str = None
    ) -> dict:
        """
        Get the intensity map from a consensusxml file. The intensity map is a dictionary with the following structure:
        - key: peptide sequence + ":_:" + charge + ":_:" + reference file
        - value: dictionary with the following structure:
          - rt: retention time
          - mz: mass to charge ratio
          - intensity: intensity
        :param consensusxml_path: path to the consensusxml file
        :param experiment_type: experiment type (e.g. lfq, tmt, etc.)
        :return: intensity map
        """
        self._consensus_xml_path = consensusxml_path
        consensus_map = oms.ConsensusMap()
        oms.ConsensusXMLFile().load(self._consensus_xml_path, consensus_map)
        df = consensus_map.get_df()
        df = df[df.sequence != "None"]

        if experiment_type is not None and "LABEL FREE" in experiment_type.upper():
            return self._get_intensity_map_lfq(df)
        elif experiment_type is not None and "TMT" in experiment_type.upper():
            return self._get_intensity_map_tmt_or_itraq(df, experiment_type)
        elif experiment_type is not None and "ITRAQ" in experiment_type.upper():
            return self._get_intensity_map_tmt_or_itraq(df, experiment_type)
        return self._get_intensity_map_lfq(
            df
        )  # If not experiment type is provided, we assume it is label free

    @staticmethod
    def _get_intensity_map_lfq(df):
        """
        Get the intensity map for label free experiments
        :param df: pandas dataframe with the consensusxml data
        :return: intensity map
        """
        peptide_columns = ["sequence", "charge", "RT", "mz", "quality"]
        intensity_columns = [
            column for column in df.columns if column not in peptide_columns
        ]
        intensity_map = {}
        for _, row in df.iterrows():
            for column in intensity_columns:
                if np.float64(row[f"{column}"]) > 0.0:
                    reference_file = column.split(".")[0]
                    key = (
                        row.sequence + ":_:" + str(row.charge) + ":_:" + reference_file
                    )
                    if key not in intensity_map:
                        intensity_map[key] = {
                            "rt": row.RT,
                            "mz": row.mz,
                            "intensity": row[column],
                        }
                    else:
                        if row[column] > intensity_map[key]["intensity"]:
                            intensity_map[key] = {
                                "rt": row.RT,
                                "mz": row.mz,
                                "intensity": row[column],
                            }
        return intensity_map

    @staticmethod
    def _get_intensity_map_tmt_or_itraq(df, experiment_type):
        """
        Get the intensity map for TMT experiments
        :param df: pandas dataframe with the consensusxml data
        :return: intensity map
        """
        peptide_columns = ["sequence", "charge", "RT", "mz", "quality", "file"]
        intensity_columns = [
            column for column in df.columns if column not in peptide_columns
        ]
        intensity_map = {}
        for _, row in df.iterrows():
            for column in intensity_columns:
                if np.float64(row[f"{column}"]) > 0.0:
                    reference_file = row.file.split(".")[0]
                    if "TMT" in experiment_type.upper():
                        channel = (
                            "TMT" + column.split("_")[1]
                        )  # A TMT channel has in consesusXML the following format:
                        # tmt10plex_129N -> TMT129N
                    else:
                        channel = "ITRAQ" + column.split("_")[1]
                    key = (
                        row.sequence
                        + ":_:"
                        + str(row.charge)
                        + ":_:"
                        + reference_file
                        + ":_:"
                        + channel
                    )
                    if key not in intensity_map:
                        intensity_map[key] = {
                            "rt": row.RT,
                            "mz": row.mz,
                            "intensity": row[column],
                            "channel": channel,
                        }
                    else:
                        if row[column] > intensity_map[key]["intensity"]:
                            intensity_map[key] = {
                                "rt": row.RT,
                                "mz": row.mz,
                                "intensity": row[column],
                                "channel": channel,
                            }
        return intensity_map

    def _get_run(
        self,
        protein_ids: List[oms.ProteinIdentification],
        peptide_id: oms.PeptideIdentification,
    ) -> Optional[str]:
        """
        Get run name from idXML using pyopenms.

        If the idXML file contains a merge index, use it to annotate the run name without file
        extension.
        """
        try:
            if peptide_id.metaValueExists("id_merge_index"):
                run = Path(
                    protein_ids[0]
                    .getMetaValue("spectra_data")[
                        peptide_id.getMetaValue("id_merge_index")
                    ]
                    .decode()
                ).stem
            elif protein_ids[0].metaValueExists("spectra_data"):
                run = Path(protein_ids[0].getMetaValue("spectra_data")[0].decode()).stem
            else:
                run = None

            # Convert back to None value
            if run == "None":
                run = None

            return run
        except Exception:
            return None

    def _is_decoy(self, peptide_hit: oms.PeptideHit) -> bool:
        """Check if PSM is target or decoy."""
        if peptide_hit.metaValueExists("target_decoy"):
            return peptide_hit.getMetaValue("target_decoy") == "decoy"
        else:
            # Check if any protein accession starts with DECOY_
            return any(
                acc.getAccession().startswith("DECOY_")
                for acc in peptide_hit.getPeptideEvidences()
            )

    def _extract_hit_metadata(self, peptide_hit: oms.PeptideHit) -> Dict[str, Any]:
        """Extract string metadata from peptide hit."""
        keys = []
        peptide_hit.getKeys(keys)
        metadata = {}

        for key in keys:
            key_str = key.decode() if isinstance(key, bytes) else key
            value = peptide_hit.getMetaValue(key_str)
            # Only include non-float metadata
            if not self._is_float(value):
                metadata[key_str] = value

        return metadata

    @staticmethod
    def _is_float(element: Any) -> bool:
        """Check if element can be coerced to a float."""
        if element is None:
            return False
        try:
            float(element)
            return True
        except ValueError:
            return False

    def _extract_modifications_from_sequence(
        self, peptide_sequence: Union[oms.AASequence, str]
    ) -> List[Dict[str, Any]]:
        """
        Extract modifications from peptide sequence, the sequence can be an string from
        OpenMS AASequence or an string written by OpenMS tools as peptidoforms in mzTab
        for example.

        If the sequence is a string, it is converted to an OpenMS AASequence object.
        The modifications list follows the format of QPX:
        [
            {
                "name": str,           # Name of the modification
                "accession": str,      # Optional accession (e.g. UNIMOD:35)
                "fields": [            # List of modification instances
                    {
                        "position": str,   # Format: "{AA}.{position}"
                                          # AA is amino acid or N-term/C-term
                                          # position is 0 for N-term, 1-based for AA, len+1 for C-term
                        "scores": [        # List of scores for this position
                            {
                                "score_name": str,
                                "score_value": str
                            }
                        ]
                    }
                ]
            }
        ]

        :param peptide_sequence: peptide sequence
        :return: list of modification dictionaries
        """

        if isinstance(peptide_sequence, str):
            peptide_sequence = oms.AASequence.fromString(peptide_sequence)

        modifications = {}
        sequence_length = peptide_sequence.size()

        # Process amino acid modifications
        for i in range(sequence_length):
            residue = peptide_sequence[i]
            if residue.isModified():
                mod = residue.getModification()
                position_str = (
                    f"{residue.getOneLetterCode()}.{i + 1}"  # 1-based position
                )

                if mod.getName() not in modifications:
                    modifications[mod.getName()] = {
                        "name": mod.getName(),
                        "accession": mod.getUniModAccession(),
                        "fields": [],
                    }

                # Add position with its scores
                position_entry = {
                    "position": position_str,
                    "scores": [],  # Can be populated with position-specific scores if available
                }
                modifications[mod.getName()]["fields"].append(position_entry)

        # Add C-terminal modification
        if peptide_sequence.hasCTerminalModification():
            c_term_mod = peptide_sequence.getCTerminalModification()
            c_term_pos = f"C-term.{sequence_length + 1}"
            modifications[c_term_mod.getName()] = {
                "name": c_term_mod.getName(),
                "accession": c_term_mod.getUniModAccession(),
                "fields": [{"position": c_term_pos, "scores": []}],
            }

        # Add N-terminal modification
        if peptide_sequence.hasNTerminalModification():
            n_term_mod = peptide_sequence.getNTerminalModification()
            n_term_pos = "N-term.0"
            modifications[n_term_mod.getName()] = {
                "name": n_term_mod.getName(),
                "accession": n_term_mod.getUniModAccession(),
                "fields": [{"position": n_term_pos, "scores": []}],
            }

        # Convert dictionary of modifications to list
        return list(modifications.values())

    def _extract_additional_scores(self, hit: oms.PeptideHit) -> List[Dict[str, Any]]:
        """
        Extract additional scores from PeptideHit

        :param hit: OpenMS PeptideHit object
        :return: list of score dictionaries
        """
        scores = []

        # Main score
        scores.append(
            {"score_name": "main_score", "score_value": float(hit.getScore())}
        )

        # Additional scores from meta values
        try:
            # Get all meta value keys - PyOpenMS getKeys needs a list to populate
            keys = []
            hit.getKeys(keys)
            for key in keys:
                try:
                    value = hit.getMetaValue(key)
                    if isinstance(value, (int, float)):
                        scores.append({"score_name": key, "score_value": float(value)})
                except (KeyError, TypeError, ValueError, RuntimeError):
                    # Skip if meta value key doesn't exist, value can't be converted,
                    # or there's a runtime error from OpenMS
                    continue
        except (AttributeError, RuntimeError):
            # Skip meta values if hit object doesn't support getKeys() or
            # there's a runtime error from OpenMS
            pass

        return scores


def get_openms_score_name(score_name: str) -> str:
    """
    Get the OpenMS score name from the MzTab score name
    :param score_name: MzTab score name
    :return: OpenMS score name
    """
    return OPENMS_NAMES_MAP.get(score_name, score_name)
