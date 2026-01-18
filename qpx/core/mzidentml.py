"""mzIdentML data processing module for QPX.

This module provides functionality to parse mzIdentML files and convert them
to QPX PSM format. It uses ElementTree XML parsing to handle mzIdentML files
that may contain non-standard modifications which pyopenms cannot parse.
"""

import gzip
import logging
import re
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

from qpx.core.format import PSM_SCHEMA

logging.basicConfig(format="%(asctime)s - %(message)s", level=logging.INFO)
logger = logging.getLogger(__name__)


# ============================================================================
# Constants
# ============================================================================

# mzIdentML namespace URIs for different versions
MZID_NAMESPACES = {
    "1.1": "http://psidev.info/psi/pi/mzIdentML/1.1",
    "1.2": "http://psidev.info/psi/pi/mzIdentML/1.2",
}

# Score name mappings with higher_better indication
# higher_better: True = higher is better, False = lower is better, None = unknown
SCORE_MAPPINGS = {
    "mascot:score": ("mascot_score", True),
    "mascot:expectation value": ("mascot_expect", False),
    "x!tandem:hyperscore": ("hyperscore", True),
    "x\\!tandem:hyperscore": ("hyperscore", True),
    "x!tandem:expect": ("xtandem_expect", False),
    "x\\!tandem:expect": ("xtandem_expect", False),
    "sequest:xcorr": ("xcorr", True),
    "sequest:deltacn": ("deltacn", True),
    "ms-gf:specevalue": ("msgf_specevalue", False),
    "ms-gf:evalue": ("msgf_evalue", False),
    "ms-gf:qvalue": ("msgf_qvalue", False),
    "ms-gf:pepqvalue": ("msgf_pepqvalue", False),
    "comet:xcorr": ("comet_xcorr", True),
    "comet:deltacn": ("comet_deltacn", True),
    "comet:expect": ("comet_expect", False),
    "percolator:score": ("percolator_score", True),
    "percolator:q-value": ("percolator_qvalue", False),
    "percolator:pep": ("percolator_pep", False),
    "peptideprophet:probability": ("peptideprophet_prob", True),
    "iprophet:probability": ("iprophet_prob", True),
    "ptmprophet:probability": ("ptmprophet_prob", True),
}

# ============================================================================
# MzIdentML Data Processor
# ============================================================================


class MzIdentML:
    """mzIdentML data processor for converting output files to QPX format.

    This class provides functionality to parse mzIdentML files and convert them
    to QPX PSM format. mzIdentML is a standard format for reporting peptide and
    protein identifications from mass spectrometry experiments.

    Attributes:
        mzid_path: Path to the mzIdentML file
        _mzml_path: Optional path to corresponding mzML file
        _spectral_data: Whether to include spectral data
        _psm_list: List of parsed PSM dictionaries
    """

    def __init__(
        self,
        mzid_path: Union[Path, str],
        mzml_path: Optional[Union[Path, str]] = None,
        mzml_folder: Optional[Union[Path, str]] = None,
        spectral_data: bool = False,
    ):
        """Initialize the mzIdentML parser.

        Args:
            mzid_path: Path to the mzIdentML file (can be .mzid or .mzid.gz)
            mzml_path: Optional path to a single mzML file for attaching spectra
            mzml_folder: Optional folder containing mzML files (matched by reference_file_name)
            spectral_data: Whether to include spectral data in output
        """
        self.mzid_path = Path(mzid_path)
        self._mzml_path: Optional[Path] = Path(mzml_path) if mzml_path else None
        self._mzml_folder: Optional[Path] = Path(mzml_folder) if mzml_folder else None
        self._spectral_data = spectral_data

        if self._spectral_data:
            logger.info("Loading spectra information into QPX")
        else:
            logger.info("Spectra information will not be loaded into QPX")

        self._namespace = None
        self._ns = {}
        self._peptide_map = {}  # peptide_ref -> peptide info
        self._dbsequence_map = {}  # dbsequence_ref -> protein info
        self._peptide_evidence_map = {}  # peptide_evidence_ref -> evidence info
        self._spectra_data_map = {}  # spectra_data_ref -> file info
        self._psm_list = []
        self._mzml_file_cache = {}  # Cache for loaded mzML files

        self._parse_mzidentml()

        # Attach mzML spectra if provided
        if self._spectral_data:
            if self._mzml_folder is not None:
                self._attach_mzml_spectra_from_folder()
            elif self._mzml_path is not None:
                self._attach_mzml_spectra()

    # ========================================================================
    # XML Parsing Methods
    # ========================================================================

    def _parse_mzidentml(self) -> None:
        """Parse mzIdentML file using ElementTree."""
        try:
            logger.info(f"Parsing mzIdentML file: {self.mzid_path}")

            # Handle gzipped files
            if str(self.mzid_path).endswith(".gz"):
                with gzip.open(self.mzid_path, "rt", encoding="utf-8") as f:
                    tree = ET.parse(f)
            else:
                tree = ET.parse(self.mzid_path)

            root = tree.getroot()

            # Detect namespace
            self._detect_namespace(root)

            # Parse reference data
            self._parse_db_sequences(root)
            self._parse_peptides(root)
            self._parse_peptide_evidences(root)
            self._parse_spectra_data(root)

            # Parse PSMs
            self._parse_spectrum_identifications(root)

            logger.info(f"Parsed {len(self._psm_list)} PSMs from mzIdentML")

        except Exception as e:
            logger.error(f"Error parsing mzIdentML file: {e}")
            raise

    def _detect_namespace(self, root) -> None:
        """Detect mzIdentML namespace from root element"""
        tag = root.tag
        if "}" in tag:
            ns_uri = tag.split("}")[0].strip("{")
            self._namespace = ns_uri
            self._ns = {"mzid": ns_uri}
        else:
            # Try common namespaces
            for version, ns_uri in MZID_NAMESPACES.items():
                self._ns = {"mzid": ns_uri}
                if root.find(".//mzid:Peptide", self._ns) is not None:
                    self._namespace = ns_uri
                    break

    def _parse_db_sequences(self, root) -> None:
        """Parse DBSequence elements (protein references)"""
        for dbseq in root.findall(".//mzid:DBSequence", self._ns):
            db_id = dbseq.get("id")
            accession = dbseq.get("accession")

            # Check for decoy
            is_decoy = 0
            for cv_param in dbseq.findall("mzid:cvParam", self._ns):
                if "decoy" in cv_param.get("name", "").lower():
                    is_decoy = 1
                    break

            self._dbsequence_map[db_id] = {
                "accession": accession,
                "is_decoy": is_decoy,
            }

    def _parse_peptides(self, root) -> None:
        """Parse Peptide elements"""
        for peptide in root.findall(".//mzid:Peptide", self._ns):
            pep_id = peptide.get("id")

            # Get sequence
            pep_seq_elem = peptide.find("mzid:PeptideSequence", self._ns)
            sequence = pep_seq_elem.text if pep_seq_elem is not None else ""

            # Parse modifications
            modifications = []
            for mod in peptide.findall("mzid:Modification", self._ns):
                mod_info = self._parse_modification(mod, sequence)
                if mod_info:
                    modifications.append(mod_info)

            # Build peptidoform string
            peptidoform = self._build_peptidoform(sequence, modifications)

            self._peptide_map[pep_id] = {
                "sequence": sequence,
                "peptidoform": peptidoform,
                "modifications": modifications if modifications else None,
            }

    def _parse_modification(self, mod_elem, sequence: str) -> Optional[Dict]:
        """Parse a single Modification element"""
        location = mod_elem.get("location")
        residues = mod_elem.get("residues", "")
        mono_mass_delta = mod_elem.get("monoisotopicMassDelta")

        # Get modification name from cvParam
        mod_name = None
        mod_accession = None
        for cv_param in mod_elem.findall("mzid:cvParam", self._ns):
            mod_name = cv_param.get("name")
            accession = cv_param.get("accession")
            if accession and accession.startswith("UNIMOD:"):
                mod_accession = accession
            break

        if not mod_name:
            return None

        # Determine position
        try:
            loc = int(location) if location else 0
        except ValueError:
            loc = 0

        # Format position string
        if loc == 0:
            position_str = "N-term.0"
        elif loc > len(sequence):
            position_str = "C-term." + str(len(sequence) + 1)
        else:
            aa = sequence[loc - 1] if loc <= len(sequence) else ""
            position_str = f"{aa}.{loc}"

        return {
            "name": mod_name,
            "accession": mod_accession,
            "positions": [{"position": position_str, "scores": None}],
        }

    def _build_peptidoform(self, sequence: str, modifications: List[Dict]) -> str:
        """Build peptidoform string from sequence and modifications"""
        if not modifications:
            return sequence

        # Simple format: sequence with modifications noted
        # For a more sophisticated format, would need proper ProForma encoding
        return sequence

    def _parse_peptide_evidences(self, root) -> None:
        """Parse PeptideEvidence elements (peptide to protein mapping)"""
        for pe in root.findall(".//mzid:PeptideEvidence", self._ns):
            pe_id = pe.get("id")
            peptide_ref = pe.get("peptide_ref")
            dbseq_ref = pe.get("dBSequence_ref")
            is_decoy = pe.get("isDecoy", "false").lower() == "true"

            self._peptide_evidence_map[pe_id] = {
                "peptide_ref": peptide_ref,
                "dbsequence_ref": dbseq_ref,
                "is_decoy": 1 if is_decoy else 0,
            }

    def _parse_spectra_data(self, root) -> None:
        """Parse SpectraData elements (source file references)"""
        for sd in root.findall(".//mzid:SpectraData", self._ns):
            sd_id = sd.get("id")
            location = sd.get("location", "")

            # Extract filename from location
            filename = Path(location).name if location else sd_id

            self._spectra_data_map[sd_id] = {
                "location": location,
                "filename": filename,
            }

    def _parse_spectrum_identifications(self, root) -> None:
        """Parse SpectrumIdentificationResult and SpectrumIdentificationItem elements"""
        for sir in root.findall(".//mzid:SpectrumIdentificationResult", self._ns):
            spectrum_id = sir.get("spectrumID", "")
            spectra_data_ref = sir.get("spectraData_ref", "")

            # Extract scan number from spectrum_id
            scan = self._extract_scan_number(spectrum_id)

            # Get reference file name
            ref_file = self._spectra_data_map.get(spectra_data_ref, {}).get(
                "filename", ""
            )

            # Parse each PSM in this spectrum
            for sii in sir.findall("mzid:SpectrumIdentificationItem", self._ns):
                psm = self._parse_spectrum_identification_item(sii, scan, ref_file)
                if psm:
                    self._psm_list.append(psm)

    def _parse_spectrum_identification_item(
        self, sii, scan: str, reference_file_name: str
    ) -> Optional[Dict]:
        """Parse a single SpectrumIdentificationItem (PSM)"""
        try:
            charge = int(sii.get("chargeState", 0))
            exp_mz = float(sii.get("experimentalMassToCharge", 0))
            calc_mz = (
                float(sii.get("calculatedMassToCharge", 0))
                if sii.get("calculatedMassToCharge")
                else None
            )
            rank = int(sii.get("rank", 1))
            pass_threshold = sii.get("passThreshold", "false").lower() == "true"
            peptide_ref = sii.get("peptide_ref", "")

            # Get peptide info
            peptide_info = self._peptide_map.get(peptide_ref, {})
            sequence = peptide_info.get("sequence", "")
            peptidoform = peptide_info.get("peptidoform", sequence)
            modifications = peptide_info.get("modifications")

            # Get protein accessions and decoy status from PeptideEvidenceRef
            protein_accessions = []
            is_decoy = 0
            for pe_ref in sii.findall("mzid:PeptideEvidenceRef", self._ns):
                pe_id = pe_ref.get("peptideEvidence_ref")
                pe_info = self._peptide_evidence_map.get(pe_id, {})

                dbseq_ref = pe_info.get("dbsequence_ref")
                if dbseq_ref:
                    dbseq_info = self._dbsequence_map.get(dbseq_ref, {})
                    accession = dbseq_info.get("accession")
                    if accession and accession not in protein_accessions:
                        protein_accessions.append(accession)
                    if dbseq_info.get("is_decoy") or pe_info.get("is_decoy"):
                        is_decoy = 1

            # Extract scores from cvParams
            additional_scores, q_value, pep = self._extract_scores(sii)

            # Extract RT if available
            rt = None
            for cv_param in sii.findall("mzid:cvParam", self._ns):
                if "retention time" in cv_param.get("name", "").lower():
                    try:
                        rt = float(cv_param.get("value", 0))
                    except ValueError:
                        pass

            return {
                "sequence": sequence,
                "peptidoform": peptidoform,
                "modifications": modifications,
                "precursor_charge": charge,
                "posterior_error_probability": pep,
                "is_decoy": is_decoy,
                "calculated_mz": calc_mz,
                "observed_mz": exp_mz,
                "additional_scores": additional_scores if additional_scores else None,
                "protein_accessions": (
                    protein_accessions if protein_accessions else None
                ),
                "predicted_rt": None,
                "reference_file_name": reference_file_name,
                "cv_params": None,
                "scan": scan,
                "rt": rt,
                "q_value": q_value,
                "ion_mobility": None,
                "number_peaks": None,
                "mz_array": None,
                "intensity_array": None,
                "charge_array": None,
                "ion_type_array": None,
                "ion_mobility_array": None,
            }

        except Exception as e:
            logger.warning(f"Error parsing SpectrumIdentificationItem: {e}")
            return None

    def _extract_scan_number(self, spectrum_id: str) -> str:
        """Extract scan number from spectrum ID.

        Supports multiple vendor formats:
        - Thermo: scan=1234
        - Waters/Agilent: sample=1 period=1 cycle=1234 experiment=1
        - Generic: index=1234, spectrum=1234
        """
        patterns = [
            r"scan[=:](\d+)",
            r"cycle[=:](\d+)",  # Waters/Agilent format
            r"index[=:](\d+)",
            r"spectrum[=:](\d+)",
            r"_(\d+)$",
        ]

        for pattern in patterns:
            match = re.search(pattern, spectrum_id, re.IGNORECASE)
            if match:
                return match.group(1)

        return spectrum_id

    def _extract_scores(
        self, sii
    ) -> Tuple[List[Dict], Optional[float], Optional[float]]:
        """Extract scores from SpectrumIdentificationItem.

        Args:
            sii: SpectrumIdentificationItem XML element

        Returns:
            Tuple of (additional_scores list, q_value, posterior_error_probability)
        """
        additional_scores = []
        q_value = None
        pep = None

        for cv_param in sii.findall("mzid:cvParam", self._ns):
            name = cv_param.get("name", "").lower()
            value_str = cv_param.get("value")

            if value_str is None:
                continue

            try:
                value = float(value_str)
            except ValueError:
                continue

            # Check for q-value
            if "q-value" in name or "qvalue" in name:
                q_value = value
                continue

            # Check for PEP
            if "posterior error probability" in name or name == "pep":
                pep = value
                continue

            # Map to standard score names using module-level constant
            for pattern, (score_name, higher_better) in SCORE_MAPPINGS.items():
                if pattern in name:
                    additional_scores.append(
                        {
                            "score_name": score_name,
                            "score_value": value,
                            "higher_better": higher_better,
                        }
                    )
                    break
            else:
                # Unknown score - add with unknown direction
                if any(x in name for x in ["score", "probability", "expect", "value"]):
                    additional_scores.append(
                        {
                            "score_name": name.replace(" ", "_"),
                            "score_value": value,
                            "higher_better": None,
                        }
                    )

        return additional_scores, q_value, pep

    # ========================================================================
    # mzML Spectra Attachment
    # ========================================================================

    def _attach_mzml_spectra(self) -> None:
        """Attach spectral data from mzML file to PSMs.

        Uses scan number to match spectra from mzML file to PSMs.
        Populates number_peaks, mz_array, and intensity_array fields.
        """
        try:
            from qpx.core.openms import OpenMSHandler

            logger.info(f"Attaching spectra from mzML file: {self._mzml_path}")
            handler = OpenMSHandler()

            attached_count = 0
            for psm in self._psm_list:
                scan = psm.get("scan")
                if scan is None:
                    continue

                try:
                    scan_int = int(str(scan))
                    num_peaks, mzs, intens = handler.get_spectrum_from_scan(
                        str(self._mzml_path), scan_int
                    )
                    psm["number_peaks"] = int(num_peaks)
                    psm["mz_array"] = [float(x) for x in mzs]
                    psm["intensity_array"] = [float(x) for x in intens]
                    attached_count += 1
                except Exception as e:
                    logger.debug(f"Could not attach spectrum for scan {scan}: {e}")

            logger.info(
                f"Attached spectra to {attached_count}/{len(self._psm_list)} PSMs"
            )

        except Exception as e:
            logger.warning(f"Failed to attach mzML spectra: {e}")

    def _attach_mzml_spectra_from_folder(self) -> None:
        """Attach spectral data from multiple mzML files in a folder.

        Uses reference_file_name from each PSM to find the matching mzML file
        in the specified folder. Supports both .mzML and .mzML.gz files.
        """
        try:
            from qpx.core.openms import OpenMSHandler

            logger.info(f"Attaching spectra from mzML folder: {self._mzml_folder}")
            handler = OpenMSHandler()

            # Build file lookup map (case-insensitive)
            mzml_files = {}
            for ext in ["*.mzML", "*.mzml", "*.mzML.gz", "*.mzml.gz"]:
                for f in self._mzml_folder.glob(ext):
                    # Store with multiple keys for flexible matching
                    base_name = f.name.lower()
                    mzml_files[base_name] = f
                    # Also store without .gz extension
                    if base_name.endswith(".gz"):
                        mzml_files[base_name[:-3]] = f
                    # Store stem for matching without extension
                    stem = f.stem.lower()
                    if stem.endswith(".mzml"):
                        stem = stem[:-5]
                    mzml_files[stem] = f

            logger.info(f"Found {len(set(mzml_files.values()))} unique mzML files")

            attached_count = 0
            missing_files = set()

            for psm in self._psm_list:
                scan = psm.get("scan")
                ref_file = psm.get("reference_file_name")

                if scan is None or ref_file is None:
                    continue

                # Normalize reference file name for matching
                ref_lower = ref_file.lower()
                ref_basename = Path(ref_file).name.lower()
                ref_stem = Path(ref_file).stem.lower()
                if ref_stem.endswith(".mzml"):
                    ref_stem = ref_stem[:-5]

                # Try to find matching mzML file
                mzml_path = None
                for key in [ref_basename, ref_lower, ref_stem, f"{ref_basename}.gz"]:
                    if key in mzml_files:
                        mzml_path = mzml_files[key]
                        break

                if mzml_path is None:
                    if ref_file not in missing_files:
                        missing_files.add(ref_file)
                    continue

                try:
                    scan_int = int(str(scan))
                    num_peaks, mzs, intens = handler.get_spectrum_from_scan(
                        str(mzml_path), scan_int
                    )
                    if num_peaks > 0:
                        psm["number_peaks"] = int(num_peaks)
                        psm["mz_array"] = [float(x) for x in mzs]
                        psm["intensity_array"] = [float(x) for x in intens]
                        attached_count += 1
                except Exception as e:
                    logger.debug(f"Could not attach spectrum for scan {scan}: {e}")

            if missing_files:
                logger.warning(
                    f"Could not find mzML files for: {sorted(missing_files)[:5]}..."
                )

            logger.info(
                f"Attached spectra to {attached_count}/{len(self._psm_list)} PSMs"
            )

        except Exception as e:
            logger.warning(f"Failed to attach mzML spectra from folder: {e}")

    # ========================================================================
    # Output Methods
    # ========================================================================

    def to_dataframe(self) -> pd.DataFrame:
        """Convert parsed PSMs to a pandas DataFrame.

        Returns:
            DataFrame containing PSM data
        """
        if not self._psm_list:
            return pd.DataFrame()

        df = pd.DataFrame(self._psm_list)
        return df

    def to_parquet(self, output_path: Union[Path, str]) -> None:
        """Convert mzIdentML data to parquet format and save to file.

        Args:
            output_path: Output file path for parquet file
        """
        df = self.to_dataframe()
        if df.empty:
            logger.warning("No PSMs found to convert")
            return

        def convert_numpy_to_list(obj):
            """Recursively convert numpy arrays to Python lists"""
            if hasattr(obj, "tolist"):
                return obj.tolist()
            elif isinstance(obj, list):
                return [convert_numpy_to_list(item) for item in obj]
            elif isinstance(obj, dict):
                return {key: convert_numpy_to_list(value) for key, value in obj.items()}
            else:
                return obj

        for col in df.columns:
            if col in [
                "modifications",
                "protein_accessions",
                "additional_scores",
                "cv_params",
                "mz_array",
                "intensity_array",
            ]:
                df[col] = df[col].apply(convert_numpy_to_list)

        table = pa.Table.from_pandas(df, schema=PSM_SCHEMA)
        pq.write_table(table, output_path)
        logger.info(f"Successfully converted mzIdentML to parquet: {output_path}")

    def get_psm_count(self) -> int:
        """Get the total number of PSMs found."""
        return len(self._psm_list)
