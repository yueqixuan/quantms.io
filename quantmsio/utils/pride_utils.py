"""
This file contains utility functions for parsing PRIDE JSON files
"""

import itertools
import logging
import re
import time
from collections import defaultdict
from typing import Optional

import pandas as pd
from Bio import SeqIO

from quantmsio.utils.mztab_utils import get_modifications_object_from_mztab_line

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def get_unanimous_name(protein_accessions: str | list, map_dict: dict) -> list | None:
    if isinstance(protein_accessions, str):
        if ";" in protein_accessions:
            protein_accessions = protein_accessions.split(";")
        else:
            protein_accessions = protein_accessions.split(",")
    unqnimous_names = []
    for accession in protein_accessions:
        if accession in map_dict:
            unqnimous_names.append(list(map_dict[accession])[0])
    if len(unqnimous_names) > 0:
        return unqnimous_names
    else:
        return None


def generate_gene_name_map(fasta: str, map_parameter: str) -> defaultdict:
    """
    according fasta database to map the proteins accessions to uniprot names.
    :param fasta: Reference fasta database
    :param map_parameter: map_protein_name or map_protein_accession
    retrun: None
    """
    map_gene_names = defaultdict(set)
    if map_parameter == "map_protein_name":
        for seq_record in SeqIO.parse(fasta, "fasta"):
            name = seq_record.id.split("|")[-1]
            gene_list = re.findall(r"GN=(\S+)", seq_record.description)
            gene_name = gene_list[0] if len(gene_list) > 0 else None
            map_gene_names[name].add(gene_name)
    else:
        for seq_record in SeqIO.parse(fasta, "fasta"):
            accession = seq_record.id.split("|")[-2]
            gene_list = re.findall(r"GN=(\S+)", seq_record.description)
            gene_name = gene_list[0] if len(gene_list) > 0 else None
            map_gene_names[accession].add(gene_name)
    return map_gene_names


def get_gene_accessions(gene_list: list | None, map_dict: dict) -> list | None:
    if gene_list is None or len(gene_list) == 0:
        return None
    else:
        accessions = []
        for gene in gene_list:
            if gene in map_dict:
                accession = map_dict[gene]
                if accession:
                    accessions.append(accession)
        if len(accessions) > 0:
            return accessions
        else:
            return None


def generate_scan_number(spectra_ref: str) -> str:
    is_scan = re.findall(r"controllerType=0 controllerNumber=1", spectra_ref)
    if len(is_scan) != 0:
        return re.findall(r"scan=(\d+)", spectra_ref)[0]
    else:
        return ",".join(re.findall(r"=(\d+)", spectra_ref))


def get_pubmed_id_pride_json(pride_json: dict) -> Optional[str]:
    """
    Parse the PubMed ID from the PRIDE JSON file
    :param pride_json: PRIDE JSON file
    :return: PubMed ID
    """
    pubmed_id = None
    if "references" in pride_json:
        for reference in pride_json["references"]:
            if "pubmedId" in reference:
                pubmed_id = reference["pubmedId"]
    return pubmed_id


def get_set_of_experiment_keywords(pride_json: dict) -> list:
    """
    Parse the experiment type from the PRIDE JSON file
    :param pride_json: PRIDE JSON file
    :return: Set of experiment types
    """
    experiment_types: set = set()
    if "projectTags" in pride_json:
        for experiment_type in pride_json["projectTags"]:
            experiment_types.add(experiment_type)
    if "keywords" in pride_json:
        for experiment_type in pride_json["keywords"]:
            experiment_types.add(experiment_type)
    return list(experiment_types)


def get_key_peptide_combination(
    msstats_peptidoform: str, charge: str, reference_file: Optional[str] = None
) -> str:
    """
    Get the key for a peptide or psm in mztab. For peptides the key is:
     - sequence peptidoform in msstats notation (e.g. PEPTIDE(Oxidation))
     - charge (e.g. 2)
    For psms the key is:
     - sequence peptidoform in msstats notation (e.g. PEPTIDE(Oxidation))
     - charge (e.g. 2)
     - reference_file (e.g. Run1-RAW)
    :param msstats_peptidoform: sequence peptidoform in msstats notation (e.g. PEPTIDE(Oxidation))
    :param charge: charge (e.g. 2)
    :param reference_file: reference file (e.g. Run1-RAW)
    :return: key
    """
    if reference_file is not None:
        # We use a different separator for the key to avoid conflicts with the separator of the msruns
        # which can sometimes use the separator "_".
        return msstats_peptidoform + ":_:" + charge + ":_:" + reference_file
    return msstats_peptidoform + ":_:" + charge


def decompose_key_peptide_combination(key: str) -> list:
    """
    decompose the key for a peptide in mztab. The key is a combination of msstats_peptidoform,
    charge and in some cases the reference file.
    :param key: key
    :return: sequence, modification, charge and retention time
    """
    if ":_:" in key:
        return key.split(":_:")
    return key.split("_")


def clean_peptidoform_sequence(sequence: str) -> str:
    """
    clean any peptidoform a proforma sequence or a msstats sequence. This function removes any modification from the
    sequence.
    :param sequence: sequence
    """
    if "(" not in sequence and "[" not in sequence:  # no modification
        return sequence

    # remove any modification if a proforma sequence is provided
    sequence = re.sub(r"[\(\[].*?[\)\]]", "", sequence)
    sequence = sequence.replace(".", "").replace(" ", "").replace("-", "")
    return sequence


def compare_protein_lists(protein_list_1: list, protein_list_2: list) -> bool:
    """
    Compare two protein lists
    :param protein_list_1: protein list 1
    :param protein_list_2: protein list 2
    :return: True if the three protein lists are the same
    """
    protein_list_1 = set(x.strip() for x in protein_list_1)
    protein_list_2 = set(x.strip() for x in protein_list_2)
    return protein_list_1 == protein_list_2


def standardize_protein_string_accession(
    protein_string: str, is_sorted: bool = False
) -> str:
    """
    standardize the protein string accession, in some cases the protein string accession is decided by commas
    instead of semicolons.
    :param protein_string: protein string
    :param is_sorted: sort the protein string
    :return: standardized protein string
    """
    protein_string = protein_string.replace(",", ";").strip()
    if is_sorted:
        accessions = protein_string.split(";")
        accessions.sort()
        return ";".join(accessions)
    return protein_string


def standardize_protein_list_accession(protein_string: str) -> list:
    """
    get the list of protein accessions from a protein string join by semicolons.
    :param protein_string: protein string
    :return: list of protein accessions
    """
    return [x.strip() for x in protein_string.split(";")]


def parse_score_name_in_mztab(score_name_mztab_line: str) -> str:
    """
    parse the score name in mztab. The score name in mztab is a combination of the score name and the score type.
    :param score_name_mztab_line: score name in mztab
    :return: score name
    """
    lines = score_name_mztab_line.split("\t")
    score_values = lines[2].replace("[", "").replace("]", "").split(",")
    score_name = score_values[2].strip()
    if ":" in score_name:
        score_name = f"'{score_name}'"  # add quotes to the score name if it contains a colon like
        # "OpenMS:Target-decoy protein q-value"
    return score_name


def fetch_peptide_spectra_ref(peptide_spectra_ref: str) -> tuple[str, str]:
    """
    Get the ms run and the scan number from a spectra ref. The spectra ref is in the format:
    ms_run[1]:controllerType=0 controllerNumber=1 scan=1
    :param peptide_spectra_ref: spectra ref
    :return: ms run and scan number
    """
    ms_run = peptide_spectra_ref.split(":")[0]
    # scan_number = peptide_spectra_ref.split(":")[1].split(" ")[-1].split("=")[-1]
    scan_number = generate_scan_number(peptide_spectra_ref)
    return ms_run, scan_number


def compare_msstats_peptidoform_with_proforma(
    msstats_peptidoform: str, proforma: str
) -> bool:
    """
    Compare a msstats peptidoform with a proforma representation.
    - [Acetyl]-EM[Oxidation]EVEES[Phospho]PEK vs .(Acetyl)EM(Oxidation)EVEES(Phospho)PEK
    :param msstats_peptidoform: msstats peptidoform
    :param proforma: proforma peptidoform
    :return: True if the peptidoforms are the same, False otherwise
    """
    proforma = proforma.replace("[", "(").replace("]", ")")
    if "-" in proforma:
        proforma_parts = proforma.split("-")
        if len(proforma_parts) > 3:
            raise Exception(
                "The proforma peptidoform is not valid has more than 3 termini modifications"
            )
        result_proforma: str = ""
        par_index = 0
        for part in proforma_parts:
            if part.startswith("(") and par_index == 0:  # n-term modification
                result_proforma = result_proforma + "." + part
            elif part.startswith("(") and par_index > 1:
                result_proforma = result_proforma + "." + part
            else:
                result_proforma = result_proforma + part
            par_index += 1
        proforma = result_proforma
    return msstats_peptidoform == proforma


def get_petidoform_msstats_notation(
    peptide_sequence: str,
    modification_string: Optional[str],
    modifications_definition: dict,
) -> str:
    """
    Get the peptidoform in msstats notation from a mztab line definition, it could be for a PEP definition
    or for a PSM. The peptidoform in msstats notation is defined as:
        - EM(Oxidation)EVEES(Phospho)PEK
        - .(iTRAQ4plex)EMEVNESPEK.(Methyl)
    :param peptide_sequence: peptide sequence
    :param modification_string: modification string
    :param modifications_definition: dictionary modifications definition
    :return: peptidoform in msstats
    """
    if (
        modification_string is None
        or modification_string == "null"
        or modification_string == ""
    ):
        return peptide_sequence

    modifications = get_modifications_object_from_mztab_line(
        modification_string=modification_string,
        modifications_definition=modifications_definition,
    )

    aa_index = 0
    result_peptide: str = ""
    peptide_sequence = list(peptide_sequence)
    # Add n-term modification if it is present
    for key_index, value_index in modifications.items():
        if aa_index in value_index["position"]:
            result_peptide = "." + "(" + key_index + ")" + result_peptide

    aa_index += 1
    for aa in peptide_sequence:
        add_mod: str = ""
        for key_index, value_index in modifications.items():
            if aa_index in value_index["position"]:
                add_mod = add_mod + "(" + key_index + ")"
        result_peptide = result_peptide + aa + add_mod
        aa_index += 1
    # Add c-term modification if it is present
    for key_index, value_index in modifications.items():
        if aa_index in value_index["position"]:
            result_peptide = result_peptide + ".(" + key_index + ")"
    return result_peptide


def fetch_peptide_from_mztab_line(
    peptide_dict: dict,
    ms_runs: dict | None = None,
    modification_definition: dict | None = None,
) -> dict:
    """
    Get the peptide from a mztab line include the post.
    :param peptide_dict: dictionary with the peptide information
    :param ms_runs: ms runs dictionary
    :param modification_definition: modification definition
    :return: peptide dictionary
    """
    keys = ["sequence", "modifications", "charge", "retention_time", "accession"]
    peptide: dict = dict(zip(keys, [peptide_dict[k] for k in keys]))

    peptide["accession"] = standardize_protein_string_accession(peptide["accession"])
    peptide["score"] = peptide_dict["best_search_engine_score[1]"]

    if (
        "opt_global_cv_MS:1002217_decoy_peptide" in peptide_dict
    ):  # check if the peptide is a decoy
        peptide["is_decoy"] = peptide_dict["opt_global_cv_MS:1002217_decoy_peptide"]
    else:
        peptide["is_decoy"] = None  # if the information is not available

    if "opt_global_cv_MS:1000889_peptidoform_sequence" in peptide_dict:
        peptide["peptidoform"] = peptide_dict[
            "opt_global_cv_MS:1000889_peptidoform_sequence"
        ]
    else:
        peptide["peptidoform"] = get_petidoform_msstats_notation(
            peptide["sequence"], peptide["modifications"], modification_definition
        )  # if the information is not available

    peptide["proforma_peptidoform"] = get_peptidoform_proforma_version_in_mztab(
        peptide_sequence=peptide["sequence"],
        modification_string=peptide["modifications"],
        modifications_definition=modification_definition,
    )

    if "spectra_ref" in peptide_dict:
        ms_run, scan_number = fetch_peptide_spectra_ref(peptide_dict["spectra_ref"])
        if ms_runs is not None and ms_run in ms_runs:
            peptide["ms_run"] = ms_runs[ms_run]
        elif ms_runs is not None and ms_run not in ms_runs:
            raise Exception(f"The ms run {ms_run} is not in the ms runs index")
        else:
            peptide["ms_run"] = ms_run
        peptide["scan_number"] = scan_number

    else:
        peptide["ms_run"] = None  # if the information is not available
        peptide["scan_number"] = None  # if the information is not available

    comparison_representation = compare_msstats_peptidoform_with_proforma(
        msstats_peptidoform=peptide["peptidoform"],
        proforma=peptide["proforma_peptidoform"],
    )

    if not comparison_representation:
        raise Exception(
            f"The proforma peptidoform {peptide['proforma_peptidoform']} & the msstats peptidoform {peptide['peptidoform']} are not equal"
        )
    return peptide


def fetch_protein_from_mztab_line(protein_dict: dict):
    """
    get the protein from a mztab line.
    :param protein_dict: dictionary with the protein information
    :return: protein dictionary
    """
    accession_string = standardize_protein_string_accession(
        protein_dict["ambiguity_members"]
    )
    return {
        "accession": accession_string,
        "score": protein_dict["best_search_engine_score[1]"],
    }


def _get_required_keys(es: dict) -> list:
    """Get the list of required keys for PSM data."""
    keys = [
        "sequence",
        "modifications",
        "charge",
        "retention_time",
        "accession",
        "start",
        "end",
        "calc_mass_to_charge",
        "exp_mass_to_charge",
        "opt_global_Posterior_Error_Probability_score",
        "opt_global_q-value",
        "opt_global_consensus_support",
    ]

    # Remove optional keys that don't exist in the data
    optional_keys = [
        "opt_global_q-value",
        "opt_global_consensus_support",
        "opt_global_Posterior_Error_Probability_score",
    ]

    for key in optional_keys:
        if key not in es:
            keys.remove(key)

    return keys


def _set_optional_psm_fields(psm: dict, es: dict) -> None:
    """Set optional PSM fields with proper fallbacks."""
    # Handle decoy peptide information
    if "opt_global_cv_MS:1002217_decoy_peptide" in es:
        psm["is_decoy"] = es["opt_global_cv_MS:1002217_decoy_peptide"]
    else:
        psm["is_decoy"] = None

    # Handle posterior error probability
    if "opt_global_Posterior_Error_Probability_score" in es:
        psm["posterior_error_probability"] = es[
            "opt_global_Posterior_Error_Probability_score"
        ]
    else:
        psm["posterior_error_probability"] = None

    # Handle global q-value
    if "opt_global_q-value" in es:
        psm["global_qvalue"] = es["opt_global_q-value"]
    else:
        psm["global_qvalue"] = None

    # Handle consensus support
    if "opt_global_consensus_support" in es:
        psm["consensus_support"] = es["opt_global_consensus_support"]
    else:
        psm["consensus_support"] = None


def _set_peptidoform_fields(
    psm: dict, es: dict, modifications_definition: dict
) -> None:
    """Set peptidoform-related fields."""
    if "opt_global_cv_MS:1000889_peptidoform_sequence" in es:
        psm["peptidoform"] = es["opt_global_cv_MS:1000889_peptidoform_sequence"]
    else:
        psm["peptidoform"] = get_petidoform_msstats_notation(
            psm["sequence"], psm["modifications"], modifications_definition
        )

    psm["proforma_peptidoform"] = get_peptidoform_proforma_version_in_mztab(
        psm["sequence"], psm["modifications"], modifications_definition
    )


def _set_spectra_fields(psm: dict, es: dict, ms_runs: dict | None) -> None:
    """Set spectra reference and scan number fields."""
    if "spectra_ref" in es:
        ms_run, scan_number = fetch_peptide_spectra_ref(es["spectra_ref"])

        if ms_runs is not None:
            if ms_run in ms_runs:
                psm["ms_run"] = ms_runs[ms_run]
            else:
                raise Exception(f"The ms run {ms_run} is not in the ms runs index")
        else:
            psm["ms_run"] = ms_run

        psm["scan_number"] = scan_number
    else:
        psm["ms_run"] = None
        psm["scan_number"] = None


def fetch_psm_from_mztab_line(
    es: dict, ms_runs: dict | None = None, modifications_definition: dict | None = None
) -> dict:
    """
    Get the psm from a mztab line include the post.
    :param es: dictionary with the psm information
    :param ms_runs: ms runs dictionary
    :param modifications_definition: modifications definition
    :return: psm dictionary
    """
    # Get required keys and create initial PSM dictionary
    keys = _get_required_keys(es)
    psm: dict = dict(zip(keys, [es[k] for k in keys]))

    # Set core fields
    psm["accession"] = standardize_protein_string_accession(psm["accession"])
    psm["score"] = es["search_engine_score[1]"]

    # Set optional fields
    _set_optional_psm_fields(psm, es)

    # Set peptidoform fields
    _set_peptidoform_fields(psm, es, modifications_definition)

    # Set spectra fields
    _set_spectra_fields(psm, es, ms_runs)

    return psm


def fetch_modifications_from_mztab_line(line: str, _modifications: dict) -> dict:
    """
    get the modifications from a mztab line. An mzTab modification could be a fixed or variable modification.
    The structure of a fixed is the following:
      MTD	fixed_mod[1]	[UNIMOD, UNIMOD:4, Carbamidomethyl, ]
      MTD	fixed_mod[1]-site	C
      MTD	fixed_mod[1]-position	Anywhere
    while the structure of a variable modification is the following:
      MTD	var_mod[1]	[UNIMOD, UNIMOD:21, Phospho, ]
      MTD	var_mod[1]-site	S
      MTD   var_mod[1]-position	Anywhere

    :param line: mztab line
    :param _modifications: modifications dictionary
    :return: modification dictionary
    """
    line = line.strip()
    line_parts = line.split("\t")
    if line_parts[0] == "MTD" and "_mod[" in line_parts[1]:
        if "site" not in line_parts[1] and "position" not in line_parts[1]:
            values = line_parts[2].replace("[", "").replace("]", "").split(",")
            accession = values[1].strip()
            name = values[2].strip()
            index = line_parts[1].split("[")[1].split("]")[0]
            _modifications[accession] = [name, index, None, None]
        elif "site" in line_parts[1]:
            index = line_parts[1].split("[")[1].split("]")[0]
            accession = None
            for key, value in _modifications.items():
                if value[1] == index:
                    accession = key
            if accession is None:
                raise Exception("The accession for the modification is None")
            _modifications[accession][2] = line_parts[2]
        elif "position" in line_parts[1]:
            index = line_parts[1].split("[")[1].split("]")[0]
            accession = None
            for key, value in _modifications.items():
                if value[1] == index:
                    accession = key
            if accession is None:
                raise Exception("The accession for the modification is None")
            _modifications[accession][3] = line_parts[2]
    return _modifications


def get_peptidoform_proforma_version_in_mztab(
    peptide_sequence: str,
    modification_string: Optional[str],
    modifications_definition: dict,
) -> str:
    """
    Get the peptidoform in propoforma notation from a mztab line definition, it could be for a PEP definition
    or for a PSM. The peptidoform in proforma notation is defined as:
     - EM[Oxidation]EVEES[Phospho]PEK
     - [iTRAQ4plex]-EMEVNESPEK-[Methyl]
    The modification string is defined as:
     - 3|4|8-UNIMOD:21, 2-UNIMO:35
     - 2[MS,MS:1001876, modification probability, 0.8]-UNIMOD:21, 3[MS,MS:1001876, modification probability, 0.8]-UNIMOD:31
     This modification means that Phospho can be in position 3, 4 or 8 and Oxidation is in position 2.
    :param peptide_sequence: Peptide sequence
    :param modification_string: modification string
    :param modifications_definition: dictionary modifications definition
    :return: peptidoform in proforma
    """
    if (
        modification_string is None
        or modification_string == "null"
        or modification_string == ""
        or pd.isna(modification_string)
    ):
        return peptide_sequence

        modifications = get_modifications_object_from_mztab_line(
            modification_string=modification_string,
            modifications_definition=modifications_definition,
        )

        aa_index = 0
        result_peptide: str = ""
        peptide_sequence = list(peptide_sequence)
        # Add n-term modification if it is present
        for key_index, value_index in modifications.items():
            if aa_index in value_index["position"]:
                result_peptide = "[" + key_index + "]" + result_peptide
        if len(result_peptide) > 0:
            result_peptide = result_peptide + "-"
    aa_index += 1
    for aa in peptide_sequence:
        add_mod: str = ""
        for key_index, value_index in modifications.items():
            if aa_index in value_index["position"]:
                add_mod = add_mod + "[" + key_index + "]"
        result_peptide = result_peptide + aa + add_mod
        aa_index += 1
    # Add c-term modification if it is present
    for key_index, value_index in modifications.items():
        if aa_index in value_index["position"]:
            result_peptide = result_peptide + "-[" + key_index + "]"
    return result_peptide


def get_permutations_of_original_list(original_elems: list) -> list:
    """
    get all the posible list by permutating elements of an original list.
    :param original_elems
    :return: list of elements
    """
    for permutation in itertools.permutations(original_elems):
        if any(left == right for left, right in zip(permutation, original_elems)):
            continue
        else:
            yield permutation


def print_estimated_time(original_time: float, step: str):
    """
    Print the estimated time of a step
    :param original_time: original time
    :param step: step
    """
    end = time.time() - original_time
    logger.info("Estimated time for {} is {} seconds".format(step, str(end)))
