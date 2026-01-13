import re
from typing import Optional, Dict, List


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
            for (
                key,
                value,
            ) in (
                _modifications.items()
            ):  # for name, age in dictionary.iteritems():  (for Python 2.x)
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


def get_modifications_object_from_mztab_line(
    modification_string: str, modifications_definition: dict
) -> dict:
    """
    get the modifications from a mztab line. This method is used to transform peptide + modification strings to
    proteoform notations, for msstats notation and for proforma notation.
    :param modification_string: modification string
    :param modifications_definition: dictionary modifications definition
    :return: modifications dictionary
    """
    if modification_string is None or modification_string == "null":
        return {}

    modifications: dict = {}
    modification_values = re.split(r",(?![^\[]*\])", modification_string)
    for modification in modification_values:
        modification = modification.strip()
        if modification == "":
            return {}
        accession = modification.split("-")[1]
        unimod_accession = accession
        if accession not in modifications_definition:
            raise Exception(
                f"The modification {accession} is not in the modifications definition"
            )
        accession = modifications_definition[accession][
            0
        ]  # get the name of the modification
        position: list = []
        position_probability_string = modification.split("-")[0]
        if (
            "[" not in position_probability_string
            and "|" not in position_probability_string
        ):  # only one position
            position = [position_probability_string]
        elif (
            "[" not in position_probability_string
            and "|" in position_probability_string
        ):  # multiple positions not probability
            position = position_probability_string.split("|")
        else:
            positions_probabilities = position_probability_string.split("|")
            for position_probability in positions_probabilities:
                if "[" not in position_probability:
                    position.append(position_probability)
                else:
                    position_with_probability = position_probability.split("[")[0]
                    position.append(position_with_probability)
        position = [int(i) for i in position]
        if accession in modifications:
            position = modifications[accession]["position"] + position
        modifications[accession] = {
            "position": position,
            "unimod_accession": unimod_accession,
        }
    return modifications


def generate_modification_list(modification_str: str, modifications):
    modifications = get_modifications_object_from_mztab_line(
        modification_str, modifications
    )
    modifications_string = ""
    for _, value in modifications.items():
        modifications_string += "|".join(map(str, value["position"]))
        modifications_string = (
            modifications_string + "-" + value["unimod_accession"] + ","
        )
    modifications_string = modifications_string[:-1]  # Remove last comma
    modification_list = modifications_string.split(",")
    return modification_list


def parse_metadata_line(line: str) -> Dict[str, Optional[str]]:
    """
    Parse metadata line into dictionary.

    Parses mzTab metadata lines (starting with 'MTD') into key-value pairs.

    Args:
        line: Tab-separated metadata line from mzTab file

    Returns:
        Dictionary with 'key' and 'value' fields

    Example:
        >>> parse_metadata_line("MTD\tfixed_mod[1]\tCarbamidomethyl")
        {'key': 'fixed_mod[1]', 'value': 'Carbamidomethyl'}
    """
    parts = line.split("\t")
    return {"key": parts[1], "value": parts[2] if len(parts) > 2 else None}


def parse_header_line(line: str) -> List[str]:
    """
    Parse header line into column names.

    Parses mzTab header lines (starting with 'PRH' or 'PSH') into column names.

    Args:
        line: Tab-separated header line from mzTab file

    Returns:
        List of column names (excluding the section identifier)

    Example:
        >>> parse_header_line("PRH\taccession\tdescription\tscore")
        ['accession', 'description', 'score']
    """
    return line.split("\t")[1:]


def parse_protein_line(line: str, protein_header: List[str]) -> Dict[str, str]:
    """
    Parse protein line into dictionary.

    Parses mzTab protein lines (starting with 'PRT') into a dictionary with
    column names as keys.

    Args:
        line: Tab-separated protein line from mzTab file
        protein_header: List of column names from the protein header

    Returns:
        Dictionary with protein data

    Example:
        >>> header = ["accession", "description", "score"]
        >>> parse_protein_line("PRT\tP12345\tTest protein\t0.95", header)
        {'accession': 'P12345', 'description': 'Test protein', 'score': '0.95'}
    """
    values = line.split("\t")
    protein_dict = dict(zip(protein_header, values[1:]))
    return protein_dict


def parse_psm_line(line: str, psm_header: List[str]) -> Dict[str, str]:
    """
    Parse PSM line into dictionary.

    Parses mzTab PSM lines (starting with 'PSM') into a dictionary with
    column names as keys.

    Args:
        line: Tab-separated PSM line from mzTab file
        psm_header: List of column names from the PSM header

    Returns:
        Dictionary with PSM data

    Example:
        >>> header = ["accession", "sequence", "charge"]
        >>> parse_psm_line("PSM\tP12345\tPEPTIDE\t2", header)
        {'accession': 'P12345', 'sequence': 'PEPTIDE', 'charge': '2'}
    """
    values = line.split("\t")
    psm_dict = dict(zip(psm_header, values[1:]))
    return psm_dict


def is_mztab_line_type(line: str, line_type: str) -> bool:
    """
    Check if an mzTab line is of a specific type.

    Args:
        line: mzTab line to check
        line_type: Type identifier (e.g., 'MTD', 'PRH', 'PRT', 'PSH', 'PSM')

    Returns:
        True if the line is of the specified type, False otherwise

    Example:
        >>> is_mztab_line_type("MTD\tfixed_mod[1]\tCarbamidomethyl", "MTD")
        True
        >>> is_mztab_line_type("PRT\tP12345\tTest protein", "PSM")
        False
    """
    return line.strip().startswith(line_type)


def get_mztab_line_type(line: str) -> Optional[str]:
    """
    Get the type of an mzTab line.

    Args:
        line: mzTab line to check

    Returns:
        Line type identifier (e.g., 'MTD', 'PRH', 'PRT', 'PSH', 'PSM') or None

    Example:
        >>> get_mztab_line_type("MTD\tfixed_mod[1]\tCarbamidomethyl")
        'MTD'
        >>> get_mztab_line_type("PRT\tP12345\tTest protein")
        'PRT'
    """
    line = line.strip()
    if not line:
        return None
    parts = line.split("\t")
    return parts[0] if parts else None


def _clean_file_path(file_path: str) -> str:
    """
    Clean a file path to extract just the filename without extension.

    This helper function removes the file:// protocol and extracts the filename
    without the file extension from various file path formats.

    Args:
        file_path: File path that may include file:// protocol and extensions

    Returns:
        Clean filename without protocol and extension

    Example:
        >>> _clean_file_path("file:///path/to/file.mzML")
        'file'
        >>> _clean_file_path("/data/sample_001.mzML")
        'sample_001'
        >>> _clean_file_path("relative/path/file.txt")
        'file'
    """
    if file_path.startswith("file://"):
        # Remove file:// protocol
        clean_path = file_path[7:]
    else:
        clean_path = file_path

    # Extract filename without extension
    filename = clean_path.split("/")[-1].split(".")[0]
    return filename


def _is_ms_run_location_line(line_parts: List[str]) -> bool:
    """
    Check if a line represents an MS run location entry.

    Args:
        line_parts: List of tab-separated parts from an mzTab line

    Returns:
        True if the line is an MS run location entry, False otherwise
    """
    return (
        len(line_parts) >= 3
        and line_parts[0] == "MTD"
        and line_parts[1].split("-")[-1] == "location"
    )


def _extract_ms_run_id_from_key(key: str) -> int:
    """
    Extract MS run numeric index from a metadata key.

    Args:
        key: Metadata key in format 'ms_run[X]-location'

    Returns:
        MS run numeric index (X)

    Example:
        >>> _extract_ms_run_id_from_key("ms_run[1]-location")
        1
        >>> _extract_ms_run_id_from_key("ms_run[10]-location")
        10
    """
    return int(key.split("[")[1].split("]")[0])


def fetch_ms_runs_from_mztab_line(mztab_line: str, ms_runs: dict) -> dict:
    """
    Extract MS run information from a single mzTab metadata line.

    This function parses mzTab metadata lines that contain MS run location information.
    It extracts the MS run numeric index and cleans the file path to return just the filename
    without the file:// protocol and file extension.

    Args:
        mztab_line: A single mzTab metadata line (tab-separated)
        ms_runs: Dictionary to store the extracted MS run information

    Returns:
        Updated ms_runs dictionary with new MS run entries using numeric keys

    Example:
        >>> ms_runs = {}
        >>> line = "MTD\tms_run[1]-location\tfile:///path/to/file.mzML"
        >>> fetch_ms_runs_from_mztab_line(line, ms_runs)
        {1: 'file'}

        >>> line = "MTD\tms_run[2]-location\tfile:///data/sample_001.mzML"
        >>> fetch_ms_runs_from_mztab_line(line, ms_runs)
        {1: 'file', 2: 'sample_001'}

    Note:
        - Only processes lines that start with "MTD" and end with "-location"
        - Removes "file://" protocol from file paths
        - Extracts filename without extension
        - Returns numeric indices as keys for better performance
        - Handles various file path formats gracefully
    """
    mztab_line = mztab_line.strip()
    line_parts = mztab_line.split("\t")

    # Check if this is an MS run location line
    if _is_ms_run_location_line(line_parts):
        ms_run_id = _extract_ms_run_id_from_key(line_parts[1])
        file_path = line_parts[2]
        filename = _clean_file_path(file_path)
        ms_runs[ms_run_id] = filename

    return ms_runs


def extract_ms_runs_from_metadata(metadata_df) -> Dict[int, str]:
    """
    Extract MS runs from metadata DataFrame.

    This function processes a metadata DataFrame to extract MS run information.
    It looks for metadata entries that follow the pattern 'ms_run[X]-location'
    and extracts clean filenames without file:// protocol and extensions.

    Args:
        metadata_df: DataFrame containing metadata with 'key' and 'value' columns

    Returns:
        Dictionary mapping MS run numeric indices to clean filenames

    Example:
        >>> metadata_df = pd.DataFrame({
        ...     'key': ['ms_run[1]-location', 'ms_run[2]-location'],
        ...     'value': ['file:///path/to/file1.mzML', 'file:///path/to/file2.mzML']
        ... })
        >>> extract_ms_runs_from_metadata(metadata_df)
        {1: 'file1', 2: 'file2'}

    Note:
        - Removes "file://" protocol from file paths
        - Extracts filename without extension
        - Returns numeric indices as keys for better performance
        - Handles various file path formats gracefully
    """
    ms_runs = {}

    if metadata_df is None or metadata_df.empty:
        return ms_runs

    # Filter for MS run location entries
    ms_run_mask = metadata_df["key"].str.startswith("ms_run[") & metadata_df[
        "key"
    ].str.endswith("-location")

    ms_run_entries = metadata_df[ms_run_mask]

    for _, row in ms_run_entries.iterrows():
        ms_run_id = _extract_ms_run_id_from_key(row["key"])
        file_path = row["value"]
        filename = _clean_file_path(file_path)
        ms_runs[ms_run_id] = filename

    return ms_runs


def extract_ms_runs_from_lines(mztab_lines: List[str]) -> Dict[int, str]:
    """
    Extract MS runs from a list of mzTab lines.

    This function processes a list of mzTab lines to extract MS run information.
    It's useful when processing mzTab files line by line.

    Args:
        mztab_lines: List of mzTab lines (strings)

    Returns:
        Dictionary mapping MS run numeric indices to clean filenames

    Example:
        >>> lines = [
        ...     "MTD\tms_run[1]-location\tfile:///path/to/file1.mzML",
        ...     "MTD\tms_run[2]-location\tfile:///path/to/file2.mzML"
        ... ]
        >>> extract_ms_runs_from_lines(lines)
        {1: 'file1', 2: 'file2'}

    Note:
        - Uses the same cleaning logic as other MS run extraction functions
        - Processes each line individually using fetch_ms_runs_from_mztab_line
        - Returns numeric indices as keys for better performance
    """
    ms_runs = {}

    for line in mztab_lines:
        ms_runs = fetch_ms_runs_from_mztab_line(line, ms_runs)

    return ms_runs


def parse_pepidoform_with_modifications(
    openms_peptidoform: str,
    openms_modifications: str,
    reference_modifications: dict,
) -> tuple[str, list]:
    """
    Parses a peptide sequence in mzTab format with modifications into a ProForma-like string and a structured list.

    The Peptidoform is a string in the following format: PEPTIDE(Oxidation)R
    The Modification string is a string in the following format:

    Args:
        openms_peptidoform (str): The input peptide sequence string, with modifications enclosed
            in parentheses, e.g.,
            "PEPTIDE(Oxidation)R".
            ".(Nterm mod)PEPTIDE(Oxidation)R"
            "PEPTIDE(Oxidation)R.(Cterm mod)"

        openms_modifications (str): A string in the following format:
            1-UNIMOD:35,2-UNIMOD:35,3-UNIMOD:35
            more complex examples:
            3[MS,MS:1001876, modification probability, 0.8]|4[MS,MS:1001876, modification probability, 0.2] MOD:00412, 8-MOD:00412

        reference_modifications (dict): A dictionary with the following format:
            {
                "Oxidation":
                   ["UNIMOD:35",
                   [index 1, 2, 3],
                   [site 1, site 2, site 3],
                   [position 1, position 2, position 3]
                   ]
            }

    Returns:
        tuple[str, list]: A tuple containing:
            - peptidoform (str): The peptide sequence formatted in ProForma-like
              notation, e.g., "PEPTIDE[UNIMOD:35]R".
            - modification_details (list): A list of dictionaries containing modification details.
              Example:
              [
                  {
                      "name": "Oxidation",
                      "accession": "UNIMOD:35",
                      "positions": [
                          {
                              "position": "M.4",
                              "scores": [
                                  {
                                      "score_name": "localization_probability",
                                      "score_value": 0.99
                                  }
                              ]
                          }
                      ]
                  }
              ]
    """
    # If there is no modification, return the peptidoform as is
    if "(" not in openms_peptidoform:
        return openms_peptidoform, []

    # Converting peptidoform from mzTab OpenMS to ProForma
    peptidoform = parse_peptidoform_openms(openms_peptidoform)

    # Get the pure sequence without modifications for position mapping
    pure_sequence = re.sub(r"\[.*?\]", "", peptidoform)

    # Parsing the modifications
    modification_details = parse_modifications_openms(
        openms_peptidoform=openms_peptidoform,
        openms_modifications=openms_modifications,
        reference_modifications=reference_modifications,
        pure_sequence=pure_sequence,
    )

    return peptidoform, modification_details


def _create_mod_name_to_accession_map(reference_modifications: dict) -> dict:
    """Create a mapping from modification name to accession."""
    return {val[0]: key for key, val in reference_modifications.items()}


def _add_modification_to_parsed_mods(
    parsed_mods: dict, accession: str, mod_name: str, position: str
):
    """Add a modification position to the parsed modifications dictionary."""
    if accession not in parsed_mods:
        parsed_mods[accession] = {"name": mod_name, "positions": []}
    parsed_mods[accession]["positions"].append({"position": position, "scores": None})


def _process_n_terminal_modification(
    peptidoform: str, parsed_mods: dict, mod_name_to_accession: dict
) -> str:
    """Process N-terminal modification and return the remaining peptidoform."""
    if not peptidoform.startswith(".("):
        return peptidoform

    match = re.match(r"\.\(([^)]+)\)", peptidoform)
    if match:
        mod_name = match.group(1)
        accession = mod_name_to_accession.get(mod_name, mod_name)
        _add_modification_to_parsed_mods(parsed_mods, accession, mod_name, "N-term.0")
        return peptidoform[len(match.group(0)) :]

    return peptidoform


def _extract_c_terminal_modification(peptidoform: str) -> tuple[str, str]:
    """Extract C-terminal modification name and return cleaned peptidoform."""
    if not (peptidoform.endswith(")") and peptidoform.rfind(".(") > 0):
        return peptidoform, None

    match = re.search(r"\.\(([^)]+)\)$", peptidoform)
    if match:
        c_term_mod_name = match.group(1)
        cleaned_peptidoform = peptidoform[: -len(match.group(0))]
        return cleaned_peptidoform, c_term_mod_name

    return peptidoform, None


def _process_inline_modifications(
    peptidoform: str, parsed_mods: dict, mod_name_to_accession: dict, pure_sequence: str
):
    """Process inline modifications in the peptidoform."""
    in_mod = False
    mod_buffer = ""
    seq_len = 0

    for char in peptidoform:
        if char == "(":
            in_mod = True
        elif char == ")":
            in_mod = False
            mod_name = mod_buffer
            accession = mod_name_to_accession.get(mod_name, mod_name)
            # Get the amino acid at this position
            aa = pure_sequence[seq_len - 1] if seq_len > 0 else "?"
            position = f"{aa}.{seq_len}"
            _add_modification_to_parsed_mods(parsed_mods, accession, mod_name, position)
            mod_buffer = ""
        elif in_mod:
            mod_buffer += char
        else:
            seq_len += 1


def _process_c_terminal_modification(
    c_term_mod_name: str,
    parsed_mods: dict,
    mod_name_to_accession: dict,
    pure_sequence: str,
):
    """Process C-terminal modification if present."""
    if not c_term_mod_name:
        return

    accession = mod_name_to_accession.get(c_term_mod_name, c_term_mod_name)
    position = f"C-term.{len(pure_sequence) + 1}"
    _add_modification_to_parsed_mods(parsed_mods, accession, c_term_mod_name, position)


def _parse_from_peptidoform(
    openms_peptidoform: str,
    parsed_mods: dict,
    mod_name_to_accession: dict,
    pure_sequence: str,
):
    """Parse modifications from peptidoform string when no explicit modifications are provided."""
    peptidoform_to_process = openms_peptidoform

    # Handle N-terminal modification
    peptidoform_to_process = _process_n_terminal_modification(
        peptidoform_to_process, parsed_mods, mod_name_to_accession
    )

    # Handle C-terminal modification
    peptidoform_to_process, c_term_mod_name = _extract_c_terminal_modification(
        peptidoform_to_process
    )

    # Process inline modifications
    _process_inline_modifications(
        peptidoform_to_process, parsed_mods, mod_name_to_accession, pure_sequence
    )

    # Add C-terminal modification if found
    _process_c_terminal_modification(
        c_term_mod_name, parsed_mods, mod_name_to_accession, pure_sequence
    )


def _format_position_string(pos: int, pure_sequence: str) -> str:
    """Format position string based on position value."""
    if pos == 0:
        return "N-term.0"
    elif pos == len(pure_sequence) + 1:
        return f"C-term.{len(pure_sequence) + 1}"
    else:
        aa = pure_sequence[pos - 1] if 0 < pos <= len(pure_sequence) else "?"
        return f"{aa}.{pos}"


def _parse_from_modifications_string(
    openms_modifications: str,
    parsed_mods: dict,
    reference_modifications: dict,
    pure_sequence: str,
):
    """Parse modifications from explicit modifications string."""
    mods = openms_modifications.split(",")
    for mod in mods:
        parts = mod.split("-")
        accession = parts[-1]
        position_str = "-".join(parts[:-1])

        mod_name = reference_modifications.get(accession, [accession])[0]

        positions_with_probs = re.split(r"\|", position_str)
        for p_w_p in positions_with_probs:
            pos_match = re.match(r"(\d+)", p_w_p)
            if pos_match:
                pos = int(pos_match.group(1))
                formatted_position = _format_position_string(pos, pure_sequence)
                _add_modification_to_parsed_mods(
                    parsed_mods, accession, mod_name, formatted_position
                )


def _format_output_list(parsed_mods: dict) -> list:
    """Format the parsed modifications dictionary into the output list format."""
    output_list = []
    for accession, mod_data in parsed_mods.items():
        output_list.append(
            {
                "name": mod_data["name"],
                "accession": accession,
                "positions": mod_data["positions"],
            }
        )
    return output_list


def parse_modifications_openms(
    openms_peptidoform: str,
    openms_modifications: str,
    reference_modifications: dict,
    pure_sequence: str,
) -> list:
    """
    Parse modifications from an OpenMS peptidoform string.

    Args:
        openms_peptidoform: The peptidoform string containing modifications
        openms_modifications: Optional explicit modifications string
        reference_modifications: Dictionary mapping accessions to modification info
        pure_sequence: The clean peptide sequence without modifications

    Returns:
        List of modification dictionaries with name, accession, and positions
    """
    if "(" not in openms_peptidoform:
        return []

    parsed_mods = {}
    mod_name_to_accession = _create_mod_name_to_accession_map(reference_modifications)

    if not openms_modifications or openms_modifications == "null":
        # Case 1: Parse from peptidoform
        _parse_from_peptidoform(
            openms_peptidoform, parsed_mods, mod_name_to_accession, pure_sequence
        )
    else:
        # Case 2: Parse from explicit modifications string
        _parse_from_modifications_string(
            openms_modifications, parsed_mods, reference_modifications, pure_sequence
        )

    return _format_output_list(parsed_mods)


def parse_peptidoform_openms(sequence: str) -> str:
    """
    Notes
    -----
    Implemented according to the documentation on
    `github.com/OpenMS/OpenMS <https://github.com/OpenMS/OpenMS/blob/8cb90/src/openms/include/OpenMS/CHEMISTRY/AASequence.h>`_
    . The differentiation between square- and round bracket notation is removed after parsing.
    """
    MOD_PATTERN = re.compile(r"\(((?:[^)(]+|\((?:[^)(]+|\([^)(]*\))*\))*)\)")
    MOD_PATTERN_NTERM = re.compile(r"^\.\[((?:[^][]+|\[(?:[^][]+|\[[^][]*\])*\])*)\]")
    MOD_PATTERN_CTERM = re.compile(r"\.\[((?:[^][]+|\[(?:[^][]+|\[[^][]*\])*\])*)\]$")

    sequence = MOD_PATTERN.sub(r"[\1]", sequence)
    if sequence[:2] == ".[":
        sequence = MOD_PATTERN_NTERM.sub(r"[\1]-", sequence)
    if sequence[-1] == "]":
        sequence = MOD_PATTERN_CTERM.sub(r"-[\1]", sequence)
    sequence = sequence.strip(".")
    return sequence
