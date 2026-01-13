"""
General test file for individual functions.

This module contains tests for standalone functions that don't require
full dataset processing or complex setup.
"""

from qpx.utils.mztab_utils import (
    fetch_modifications_from_mztab_line,
    generate_modification_list,
    get_modifications_object_from_mztab_line,
)

metadata_block = (
    "MTD\tfixed_mod[1]\t[UNIMOD, UNIMOD:4, Carbamidomethyl, ]\n"
    "MTD\tfixed_mod[1]-site\tC\n"
    "MTD\tfixed_mod[1]-position\tAnywhere\n"
    "MTD\tvar_mod[1]\t[UNIMOD, UNIMOD:21, Phospho, ]\n"
    "MTD\tvar_mod[1]-site\tS\n"
    "MTD\tvar_mod[1]-position\tAnywhere"
)


def test_fetch_modifications_metadata_block():
    """Test parsing a full metadata block with multiple modifications."""
    modifications = {}
    for line in metadata_block.strip().split("\n"):
        modifications = fetch_modifications_from_mztab_line(line, modifications)
    expected = {
        "UNIMOD:4": ["Carbamidomethyl", "1", "C", "Anywhere"],
        "UNIMOD:21": ["Phospho", "1", "S", "Anywhere"],
    }
    assert modifications == expected


def test_get_modifications_object_from_mztab_line():
    modifications = {}
    for line in metadata_block.strip().split("\n"):
        modifications = fetch_modifications_from_mztab_line(line, modifications)
    expected = {
        "Carbamidomethyl": {"position": [2, 5], "unimod_accession": "UNIMOD:4"},
        "Phospho": {"position": [1, 3, 7], "unimod_accession": "UNIMOD:21"},
    }
    modification_string = "2|5-UNIMOD:4,1|3|7-UNIMOD:21"
    modification_list = get_modifications_object_from_mztab_line(
        modification_string, modifications
    )
    assert modification_list == expected


def test_generate_modification_list():
    modifications = {}
    for line in metadata_block.strip().split("\n"):
        modifications = fetch_modifications_from_mztab_line(line, modifications)
    modification_string = "2|5-UNIMOD:4,1|3|7-UNIMOD:21"
    modification_list = generate_modification_list(modification_string, modifications)
    assert modification_list == ["2|5-UNIMOD:4", "1|3|7-UNIMOD:21"]
