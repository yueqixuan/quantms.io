from pathlib import Path

from qpx.core.sdrf import SDRFHandler

TEST_DATA_ROOT = Path(__file__).parents[1] / "examples"


def test_load_sdrf_info():
    """Test loading SDRF information."""
    # Resolve file path
    file = TEST_DATA_ROOT / "quantms/dda-lfq-full/PXD007683-LFQ.sdrf.tsv"

    # Initialize SDRF handler
    sdrf_handler = SDRFHandler(file)
    sdrf_handler._load_sdrf_info(file)

    # Assert that the SDRF information has been loaded correctly
    assert sdrf_handler.get_organisms() == ["Homo sapiens"]
    assert sdrf_handler.get_instruments() == ["Orbitrap Fusion Lumos"]
    assert sdrf_handler.get_diseases() == ["Neuroblastoma"]
    assert sdrf_handler.get_enzymes() == ["Trypsin/P"]
    assert sdrf_handler.get_cell_lines() == ["SH-SY5Y"]
    assert sdrf_handler.get_acquisition_properties() == [
        {"proteomics data acquisition method": "Label free"},
        {"proteomics data acquisition method": "Data-dependent acquisition"},
        {"dissociation method": "HCD"},
        {"precursor mass tolerance": "20 ppm"},
        {"fragment mass tolerance": "0.4 Da"},
    ]

    # Get experiment type
    experiment_type = sdrf_handler.get_experiment_type_from_sdrf()
    assert experiment_type is not None
    assert experiment_type in [
        "LFQ",
        "SILAC",
        "TMT6",
        "TMT10",
        "TMT11",
        "TMT16",
        "ITRAQ4",
        "ITRAQ8",
    ]


def test_get_labels():
    """Test getting labels from SDRF."""
    # Resolve file path
    file = TEST_DATA_ROOT / "quantms/dda-plex-full/PXD007683-TMT.sdrf.tsv"

    # Initialize SDRF handler
    sdrf_handler = SDRFHandler(file)

    # Get sample labels
    labels = sdrf_handler.get_sample_labels()
    assert labels is not None
    assert len(labels) == 11

    # Get experiment type
    experiment_type = sdrf_handler.get_experiment_type_from_sdrf()
    assert experiment_type is not None
    assert experiment_type in [
        "LFQ",
        "SILAC",
        "TMT6",
        "TMT10",
        "TMT11",
        "TMT16",
        "ITRAQ4",
        "ITRAQ8",
    ]
    assert experiment_type == "TMT11"


def test_get_sample_map():
    """Test getting sample map from SDRF."""
    # Resolve file path
    file = TEST_DATA_ROOT / "quantms/dda-lfq-full/PXD007683-LFQ.sdrf.tsv"

    # Initialize SDRF handler
    sdrf_handler = SDRFHandler(file)

    # Get sample map
    sample_map = sdrf_handler.get_sample_map()
    assert sample_map is not None
    assert isinstance(sample_map, dict)
    assert len(sample_map) > 0

    # Get sample map run
    sample_map_run = sdrf_handler.get_sample_map_run()
    assert sample_map_run is not None
    assert isinstance(sample_map_run, dict)
    assert len(sample_map_run) > 0


def test_get_mods_dict():
    """Test getting modifications dictionary from SDRF."""
    # Resolve file path
    file = TEST_DATA_ROOT / "quantms/dda-lfq-full/PXD007683-LFQ.sdrf.tsv"

    # Initialize SDRF handler
    sdrf_handler = SDRFHandler(file)

    # Get modifications dictionary
    mods_dict = sdrf_handler.get_mods_dict()
    assert mods_dict is not None
    assert isinstance(mods_dict, dict)

    # Check if common modifications are included
    expected_mods = ["Carbamidomethyl", "Oxidation"]
    for mod in expected_mods:
        assert any(
            mod in key for key in mods_dict.keys()
        ), f"Expected modification {mod} not found in mods_dict"
