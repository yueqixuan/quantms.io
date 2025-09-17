from pathlib import Path

from quantmsio.core.ae import AbsoluteExpressionHander

TEST_DATA_ROOT = Path(__file__).parents[2] / "examples"


def test_ae():
    """Test absolute expression handler."""
    # Resolve file paths
    project_path, ibaq_path, sdrf_path = (
        TEST_DATA_ROOT / "AE/project.json",
        TEST_DATA_ROOT / "AE/PXD016999.1-ibaq.tsv",
        TEST_DATA_ROOT / "AE/PXD016999-first-instrument.sdrf.tsv",
    )

    # Initialize AbsoluteExpressionHander
    ae_handler = AbsoluteExpressionHander()

    # Load a project file
    ae_handler.load_project_file(str(project_path))
    assert ae_handler.project_manager is not None
    assert "project_accession" in ae_handler.project_manager.project.project_info

    # Load ibaq file
    ae_handler.load_ibaq_file(ibaq_path)
    assert ae_handler.ibaq_df is not None
    assert len(ae_handler.ibaq_df) > 0
    assert "protein" in ae_handler.ibaq_df.columns
    assert "ibaq" in ae_handler.ibaq_df.columns

    # Load sdrf file
    ae_handler.load_sdrf_file(str(sdrf_path))
    assert ae_handler.sdrf_manager is not None
    assert hasattr(ae_handler.sdrf_manager, "get_sample_map")

    # Test generating absolute expression
    assert ae_handler.ibaq_df is not None
    assert len(ae_handler.ibaq_df) > 0
    assert "protein" in ae_handler.ibaq_df
