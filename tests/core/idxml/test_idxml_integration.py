import tempfile
import os
from pathlib import Path
import pytest
from click.testing import CliRunner

from quantmsio.commands.convert.quantms import convert_idxml_psm_cmd

# Test data path
TEST_DATA_ROOT = Path(__file__).parent / ".." / ".." / "examples"
IDXML_TEST_PATH = (
    TEST_DATA_ROOT
    / "idxml/20111219_EXQ5_KiSh_SA_LabelFree_HeLa_Proteome_Control_rep1_pH3_consensus_fdr_filter_pep_luciphor.idXML"
)


def test_idxml_psm_command():
    """Test the idXML PSM conversion command."""
    # Skip if test data doesn't exist
    if not IDXML_TEST_PATH.exists():
        pytest.skip(f"Test data not found: {IDXML_TEST_PATH}")

    runner = CliRunner()

    with tempfile.TemporaryDirectory() as temp_dir:
        # Run the command
        result = runner.invoke(
            convert_idxml_psm_cmd,
            [
                "--idxml-path",
                str(IDXML_TEST_PATH),
                "--output-folder",
                temp_dir,
                "--output-prefix",
                "test-idxml",
                "--verbose",
            ],
        )

        # Check that command succeeded
        assert result.exit_code == 0, f"Command failed with output: {result.output}"

        # Check that output file was created
        output_files = list(Path(temp_dir).glob("*.parquet"))
        assert (
            len(output_files) == 1
        ), f"Expected 1 output file, found {len(output_files)}"

        output_file = output_files[0]
        assert output_file.suffix == ".parquet"
        assert "test-idxml" in output_file.name
        assert output_file.stat().st_size > 0

        print(f"Command succeeded, created file: {output_file}")


def test_idxml_psm_command_missing_file():
    """Test the idXML PSM conversion command with missing input file."""
    runner = CliRunner()

    with tempfile.TemporaryDirectory() as temp_dir:
        # Run the command with non-existent file
        result = runner.invoke(
            convert_idxml_psm_cmd,
            [
                "--idxml-path",
                "non_existent_file.idXML",
                "--output-folder",
                temp_dir,
            ],
        )

        # Command should fail due to missing file
        assert result.exit_code != 0
        assert "does not exist" in result.output or "Error" in result.output


def test_idxml_psm_command_help():
    """Test the idXML PSM conversion command help."""
    runner = CliRunner()

    result = runner.invoke(convert_idxml_psm_cmd, ["--help"])

    assert result.exit_code == 0
    assert "Convert PSM data from idXML to quantms.io format" in result.output
    assert "--idxml-path" in result.output
    assert "--output-folder" in result.output


if __name__ == "__main__":
    # Run tests manually if needed
    test_idxml_psm_command()
    test_idxml_psm_command_missing_file()
    test_idxml_psm_command_help()
    print("All integration tests passed!")
