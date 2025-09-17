from pathlib import Path
from unittest import TestCase
from unittest.mock import patch

import requests
from ddt import ddt

from quantmsio.core.project import ProjectHandler

TEST_DATA_ROOT = Path(__file__).parents[1] / "examples"
test_dataset = (
    "PXD007683",
    TEST_DATA_ROOT / "quantms/dda-plex-full/PXD007683-TMT.sdrf.tsv",
)


@ddt
class TestProjectHandler(TestCase):

    @patch("requests.get")
    def test_populate_from_pride_archive_successful(self, mock_get):
        # Mock the API response for a successful request
        mock_response = requests.models.Response()
        mock_response.status_code = 200
        mock_response.json = lambda: {
            "title": "Test Project",
            "projectDescription": "Test description",
            "sampleProcessingProtocol": "Test sample processing protocol",
            "dataProcessingProtocol": "Test data processing protocol",
            "references": [{"pubmedId": "12345678"}],
            "keywords": ["keyword1", "keyword2"],
            "projectTags": ["tag1", "tag2"],
            # Add other mock data as needed
        }
        mock_get.return_value = mock_response

        project_accession = "PXD123456"
        project_manager = ProjectHandler(project_accession)
        project_manager.populate_from_pride_archive()

        # Assert that the project_info has been updated
        self.assertEqual(
            project_manager.project.project_info["project_title"], "Test Project"
        )
        self.assertEqual(
            project_manager.project.project_info["project_description"],
            "Test description",
        )

    def test_populate_from_pride_archive_api(self):
        project_accession = "PXD020453"
        project_manager = ProjectHandler(project_accession)
        project_manager.populate_from_pride_archive()

        # Assert that the project_info has been updated
        self.assertEqual(
            project_manager.project.project_info["project_title"],
            "Structural insights into Cullin4-RING ubiquitin ligase remodelling by Vpr from simian immunodeficiency viruses",
        )
        print(project_manager.project.project_info)

    def test_save_project_info(self):
        project_accession = test_dataset[0]
        sdrf_file = test_dataset[1]

        project_manager = ProjectHandler(project_accession)
        project_manager.populate_from_pride_archive()
        project_manager.populate_from_sdrf(sdrf_file)
