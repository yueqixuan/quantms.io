"""Tests for PRIDE API client."""

import json
import pytest
from unittest.mock import patch, MagicMock

from qpx.core.pride import fetch_pride_metadata


# Sample PRIDE API v3 response
MOCK_PRIDE_RESPONSE = {
    "accession": "PXD014414",
    "title": "Quantitative proteomic landscape of metaplastic breast carcinoma",
    "projectDescription": "This study provides a comprehensive proteomic analysis...",
    "doi": "10.6019/PXD014414",
    "submissionDate": "2019-06-25",
    "publicationDate": "2020-04-08",
    "keywords": ["Proteomics", "Breast cancer"],
    "organisms": [
        {"@type": "CvParam", "name": "Homo sapiens (human)", "accession": "NEWT:9606"}
    ],
    "instruments": [
        {"@type": "CvParam", "name": "Q Exactive HF", "accession": "MS:1002523"}
    ],
    "references": [
        {"pubmedID": 32265444, "doi": "10.1038/s41467-020-15283-z"}
    ],
}


class TestFetchPrideMetadata:
    """Test PRIDE API client with mocked HTTP responses."""

    @patch("qpx.core.pride.urlopen")
    def test_basic_fetch(self, mock_urlopen):
        """Should parse all fields from PRIDE API response."""
        mock_resp = MagicMock()
        mock_resp.read.return_value = json.dumps(MOCK_PRIDE_RESPONSE).encode()
        mock_resp.__enter__ = lambda s: s
        mock_resp.__exit__ = MagicMock(return_value=False)
        mock_urlopen.return_value = mock_resp

        result = fetch_pride_metadata("PXD014414")

        assert result["project_title"] == "Quantitative proteomic landscape of metaplastic breast carcinoma"
        assert "proteomic analysis" in result["project_description"]
        assert result["pubmed_id"] == "32265444"
        assert result["doi"] == "10.6019/PXD014414"
        assert result["organisms"] == ["Homo sapiens (human)"]
        assert result["instruments"] == ["Q Exactive HF"]
        assert result["keywords"] == ["Proteomics", "Breast cancer"]
        assert result["submission_date"] == "2019-06-25"
        assert result["publication_date"] == "2020-04-08"

    @patch("qpx.core.pride.urlopen")
    def test_no_references(self, mock_urlopen):
        """Should handle missing references gracefully."""
        response = {**MOCK_PRIDE_RESPONSE, "references": []}
        mock_resp = MagicMock()
        mock_resp.read.return_value = json.dumps(response).encode()
        mock_resp.__enter__ = lambda s: s
        mock_resp.__exit__ = MagicMock(return_value=False)
        mock_urlopen.return_value = mock_resp

        result = fetch_pride_metadata("PXD014414")
        assert result["pubmed_id"] is None

    @patch("qpx.core.pride.urlopen")
    def test_no_organisms(self, mock_urlopen):
        """Should handle missing organisms gracefully."""
        response = {**MOCK_PRIDE_RESPONSE, "organisms": None}
        mock_resp = MagicMock()
        mock_resp.read.return_value = json.dumps(response).encode()
        mock_resp.__enter__ = lambda s: s
        mock_resp.__exit__ = MagicMock(return_value=False)
        mock_urlopen.return_value = mock_resp

        result = fetch_pride_metadata("PXD014414")
        assert result["organisms"] == []

    @patch("qpx.core.pride.urlopen")
    def test_404_raises_value_error(self, mock_urlopen):
        """Should raise ValueError for unknown accession."""
        from urllib.error import HTTPError
        mock_urlopen.side_effect = HTTPError(
            url="", code=404, msg="Not Found", hdrs=None, fp=None
        )
        with pytest.raises(ValueError, match="not found"):
            fetch_pride_metadata("PXD999999")

    @patch("qpx.core.pride.urlopen")
    def test_network_error_raises_connection_error(self, mock_urlopen):
        """Should raise ConnectionError for network failures."""
        from urllib.error import URLError
        mock_urlopen.side_effect = URLError("Connection refused")
        with pytest.raises(ConnectionError, match="Cannot reach"):
            fetch_pride_metadata("PXD014414")


class TestDatasetEnrichFromPride:
    """Test Dataset.enrich_from_pride() with mocked API."""

    @patch("qpx.core.pride.urlopen")
    def test_enrich_updates_dataset_parquet(self, mock_urlopen, tmp_path):
        """enrich_from_pride should update dataset.parquet fields."""
        mock_resp = MagicMock()
        mock_resp.read.return_value = json.dumps(MOCK_PRIDE_RESPONSE).encode()
        mock_resp.__enter__ = lambda s: s
        mock_resp.__exit__ = MagicMock(return_value=False)
        mock_urlopen.return_value = mock_resp

        # Create a minimal dataset with just dataset.parquet
        from qpx.writers.dataset import DatasetWriter
        from qpx.dataset import Dataset
        from datetime import datetime

        dataset_path = tmp_path / "test.dataset.parquet"
        with DatasetWriter(dataset_path) as w:
            w.write_batch([{
                "project_accession": "PXD014414",
                "creation_date": datetime.now().isoformat(),
            }])

        ds = Dataset(tmp_path)
        result = ds.enrich_from_pride()

        assert result["project_title"] is not None

        # Re-open and check the updated dataset.parquet
        ds2 = Dataset(tmp_path)
        meta_df = ds2.dataset_meta.to_df()
        assert meta_df["project_title"].iloc[0] == "Quantitative proteomic landscape of metaplastic breast carcinoma"
        assert meta_df["pubmed_id"].iloc[0] == "32265444"
        ds.close()
        ds2.close()
