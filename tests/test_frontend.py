from unittest.mock import patch, MagicMock
import json

@patch("requests.get")
def test_fingerprint_api_mocked_success(mock_get):
    mock_response = MagicMock()
    mock_response.status_code = 200
    mock_response.json.return_value = {
        "standardized_smiles": "CCO",
        "fingerprint": [0, 1, 0, 0, 1] + [0] * 1019,
    }
    mock_get.return_value = mock_response

    from frontend import API_URL  # wherever your Streamlit function is

    smiles = "CCO"
    response = mock_get(f"{API_URL}/getFingerPrint_bySMILES/{smiles}")
    assert response.status_code == 200
    result = response.json()
    assert "fingerprint" in result
    assert result["standardized_smiles"] == "CCO"

@patch("requests.get")
def test_fingerprint_api_mocked_failure(mock_get):
    mock_response = MagicMock()
    mock_response.status_code = 500
    mock_response.json.return_value = {"detail": "Kekulization failed"}
    mock_get.return_value = mock_response

    from frontend import API_URL

    smiles = "INVALID_SMILES"
    response = mock_get(f"{API_URL}/getFingerPrint_bySMILES/{smiles}")
    assert response.status_code == 500
    assert response.json()["detail"] == "Kekulization failed"
