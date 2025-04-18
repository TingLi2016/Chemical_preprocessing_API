import pytest
from fastapi.testclient import TestClient
from app.main import app  

client = TestClient(app)

def test_fingerprint_api_valid_smiles():
    response = client.get("/getFingerPrint_bySMILES/CCO")
    assert response.status_code == 200
    data = response.json()
    assert "standardized_smiles" in data
    assert "fingerprint" in data
    assert isinstance(data["fingerprint"], list)
    assert all(bit in [0, 1] for bit in data["fingerprint"])

def test_fingerprint_api_invalid_smiles():
    response = client.get("/getFingerPrint_bySMILES/INVALID_SMILES")
    assert response.status_code == 500 or response.status_code == 422  # either FastAPI validation or handled error
    data = response.json()
    assert "detail" in data
