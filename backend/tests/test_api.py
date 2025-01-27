import pytest
from fastapi.testclient import TestClient
from unittest.mock import patch, MagicMock
import numpy as np

from app.utils.processing import smiles_to_features

# Mock the ONNX session before importing app
mock_session = MagicMock()
mock_session.get_inputs.return_value = [MagicMock(name="float_input")]
mock_session.run.return_value = [np.array([[0.3, 0.7]])]

with patch("onnxruntime.InferenceSession", return_value=mock_session):
    from app.main import app

client = TestClient(app)

# Test data
VALID_SMILES = "CC(=O)OC1=CC=CC=C1C(=O)O"  # Aspirin
INVALID_SMILES = "NOT_A_SMILES"


def test_health_check():
    response = client.get("/health")
    assert response.status_code == 200
    assert "status" in response.json()
    assert "model_loaded" in response.json()


def test_smiles_to_features_valid():
    features = smiles_to_features(VALID_SMILES)
    assert isinstance(features, np.ndarray)
    assert features.dtype == np.float32
    assert features.shape == (4200,)  # Based on your FEATURE_COUNT


def test_smiles_to_features_invalid():
    with pytest.raises(ValueError):
        smiles_to_features(INVALID_SMILES)


def test_predict_single_valid():
    response = client.post("/predict/single", json={"smiles": VALID_SMILES})
    assert response.status_code == 200
    data = response.json()
    assert "prediction" in data
    assert "probability" in data
    assert data["smiles"] == VALID_SMILES
    # Test the probability is between 0 and 1
    assert 0 <= data["probability"] <= 1


def test_predict_single_invalid():
    response = client.post("/predict/single", json={"smiles": INVALID_SMILES})
    assert response.status_code == 200
    data = response.json()
    assert data["error"] is not None


@pytest.mark.asyncio
async def test_predict_batch():
    csv_content = "smiles\n" + VALID_SMILES + "\n" + VALID_SMILES
    response = client.post("/predict/batch", files={"file": ("test.csv", csv_content.encode(), "text/csv")})
    assert response.status_code == 200
    data = response.json()
    assert len(data) == 2
    assert all("prediction" in item for item in data)
    assert all("probability" in item for item in data)
    assert all(0 <= item["probability"] <= 1 for item in data)


@pytest.mark.asyncio
async def test_predict_batch_invalid_file():
    response = client.post("/predict/batch", files={"file": ("test.txt", b"not a csv", "text/plain")})
    assert response.status_code == 400
    assert "Only CSV files are supported" in response.json()["detail"]
