import pytest
from fastapi.testclient import TestClient
from unittest.mock import patch
import numpy as np

from app.main import app
from app.utils.processing import smiles_to_features
from app.models import SMILESInput

client = TestClient(app)

# Test data
VALID_SMILES = "CC(=O)OC1=CC=CC=C1C(=O)O"  # Aspirin
INVALID_SMILES = "NOT_A_SMILES"


@pytest.fixture
def mock_onnx_session():
    with patch("onnxruntime.InferenceSession") as mock_session:
        # Mock prediction output shape: [batch_size, 2] for binary classification
        mock_session.return_value.run.return_value = [np.array([[0.3, 0.7]])]
        yield mock_session


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


def test_predict_single_valid(mock_onnx_session):
    response = client.post("/predict/single", json={"smiles": VALID_SMILES})
    assert response.status_code == 200
    data = response.json()
    assert "prediction" in data
    assert "probability" in data
    assert data["smiles"] == VALID_SMILES


def test_predict_single_invalid(mock_onnx_session):
    response = client.post("/predict/single", json={"smiles": INVALID_SMILES})
    assert response.status_code == 200
    data = response.json()
    assert data["error"] is not None


def test_predict_batch_valid(mock_onnx_session):
    csv_content = "smiles\n" + VALID_SMILES + "\n" + VALID_SMILES
    response = client.post("/predict/batch", files={"file": ("test.csv", csv_content.encode(), "text/csv")})
    assert response.status_code == 200
    data = response.json()
    assert len(data) == 2
    assert all("prediction" in item for item in data)
