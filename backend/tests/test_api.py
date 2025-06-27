import pytest
from fastapi.testclient import TestClient
from unittest.mock import patch, MagicMock, AsyncMock
import numpy as np
import uuid
from datetime import datetime

from app.utils.processing import smiles_to_features

# Mock model validation and Celery tasks before importing app
mock_models = {
    "classifier": MagicMock(),
    "regressor": MagicMock()
}
mock_models["classifier"].predict.return_value = np.array([1])
mock_models["classifier"].predict_proba.return_value = np.array([[0.3, 0.7]])
mock_models["regressor"].predict.return_value = np.array([0.85])

mock_celery_task = MagicMock()
mock_celery_task.id = str(uuid.uuid4())
mock_celery_task.state = 'SUCCESS'
mock_celery_task.result = {
    'results': [{
        'smiles': 'CC(=O)OC1=CC=CC=C1C(=O)O',
        'prediction': 0.85,
        'confidence': 0.7,
        'classifier_prediction': 1,
        'features': {
            'morgan_fingerprint': [1, 0, 1, 0] * 1050,  # 4200 features
            'descriptors': {
                'MolWt': 180.16,
                'LogP': 1.19,
                'TPSA': 63.6,
                'NumHDonors': 1,
                'NumHAcceptors': 4,
                'NumRotatableBonds': 3,
                'NumAromaticRings': 1
            }
        }
    }],
    'total_processed': 1,
    'successful': 1,
    'failed': 0
}

with patch("app.utils.validation.validate_models", return_value=mock_models), \
     patch("app.worker.celery_app.send_task", return_value=mock_celery_task), \
     patch("app.worker.celery_app.AsyncResult", return_value=mock_celery_task):
    from app.main import app

client = TestClient(app)

# Test data
VALID_SMILES = "CC(=O)OC1=CC=CC=C1C(=O)O"  # Aspirin
INVALID_SMILES = "NOT_A_SMILES"


def test_smiles_to_features_valid():
    """Test SMILES to features conversion with valid input."""
    features = smiles_to_features(VALID_SMILES)
    assert isinstance(features, np.ndarray)
    assert features.dtype == np.int8
    assert features.shape == (4200,)


def test_smiles_to_features_invalid():
    """Test SMILES to features conversion with invalid input."""
    with pytest.raises(ValueError):
        smiles_to_features(INVALID_SMILES)


def test_submit_prediction_job_single():
    """Test GraphQL mutation to submit a single prediction job."""
    query = """
    mutation SubmitPredictionJob($input: PredictionJobInput!) {
        submitPredictionJob(jobInput: $input) {
            jobId
            status
            createdAt
            progress
            error
        }
    }
    """
    variables = {
        "input": {
            "smilesList": [VALID_SMILES],
            "jobName": "test_single_prediction"
        }
    }
    
    response = client.post("/graphql", json={"query": query, "variables": variables})
    assert response.status_code == 200
    
    data = response.json()
    assert "data" in data
    assert "submitPredictionJob" in data["data"]
    
    job_data = data["data"]["submitPredictionJob"]
    assert job_data["jobId"] is not None
    assert job_data["status"] == "submitted"
    assert job_data["error"] is None


def test_submit_prediction_job_batch():
    """Test GraphQL mutation to submit a batch prediction job."""
    query = """
    mutation SubmitPredictionJob($input: PredictionJobInput!) {
        submitPredictionJob(jobInput: $input) {
            jobId
            status
            createdAt
            progress
            error
        }
    }
    """
    variables = {
        "input": {
            "smilesList": [VALID_SMILES, "CCO", "CCC"],  # Multiple SMILES
            "jobName": "test_batch_prediction"
        }
    }
    
    response = client.post("/graphql", json={"query": query, "variables": variables})
    assert response.status_code == 200
    
    data = response.json()
    job_data = data["data"]["submitPredictionJob"]
    assert job_data["jobId"] is not None
    assert job_data["status"] == "submitted"
    assert "3 compounds" in job_data["progress"]


def test_submit_prediction_job_empty_list():
    """Test GraphQL mutation with empty SMILES list should return error."""
    query = """
    mutation SubmitPredictionJob($input: PredictionJobInput!) {
        submitPredictionJob(jobInput: $input) {
            jobId
            status
            error
        }
    }
    """
    variables = {
        "input": {
            "smilesList": [],
            "jobName": "test_empty"
        }
    }
    
    response = client.post("/graphql", json={"query": query, "variables": variables})
    assert response.status_code == 200
    
    data = response.json()
    job_data = data["data"]["submitPredictionJob"]
    assert job_data["status"] == "error"
    assert "cannot be empty" in job_data["error"]


def test_get_prediction_result():
    """Test GraphQL query to get prediction results."""
    query = """
    query GetPredictionResult($jobId: String!) {
        getPredictionResult(jobId: $jobId) {
            ... on JobResult {
                status
                results {
                    smiles
                    prediction
                    confidence
                    classifierPrediction
                    features {
                        morganFingerprint
                        descriptors {
                            molWt
                            logP
                            tpsa
                            numHDonors
                            numHAcceptors
                            numRotatableBonds
                            numAromaticRings
                        }
                    }
                    error
                }
                totalProcessed
                successful
                failed
                jobId
                createdAt
                completedAt
            }
            ... on JobStatus {
                jobId
                status
                createdAt
                progress
                error
            }
        }
    }
    """
    variables = {"jobId": mock_celery_task.id}
    
    response = client.post("/graphql", json={"query": query, "variables": variables})
    assert response.status_code == 200
    
    data = response.json()
    assert "data" in data
    result_data = data["data"]["getPredictionResult"]
    
    assert result_data["status"] == "completed"
    assert result_data["totalProcessed"] == 1
    assert result_data["successful"] == 1
    assert result_data["failed"] == 0
    assert len(result_data["results"]) == 1
    
    prediction = result_data["results"][0]
    assert prediction["smiles"] == VALID_SMILES
    assert prediction["prediction"] == 0.85
    assert prediction["confidence"] == 0.7
    assert prediction["classifierPrediction"] == 1
    assert prediction["features"] is not None


def test_get_job_status():
    """Test GraphQL query to get job status."""
    query = """
    query GetJobStatus($jobId: String!) {
        getJobStatus(jobId: $jobId) {
            jobId
            status
            createdAt
            progress
            error
        }
    }
    """
    variables = {"jobId": mock_celery_task.id}
    
    response = client.post("/graphql", json={"query": query, "variables": variables})
    assert response.status_code == 200
    
    data = response.json()
    status_data = data["data"]["getJobStatus"]
    
    assert status_data["jobId"] == mock_celery_task.id
    assert status_data["status"] == "completed"
    assert status_data["error"] is None


def test_get_nonexistent_job():
    """Test GraphQL query with non-existent job ID."""
    fake_job_id = str(uuid.uuid4())
    
    # Create a mock for non-existent job
    mock_nonexistent_task = MagicMock()
    mock_nonexistent_task.state = 'PENDING'
    
    query = """
    query GetJobStatus($jobId: String!) {
        getJobStatus(jobId: $jobId) {
            jobId
            status
            error
        }
    }
    """
    variables = {"jobId": fake_job_id}
    
    with patch("app.worker.celery_app.AsyncResult", return_value=mock_nonexistent_task):
        response = client.post("/graphql", json={"query": query, "variables": variables})
        assert response.status_code == 200
        
        data = response.json()
        status_data = data["data"]["getJobStatus"]
        assert status_data["status"] == "pending"
