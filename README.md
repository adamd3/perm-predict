# Perm-Predict

A web application for AI-based prediction of chemical accumulation in bacteria. Perm-Predict takes chemical compounds (SMILES strings) as input and can predict molecule permeability, model confidence scores and perform molecular feature analysis.

## Overview

Perm-Predict uses a sophisticated two-stage machine learning approach:

1. **Binary Classifier**: First determines if a compound has "near-zero" accumulation
2. **Regression Model**: If not near-zero, predicts the specific permeability level using an ensemble of complementary models (XGBoost, AttentiveFP, DimeNet++)

The application converts SMILES strings into rich molecular representations including learned graph embeddings, molecular descriptors, and Morgan fingerprints for accurate predictions.

## Architecture

The repo consists of two main packages:

```
perm-predict/
├── backend/           # FastAPI + GraphQL + Celery
│   ├── app/
│   │   ├── main.py           # GraphQL API server (FastAPI + Strawberry)
│   │   ├── worker.py         # Celery worker for ML predictions
│   │   ├── schema.py         # GraphQL types, queries, mutations
│   │   ├── models.py         # Pydantic models
│   │   ├── config.py         # Application settings
│   │   ├── ml_models/        # Serialized ML models (.pkl)
│   │   └── utils/
│   │       ├── processing.py # Feature extraction
│   │       ├── validation.py # Model validation
│   │       └── logger.py     # Logging utilities
│   ├── tests/                # Pytest test suite
│   └── requirements.txt      # Python dependencies
├── frontend/          # Next.js + TypeScript + Tailwind
│   ├── src/
│   │   ├── app/             # Next.js app router
│   │   └── components/      # React components
│   └── package.json         # Node dependencies
└── docker-compose.yml       # Multi-service orchestration
```

### Backend (FastAPI + Celery + GraphQL)

**Asynchronous Architecture**: The backend uses a decoupled design to handle computationally intensive ML predictions without blocking the API.

- **GraphQL API Server** (`main.py`): FastAPI application serving GraphQL via Strawberry
- **Celery Worker** (`worker.py`): Handles prediction jobs asynchronously
- **Redis**: Message broker and result store for Celery
- **ML Pipeline**: Two-stage classifier → regressor with rich feature extraction

### Frontend (Next.js + TypeScript + Tailwind)

- **Modern React**: Next.js 14+ with App Router
- **GraphQL Client**: Communicates with backend via GraphQL queries/mutations
- **UI Components**: Built with Tailwind CSS and shadcn/ui
- **Molecular Visualization**: Interactive display of results and molecular features

### Data Flow

1. **Submit**: User submits SMILES string(s) via frontend
2. **Queue**: Frontend sends `submitPredictionJob` GraphQL mutation → FastAPI adds job to Celery queue
3. **Process**: Celery worker processes ML pipeline (feature extraction → classifier → regressor)
4. **Poll**: Frontend periodically queries `getPredictionResult` using job ID
5. **Display**: Results shown with predictions, confidence scores, and molecular features

## Quick Start

### Using Docker (Recommended)

```bash
# Start all services (API, worker, Redis, frontend)
docker-compose up --build
```

Access the application at `http://localhost:3000`

### Manual Development Setup

**Backend:**

```bash
cd backend
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate
pip install -r requirements.txt

# Terminal 1: Start API server
uvicorn main:app --reload

# Terminal 2: Start Celery worker
celery -A worker.celery_app worker --loglevel=info
```

**Frontend:**

```bash
cd frontend
npm install
npm run dev
```

## Configuration

### Backend Environment Variables (`backend/.env`)

```
CELERY_BROKER_URL=redis://redis:6379/0
CELERY_RESULT_BACKEND=redis://redis:6379/0
MODEL_CLASSIFIER_PATH=app/ml_models/classifier_model.pkl
MODEL_REGRESSOR_PATH=app/ml_models/regressor_model.pkl
```

### Frontend Environment Variables (`frontend/.env.local`)

```
NEXT_PUBLIC_GRAPHQL_ENDPOINT=http://localhost:8000/graphql
```

## API

The application provides a GraphQL API with the following main operations:

- **`submitPredictionJob`**: Submit SMILES strings for prediction
- **`getPredictionResult`**: Get complete prediction results by job ID
- **`getJobStatus`**: Check current status of a prediction job

Example GraphQL query:

```graphql
mutation {
  submitPredictionJob(
    jobInput: {
      smilesList: ["CC(=O)OC1=CC=CC=C1C(=O)O"]
      jobName: "aspirin_test"
    }
  ) {
    jobId
    status
  }
}
```

## Testing

```bash
cd backend
pytest tests/
```
