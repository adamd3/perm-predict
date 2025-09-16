# CLAUDE.md

## Project Overview

Perm-Predict is a web application for advanced, machine learning-based prediction of chemical accumulation in bacteria. It takes a single chemical compound (represented as a SMILES string) or list of chemical compounds (SMILES strings) as input and uses a machine learning-based pipeline to predict its permeance (accumulation in bacterial cells). Additionally, the app is able to suggest potential modifications or alternative compounds, when a user uploads their compound(s). It suggests improvements by finding analogs of the compound in the Enamine library, querying the classification model with those and seeing if the analogs are more likely permeant.

## Architecture

The project is a **monorepo** with two main packages: `/backend` and `/frontend`.

Performing the two-step prediction and the calibrated confidence intervals depends on a number of files. The ml_models directory should contain:

- **The Trained Classifier Model**: - `./backend/app/ml_models/classifier_model.pkl`. When the Celery worker processes a job, it will load this file.
- **The Script used to prepare inputs for the classifier model from SMILES strings**: Initially `./backend/app/ml_models/input_generation.py`, this has been updated to use `./backend/app/ml_models/alvadesc_feature_generation.py` for generating chemical descriptors from SMILES strings. This script will integrate with the `alvadesccliwrapper` package.

## Model Output Format

Model Type: XGBoost Classifier (binary classification)

predict() method:

- Returns: numpy.ndarray with shape (n_samples,)
- Values: 0 (impermeant) or 1 (permeant)
- Data type: int64

predict_proba() method:

- Returns: numpy.ndarray with shape (n_samples, 2)
- Values: [[prob_class_0, prob_class_1]]
- Data type: float32
- Example: [[0.994, 0.006]] = 99.4% confidence impermeant, 0.6% confidence permeant

How This Fits Our Code:

Our current worker.py implementation:

classifier_pred = classifier_model.predict(feature_vector)[0] # Returns 0 or 1
classifier_prob = classifier_model.predict_proba(feature_vector)[0] # Returns [prob_0, prob_1]

The confidence calculation takes the maximum probability:
confidence = max(classifier_prob) # Highest probability = confidence

Key Requirements:

- Feature Count: Model expects exactly 10,160 features (your feature extraction must match this)
- Classes: 0 = impermeant, 1 = permeant
- XGBoost Dependency: Need to ensure xgboost is in requirements.txt

The current implementation in worker.py correctly handles this XGBoost classifier's output format!

### Backend (FastAPI + Celery + GraphQL)

The backend uses a decoupled, asynchronous architecture to handle computationally intensive ML predictions.

- **GraphQL API Server**: `backend/main.py` - A FastAPI application that serves the GraphQL API using Strawberry. It handles incoming mutations and queries.
- **Asynchronous Worker**: `backend/worker.py` - A Celery worker that consumes prediction jobs from a queue. This is where the heavy ML model inference happens.
- **GraphQL Schema**: `backend/schema.py` - Defines all GraphQL types, queries, and mutations.
- **Message Broker / Result Store**: Redis is used by Celery to manage the task queue and store prediction results.
- **ML Models**: The ML pipeline is comprised of an ML model in `backend/ml_models/`.

### Frontend (Next.js + TypeScript + Tailwind)

- **Main Page**: `frontend/src/app/page.tsx` - The application's entry point, containing UI components for submitting prediction jobs and a menu that allows access to a graphical chemical interface (something like Ketcher, not sure if there's a better alternative?) which allows the user to design chemical compounds on the web app and export them in SMILES and potentially other formats.
- **GraphQL Client**: The frontend uses a GraphQL client (like Apollo Client) to communicate with the backend API.
- **Components**: `frontend/src/components/` - React components for the input form, results display, loading states, and molecular viewers.
- **Styling**: Use Tailwind CSS with custom theme configuration and shadcn/ui components.

### Key Data Flow (Asynchronous Polling)

1.  User submits a SMILES string via the Next.js frontend.
2.  The frontend sends a `submitPredictionJob` **GraphQL Mutation** to the backend.
3.  The FastAPI server adds the job to the Celery queue and immediately returns a unique `jobId`.
4.  The frontend then periodically sends a `getPredictionResult` **GraphQL Query** using the `jobId`.
5.  A Celery worker picks up the job, runs the full ML prediction pipeline, and saves the result to Redis.
6.  Once the result is available, the `getPredictionResult` query returns the complete JSON data, and the frontend displays it.

## Development Commands

### Backend

```bash
# In one terminal, from the root directory
cd backend && source venv/bin/activate && uvicorn main:app --reload

# In a second terminal, from the root directory
cd backend && source venv/bin/activate && celery -A worker.celery_app worker --loglevel=info
```

### Frontend

```bash
cd frontend
npm run dev
```

### Docker (recommended)

```bash
# This single command starts the API, worker, and Redis services
docker-compose up --build
```

## Key Configuration

### Backend Environment Variables (`backend/.env`)

- `CELERY_BROKER_URL`: URL for the Celery message broker (e.g., `redis://redis:6379/0`)
- `CELERY_RESULT_BACKEND`: URL for the Celery result backend (e.g., `redis://redis:6379/0`)

**ML Model Paths:**

- `MODEL_CLASSIFIER_PATH`: (`./backend/app/ml_models/classifier_model.pkl`).

### Frontend Environment Variables (`frontend/.env.local`)

- `NEXT_PUBLIC_GRAPHQL_ENDPOINT`: The full URL for the backend GraphQL API (e.g., `http://localhost:8000/graphql`)

## Testing

Backend tests are located in `backend/tests/` and use pytest. The test configuration is in `backend/pyproject.toml`.

**Run Tests:**

```bash
cd backend && source venv/bin/activate && python -m pytest tests/ -v
```

## Rate Limiting

The GraphQL API should implement rate limiting to prevent abuse. Suggested limits:

- `submitPredictionJob` mutation: 20/minute
- `getPredictionResult` query: 100/minute

## Model Details

The application uses a machine learning pipeline to predict the permeability of a chemical compound (or multiple compounds).

- **Input Features**: The pipeline begins by converting SMILES strings into rich representations. This now primarily involves generating molecular descriptors using the `alvaDesc` Python package (via `alvadesc_feature_generation.py`), which will eventually replace or complement other feature types like Morgan fingerprints.
- **Classification / prediction**: the classifier in ./backend/app/ml_models/classifier_model.pkl is used to make a prediction on the permeance of a compound (i.e. its accumulation in bacterial cells).

## Plans for development

To start with, we will focus only using our classification model to obtain the probability of classify a user's compound (or list of compounds) as permeant/impermeant.

The plan is to build the application in logical stages, following a standard "crawl, walk, run" development approach:

- **Crawl (Phase 1-3)**: Build the core prediction engine. This involves creating the robust, asynchronous backend and a functional frontend that can take a user's molecule (SMILES format) and return a high-quality prediction with a confidence interval. This is our current focus. A reliable predictor is the essential foundation upon which everything else is built.
- **Walk (Phase 4)**: Implement the generative suggestions. Once the core predictor is built, tested, and deployed, the next major step is to build out this feature. This involves:

1.  Adding the logic to the Celery worker to propose modifications or alternative molecules (Enamine library).
2.  Creating new GraphQL endpoints to serve these suggestions.
3.  Building the UI components in the Next.js frontend to display the suggested molecular modifications to the user.

## Current Status & Immediate Next Steps

**✅ Phase 1 Backend Core - COMPLETED:**

- ✅ Comprehensive feature extraction with validation and normalization (now using `alvadesc_feature_generation.py` for descriptor generation)
- ✅ GraphQL API with async job processing via Celery
- ✅ Core backend tests for feature processing functionality

**🚧 TODO:**

- Test full pipeline with actual model files

**📋 Phase 1 Remaining Tasks:**

- Complete GraphQL API testing with model integration
- Performance optimization and error handling refinement
- **Phase 2**: Develop the core Next.js frontend UI and data-fetching logic
- **Phase 3**: Deploy the application to a cloud provider with GPU support
- **Phase 4**: Implement generative suggestions feature

Steps of Phase 4:

1.  Develop Generative Backend Logic:

- Create a new Celery task in worker.py that suggests alternative molecules that are more likely to be permeant.

2.  Extend the GraphQL Schema:

- Add a new mutation (e.g., submitImprovementJob) and a corresponding query to our schema.py to expose this new functionality.

3.  Build the "Suggestions" UI:

- Create new components in our Next.js frontend that allow users to request and visualize the suggested molecular improvements.

### What We Can Build Now:

The next logical step is to begin Phase 2: Develop the core Next.js frontend UI and data-fetching logic.

This involves building the user interface for submitting SMILES strings and displaying the prediction results. While we proceed with the frontend, it will be crucial to thoroughly test the integrated
backend pipeline (GraphQL API, Celery worker, and the feature generation script) to ensure it's robust and reliable. We can use the existing backend tests and potentially write new ones to cover the full
data flow.

We can work on the following aspects of the frontend:

1. GraphQL Client with Mocked Data - Replace the current REST API with Apollo Client and create mock resolvers that simulate the async job submission/polling pattern described in your architecture.
2. Chemical Structure Editor - For the graphical interface, Ketcher is an excellent choice. It's the gold standard for web-based chemical structure editing, supports SMILES export, and has good
   React integration options.
3. Enhanced UI Components - Your current foundation is solid, but we can add:

- Better results visualization with confidence intervals
- Loading states that show polling progress
- Molecular structure viewers for results
- Batch processing improvements
