from fastapi import FastAPI, HTTPException, UploadFile, Request
from fastapi.responses import JSONResponse
import pandas as pd
from typing import Optional
from slowapi import Limiter
from slowapi.util import get_remote_address
from slowapi.errors import RateLimitExceeded
import datetime
import onnxruntime as ort

from app.utils.validation import validate_model
from app.utils.processing import process_batch
from app.utils.logger import setup_logging
from app.config import settings
from app.models import SMILESInput, PredictionResponse

limiter = Limiter(key_func=get_remote_address)
app = FastAPI()
app.state.limiter = limiter

# Load ONNX model
model_path = settings.MODEL_PATH
session = ort.InferenceSession(model_path)


@app.exception_handler(RateLimitExceeded)
async def rate_limit_handler(request, exc):
    return JSONResponse({"error": "Rate limit exceeded"}, status_code=429)


@app.on_event("startup")
async def startup_event():
    global session
    session = await validate_model()


@app.post("/predict/single", response_model=PredictionResponse)
@limiter.limit("100/minute")
async def predict_single(request: Request, input_data: SMILESInput):
    results = await process_batch([input_data.smiles], session)
    return results[0]


@app.post("/predict/batch")
@limiter.limit("20/minute")
async def predict_batch(request: Request, file: UploadFile):
    if not file.filename.endswith(".csv"):
        raise HTTPException(status_code=400, detail="Only CSV files are supported")
    try:
        df = pd.read_csv(file.file)
        if "smiles" not in df.columns:
            raise HTTPException(status_code=400, detail="CSV must contain 'smiles' column")

        results = []
        for i in range(0, len(df), settings.BATCH_SIZE):
            batch = df["smiles"][i : i + settings.BATCH_SIZE].tolist()
            batch_results = await process_batch(batch, session)
            results.extend(batch_results)
        return results
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))


@app.get("/health")
async def health_check():
    """Health check endpoint."""
    return {"status": "healthy", "model_loaded": bool(session), "timestamp": datetime.datetime.now().isoformat()}
