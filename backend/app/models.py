from pydantic import BaseModel
from typing import Optional


class SMILESInput(BaseModel):
    smiles: str


class PredictionResponse(BaseModel):
    smiles: str
    prediction: float
    probability: float
    error: Optional[str] = None
