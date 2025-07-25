from pydantic import BaseModel, Field
from typing import Optional, List, Dict, Any
from datetime import datetime


class MolecularDescriptors(BaseModel):
    mol_wt: float = Field(..., description="Molecular weight")
    log_p: float = Field(..., description="LogP (octanol-water partition coefficient)")
    tpsa: float = Field(..., description="Topological polar surface area")
    num_h_donors: int = Field(..., description="Number of hydrogen bond donors")
    num_h_acceptors: int = Field(..., description="Number of hydrogen bond acceptors")
    num_rotatable_bonds: int = Field(..., description="Number of rotatable bonds")
    num_aromatic_rings: int = Field(..., description="Number of aromatic rings")


class PredictionFeatures(BaseModel):
    morgan_fingerprint: List[int] = Field(..., description="Morgan fingerprint bit vector")
    descriptors: MolecularDescriptors = Field(..., description="Molecular descriptors")


class PredictionResult(BaseModel):
    smiles: str = Field(..., description="Input SMILES string")
    prediction: float = Field(..., description="Predicted permeability value")
    confidence: float = Field(..., description="Model confidence score")
    uncertainty: Optional[float] = Field(None, description="Prediction uncertainty from ensemble variance")
    ensemble_std: Optional[float] = Field(None, description="Standard deviation of ensemble predictions")
    classifier_prediction: int = Field(..., description="Binary classifier prediction (0 or 1)")
    ensemble_predictions: Optional[List[float]] = Field(None, description="Individual ensemble model predictions")
    features: Optional[PredictionFeatures] = Field(None, description="Extracted molecular features")
    error: Optional[str] = Field(None, description="Error message if prediction failed")


class JobResult(BaseModel):
    status: str = Field(..., description="Job completion status")
    results: List[PredictionResult] = Field(..., description="List of prediction results")
    total_processed: int = Field(..., description="Total number of compounds processed")
    successful: int = Field(..., description="Number of successful predictions")
    failed: int = Field(..., description="Number of failed predictions")
    job_id: str = Field(..., description="Unique job identifier")
    created_at: str = Field(..., description="Job creation timestamp")
    completed_at: Optional[str] = Field(None, description="Job completion timestamp")


class JobStatus(BaseModel):
    job_id: str = Field(..., description="Unique job identifier")
    status: str = Field(..., description="Current job status")
    created_at: str = Field(..., description="Job creation timestamp")
    progress: Optional[str] = Field(None, description="Current progress information")
    error: Optional[str] = Field(None, description="Error message if job failed")


class PredictionJobInput(BaseModel):
    smiles_list: List[str] = Field(..., min_items=1, description="List of SMILES strings to process")
    job_name: Optional[str] = Field(None, description="Optional job name for identification")


class SMILESInput(BaseModel):
    smiles: str = Field(..., description="SMILES string representation of molecule")


class PredictionResponse(BaseModel):
    smiles: str = Field(..., description="Input SMILES string")
    prediction: float = Field(..., description="Predicted permeability value")
    confidence: float = Field(..., description="Model confidence score")
    error: Optional[str] = Field(None, description="Error message if prediction failed")
