import strawberry
from typing import List, Optional, Dict, Any
from datetime import datetime
import uuid

from app.worker import celery_app

@strawberry.type
class MolecularDescriptors:
    mol_wt: float
    log_p: float
    tpsa: float
    num_h_donors: int
    num_h_acceptors: int
    num_rotatable_bonds: int
    num_aromatic_rings: int

@strawberry.type
class PredictionFeatures:
    morgan_fingerprint: List[int]
    descriptors: MolecularDescriptors

@strawberry.type
class PredictionResult:
    smiles: str
    prediction: float
    confidence: float
    classifier_prediction: int
    features: Optional[PredictionFeatures] = None
    error: Optional[str] = None

@strawberry.type
class JobResult:
    status: str
    results: List[PredictionResult]
    total_processed: int
    successful: int
    failed: int
    job_id: str
    created_at: str
    completed_at: Optional[str] = None

@strawberry.type
class JobStatus:
    job_id: str
    status: str
    created_at: str
    progress: Optional[str] = None
    error: Optional[str] = None

@strawberry.input
class PredictionJobInput:
    smiles_list: List[str]
    job_name: Optional[str] = None

def convert_features_to_graphql(features_dict: Dict[str, Any]) -> Optional[PredictionFeatures]:
    """Convert dictionary features to GraphQL type."""
    if not features_dict:
        return None
    
    descriptors = MolecularDescriptors(
        mol_wt=features_dict['descriptors']['MolWt'],
        log_p=features_dict['descriptors']['LogP'],
        tpsa=features_dict['descriptors']['TPSA'],
        num_h_donors=features_dict['descriptors']['NumHDonors'],
        num_h_acceptors=features_dict['descriptors']['NumHAcceptors'],
        num_rotatable_bonds=features_dict['descriptors']['NumRotatableBonds'],
        num_aromatic_rings=features_dict['descriptors']['NumAromaticRings']
    )
    
    return PredictionFeatures(
        morgan_fingerprint=features_dict['morgan_fingerprint'],
        descriptors=descriptors
    )

@strawberry.type
class Query:
    @strawberry.field
    def get_prediction_result(self, job_id: str) -> Optional[JobResult]:
        """Get the result of a prediction job by job ID."""
        try:
            # Get task result from Celery
            task = celery_app.AsyncResult(job_id)
            
            if not task:
                return None
            
            if task.state == 'PENDING':
                return JobStatus(
                    job_id=job_id,
                    status='pending',
                    created_at=datetime.now().isoformat(),
                    progress="Job is queued and waiting to be processed"
                )
            
            elif task.state == 'PROGRESS':
                return JobStatus(
                    job_id=job_id,
                    status='processing',
                    created_at=datetime.now().isoformat(),
                    progress="Job is currently being processed"
                )
            
            elif task.state == 'SUCCESS':
                result = task.result
                
                # Convert results to GraphQL types
                prediction_results = []
                for r in result.get('results', []):
                    features = convert_features_to_graphql(r.get('features'))
                    prediction_results.append(
                        PredictionResult(
                            smiles=r['smiles'],
                            prediction=r['prediction'],
                            confidence=r['confidence'],
                            classifier_prediction=r['classifier_prediction'],
                            features=features,
                            error=r.get('error')
                        )
                    )
                
                return JobResult(
                    status='completed',
                    results=prediction_results,
                    total_processed=result.get('total_processed', 0),
                    successful=result.get('successful', 0),
                    failed=result.get('failed', 0),
                    job_id=job_id,
                    created_at=datetime.now().isoformat(),
                    completed_at=datetime.now().isoformat()
                )
            
            else:  # FAILURE or other error states
                return JobStatus(
                    job_id=job_id,
                    status='failed',
                    created_at=datetime.now().isoformat(),
                    error=str(task.info) if task.info else "Unknown error occurred"
                )
                
        except Exception as e:
            return JobStatus(
                job_id=job_id,
                status='error',
                created_at=datetime.now().isoformat(),
                error=f"Failed to retrieve job result: {str(e)}"
            )
    
    @strawberry.field
    def get_job_status(self, job_id: str) -> Optional[JobStatus]:
        """Get the current status of a prediction job."""
        try:
            task = celery_app.AsyncResult(job_id)
            
            if not task:
                return None
            
            status_map = {
                'PENDING': 'pending',
                'PROGRESS': 'processing',
                'SUCCESS': 'completed',
                'FAILURE': 'failed',
                'RETRY': 'retrying',
                'REVOKED': 'cancelled'
            }
            
            return JobStatus(
                job_id=job_id,
                status=status_map.get(task.state, 'unknown'),
                created_at=datetime.now().isoformat(),
                progress=f"Task state: {task.state}",
                error=str(task.info) if task.state == 'FAILURE' and task.info else None
            )
            
        except Exception as e:
            return JobStatus(
                job_id=job_id,
                status='error',
                created_at=datetime.now().isoformat(),
                error=f"Failed to get job status: {str(e)}"
            )

@strawberry.type
class Mutation:
    @strawberry.field
    def submit_prediction_job(self, job_input: PredictionJobInput) -> JobStatus:
        """Submit a new prediction job and return the job ID."""
        try:
            # Validate input
            if not job_input.smiles_list:
                return JobStatus(
                    job_id="",
                    status='error',
                    created_at=datetime.now().isoformat(),
                    error="SMILES list cannot be empty"
                )
            
            # Submit job to Celery
            task = celery_app.send_task(
                'predict_permeability',
                args=[job_input.smiles_list]
            )
            
            return JobStatus(
                job_id=task.id,
                status='submitted',
                created_at=datetime.now().isoformat(),
                progress=f"Job submitted with {len(job_input.smiles_list)} compounds"
            )
            
        except Exception as e:
            return JobStatus(
                job_id="",
                status='error',
                created_at=datetime.now().isoformat(),
                error=f"Failed to submit job: {str(e)}"
            )

schema = strawberry.Schema(
    query=Query,
    mutation=Mutation
)