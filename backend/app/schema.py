import strawberry
from strawberry.experimental import pydantic
from typing import List, Optional, Dict, Any
from datetime import datetime
import uuid

from app.worker import celery_app
from app.models import (
    MolecularDescriptors as MolecularDescriptorsModel,
    PredictionFeatures as PredictionFeaturesModel,
    PredictionResult as PredictionResultModel,
    JobResult as JobResultModel,
    JobStatus as JobStatusModel,
    PredictionJobInput as PredictionJobInputModel
)

@pydantic.type(model=MolecularDescriptorsModel)
class MolecularDescriptors:
    pass

@pydantic.type(model=PredictionFeaturesModel)
class PredictionFeatures:
    pass

@pydantic.type(model=PredictionResultModel)
class PredictionResult:
    pass

@pydantic.type(model=JobResultModel)
class JobResult:
    pass

@pydantic.type(model=JobStatusModel)
class JobStatus:
    pass

@pydantic.input(model=PredictionJobInputModel)
class PredictionJobInput:
    pass

def convert_features_to_model(features_dict: Dict[str, Any]) -> Optional[PredictionFeaturesModel]:
    """Convert dictionary features to Pydantic model."""
    if not features_dict:
        return None
    
    descriptors = MolecularDescriptorsModel(
        mol_wt=features_dict['descriptors']['MolWt'],
        log_p=features_dict['descriptors']['LogP'],
        tpsa=features_dict['descriptors']['TPSA'],
        num_h_donors=features_dict['descriptors']['NumHDonors'],
        num_h_acceptors=features_dict['descriptors']['NumHAcceptors'],
        num_rotatable_bonds=features_dict['descriptors']['NumRotatableBonds'],
        num_aromatic_rings=features_dict['descriptors']['NumAromaticRings']
    )
    
    return PredictionFeaturesModel(
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
                return JobStatusModel(
                    job_id=job_id,
                    status='pending',
                    created_at=datetime.now().isoformat(),
                    progress="Job is queued and waiting to be processed"
                )
            
            elif task.state == 'PROGRESS':
                return JobStatusModel(
                    job_id=job_id,
                    status='processing',
                    created_at=datetime.now().isoformat(),
                    progress="Job is currently being processed"
                )
            
            elif task.state == 'SUCCESS':
                result = task.result
                
                # Convert results to Pydantic models
                prediction_results = []
                for r in result.get('results', []):
                    features = convert_features_to_model(r.get('features'))
                    prediction_results.append(
                        PredictionResultModel(
                            smiles=r['smiles'],
                            prediction=r['prediction'],
                            confidence=r['confidence'],
                            classifier_prediction=r['classifier_prediction'],
                            features=features,
                            error=r.get('error')
                        )
                    )
                
                return JobResultModel(
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
                return JobStatusModel(
                    job_id=job_id,
                    status='failed',
                    created_at=datetime.now().isoformat(),
                    error=str(task.info) if task.info else "Unknown error occurred"
                )
                
        except Exception as e:
            return JobStatusModel(
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
            
            return JobStatusModel(
                job_id=job_id,
                status=status_map.get(task.state, 'unknown'),
                created_at=datetime.now().isoformat(),
                progress=f"Task state: {task.state}",
                error=str(task.info) if task.state == 'FAILURE' and task.info else None
            )
            
        except Exception as e:
            return JobStatusModel(
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
                return JobStatusModel(
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
            
            return JobStatusModel(
                job_id=task.id,
                status='submitted',
                created_at=datetime.now().isoformat(),
                progress=f"Job submitted with {len(job_input.smiles_list)} compounds"
            )
            
        except Exception as e:
            return JobStatusModel(
                job_id="",
                status='error',
                created_at=datetime.now().isoformat(),
                error=f"Failed to submit job: {str(e)}"
            )

schema = strawberry.Schema(
    query=Query,
    mutation=Mutation
)