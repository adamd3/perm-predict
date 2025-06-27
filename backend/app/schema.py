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

@pydantic.type(model=MolecularDescriptorsModel, all_fields=True)
class MolecularDescriptors:
    pass

@pydantic.type(model=PredictionFeaturesModel, all_fields=True)
class PredictionFeatures:
    pass

@pydantic.type(model=PredictionResultModel, all_fields=True)
class PredictionResult:
    pass

@pydantic.type(model=JobResultModel, all_fields=True)
class JobResult:
    pass

@pydantic.type(model=JobStatusModel, all_fields=True)
class JobStatus:
    pass

@pydantic.input(model=PredictionJobInputModel, all_fields=True)
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
        """Get the result of a prediction job by job ID. Only returns results for completed jobs."""
        try:
            # Get task result from Celery
            task = celery_app.AsyncResult(job_id)
            
            if not task:
                return None
            
            # Only return JobResult for successfully completed tasks
            if task.state == 'SUCCESS':
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
                    created_at=result.get('created_at', datetime.now().isoformat()),
                    completed_at=result.get('completed_at', datetime.now().isoformat())
                )
            
            # For any other state (PENDING, PROGRESS, FAILURE), return None
            # Clients should use get_job_status to check status
            return None
                
        except Exception as e:
            # Log the error but return None to maintain consistent return type
            logger.error(f"Failed to retrieve job result for {job_id}: {str(e)}")
            return None
    
    @strawberry.field
    def get_job_status(self, job_id: str) -> Optional[JobStatus]:
        """Get the current status of a prediction job."""
        try:
            task = celery_app.AsyncResult(job_id)
            
            if not task:
                return None
            
            # Retrieve job metadata from Redis
            try:
                metadata = celery_app.backend.get(f"job_metadata:{job_id}")
                created_at = metadata.get('created_at') if metadata else datetime.now().isoformat()
            except:
                created_at = datetime.now().isoformat()
            
            status_map = {
                'PENDING': 'pending',
                'PROGRESS': 'processing',
                'SUCCESS': 'completed',
                'FAILURE': 'failed',
                'RETRY': 'retrying',
                'REVOKED': 'cancelled'
            }
            
            # Generate appropriate progress message
            if task.state == 'PENDING':
                progress = "Job is queued and waiting to be processed"
            elif task.state == 'PROGRESS':
                progress = "Job is currently being processed"
            elif task.state == 'SUCCESS':
                progress = "Job completed successfully"
            elif task.state == 'FAILURE':
                progress = "Job failed during processing"
            else:
                progress = f"Task state: {task.state}"
            
            return JobStatusModel(
                job_id=job_id,
                status=status_map.get(task.state, 'unknown'),
                created_at=created_at,
                progress=progress,
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
            
            # Capture creation timestamp
            created_at = datetime.now().isoformat()
            
            # Submit job to Celery with metadata
            task = celery_app.send_task(
                'predict_permeability',
                args=[job_input.smiles_list],
                kwargs={
                    'created_at': created_at,
                    'job_name': job_input.job_name
                }
            )
            
            # Store job metadata in Redis for timestamp tracking
            celery_app.backend.set(
                f"job_metadata:{task.id}",
                {
                    'created_at': created_at,
                    'job_name': job_input.job_name,
                    'smiles_count': len(job_input.smiles_list)
                }
            )
            
            return JobStatusModel(
                job_id=task.id,
                status='submitted',
                created_at=created_at,
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