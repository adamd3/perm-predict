from celery import Celery
from celery.signals import worker_process_init
from typing import List, Dict, Any, Union
import logging
import pickle
import joblib
import numpy as np
import os
from datetime import datetime
from pathlib import Path

from app.config import settings
from app.utils.logger import setup_logging
from app.utils.processing import smiles_to_comprehensive_features, combine_features

setup_logging()
logger = logging.getLogger(__name__)

celery_app = Celery(
    "perm_predict_worker",
    broker=settings.CELERY_BROKER_URL,
    backend=settings.CELERY_RESULT_BACKEND,
    include=["app.worker"]
)

celery_app.conf.update(
    task_serializer="json",
    accept_content=["json"],
    result_serializer="json",
    timezone="UTC",
    enable_utc=True,
)

# Global model variables - initialized by worker_process_init
classifier_model = None
regressor_models = None  # Will hold ensemble of models

def _load_model(model_path: Union[str, Path]):
    """Load a model from file, trying different formats."""
    model_path = Path(model_path)
    try:
        # Try joblib first (scikit-learn models)
        return joblib.load(model_path)
    except Exception:
        try:
            # Try pickle as fallback
            with open(model_path, 'rb') as f:
                return pickle.load(f)
        except Exception as e:
            raise RuntimeError(f"Could not load model from {model_path}: {str(e)}")

def _load_ensemble_regressors():
    """
    Load ensemble of regressor models. 
    Currently loads single model but structured for future ensemble expansion.
    """
    # TODO: In future, load multiple models (XGBoost, AttentiveFP, DimeNet++)
    # For now, load the single regressor model
    base_regressor = _load_model(settings.MODEL_REGRESSOR_PATH)
    
    # Structure as list for future ensemble expansion
    return {
        'models': [base_regressor],
        'model_names': ['base_regressor'],
        'weights': [1.0]  # Equal weighting for future ensemble
    }

@worker_process_init.connect
def init_worker(**kwargs):
    """
    Initialize models when worker process starts.
    This is called once per worker process and is more robust than checking globals.
    """
    global classifier_model, regressor_models
    
    try:
        logger.info("Initializing worker process - loading ML models...")
        
        # Load binary classifier
        classifier_model = _load_model(settings.MODEL_CLASSIFIER_PATH)
        logger.info(f"Classifier model loaded: {settings.MODEL_CLASSIFIER_PATH}")
        
        # Load ensemble regressors
        regressor_models = _load_ensemble_regressors()
        logger.info(f"Regressor ensemble loaded with {len(regressor_models['models'])} model(s)")
        
        logger.info("Worker process initialization complete")
        
    except Exception as e:
        logger.error(f"Failed to initialize worker process: {e}")
        raise


def _predict_with_ensemble(feature_vector: np.ndarray) -> Dict[str, float]:
    """
    Predict using ensemble of regressor models and calculate confidence interval.
    Currently uses single model but structured for future ensemble expansion.
    """
    global regressor_models
    
    predictions = []
    
    # Get predictions from all models in ensemble
    for i, model in enumerate(regressor_models['models']):
        weight = regressor_models['weights'][i]
        pred = model.predict(feature_vector)[0]
        predictions.append(pred * weight)
    
    # Calculate ensemble statistics
    ensemble_prediction = np.mean(predictions)
    ensemble_std = np.std(predictions) if len(predictions) > 1 else 0.0
    
    # Calculate confidence based on ensemble variance
    # Lower variance = higher confidence (inverse relationship)
    # TODO: Calibrate this confidence calculation based on validation data
    confidence_from_variance = max(0.1, 1.0 - min(ensemble_std / ensemble_prediction, 0.9)) if ensemble_prediction > 0 else 0.1
    
    return {
        'prediction': float(ensemble_prediction),
        'confidence_from_ensemble': float(confidence_from_variance),
        'ensemble_std': float(ensemble_std),
        'individual_predictions': [float(p) for p in predictions]
    }

def _ensure_feature_consistency(features: Dict[str, Any]) -> np.ndarray:
    """
    Ensure deterministic feature vector ordering for model input.
    Combines molecular descriptors and Morgan fingerprint in consistent order.
    """
    # Get molecular descriptors in sorted key order for consistency
    descriptor_values = [features['descriptors'][key] for key in sorted(features['descriptors'].keys())]
    
    # Combine with Morgan fingerprint
    feature_vector = np.array(descriptor_values + features['morgan_fingerprint'])
    
    # Reshape for model input (models expect 2D array)
    return feature_vector.reshape(1, -1)

@celery_app.task(bind=True, name="predict_permeability")
def predict_permeability(self, smiles_list: List[str], created_at: str = None, job_name: str = None) -> Dict[str, Any]:
    """
    Predict permeability for a list of SMILES strings using the two-stage pipeline.
    """
    try:
        # Models should be loaded by worker_process_init signal
        if classifier_model is None or regressor_models is None:
            raise RuntimeError("Models not initialized. Worker process initialization may have failed.")
        
        results = []
        
        for smiles in smiles_list:
            try:
                # Extract features using processing.py logic
                features = smiles_to_comprehensive_features(smiles)
                feature_vector = _ensure_feature_consistency(features)
                
                # Stage 1: Binary classification (near-zero vs non-zero accumulation)
                classifier_pred = classifier_model.predict(feature_vector)[0]
                classifier_prob = classifier_model.predict_proba(feature_vector)[0]
                
                if classifier_pred == 0:  # Near-zero accumulation
                    prediction = 0.0
                    confidence = float(classifier_prob[0])  # Confidence in "near-zero" prediction
                    ensemble_info = None
                else:
                    # Stage 2: Ensemble regression for specific permeability level
                    ensemble_result = _predict_with_ensemble(feature_vector)
                    prediction = ensemble_result['prediction']
                    
                    # Combine classifier confidence with ensemble confidence
                    classifier_confidence = float(classifier_prob[1])  # Confidence in "non-zero" prediction
                    ensemble_confidence = ensemble_result['confidence_from_ensemble']
                    
                    # Overall confidence is combination of both stages
                    # TODO: Calibrate this combination based on validation data
                    confidence = (classifier_confidence + ensemble_confidence) / 2.0
                    
                    ensemble_info = {
                        'ensemble_std': ensemble_result['ensemble_std'],
                        'individual_predictions': ensemble_result['individual_predictions'],
                        'classifier_confidence': classifier_confidence,
                        'ensemble_confidence': ensemble_confidence
                    }
                
                result = {
                    'smiles': smiles,
                    'prediction': prediction,
                    'confidence': confidence,
                    'classifier_prediction': int(classifier_pred),
                    'features': features,
                    'ensemble_info': ensemble_info,  # Additional debugging info
                    'error': None
                }
                
            except Exception as e:
                logger.error(f"Error processing SMILES {smiles}: {e}")
                result = {
                    'smiles': smiles,
                    'prediction': 0.0,
                    'confidence': 0.0,
                    'classifier_prediction': 0,
                    'features': None,
                    'ensemble_info': None,
                    'error': str(e)
                }
            
            results.append(result)
        
        # Capture completion timestamp
        completed_at = datetime.now().isoformat()
        
        return {
            'status': 'completed',
            'results': results,
            'total_processed': len(results),
            'successful': len([r for r in results if r['error'] is None]),
            'failed': len([r for r in results if r['error'] is not None]),
            'created_at': created_at or datetime.now().isoformat(),
            'completed_at': completed_at,
            'job_name': job_name
        }
        
    except Exception as e:
        logger.error(f"Task failed: {e}")
        self.retry(countdown=60, max_retries=3)
        return {
            'status': 'failed',
            'error': str(e),
            'results': []
        }

# Models are now initialized by @worker_process_init signal
# This ensures proper initialization per worker process