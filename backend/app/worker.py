from celery import Celery
from typing import List, Dict, Any
import logging
import pickle
import numpy as np
import os

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

# Load ML models at startup
classifier_model = None
regressor_model = None

def load_models():
    global classifier_model, regressor_model
    try:
        with open(settings.MODEL_CLASSIFIER_PATH, 'rb') as f:
            classifier_model = pickle.load(f)
        with open(settings.MODEL_REGRESSOR_PATH, 'rb') as f:
            regressor_model = pickle.load(f)
        logger.info("ML models loaded successfully")
    except Exception as e:
        logger.error(f"Failed to load models: {e}")
        raise


@celery_app.task(bind=True, name="predict_permeability")
def predict_permeability(self, smiles_list: List[str]) -> Dict[str, Any]:
    """
    Predict permeability for a list of SMILES strings using the two-stage pipeline.
    """
    try:
        # Ensure models are loaded
        if classifier_model is None or regressor_model is None:
            load_models()
        
        results = []
        
        for smiles in smiles_list:
            try:
                # Extract features using processing.py logic
                features = smiles_to_comprehensive_features(smiles)
                feature_vector = combine_features(features)
                
                # Stage 1: Binary classification (near-zero vs non-zero accumulation)
                classifier_pred = classifier_model.predict(feature_vector)[0]
                classifier_prob = classifier_model.predict_proba(feature_vector)[0]
                
                if classifier_pred == 0:  # Near-zero accumulation
                    prediction = 0.0
                    confidence = float(classifier_prob[0])
                else:
                    # Stage 2: Regression for specific permeability level
                    regressor_pred = regressor_model.predict(feature_vector)[0]
                    prediction = float(regressor_pred)
                    confidence = float(classifier_prob[1])
                
                result = {
                    'smiles': smiles,
                    'prediction': prediction,
                    'confidence': confidence,
                    'classifier_prediction': int(classifier_pred),
                    'features': features,
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
                    'error': str(e)
                }
            
            results.append(result)
        
        return {
            'status': 'completed',
            'results': results,
            'total_processed': len(results),
            'successful': len([r for r in results if r['error'] is None]),
            'failed': len([r for r in results if r['error'] is not None])
        }
        
    except Exception as e:
        logger.error(f"Task failed: {e}")
        self.retry(countdown=60, max_retries=3)
        return {
            'status': 'failed',
            'error': str(e),
            'results': []
        }

# Initialize models when worker starts
try:
    load_models()
except Exception as e:
    logger.warning(f"Models not loaded at startup: {e}")