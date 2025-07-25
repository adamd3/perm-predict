from celery import Celery
from typing import List, Dict, Any, Optional
import logging
import pickle
import numpy as np
import os
import torch

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

def load_models():
    """Load the classification model (simplified for classification-only mode)."""
    global classifier_model
    
    try:
        # Load binary classifier
        classifier_path = settings.MODEL_CLASSIFIER_PATH
        if os.path.exists(classifier_path):
            with open(classifier_path, 'rb') as f:
                classifier_model = pickle.load(f)
            logger.info("Classifier model loaded successfully")
        else:
            logger.warning(f"Classifier model not found at {classifier_path}")
            
        if classifier_model is None:
            logger.error("Classifier model could not be loaded - check model file path")
            
    except Exception as e:
        logger.error(f"Failed to load classifier model: {e}")
        raise


# COMMENTED OUT - Regression/Ensemble functionality (for future use)
# ensemble_regressors = {}
# blender_model = None

# def load_ensemble_models():
#     """Load all models for the two-step ensemble prediction pipeline."""
#     global ensemble_regressors, blender_model
#     
#     try:
#         # Load ensemble regressors
#         regressor_paths = {
#             'xgboost': os.path.join(os.path.dirname(settings.MODEL_CLASSIFIER_PATH), 'xgboost_regressor.pkl'),
#             'attentivefp': os.path.join(os.path.dirname(settings.MODEL_CLASSIFIER_PATH), 'attentivefp_regressor.pt'),
#             'dimenet': os.path.join(os.path.dirname(settings.MODEL_CLASSIFIER_PATH), 'dimenet_regressor.pt')
#         }
#         
#         for name, path in regressor_paths.items():
#             if os.path.exists(path):
#                 if path.endswith('.pkl'):
#                     with open(path, 'rb') as f:
#                         ensemble_regressors[name] = pickle.load(f)
#                 elif path.endswith('.pt'):
#                     # For PyTorch models, we'll need the model architecture loaded separately
#                     # For now, just log that we found the file
#                     logger.info(f"Found {name} model at {path} (PyTorch loading not implemented yet)")
#                 logger.info(f"Regressor {name} loaded successfully")
#             else:
#                 logger.warning(f"Regressor {name} not found at {path}")
#         
#         # Load blender model
#         blender_path = os.path.join(os.path.dirname(settings.MODEL_CLASSIFIER_PATH), 'blender_model.pkl')
#         if os.path.exists(blender_path):
#             with open(blender_path, 'rb') as f:
#                 blender_model = pickle.load(f)
#             logger.info("Blender model loaded successfully")
#         else:
#             logger.warning(f"Blender model not found at {blender_path}")
#             
#     except Exception as e:
#         logger.error(f"Failed to load ensemble models: {e}")
#         raise


# COMMENTED OUT - Ensemble prediction functions (for future use)
# def get_ensemble_predictions(feature_vector: np.ndarray) -> List[float]:
#     """Get predictions from all available ensemble regressors."""
#     predictions = []
#     
#     # XGBoost regressor
#     if 'xgboost' in ensemble_regressors:
#         try:
#             xgb_pred = ensemble_regressors['xgboost'].predict(feature_vector)[0]
#             predictions.append(float(xgb_pred))
#         except Exception as e:
#             logger.warning(f"XGBoost regressor failed: {e}")
#     
#     # PyTorch models (AttentiveFP, DimeNet++) - placeholder for now
#     # TODO: Implement when PyTorch model architectures are available
#     for model_name in ['attentivefp', 'dimenet']:
#         if model_name in ensemble_regressors:
#             logger.warning(f"{model_name} prediction not yet implemented")
#     
#     return predictions


# def calculate_confidence_interval(predictions: List[float], classifier_confidence: float) -> Dict[str, float]:
#     """Calculate calibrated confidence interval from ensemble variance."""
#     if len(predictions) == 0:
#         return {'confidence': 0.0, 'uncertainty': 1.0, 'ensemble_std': 0.0}
#     
#     if len(predictions) == 1:
#         # Single model - use classifier confidence
#         return {
#             'confidence': classifier_confidence,
#             'uncertainty': 1.0 - classifier_confidence,
#             'ensemble_std': 0.0
#         }
#     
#     # Multiple models - calculate ensemble statistics
#     ensemble_mean = np.mean(predictions)
#     ensemble_std = np.std(predictions)
#     
#     # Combine classifier confidence with ensemble uncertainty
#     # Higher std = lower confidence
#     ensemble_confidence = classifier_confidence * np.exp(-ensemble_std)
#     
#     return {
#         'confidence': float(ensemble_confidence),
#         'uncertainty': float(ensemble_std),
#         'ensemble_std': float(ensemble_std)
#     }


def calculate_classification_confidence(classifier_proba: np.ndarray) -> Dict[str, float]:
    """Calculate confidence metrics for classification predictions."""
    max_proba = float(np.max(classifier_proba))
    return {
        'confidence': max_proba,
        'uncertainty': 1.0 - max_proba,
        'class_probabilities': classifier_proba.tolist()
    }


@celery_app.task(bind=True, name="predict_permeability")
def predict_permeability(self, smiles_list: List[str]) -> Dict[str, Any]:
    """
    Predict permeability for a list of SMILES strings using classification model only.
    """
    try:
        # Ensure models are loaded
        if classifier_model is None:
            load_models()
        
        results = []
        
        for smiles in smiles_list:
            try:
                # Extract features using processing.py logic
                features = smiles_to_comprehensive_features(smiles)
                feature_vector = combine_features(features)
                
                # Classification prediction
                if classifier_model is not None:
                    classifier_pred = classifier_model.predict(feature_vector)[0]
                    classifier_prob = classifier_model.predict_proba(feature_vector)[0]
                    confidence_stats = calculate_classification_confidence(classifier_prob)
                else:
                    # Fallback if no classifier
                    classifier_pred = 0
                    confidence_stats = {'confidence': 0.0, 'uncertainty': 1.0, 'class_probabilities': [0.5, 0.5]}
                    logger.warning("No classifier model available - defaulting to non-permeant")
                
                # Convert classification to binary prediction
                prediction = 1 if classifier_pred == 1 else 0  # 1 = permeant, 0 = non-permeant
                
                result = {
                    'smiles': smiles,
                    'prediction': prediction,
                    'confidence': confidence_stats['confidence'],
                    'uncertainty': confidence_stats['uncertainty'],
                    'class_probabilities': confidence_stats['class_probabilities'],
                    'classifier_prediction': int(classifier_pred),
                    'features': features,
                    'error': None
                }
                
            except Exception as e:
                logger.error(f"Error processing SMILES {smiles}: {e}")
                result = {
                    'smiles': smiles,
                    'prediction': 0,
                    'confidence': 0.0,
                    'uncertainty': 1.0,
                    'class_probabilities': [0.5, 0.5],
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