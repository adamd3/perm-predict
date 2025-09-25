from typing import List, Dict, Any, Optional
import logging
import pickle
import numpy as np
import os
import torch
import traceback  # Import traceback
import xgboost as xgb

from app.config import settings
from app.utils.logger import logger
from app.utils.processing import smiles_to_comprehensive_features, combine_features
from app.ml_models.alvadesc_feature_generation import generate_all_features, FULL_FEATURE_NAMES
from app.models import PredictionFeatures, MolecularDescriptors  # Import Pydantic models
from app.celery_instance import celery_app  # Import celery_app from the new instance file


# Load ML models at startup
classifier_model = None


def load_models():
    """Load the classification model (simplified for classification-only mode)."""
    global classifier_model

    try:
        classifier_path = settings.MODEL_CLASSIFIER_PATH
        logger.info(f"Attempting to load classifier model from: {classifier_path}")
        if os.path.exists(classifier_path):
            with open(classifier_path, "rb") as f:
                classifier_model = xgb.Booster()
                classifier_model.load_model(classifier_path)
            logger.info("Classifier model loaded successfully")
            if classifier_model.feature_names:
                logger.info(f"Classifier model feature names (first 10): {classifier_model.feature_names[:10]}")
                logger.info(f"Classifier model feature names (last 10): {classifier_model.feature_names[-10:]}")
                logger.info(f"Classifier model feature count: {len(classifier_model.feature_names)}")
            else:
                logger.info("Classifier model has no feature names.")

        else:
            logger.error("Classifier model could not be loaded - check model file path")

    except Exception as e:
        raise


load_models()


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
    return {"confidence": max_proba, "uncertainty": 1.0 - max_proba, "class_probabilities": classifier_proba.tolist()}


@celery_app.task(bind=True, name="predict_permeability")
def predict_permeability(
    self, smiles_list: List[str], created_at: Optional[str] = None, job_name: Optional[str] = None
) -> Dict[str, Any]:
    """
    Predict permeability for a list of SMILES strings using classification model only.
    """
    logger.info(f"predict_permeability received smiles_list: type={type(smiles_list)}, content={smiles_list}")
    if not isinstance(smiles_list, list):
        logger.error(f"smiles_list is not a list! Type: {type(smiles_list)}")
        raise TypeError("smiles_list must be a list of strings")
    if not all(isinstance(s, str) for s in smiles_list):
        logger.error(f"smiles_list contains non-string elements: {smiles_list}")
        raise TypeError("All elements in smiles_list must be strings")

    # --- Dummy Result for Testing ---
    # dummy_results = []
    # for smiles in smiles_list:
    #     dummy_results.append({
    #         'smiles': smiles,
    #         'prediction': 1, # Dummy permeant
    #         'confidence': 0.99,
    #         'uncertainty': 0.01,
    #         'class_probabilities': [0.01, 0.99],
    #         'classifier_prediction': 1,
    #         'features': None,
    #         'error': None
    #     })

    # return {
    #     'status': 'completed',
    #     'results': dummy_results,
    #     'total_processed': len(dummy_results),
    #     'successful': len(dummy_results),
    #     'failed': 0
    # }
    # --- End Dummy Result ---

    try:
        # Ensure models are loaded
        if classifier_model is None:
            load_models()

        results = []

        for smiles in smiles_list:
            try:
                # Extract features using alvaDesc CLI wrapper and Morgan fingerprints
                logger.info(f"Generating features for SMILES: {smiles}")
                all_features_df = generate_all_features(smiles)
                logger.info(f"Features generated for SMILES: {smiles}")

                # Convert DataFrame to numpy array for model input
                feature_vector = all_features_df.values.astype(np.float32)

                # Ensure feature_vector has the correct shape (1, 10160) for a single sample
                if feature_vector.ndim == 1:
                    feature_vector = feature_vector.reshape(1, -1)
                elif feature_vector.ndim > 2:
                    # For now, we\'ll take the first row if multiple are returned for a single SMILES
                    feature_vector = feature_vector[0].reshape(1, -1)

                # Convert numpy array to DMatrix for XGBoost Booster model
                dmatrix_feature_vector = xgb.DMatrix(feature_vector, feature_names=FULL_FEATURE_NAMES)

                # Store the generated features for the result
                processed_features_for_result = {
                    "feature_vector": feature_vector.tolist()  # Store as list for JSON serialization
                }

                logger.info(f"Feature vector shape: {feature_vector.shape}")
                logger.info("Attempting classification prediction...")
                # Classification prediction
                if classifier_model is not None:
                    logger.info(f"Classifier model type: {type(classifier_model)}")
                    logger.info(f"Feature vector dtype: {feature_vector.dtype}, shape: {feature_vector.shape}")
                    logger.info(f"Feature vector min: {np.min(feature_vector)}")
                    logger.info(f"Feature vector max: {np.max(feature_vector)}")
                    logger.info(f"Feature vector mean: {np.mean(feature_vector)}")
                    logger.info(f"Feature vector std: {np.std(feature_vector)}")
                    logger.info(f"Feature vector sample (first 10): {feature_vector.flatten()[:10]}")
                    np.save("/app/feature_vector_debug.npy", feature_vector)
                    logger.info(f"Feature vector saved to /app/feature_vector_debug.npy")
                    logger.info(f"FULL_FEATURE_NAMES length: {len(FULL_FEATURE_NAMES)}")
                    logger.info(f"FULL_FEATURE_NAMES sample (first 10): {FULL_FEATURE_NAMES[:10]}")
                    logger.info(f"FULL_FEATURE_NAMES sample (last 10): {FULL_FEATURE_NAMES[-10:]}")
                    if np.isnan(feature_vector).any():
                        logger.error("Feature vectors contain NaN values before prediction.")
                        raise ValueError("Feature vectors contain NaN values.")
                    if np.isinf(feature_vector).any():
                        logger.error("Feature vectors contain Inf values before prediction.")
                        raise ValueError("Feature vectors contain Inf values.")

                    try:
                        classifier_pred = classifier_model.predict(dmatrix_feature_vector)[0]
                        logger.info(f"Classifier prediction: {classifier_pred}")
                        classifier_prob = classifier_model.predict_proba(dmatrix_feature_vector)[0]
                        logger.info(f"Classifier probabilities: {classifier_prob}")
                        confidence_stats = calculate_classification_confidence(classifier_prob)
                    except Exception as xgb_e:
                        logger.exception(f"XGBoost prediction failed: {xgb_e}")
                        raise  # Re-raise to be caught by the outer try-except
                    logger.info(f"Confidence stats: {confidence_stats}")
                else:
                    # Fallback if no classifier
                    classifier_pred = 0
                    confidence_stats = {"confidence": 0.0, "uncertainty": 1.0, "class_probabilities": [0.5, 0.5]}
                    logger.warning("No classifier model available - defaulting to non-permeant")

                logger.info("Classification prediction completed.")

                # Convert classification to binary prediction
                prediction = 1 if classifier_pred == 1 else 0  # 1 = permeant, 0 = non-permeant

                result = {
                    "smiles": smiles,
                    "prediction": prediction,
                    "confidence": confidence_stats["confidence"],
                    "uncertainty": confidence_stats["uncertainty"],
                    "class_probabilities": confidence_stats["class_probabilities"],
                    "classifier_prediction": int(classifier_pred),
                    "features": processed_features_for_result,  # Use the updated features
                    "error": None,
                }

            except Exception as e:
                logger.exception(f"Error processing SMILES {smiles}: {e}")
                result = {
                    "smiles": smiles,
                    "prediction": 0,
                    "confidence": 0.0,
                    "uncertainty": 1.0,
                    "class_probabilities": [0.5, 0.5],
                    "classifier_prediction": 0,
                    "features": processed_features_for_result,
                    "error": str(e),
                }
            results.append(result)

        return {
            "status": "completed",
            "results": results,
            "total_processed": len(results),
            "successful": len([r for r in results if r["error"] is None]),
            "failed": len([r for r in results if r["error"] is not None]),
        }

    except Exception as e:
        logger.error(f"Task failed: {e}")
        traceback.print_exc()  # Print full traceback
        self.retry(countdown=60, max_retries=3)
        return {"status": "failed", "error": str(e), "results": []}
