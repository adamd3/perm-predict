from fastapi import HTTPException
import joblib
import pickle
import logging
import os
from pathlib import Path
from app.config import settings

logger = logging.getLogger(__name__)


async def validate_models():
    """Validate classifier and regressor models at startup."""
    classifier_path = Path(settings.MODEL_CLASSIFIER_PATH)
    regressor_path = Path(settings.MODEL_REGRESSOR_PATH)
    
    # Check if model files exist
    if not classifier_path.exists():
        raise RuntimeError(f"Classifier model not found at {classifier_path}")
    if not regressor_path.exists():
        raise RuntimeError(f"Regressor model not found at {regressor_path}")

    try:
        # Load classifier model
        classifier = _load_model(classifier_path)
        _validate_classifier(classifier)
        logger.info(f"Classifier model validated successfully: {classifier_path}")
        
        # Load regressor model
        regressor = _load_model(regressor_path)
        _validate_regressor(regressor)
        logger.info(f"Regressor model validated successfully: {regressor_path}")
        
        return {"classifier": classifier, "regressor": regressor}
        
    except Exception as e:
        logger.error(f"Model validation failed: {str(e)}")
        raise


def _load_model(model_path: Path):
    """Load a model from file, trying different formats."""
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


def _validate_classifier(model):
    """Validate that the classifier has required methods."""
    if not hasattr(model, 'predict'):
        raise ValueError("Classifier model must have 'predict' method")
    if not hasattr(model, 'predict_proba'):
        raise ValueError("Classifier model must have 'predict_proba' method")
    
    # Test with dummy data to ensure it works
    import numpy as np
    dummy_features = np.zeros((1, settings.FEATURE_COUNT))
    try:
        prediction = model.predict(dummy_features)
        probabilities = model.predict_proba(dummy_features)
        logger.info(f"Classifier test passed - prediction shape: {prediction.shape}, probabilities shape: {probabilities.shape}")
    except Exception as e:
        raise ValueError(f"Classifier model failed test prediction: {str(e)}")


def _validate_regressor(model):
    """Validate that the regressor has required methods."""
    if not hasattr(model, 'predict'):
        raise ValueError("Regressor model must have 'predict' method")
    
    # Test with dummy data to ensure it works
    import numpy as np
    dummy_features = np.zeros((1, settings.FEATURE_COUNT))
    try:
        prediction = model.predict(dummy_features)
        logger.info(f"Regressor test passed - prediction shape: {prediction.shape}")
    except Exception as e:
        raise ValueError(f"Regressor model failed test prediction: {str(e)}")


# Backward compatibility alias
async def validate_model():
    """Backward compatibility wrapper."""
    models = await validate_models()
    return models["classifier"]  # Return classifier for backward compatibility
