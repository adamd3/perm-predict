from fastapi import HTTPException
import onnxruntime as ort
import logging
import os
from app.config import settings

logger = logging.getLogger(__name__)


async def validate_model():
    """Validate ONNX model at startup."""
    if not os.path.exists(settings.MODEL_PATH):
        raise RuntimeError(f"Model not found at {settings.MODEL_PATH}")

    try:
        session = ort.InferenceSession(settings.MODEL_PATH)
        input_shape = session.get_inputs()[0].shape
        if input_shape[1] != settings.FEATURE_COUNT:
            raise ValueError(f"Model expects {input_shape[1]} features, configured for {settings.FEATURE_COUNT}")
        logger.info(f"Model validated successfully. Input shape: {input_shape}")
        return session
    except Exception as e:
        logger.error(f"Model validation failed: {str(e)}")
        raise
