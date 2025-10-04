import pickle
import numpy as np
import os
import sys
import xgboost as xgb

# Add app directory to path to resolve imports
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from config import settings
from app.ml_models.alvadesc_feature_generation import FULL_FEATURE_NAMES

# Load the classifier model
try:
    with open(settings.MODEL_CLASSIFIER_PATH, 'rb') as f:
        classifier_model = xgb.Booster()
        classifier_model.load_model(settings.MODEL_CLASSIFIER_PATH)
    print(f"Classifier model loaded from {settings.MODEL_CLASSIFIER_PATH}")
    if classifier_model.feature_names:
        print(f"Classifier model feature names (first 10): {classifier_model.feature_names[:10]}")
        print(f"Classifier model feature names (last 10): {classifier_model.feature_names[-10:]}")
        print(f"Classifier model feature count: {len(classifier_model.feature_names)}")
    else:
        print("Classifier model has no feature names.")
except Exception as e:
    print(f"Error loading classifier model: {e}")
    sys.exit(1)

# Create a dummy feature vector
# It should have the same shape and dtype as the actual feature vectors
dummy_feature_vector = np.random.rand(1, len(FULL_FEATURE_NAMES)).astype(np.float32)
print(f"Dummy feature vector shape: {dummy_feature_vector.shape}, dtype: {dummy_feature_vector.dtype}")
dummy_dmatrix = xgb.DMatrix(dummy_feature_vector, feature_names=FULL_FEATURE_NAMES)

try:
    print("Attempting XGBoost prediction...")
    pred = classifier_model.predict(dummy_dmatrix)
    print(f"Prediction: {pred}")
    prob = classifier_model.predict_proba(dummy_dmatrix)
    print(f"Probabilities: {prob}")
    print("XGBoost prediction successful.")
except Exception as e:
    print(f"XGBoost prediction failed with Python exception: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)
except: # Catch all other exceptions, including segfaults if possible
    print("XGBoost prediction failed with an unknown error (possibly a low-level crash).")
    sys.exit(1)

sys.exit(0)