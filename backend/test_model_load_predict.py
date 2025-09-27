import xgboost as xgb
import numpy as np
import sys
import os

# Assuming model_mtb.json is the correct model file
MODEL_PATH = "/app/app/ml_models/model_mtb.json"

try:
    print(f"XGBoost version: {xgb.__version__}")
    xgb.set_config(verbosity=3) # Enable verbose XGBoost logging

    # 1. Load the model
    print(f"Attempting to load model from: {MODEL_PATH}")
    classifier_model = xgb.Booster()
    classifier_model.load_model(MODEL_PATH)
    print("Model loaded successfully.")

    # 2. Create a dummy DMatrix
    dummy_feature_vector = np.zeros((1, 9904), dtype=np.float32)
    dummy_feature_names = [f"f{i}" for i in range(9904)] # Use generic feature names
    dummy_dmatrix = xgb.DMatrix(dummy_feature_vector, feature_names=dummy_feature_names)
    print("Dummy DMatrix created successfully.")

    # 3. Attempt prediction
    print("Attempting prediction with dummy DMatrix...")
    pred = classifier_model.predict(dummy_dmatrix, validate_features=False)
    print(f"Prediction successful. Result: {pred}")

    # 4. Attempt predict_proba
    print("Attempting predict_proba with dummy DMatrix...")
    prob = classifier_model.predict_proba(dummy_dmatrix, validate_features=False)
    print(f"Predict_proba successful. Result: {prob}")

except Exception as e:
    print(f"Error in test_model_load_predict.py: {e}", file=sys.stderr)
    sys.exit(1)
