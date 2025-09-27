import numpy as np
import xgboost as xgb
import sys
import os

# Add app directory to path for imports
sys.path.append(os.path.join(os.path.dirname(__file__), "app"))

from app.ml_models.alvadesc_feature_generation import generate_all_features, FULL_FEATURE_NAMES

try:
    print(f"XGBoost version: {xgb.__version__}")

    # Use a valid SMILES string that is present in desc_smiles.tsv
    valid_smiles = "CCOc(cccc1C=NN2CCN(Cc3ccccc3)CC2)c1O"
    print(f"Attempting to generate features for valid SMILES: {valid_smiles}")

    # Ensure MOCK_ALVADESC is True for this test
    # In a real scenario, this would come from settings, but for a standalone test, we pass it directly
    all_features_df = generate_all_features(valid_smiles, mock_alvadesc=True)
    print(f"Features generated successfully. Shape: {all_features_df.shape}")

    if all_features_df.empty:
        raise ValueError("generate_all_features returned an empty DataFrame.")

    feature_vector = all_features_df.values.astype(np.float32)
    print(f"Numpy array created. Shape: {feature_vector.shape}, Dtype: {feature_vector.dtype}")

    if feature_vector.ndim == 1:
        feature_vector = feature_vector.reshape(1, -1)
        print(f"Reshaped 1D feature_vector to {feature_vector.shape}")
    elif feature_vector.ndim > 2:
        feature_vector = feature_vector[0].reshape(1, -1)
        print(f"Reshaped multi-dimensional feature_vector to {feature_vector.shape}")
    print(f"Final feature_vector shape: {feature_vector.shape}")

    if np.isnan(feature_vector).any():
        raise ValueError("Feature vectors contain NaN values before DMatrix creation.")
    if np.isinf(feature_vector).any():
        raise ValueError("Feature vectors contain Inf values before DMatrix creation.")

    print(f"Preparing to create actual DMatrix. Feature vector shape: {feature_vector.shape}, FULL_FEATURE_NAMES length: {len(FULL_FEATURE_NAMES)}")
    
    # This is the line that was crashing in the worker
    dmatrix_feature_vector = xgb.DMatrix(feature_vector, feature_names=FULL_FEATURE_NAMES)
    print("Actual DMatrix created successfully.")
    print(f"DMatrix has {dmatrix_feature_vector.num_row()} rows and {dmatrix_feature_vector.num_col()} columns.")

    # Optional: Try a dummy prediction to see if the model loads and predicts
    # This requires classifier_model to be loaded, which is not in this script's scope.
    # If the DMatrix creation passes, the issue is likely in model loading/prediction.

except ValueError as e:
    print(f"Caught ValueError: {e}", file=sys.stderr)
    sys.exit(1)
except Exception as e:
    print(f"Caught unexpected error: {e}", file=sys.stderr)
    sys.exit(1)