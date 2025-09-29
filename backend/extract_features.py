import xgboost as xgb
import json
import os
import sys

# Add app directory to path to resolve imports
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

model_path = "app/ml_models/model_mtb.json"
output_path = "../full_feature_names.json" # Relative to backend directory

try:
    # Load the XGBoost model
    bst = xgb.Booster()
    bst.load_model(model_path)

    # Get feature names
    feature_names = bst.feature_names

    if feature_names:
        # Write feature names to JSON file
        # Ensure the output directory exists if it's not the current one
        output_dir = os.path.dirname(output_path)
        if output_dir and not os.path.exists(output_dir):
            os.makedirs(output_dir)

        with open(output_path, 'w') as f:
            json.dump(feature_names, f, indent=4)
        print(f"Successfully extracted {len(feature_names)} feature names and saved to {output_path}")
    else:
        print("No feature names found in the model.")

except Exception as e:
    print(f"Error: {e}")
    sys.exit(1)
