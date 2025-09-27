import sys
from app.ml_models.alvadesc_feature_generation import generate_all_features

try:
    # Use a SMILES string that is highly unlikely to be in desc_smiles.tsv
    test_smiles = "C"
    print(f"Attempting to generate features for SMILES: {test_smiles}")
    features = generate_all_features(test_smiles, mock_alvadesc=True)
    print("Features generated successfully (unexpected for 'C').")
    print(features)
except ValueError as e:
    print(f"Caught expected ValueError: {e}", file=sys.stderr)
    sys.exit(1)
except Exception as e:
    print(f"Caught unexpected error: {e}", file=sys.stderr)
    sys.exit(1)
