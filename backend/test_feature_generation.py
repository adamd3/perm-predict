import sys
import os
import pandas as pd

# Add the app directory to the Python path to allow absolute imports
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), 'app')))

from app.ml_models.alvadesc_feature_generation import generate_all_features

if __name__ == '__main__':
    test_smiles = "CCO"
    print(f"\nProcessing single SMILES: {test_smiles}")
    features = generate_all_features(test_smiles)
    print("Generated Features:")
    print(features.head())
    print(f"Shape: {features.shape}")

    # You can add more test cases here if needed
    # test_smiles_list = ["CCC", "C1=CC=CN=C1", "O=C(C)Oc1ccccc1C(=O)O"]
    # print(f"\nProcessing list of SMILES: {test_smiles_list}")
    # features_list = generate_all_features(test_smiles_list)
    # print("Generated Features (List):")
    # print(features_list.head())
    # print(f"Shape: {features_list.shape}")

    # SMILES from a file (for testing, create a dummy file)
    # dummy_smiles_file = "dummy_smiles.txt"
    # with open(dummy_smiles_file, "w") as f:
    #     f.write("CCC\n")
    #     f.write("C1=CC=CN=C1\n")
    #     f.write("CCO\n")
    # print(f"\nProcessing SMILES from file: {dummy_smiles_file}")
    # features_file = generate_all_features(dummy_smiles_file)
    # print(features_file.head())
    # print(f"Shape: {features_file.shape}")
    # os.remove(dummy_smiles_file) # Clean up dummy file
