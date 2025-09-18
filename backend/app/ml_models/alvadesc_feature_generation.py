import pandas as pd
import numpy as np
import os
import json
from alvadesccliwrapper.alvadesc import AlvaDesc
try:
    from .input_generation import smiles_to_morgan_fp
except ImportError:
    # When run as a script, adjust path and import directly
    import sys
    import os
    sys.path.append(os.path.join(os.path.dirname(__file__), '.'))
    from input_generation import smiles_to_morgan_fp

# --- Configuration ---
ALVADESC_CLI_PATH = "/app/backend/app/ml_models/alvaDescCLI"
ALVADESC_DESCRIPTOR_NAMES_PATH = "/app/alvadesc_descriptor_names.json"
FULL_FEATURE_NAMES_PATH = "/app/full_feature_names.json"

# --- Load feature names ---
try:
    with open(ALVADESC_DESCRIPTOR_NAMES_PATH, 'r') as f:
        ALVADESC_DESCRIPTOR_NAMES = json.load(f)
    with open(FULL_FEATURE_NAMES_PATH, 'r') as f:
        FULL_FEATURE_NAMES = json.load(f)
except FileNotFoundError as e:
    print(f"Error loading feature names JSON: {e}")
    ALVADESC_DESCRIPTOR_NAMES = []
    FULL_FEATURE_NAMES = []
    # Exit or handle this error appropriately in a real application

def generate_all_features(smiles_input):
    """
    Generates alvaDesc chemical descriptors and Morgan fingerprints for a given
    SMILES string(s) and combines them into a single DataFrame with the
    expected feature order.

    Args:
        smiles_input (str or list[str]): A single SMILES string or a list of SMILES strings.

    Returns:
        pd.DataFrame: A DataFrame where each row corresponds to a SMILES string
                      and columns are the combined and ordered features.
                      Returns an empty DataFrame if no SMILES are processed or on error.
    """
    smiles_list = []

    if isinstance(smiles_input, str):
        if os.path.isfile(smiles_input):
            with open(smiles_input, 'r') as f:
                smiles_list = [line.strip() for line in f if line.strip()]
        else:
            smiles_list = [smiles_input]
    elif isinstance(smiles_input, list):
        smiles_list = smiles_input
    else:
        print("Invalid smiles_input type. Must be str (SMILES or file path) or list of str.")
        return pd.DataFrame()

    if not smiles_list:
        print("No SMILES strings provided for feature generation.")
        return pd.DataFrame()

    # --- Generate Morgan Fingerprints ---
    morgan_fps = [smiles_to_morgan_fp(s) for s in smiles_list]
    morgan_df = pd.DataFrame(morgan_fps, columns=[f'fp_{i}' for i in range(len(morgan_fps[0]))], index=smiles_list)

    # --- Generate alvaDesc Descriptors ---
    alva_desc_df = pd.DataFrame(index=smiles_list)
    if ALVADESC_DESCRIPTOR_NAMES:
        try:
            alva_desc_client = AlvaDesc(ALVADESC_CLI_PATH)
            alva_desc_client.set_input_SMILES(smiles_list)
            
            if not alva_desc_client.calculate_descriptors(ALVADESC_DESCRIPTOR_NAMES):
                error_msg = alva_desc_client.get_error()
                print(f"Error calculating alvaDesc descriptors: {error_msg}")
                # Fill with NaNs or zeros if alvaDesc fails
                alva_desc_data = np.full((len(smiles_list), len(ALVADESC_DESCRIPTOR_NAMES)), np.nan)
                alva_desc_df = pd.DataFrame(alva_desc_data, columns=ALVADESC_DESCRIPTOR_NAMES, index=smiles_list)
            else:
                alva_desc_output = alva_desc_client.get_output()
                alva_desc_output_names = alva_desc_client.get_output_descriptors()
                
                if alva_desc_output and alva_desc_output_names:
                    alva_desc_df = pd.DataFrame(alva_desc_output, columns=alva_desc_output_names, index=smiles_list)
                    # Ensure all expected descriptors are present, fill missing with NaN
                    for desc_name in ALVADESC_DESCRIPTOR_NAMES:
                        if desc_name not in alva_desc_df.columns:
                            alva_desc_df[desc_name] = np.nan
                    alva_desc_df = alva_desc_df[ALVADESC_DESCRIPTOR_NAMES] # Ensure order
                else:
                    print("alvaDesc returned no output or descriptor names.")
                    alva_desc_data = np.full((len(smiles_list), len(ALVADESC_DESCRIPTOR_NAMES)), np.nan)
                    alva_desc_df = pd.DataFrame(alva_desc_data, columns=ALVADESC_DESCRIPTOR_NAMES, index=smiles_list)

        except Exception as e:
            print(f"An unexpected error occurred with AlvaDesc: {e}")
            alva_desc_data = np.full((len(smiles_list), len(ALVADESC_DESCRIPTOR_NAMES)), np.nan)
            alva_desc_df = pd.DataFrame(alva_desc_data, columns=ALVADESC_DESCRIPTOR_NAMES, index=smiles_list)
    else:
        print("alvaDesc descriptor names not loaded. Skipping alvaDesc generation.")
        alva_desc_data = np.full((len(smiles_list), 5960), np.nan) # Fallback if names not loaded
        alva_desc_df = pd.DataFrame(alva_desc_data, columns=[f'alvadesc_dummy_{i}' for i in range(5960)], index=smiles_list)


    # --- Combine Features ---
    combined_df = pd.concat([morgan_df, alva_desc_df], axis=1)

    # --- Reorder columns to match model's expected input ---
    if FULL_FEATURE_NAMES and len(FULL_FEATURE_NAMES) == combined_df.shape[1]:
        # Ensure all columns in FULL_FEATURE_NAMES are in combined_df
        missing_cols = [col for col in FULL_FEATURE_NAMES if col not in combined_df.columns]
        if missing_cols:
            print(f"Warning: Missing columns in combined features: {missing_cols}. Filling with NaN.")
            for col in missing_cols:
                combined_df[col] = np.nan
        
        # Reorder and select only the columns expected by the model
        final_features_df = combined_df[FULL_FEATURE_NAMES]
    else:
        print("Warning: Full feature names not loaded or mismatch in feature count. Returning combined_df as is.")
        final_features_df = combined_df

    print(f"Generated combined features for {len(smiles_list)} SMILES strings.")
    print(f"Each with {final_features_df.shape[1]} features.")

    return final_features_df

if __name__ == '__main__':
    # Example usage:
    # Single SMILES string
    test_smiles = "CCO"
    print(f"\nProcessing single SMILES: {test_smiles}")
    features_single = generate_all_features(test_smiles)
    print(features_single.head())
    print(f"Shape: {features_single.shape}")

    # List of SMILES strings
    test_smiles_list = ["CCC", "C1=CC=CN=C1", "O=C(C)Oc1ccccc1C(=O)O"]
    print(f"\nProcessing list of SMILES: {test_smiles_list}")
    features_list = generate_all_features(test_smiles_list)
    print(features_list.head())
    print(f"Shape: {features_list.shape}")

    # SMILES from a file (for testing, create a dummy file)
    dummy_smiles_file = "dummy_smiles.txt"
    with open(dummy_smiles_file, "w") as f:
        f.write("CCC\n")
        f.write("C1=CC=CN=C1\n")
        f.write("CCO\n")
    print(f"\nProcessing SMILES from file: {dummy_smiles_file}")
    features_file = generate_all_features(dummy_smiles_file)
    print(features_file.head())
    print(f"Shape: {features_file.shape}")
    os.remove(dummy_smiles_file) # Clean up dummy file

