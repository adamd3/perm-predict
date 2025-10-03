import pandas as pd
import numpy as np
import os
import json
import random
import redis
from alvadesccliwrapper.alvadesc import AlvaDesc

try:
    from .input_generation import smiles_to_morgan_fp
except ImportError:
    # When run as a script, adjust path and import directly
    import sys
    import os

    sys.path.append(os.path.join(os.path.dirname(__file__), "."))
    from input_generation import smiles_to_morgan_fp

# --- Configuration ---
ALVADESC_CLI_PATH = "/app/backend/app/ml_models/alvaDescCLI"
ALVADESC_DESCRIPTOR_NAMES_PATH = "/app/alvadesc_descriptor_names.json"
FULL_FEATURE_NAMES_PATH = "/app/full_feature_names.json"
TSV_MOCK_FILE_PATH = "/app/desc_smiles.tsv"  # Path to the TSV file for mock descriptors
SMILES_COLUMN_NAME = "Smiles"  # Name of the SMILES column in the TSV file


# --- Load feature names ---
try:
    with open(ALVADESC_DESCRIPTOR_NAMES_PATH, "r") as f:
        ALVADESC_ALL_AVAILABLE_NAMES = json.load(f)
    with open(FULL_FEATURE_NAMES_PATH, "r") as f:
        FULL_FEATURE_NAMES = json.load(f)

    # Filter ALVADESC_DESCRIPTOR_NAMES to only include those present in FULL_FEATURE_NAMES
    # and that are not Morgan fingerprints
    ALVADESC_DESCRIPTOR_NAMES = [name for name in FULL_FEATURE_NAMES if not name.startswith("fp_")]

    # Determine the exact number of Morgan fingerprints from FULL_FEATURE_NAMES
    MORGAN_FINGERPRINT_COUNT = len([f for f in FULL_FEATURE_NAMES if f.startswith("fp_")])
except FileNotFoundError as e:
    print(f"Error loading feature names JSON: {e}")
    ALVADESC_DESCRIPTOR_NAMES = []
    FULL_FEATURE_NAMES = []
    # Exit or handle this error appropriately in a real application

# --- Redis Cache Configuration ---
# TODO: Move REDIS_HOST, REDIS_PORT, and REDIS_DB to environment variables
REDIS_HOST = os.getenv("REDIS_HOST", "redis") # Use 'redis' for Docker, 'localhost' for local
REDIS_PORT = int(os.getenv("REDIS_PORT", 6379))
REDIS_DB = int(os.getenv("REDIS_DB", 0))

try:
    redis_client = redis.StrictRedis(host=REDIS_HOST, port=REDIS_PORT, db=REDIS_DB, decode_responses=True)
    redis_client.ping()
    print(f"Successfully connected to Redis at {REDIS_HOST}:{REDIS_PORT}")
except redis.exceptions.ConnectionError as e:
    print(f"Could not connect to Redis: {e}. Caching will be disabled.")
    redis_client = None

def get_cached_features(smiles_string: str) -> pd.DataFrame | None:
    if redis_client is None:
        return None
    try:
        cached_data = redis_client.get(smiles_string)
        if cached_data:
            print(f"Cache hit for SMILES: {smiles_string}")
            # Deserialize JSON string back to DataFrame
            df_dict = json.loads(cached_data)
            return pd.DataFrame.from_dict(df_dict)
        print(f"Cache miss for SMILES: {smiles_string}")
        return None
    except Exception as e:
        print(f"Error retrieving from Redis cache for {smiles_string}: {e}")
        return None

def set_cached_features(smiles_string: str, features_df: pd.DataFrame):
    if redis_client is None:
        return
    try:
        # Serialize DataFrame to JSON string
        # Using .to_json(orient='split') for better DataFrame reconstruction
        json_data = features_df.to_json(orient='split')
        redis_client.set(smiles_string, json_data)
        print(f"Cached features for SMILES: {smiles_string}")
    except Exception as e:
        print(f"Error storing to Redis cache for {smiles_string}: {e}")

mock_tsv_data_cache = None


def load_mock_descriptors_from_tsv(file_path):
    """
    Loads descriptors from a TSV file.
    Assumes the first column is 'smiles' and should be used as index.
    Caches the result for subsequent calls.
    """
    global mock_tsv_data_cache
    if mock_tsv_data_cache is not None:
        return mock_tsv_data_cache

    print(f"Loading mock data from TSV: {file_path}...")
    try:
        full_df = pd.read_csv(file_path, sep="\t", index_col=SMILES_COLUMN_NAME)

        # Ensure all expected ALVADESC_DESCRIPTOR_NAMES are present, fill missing with 0.0
        missing_from_tsv = [desc_name for desc_name in ALVADESC_DESCRIPTOR_NAMES if desc_name not in full_df.columns]
        if missing_from_tsv:
            print(f"DEBUG: Columns missing from desc_smiles.tsv: {missing_from_tsv}. Adding with 0.0.")
            for desc_name in missing_from_tsv:
                full_df[desc_name] = 0.0

        # Drop any extra columns that are not in ALVADESC_DESCRIPTOR_NAMES
        extra_cols_in_mock = [
            col for col in full_df.columns if col not in ALVADESC_DESCRIPTOR_NAMES and col != SMILES_COLUMN_NAME
        ]
        if extra_cols_in_mock:
            print(f"Warning: Extra columns in mock data from TSV: {extra_cols_in_mock}. Dropping them.")
            full_df = full_df.drop(columns=extra_cols_in_mock)

        # Select and reorder columns to match ALVADESC_DESCRIPTOR_NAMES
        full_df = full_df[ALVADESC_DESCRIPTOR_NAMES]

        mock_tsv_data_cache = full_df
        print(
            f"Successfully loaded {len(mock_tsv_data_cache)} rows from TSV with {len(ALVADESC_DESCRIPTOR_NAMES)} alvaDesc features."
        )
        return mock_tsv_data_cache
    except FileNotFoundError:
        print(f"Error: Mock TSV file not found at {file_path}")
        mock_tsv_data_cache = pd.DataFrame()
        return mock_tsv_data_cache
    except Exception as e:
        print(f"Error loading mock data from TSV: {e}")
        mock_tsv_data_cache = pd.DataFrame()
        return mock_tsv_data_cache


def generate_all_features(smiles_input, mock_alvadesc=False):
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
            with open(smiles_input, "r") as f:
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

    # --- Caching Logic ---
    cached_features = []
    smiles_to_process = []
    for smiles in smiles_list:
        cached_df = get_cached_features(smiles)
        if cached_df is not None:
            cached_features.append(cached_df)
        else:
            smiles_to_process.append(smiles)

    if not smiles_to_process:
        print("All SMILES strings found in cache.")
        return pd.concat(cached_features) if cached_features else pd.DataFrame()

    # --- Generate Morgan Fingerprints for SMILES to process ---
    morgan_fps = [smiles_to_morgan_fp(s, n_bits=MORGAN_FINGERPRINT_COUNT) for s in smiles_to_process]
    morgan_df = pd.DataFrame(
        morgan_fps, columns=[f"fp_{i}" for i in range(MORGAN_FINGERPRINT_COUNT)], index=smiles_to_process
    )

    # --- Generate alvaDesc Descriptors for SMILES to process ---
    alva_desc_df = pd.DataFrame(index=smiles_to_process)

    if mock_alvadesc:
        print("Mocking alvaDesc descriptors from TSV for uncached SMILES...")
        mock_data = load_mock_descriptors_from_tsv(TSV_MOCK_FILE_PATH)

        if mock_data.empty:
            print("Failed to load mock TSV data. Returning empty DataFrame.")
            return pd.DataFrame()

        alva_desc_rows = []
        for smiles in smiles_to_process:
            if smiles not in mock_data.index:
                print(f"DEBUG: SMILES '{smiles}' not found in mock_data.index. Raising ValueError.")
                raise ValueError("SMILES string must be in the example set.")
            alva_desc_rows.append(mock_data.loc[smiles])

        temp_alva_desc_df = pd.DataFrame(alva_desc_rows, index=smiles_to_process)

        # Ensure all expected ALVADESC_DESCRIPTOR_NAMES are present, fill missing with 0.0
        for desc_name in ALVADESC_DESCRIPTOR_NAMES:
            if desc_name not in temp_alva_desc_df.columns:
                temp_alva_desc_df[desc_name] = 0.0

        # Select and reorder columns to match ALVADESC_DESCRIPTOR_NAMES
        alva_desc_df = temp_alva_desc_df[ALVADESC_DESCRIPTOR_NAMES]
        print(f"Successfully mocked alvaDesc data for {len(smiles_to_process)} uncached SMILES strings from TSV.")

    elif ALVADESC_DESCRIPTOR_NAMES:
        try:
            alva_desc_client = AlvaDesc(ALVADESC_CLI_PATH)
            alva_desc_client.set_input_SMILES(smiles_to_process)

            if not alva_desc_client.calculate_descriptors(ALVADESC_DESCRIPTOR_NAMES):
                error_msg = alva_desc_client.get_error()
                print(f"Error calculating alvaDesc descriptors: {error_msg}")
                # Fill with zeros if alvaDesc fails due to licensing
                alva_desc_data = np.full((len(smiles_to_process), len(ALVADESC_DESCRIPTOR_NAMES)), 0.0)
                alva_desc_df = pd.DataFrame(alva_desc_data, columns=ALVADESC_DESCRIPTOR_NAMES, index=smiles_to_process)
            else:
                alva_desc_output = alva_desc_client.get_output()
                alva_desc_output_names = alva_desc_client.get_output_descriptors()

                if alva_desc_output and alva_desc_output_names:
                    alva_desc_df = pd.DataFrame(alva_desc_output, columns=alva_desc_output_names, index=smiles_to_process)
                    print(f"Head of alvaDesc output DataFrame:\n{alva_desc_df.head()}")
                    # Ensure all expected descriptors are present, fill missing with NaN
                    for desc_name in ALVADESC_DESCRIPTOR_NAMES:
                        if desc_name not in alva_desc_df.columns:
                            alva_desc_df[desc_name] = np.nan
                    alva_desc_df = alva_desc_df[ALVADESC_DESCRIPTOR_NAMES]  # Ensure order
                else:
                    print("alvaDesc returned no output or descriptor names.")
                    alva_desc_data = np.full((len(smiles_to_process), len(ALVADESC_DESCRIPTOR_NAMES)), 0.0)
                    alva_desc_df = pd.DataFrame(alva_desc_data, columns=ALVADESC_DESCRIPTOR_NAMES, index=smiles_to_process)

        except Exception as e:
            print(f"An unexpected error occurred with AlvaDesc: {e}")
            alva_desc_data = np.full((len(smiles_to_process), len(ALVADESC_DESCRIPTOR_NAMES)), 0.0)
            alva_desc_df = pd.DataFrame(alva_desc_data, columns=ALVADESC_DESCRIPTOR_NAMES, index=smiles_to_process)
    else:
        print("alvaDesc descriptor names not loaded. Skipping alvaDesc generation.")
        alva_desc_data = np.full((len(smiles_to_process), 5960), 0.0)  # Fallback if names not loaded
        alva_desc_df = pd.DataFrame(
            alva_desc_data, columns=[f"alvadesc_dummy_{i}" for i in range(5960)], index=smiles_to_process
        )

    # --- Combine Features ---
    combined_df = pd.concat([morgan_df, alva_desc_df], axis=1)

    # --- Reorder columns to match model's expected input ---
    if FULL_FEATURE_NAMES:
        # Add any missing columns to combined_df and fill with 0.0
        missing_cols = [col for col in FULL_FEATURE_NAMES if col not in combined_df.columns]
        if missing_cols:
            print(f"Warning: Missing columns in combined features: {missing_cols}. Filling with 0.0.")
            for col in missing_cols:
                combined_df[col] = 0.0

        # Drop any extra columns that are not in FULL_FEATURE_NAMES
        extra_cols = [col for col in combined_df.columns if col not in FULL_FEATURE_NAMES]
        if extra_cols:
            print(f"Warning: Extra columns in combined features: {extra_cols}. Dropping them.")
            combined_df = combined_df.drop(columns=extra_cols)

        # Reorder and select only the columns expected by the model
        final_features_df = combined_df[FULL_FEATURE_NAMES]
    else:
        print("Warning: Full feature names not loaded. Returning combined_df as is.")
        final_features_df = combined_df

    print(f"Generated combined features for {len(smiles_to_process)} SMILES strings.")
    print(f"Each with {final_features_df.shape[1]} features.")

    # --- Store newly generated features in cache ---
    for smiles in smiles_to_process:
        set_cached_features(smiles, final_features_df.loc[[smiles]]) # Pass a single-row DataFrame

    # --- Combine cached and newly generated features ---
    if cached_features:
        all_features_df = pd.concat(cached_features + [final_features_df])
    else:
        all_features_df = final_features_df

    # Ensure the final DataFrame has the correct order of SMILES strings as per input
    all_features_df = all_features_df.loc[smiles_list]

    return all_features_df


if __name__ == "__main__":
    # IMPORTANT: Replace 'YOUR_SMILES_FROM_TSV' with an actual SMILES string from your desc_smiles.tsv
    # For example, if 'CCC' is in your desc_smiles.tsv, use: test_smiles_in_tsv = "CCC"
    test_smiles_in_tsv = "Oc(cccc1)c1C(N/N=C/c1ccco1)=O"
    test_smiles_not_in_tsv = "CCO"

    # Example usage:
    # Single SMILES string present in TSV
    print(f"\nProcessing single SMILES (in TSV): {test_smiles_in_tsv}")
    features_single = generate_all_features(test_smiles_in_tsv, mock_alvadesc=True)
    print(features_single.head())
    print(f"Shape: {features_single.shape}")

    # Process the same SMILES again to demonstrate cache hit
    print(f"\nProcessing single SMILES (in TSV) again to test cache: {test_smiles_in_tsv}")
    features_single_cached = generate_all_features(test_smiles_in_tsv, mock_alvadesc=True)
    print(features_single_cached.head())
    print(f"Shape: {features_single_cached.shape}")

    # List of SMILES strings (one present, others likely not)
    test_smiles_list = [test_smiles_in_tsv, "C1=CC=CN=C1", test_smiles_not_in_tsv]
    print(f"\nProcessing list of SMILES: {test_smiles_list}")
    features_list = generate_all_features(test_smiles_list, mock_alvadesc=True)
    print(features_list.head())
    print(f"Shape: {features_list.shape}")

    # Test with a SMILES not in the mock TSV (should raise ValueError if mock_alvadesc=True)
    print(f"\nProcessing unknown SMILES (not in TSV, mock_alvadesc=True): {test_smiles_not_in_tsv}")
    try:
        features_unknown = generate_all_features(test_smiles_not_in_tsv, mock_alvadesc=True)
        print(features_unknown)
    except ValueError as e:
        print(f"Expected error caught: {e}")

    # Test with a SMILES not in the mock TSV (mock_alvadesc=False, will attempt alvaDesc)
    # NOTE: This will likely fail if alvaDesc license is not active or CLI not configured.
    print(f"\nProcessing unknown SMILES (not in TSV, mock_alvadesc=False): {test_smiles_not_in_tsv}")
    try:
        features_alvadesc = generate_all_features(test_smiles_not_in_tsv, mock_alvadesc=False)
        print(features_alvadesc.head())
        print(f"Shape: {features_alvadesc.shape}")
    except Exception as e:
        print(f"Error processing with alvaDesc (expected if license not active): {e}")

    # SMILES from a file (for testing, create a dummy file)
    dummy_smiles_file = "dummy_smiles.txt"
    with open(dummy_smiles_file, "w") as f:
        f.write(f"{test_smiles_in_tsv}\n")
        f.write("C1=CC=CN=C1\n")
        f.write("CCO\n")
    print(f"\nProcessing SMILES from file: {dummy_smiles_file}")
    features_file = generate_all_features(dummy_smiles_file, mock_alvadesc=True)
    print(features_file.head())
    print(f"Shape: {features_file.shape}")
    os.remove(dummy_smiles_file)  # Clean up dummy file
