import pandas as pd
import os

TSV_MOCK_FILE_PATH = "/Users/adamdinan/perm-predict/desc_smiles.tsv"
SMILES_COLUMN_NAME = "Smiles"

def check_smiles_in_tsv(smiles_to_check):
    print(f"Attempting to load TSV from: {TSV_MOCK_FILE_PATH}")
    if not os.path.exists(TSV_MOCK_FILE_PATH):
        print(f"Error: File not found at {TSV_MOCK_FILE_PATH}")
        return

    try:
        # Read only the SMILES column to save memory
        df = pd.read_csv(TSV_MOCK_FILE_PATH, sep="\t", usecols=[SMILES_COLUMN_NAME])
        
        # Ensure the SMILES column is set as index for efficient lookup
        df.set_index(SMILES_COLUMN_NAME, inplace=True)
        
        print(f"Successfully loaded {len(df)} SMILES strings from TSV.")
        print(f"First 5 SMILES in index: {df.index.tolist()[:5]}")
        print(f"Last 5 SMILES in index: {df.index.tolist()[-5:]}")

        if smiles_to_check in df.index:
            print(f"SMILES string '{smiles_to_check}' IS present in {TSV_MOCK_FILE_PATH}")
        else:
            print(f"SMILES string '{smiles_to_check}' IS NOT present in {TSV_MOCK_FILE_PATH}")

    except Exception as e:
        print(f"An error occurred while reading the TSV file: {e}")

if __name__ == "__main__":
    check_smiles_in_tsv("CCO")
