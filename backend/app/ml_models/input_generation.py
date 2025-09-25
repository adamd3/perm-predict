"""
Steps
-----
1. Load raw data from CSV.
2. Clean & filter rows / columns.
3. Create Morgan fingerprints (binary).
4. Collect additional numeric descriptors.
5. Assemble the final feature matrix.

"""

from __future__ import annotations

import pandas as pd
import numpy as np
from collections import Counter
from typing import List

# Chemistry
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem
from rdkit.Chem import rdFingerprintGenerator



CSV_PATH: str = "..."           # Raw dataset
FINGERPRINT_LEN: int = 4200                      # Bits in Morgan FP
FINGERPRINT_RADIUS: int = 2                      # Radius for FP
TRAIN_TEST_SPLIT: float = 0.20                   # Test fraction
RANDOM_STATE: int = 42                           # Reproducibility


def smiles_to_morgan_fp(smiles: str,
                        radius: int = FINGERPRINT_RADIUS,
                        n_bits: int = FINGERPRINT_LEN) -> np.ndarray:
    """Convert a SMILES string to a binary Morgan fingerprint."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return np.zeros(n_bits, dtype=np.int8)
    
    # Use rdFingerprintGenerator to avoid deprecation warning
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=radius, fpSize=n_bits)
    fp = gen.Get='MorganFingerprint'(mol)
    return np.array(fp, dtype=np.int8)

# Data Loading & Pre‑processing

def load_and_preprocess(path: str) -> tuple[pd.DataFrame, List[str], List[str], List[str]]: 
    """Load raw CSV, clean, and generate all feature blocks.

    Returns
    -------
    df : pd.DataFrame
        Full table including features and the target column.
    fp_cols and desc_cols : List[str]
        Column names for fingerprints and additional descriptors.
    """
    df = pd.read_csv(path)

    # Clean rows

    df = df[df["AVG_cells"].notna() & (df["AVG_cells"] != "FLAG")]
    df["AVG_cells"] = pd.to_numeric(df["AVG_cells"], errors="coerce")
    df = df.replace("na", np.nan).dropna(axis=1, how="any")

    df = df[(df["AVG_cells"] < 150) | (df["AVG_cells"] > 250)]

    df["target"] = (df["AVG_cells"] >= 200).astype(int)
    print("Class distribution:", Counter(df["target"]))

    # Morgan fingerprints 

    df["fp"] = df["Smiles"].astype(str).apply(smiles_to_morgan_fp)
    fp_cols = [f"fp_{i}" for i in range(FINGERPRINT_LEN)]
    df = pd.concat([df, pd.DataFrame(df["fp"].tolist(), columns=fp_cols)], axis=1)

    # Additional numeric descriptors 

    excluded = {"Smiles", "AVG_cells", "target", *fp_cols}
    desc_cols = [c for c in df.columns if c not in excluded and pd.api.types.is_numeric_dtype(df[c])]
    print(f"Found {len(desc_cols)} additional descriptors")

    df = df.drop(columns=["fp"])  # raw list columns no longer needed
    return df, fp_cols, desc_cols

# Main script

if __name__ == "__main__":
    df, fp_cols, desc_cols = load_and_preprocess(CSV_PATH)

    feature_cols = fp_cols + desc_cols
    X = df[feature_cols]
    y = df["target"]

    # ...
