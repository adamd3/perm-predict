from typing import List
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from app.models import PredictionResponse
from app.config import settings


def smiles_to_features(smiles: str) -> np.ndarray:
    """Convert SMILES to Morgan fingerprint features."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES string: {smiles}")

    return np.array(
        AllChem.GetMorganFingerprintAsBitVect(mol, settings.FINGERPRINT_RADIUS, settings.FEATURE_COUNT), dtype=np.int8
    )


async def process_batch(smiles_batch: List[str], session) -> List[PredictionResponse]:
    """Process a batch of SMILES strings."""
    results = []
    features_list = []
    valid_smiles = []

    for smiles in smiles_batch:
        try:
            # Get features as int8 for memory efficiency
            features = smiles_to_features(smiles)
            features_list.append(features)
            valid_smiles.append(smiles)
        except ValueError as e:
            error_msg = f"Invalid SMILES structure: {str(e)}. Please check for correct syntax and valid atoms."
            results.append(PredictionResponse(smiles=smiles, prediction=0.0, probability=0.0, error=error_msg))

    if features_list:
        # Stack the features and convert to float32 just before prediction
        features_array = np.vstack(features_list).astype(np.float32)
        input_name = session.get_inputs()[0].name
        predictions = session.run(None, {input_name: features_array})[0]

        for smiles, pred in zip(valid_smiles, predictions):
            prob = float(pred[1])
            results.append(PredictionResponse(smiles=smiles, prediction=1 if prob >= 0.5 else 0, probability=prob))

    return results
