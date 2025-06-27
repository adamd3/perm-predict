from typing import List, Dict, Any
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
from app.config import settings
from app.models import PredictionResponse


def smiles_to_features(smiles: str) -> np.ndarray:
    """Convert SMILES to Morgan fingerprint features."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES string: {smiles}")

    return np.array(
        AllChem.GetMorganFingerprintAsBitVect(mol, settings.FINGERPRINT_RADIUS, settings.FEATURE_COUNT), dtype=np.int8
    )


def smiles_to_comprehensive_features(smiles: str) -> Dict[str, Any]:
    """Convert SMILES to comprehensive molecular features including fingerprints and descriptors."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES string: {smiles}")
    
    # Morgan fingerprints
    morgan_fp = np.array(
        AllChem.GetMorganFingerprintAsBitVect(mol, settings.FINGERPRINT_RADIUS, settings.FEATURE_COUNT)
    )
    
    # Molecular descriptors
    descriptors = {
        'MolWt': Descriptors.MolWt(mol),
        'LogP': Descriptors.MolLogP(mol),
        'TPSA': Descriptors.TPSA(mol),
        'NumHDonors': Descriptors.NumHDonors(mol),
        'NumHAcceptors': Descriptors.NumHAcceptors(mol),
        'NumRotatableBonds': Descriptors.NumRotatableBonds(mol),
        'NumAromaticRings': Descriptors.NumAromaticRings(mol),
    }
    
    return {
        'morgan_fingerprint': morgan_fp.tolist(),
        'descriptors': descriptors
    }


def combine_features(features: Dict[str, Any]) -> np.ndarray:
    """Combine different feature types into a single feature vector."""
    morgan_fp = np.array(features['morgan_fingerprint'])
    descriptor_values = np.array(list(features['descriptors'].values()))
    
    # Combine features
    combined = np.concatenate([morgan_fp, descriptor_values])
    return combined.reshape(1, -1)


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
