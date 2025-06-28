from typing import List, Dict, Any
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
from app.config import settings


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


def validate_feature_dimensions(features: Dict[str, Any]) -> bool:
    """Validate that extracted features match expected dimensions."""
    try:
        # Check Morgan fingerprint dimensions
        if len(features['morgan_fingerprint']) != settings.FEATURE_COUNT:
            return False
        
        # Check descriptor keys
        expected_descriptors = {
            'MolWt', 'LogP', 'TPSA', 'NumHDonors', 
            'NumHAcceptors', 'NumRotatableBonds', 'NumAromaticRings'
        }
        if not expected_descriptors.issubset(features['descriptors'].keys()):
            return False
            
        return True
    except (KeyError, TypeError):
        return False


def normalize_features(feature_vector: np.ndarray) -> np.ndarray:
    """Apply z-score normalization to continuous features, leave binary fingerprints unchanged."""
    # First 4200 features are Morgan fingerprints (binary, leave unchanged)
    fingerprint_features = feature_vector[:, :settings.FEATURE_COUNT]
    
    # Remaining features are continuous descriptors (apply z-score normalization)
    if feature_vector.shape[1] > settings.FEATURE_COUNT:
        descriptor_features = feature_vector[:, settings.FEATURE_COUNT:]
        
        # Simple z-score normalization (mean=0, std=1)
        descriptor_mean = np.mean(descriptor_features, axis=0)
        descriptor_std = np.std(descriptor_features, axis=0)
        descriptor_std = np.where(descriptor_std == 0, 1, descriptor_std)  # Avoid division by zero
        
        normalized_descriptors = (descriptor_features - descriptor_mean) / descriptor_std
        
        # Combine normalized descriptors with unchanged fingerprints
        normalized_features = np.concatenate([fingerprint_features, normalized_descriptors], axis=1)
    else:
        # Only fingerprint features
        normalized_features = fingerprint_features
    
    return normalized_features
