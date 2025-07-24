import pytest
import numpy as np

from app.utils.processing import smiles_to_features, smiles_to_comprehensive_features, combine_features

# Test data
VALID_SMILES = "CC(=O)OC1=CC=CC=C1C(=O)O"  # Aspirin
INVALID_SMILES = "NOT_A_SMILES"


def test_smiles_to_features_valid():
    """Test SMILES to features conversion with valid input."""
    features = smiles_to_features(VALID_SMILES)
    assert isinstance(features, np.ndarray)
    assert features.dtype == np.int8
    assert features.shape == (4200,)


def test_smiles_to_features_invalid():
    """Test SMILES to features conversion with invalid input."""
    with pytest.raises(ValueError):
        smiles_to_features(INVALID_SMILES)


def test_smiles_to_comprehensive_features():
    """Test comprehensive feature extraction."""
    features = smiles_to_comprehensive_features(VALID_SMILES)
    
    # Check structure
    assert 'morgan_fingerprint' in features
    assert 'descriptors' in features
    
    # Check Morgan fingerprint
    assert len(features['morgan_fingerprint']) == 4200
    assert all(isinstance(x, int) for x in features['morgan_fingerprint'])
    
    # Check descriptors
    expected_descriptors = {
        'MolWt', 'LogP', 'TPSA', 'NumHDonors', 
        'NumHAcceptors', 'NumRotatableBonds', 'NumAromaticRings'
    }
    assert expected_descriptors.issubset(features['descriptors'].keys())
    
    # Check descriptor values are reasonable for aspirin
    assert features['descriptors']['MolWt'] > 100  # Should be around 180
    assert features['descriptors']['NumAromaticRings'] >= 1  # Aspirin has benzene ring


def test_combine_features():
    """Test feature combination into single vector."""
    features = smiles_to_comprehensive_features(VALID_SMILES)
    combined = combine_features(features)
    
    # Should be 2D array with shape (1, n_features)
    assert combined.shape[0] == 1
    assert combined.shape[1] == 4200 + 7  # 4200 Morgan + 7 descriptors
