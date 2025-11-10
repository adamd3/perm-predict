from typing import List, Dict, Any, Optional
import logging
import pickle
import numpy as np
import os
import torch
import traceback  # Import traceback
import xgboost as xgb
import threading

from app.config import settings
from app.utils.logger import logger
from app.utils.processing import smiles_to_comprehensive_features, combine_features
from app.ml_models.alvadesc_feature_generation import generate_all_features, FULL_FEATURE_NAMES, MORGAN_FINGERPRINT_COUNT
from app.models import PredictionFeatures, MolecularDescriptors  # Import Pydantic models
from app.celery_instance import celery_app  # Import celery_app from the new instance file
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import Crippen


# Removed global classifier_model and top-level load_models()

def sigmoid(x):
    return 1 / (1 + np.exp(-x))


def calculate_classification_confidence(classifier_proba: np.ndarray) -> Dict[str, float]:
    """Calculate confidence metrics for classification predictions."""
    max_proba = float(np.max(classifier_proba))
    return {"confidence": max_proba, "uncertainty": 1.0 - max_proba, "class_probabilities": classifier_proba.tolist()}


class DMatrixCreationThread(threading.Thread):
    def __init__(self, feature_vector, feature_names):
        super().__init__()
        self.feature_vector = feature_vector
        self.feature_names = feature_names
        self.dmatrix = None
        self.error = None

    def run(self):
        try:
            logger.info("DMatrixCreationThread: Starting DMatrix creation.")
            self.dmatrix = xgb.DMatrix(self.feature_vector, feature_names=self.feature_names)
            logger.info("DMatrixCreationThread: DMatrix created successfully.")
        except Exception as e:
            logger.exception(f"DMatrixCreationThread: Error during xgb.DMatrix creation in thread: {e}")
            self.error = e


@celery_app.task(bind=True, name="predict_permeability")
def predict_permeability(
    self, smiles_list: List[str], created_at: Optional[str] = None, job_name: Optional[str] = None
) -> Dict[str, Any]:
    logger.info(f"predict_permeability received smiles_list: type={type(smiles_list)}, content={smiles_list}")
    logger.info(f"XGBoost version: {xgb.__version__}")
    xgb.set_config(verbosity=3) # Enable verbose XGBoost logging
    if not isinstance(smiles_list, list):
        logger.error(f"smiles_list is not a list! Type: {type(smiles_list)}")
        raise TypeError("smiles_list must be a list of strings")
    if not all(isinstance(s, str) for s in smiles_list):
        logger.error(f"smiles_list contains non-string elements: {smiles_list}")
        raise TypeError("All elements in smiles_list must be strings")

    try:
        # Load classifier model inside the task for process isolation
        classifier_model = None
        classifier_path = settings.MODEL_CLASSIFIER_PATH
        logger.info(f"Attempting to load classifier model from: {classifier_path} inside task.")
        if os.path.exists(classifier_path):
            with open(classifier_path, "rb") as f:
                classifier_model = xgb.Booster()
                classifier_model.load_model(classifier_path)
            logger.info("Classifier model loaded successfully inside task.")
            if classifier_model.feature_names:
                logger.info(f"Classifier model feature names (first 10): {classifier_model.feature_names[:10]}")
                logger.info(f"Classifier model feature names (last 10): {classifier_model.feature_names[-10:]}")
                logger.info(f"Classifier model feature count: {len(classifier_model.feature_names)}")
            else:
                logger.info("Classifier model has no feature names.")
        else:
            logger.error("Classifier model could not be loaded inside task - check model file path")
            raise FileNotFoundError(f"Model file not found: {classifier_path}")

        results = []

        for smiles in smiles_list:
            processed_features_for_result = None  # Initialize here
            try:
                mol = Chem.MolFromSmiles(smiles)
                if not mol:
                    raise ValueError(f"Invalid SMILES string: {smiles}")

                # Extract features using alvaDesc CLI wrapper and Morgan fingerprints
                logger.info(f"Generating features for SMILES: {smiles}")
                logger.info(f"DEBUG: settings.MOCK_ALVADESC = {settings.MOCK_ALVADESC}")  # New line
                all_features_df = generate_all_features(smiles, mock_alvadesc=settings.MOCK_ALVADESC)
                logger.info(f"Features generated for SMILES: {smiles}")
                # logger.info(f"Head of generated features: {all_features_df.head().to_dict(orient='records')}")

                # Convert DataFrame to numpy array for model input
                logger.info("Converting DataFrame to numpy array...")
                feature_vector = all_features_df.values.astype(np.float32)
                logger.info(f"Numpy array created. Shape: {feature_vector.shape}, Dtype: {feature_vector.dtype}")

                # Ensure feature_vector has the correct shape (1, 10160) for a single sample
                logger.info("Checking feature_vector shape...")
                if feature_vector.ndim == 1:
                    feature_vector = feature_vector.reshape(1, -1)
                    logger.info(f"Reshaped 1D feature_vector to {feature_vector.shape}")
                elif feature_vector.ndim > 2:
                    # For now, we'll take the first row if multiple are returned for a single SMILES
                    feature_vector = feature_vector[0].reshape(1, -1)
                    logger.info(f"Reshaped multi-dimensional feature_vector to {feature_vector.shape}")
                logger.info(f"Final feature_vector shape: {feature_vector.shape}")

                # Convert numpy array to DMatrix for XGBoost Booster model
                try:
                    logger.info(f"Preparing to create actual DMatrix in a separate thread. Feature vector shape: {feature_vector.shape}, FULL_FEATURE_NAMES length: {len(FULL_FEATURE_NAMES)}")
                    dmatrix_thread = DMatrixCreationThread(feature_vector, FULL_FEATURE_NAMES)
                    dmatrix_thread.start()
                    dmatrix_thread.join() # Wait for the thread to complete

                    if dmatrix_thread.error:
                        raise dmatrix_thread.error

                    dmatrix_feature_vector = dmatrix_thread.dmatrix
                    logger.info("Actual DMatrix created successfully in main thread after thread join.")
                    logger.info(f"DMatrix has {dmatrix_feature_vector.num_row()} rows and {dmatrix_feature_vector.num_col()} columns.")

                except Exception as dmatrix_e:
                    logger.exception(f"Error creating DMatrix for SMILES {smiles}: {dmatrix_e}")
                    result = {
                        "smiles": smiles,
                        "prediction": 0,
                        "confidence": 0.0,
                        "uncertainty": 1.0,
                        "class_probabilities": [0.5, 0.5],
                        "classifier_prediction": 0,
                        "features": processed_features_for_result,
                        "error": f"Error creating DMatrix: {str(dmatrix_e)}",
                    }
                    results.append(result)
                    continue  # Skip to the next SMILES in the list

                # Store the generated features for the result
                processed_features_for_result = {
                    "morgan_fingerprint": feature_vector[0, :MORGAN_FINGERPRINT_COUNT].tolist(),
                    "descriptors": {"alvadesc_features": feature_vector[0, MORGAN_FINGERPRINT_COUNT:].tolist()}
                }

                # Extract key features for summary using RDKit
                features_summary = {
                    "MW": Descriptors.MolWt(mol),
                    "ALOGP": Crippen.MolLogP(mol),
                    "nHDon": Descriptors.NumHDonors(mol),
                    "TPSA": Descriptors.TPSA(mol),
                }
                logger.info(f"RDKit-derived summary features: {features_summary}")

                logger.info(f"Feature vector shape: {feature_vector.shape}")
                logger.info("Attempting classification prediction...")
                # Classification prediction
                if classifier_model is not None:
                    logger.info(f"Classifier model type: {type(classifier_model)}")
                    logger.info(f"Feature vector dtype: {feature_vector.dtype}, shape: {feature_vector.shape}")
                    logger.info(f"Feature vector min: {np.min(feature_vector)}")
                    logger.info(f"Feature vector max: {np.max(feature_vector)}")
                    logger.info(f"Feature vector mean: {np.mean(feature_vector)}")
                    logger.info(f"Feature vector std: {np.std(feature_vector)}")
                    logger.info(f"Feature vector sample (first 10): {feature_vector.flatten()[:10]}")
                    np.save("/app/feature_vector_debug.npy", feature_vector)
                    logger.info(f"Feature vector saved to /app/feature_vector_debug.npy")
                    logger.info(f"FULL_FEATURE_NAMES length: {len(FULL_FEATURE_NAMES)}")
                    logger.info(f"FULL_FEATURE_NAMES sample (first 10): {FULL_FEATURE_NAMES[:10]}")
                    logger.info(f"FULL_FEATURE_NAMES sample (last 10): {FULL_FEATURE_NAMES[-10:]}")

                    try:
                        # Get raw scores from predict
                        raw_scores = classifier_model.predict(dmatrix_feature_vector, output_margin=True)
                        logger.info(f"Raw scores from classifier.predict: {raw_scores}")

                        # Apply sigmoid to get probabilities for binary classification
                        classifier_prob_positive = sigmoid(raw_scores)
                        classifier_prob_negative = 1 - classifier_prob_positive
                        classifier_prob = np.array([[classifier_prob_negative[0], classifier_prob_positive[0]]])
                        logger.info(f"Classifier probabilities (after sigmoid): {classifier_prob}")

                        # Determine prediction based on probability (e.g., > 0.5 for positive class)
                        classifier_pred = (classifier_prob_positive > 0.5).astype(int)[0]
                        logger.info(f"Classifier prediction: {classifier_pred}")

                        confidence_stats = calculate_classification_confidence(classifier_prob[0])
                    except Exception as xgb_e:
                        logger.exception(f"XGBoost prediction failed for SMILES {smiles}: {xgb_e}")
                        # Set error for this specific SMILES and continue processing other SMILES
                        result = {
                            "smiles": smiles,
                            "prediction": 0,
                            "confidence": 0.0,
                            "uncertainty": 1.0,
                            "class_probabilities": [0.5, 0.5],
                            "classifier_prediction": 0,
                            "features": processed_features_for_result,
                            "features_summary": [], # Add empty features_summary on error
                            "error": f"XGBoost prediction failed: {str(xgb_e)}",
                        }
                        results.append(result)
                        continue  # Skip to the next SMILES in the list
                    logger.info(f"Confidence stats: {confidence_stats}")
                else:
                    # Fallback if no classifier
                    classifier_pred = 0
                    confidence_stats = {"confidence": 0.0, "uncertainty": 1.0, "class_probabilities": [0.5, 0.5]}
                    logger.warning("No classifier model available - defaulting to non-permeant")

                # Convert features_summary dictionary to a list of objects for frontend consumption
                features_summary_list = []
                for name, value in features_summary.items():
                    features_summary_list.append({"name": name, "value": value})

                logger.info(f"Features summary list being sent to frontend: {features_summary_list}")

                logger.info("Classification prediction completed.")

                # Convert classification to binary prediction
                prediction = 1 if classifier_pred == 1 else 0  # 1 = permeant, 0 = non-permeant

                result = {
                    "smiles": smiles,
                    "prediction": prediction,
                    "confidence": confidence_stats["confidence"],
                    "uncertainty": confidence_stats["uncertainty"],
                    "class_probabilities": confidence_stats["class_probabilities"],
                    "classifier_prediction": int(classifier_pred),
                    "features": processed_features_for_result,
                    "features_summary": features_summary_list,
                    "error": None,
                }

            except Exception as e:
                logger.exception(f"Error processing SMILES {smiles}: {e}")
                result = {
                    "smiles": smiles,
                    "prediction": 0,
                    "confidence": 0.0,
                    "uncertainty": 1.0,
                    "class_probabilities": [0.5, 0.5],
                    "classifier_prediction": 0,
                    "features": processed_features_for_result,
                    "features_summary": [], # Add empty features_summary on error
                    "error": str(e),
                }
            results.append(result)

        return {
            "status": "completed",
            "results": results,
            "total_processed": len(results),
            "successful": len([r for r in results if r["error"] is None]),
            "failed": len([r for r in results if r["error"] is not None]),
        }

    except Exception as e:
        logger.error(f"Task failed: {e}")
        traceback.print_exc()  # Print full traceback
        self.retry(countdown=60, max_retries=3)
        return {"status": "failed", "error": str(e), "results": []}


@celery_app.task(bind=True, name="explore_permeability")
def explore_permeability(
    self, original_smiles: str, created_at: Optional[str] = None, job_name: Optional[str] = None
) -> Dict[str, Any]:
    logger.info(f"explore_permeability received original_smiles: {original_smiles}")
    
    # Load classifier model (similar to predict_permeability)
    classifier_model = None
    classifier_path = settings.MODEL_CLASSIFIER_PATH
    if os.path.exists(classifier_path):
        with open(classifier_path, "rb") as f:
            classifier_model = xgb.Booster()
            classifier_model.load_model(classifier_path)
        logger.info("Classifier model loaded successfully inside explore_permeability task.")
    else:
        logger.error("Classifier model could not be loaded inside explore_permeability task")
        return {"status": "failed", "error": "Classifier model not found", "results": []}

    results = []
    permuted_smiles_list = []

    try:
        mol = Chem.MolFromSmiles(original_smiles)
        if not mol:
            raise ValueError("Invalid SMILES string provided.")

        # --- RDKit-based Permutation Logic (Refined) ---
        # This refined approach focuses on generating valid SMILES through simpler, more controlled modifications.
        # We will try to add small groups (e.g., methyl, hydroxyl) at various positions.
        
        # Ensure the original molecule is valid
        if not mol:
            raise ValueError("Invalid SMILES string provided.")

        # List to store unique valid permuted SMILES
        unique_permuted_smiles = set()
        
        # Add original SMILES to ensure it's always part of the results
        unique_permuted_smiles.add(original_smiles)

        # Strategy 1: Add a methyl group at various positions
        for _ in range(5): # Attempt up to 5 methyl additions
            if len(unique_permuted_smiles) >= 10: break
            try:
                temp_mol = Chem.Mol(mol)
                rw_mol = Chem.RWMol(temp_mol)
                
                # Find a random atom to attach a methyl group
                # Prioritize carbons that are not fully saturated
                eligible_atoms = [atom.GetIdx() for atom in rw_mol.GetAtoms() if atom.GetAtomicNum() == 6 and (atom.GetImplicitValence() > 0 or atom.GetExplicitValence() < atom.GetTotalValence())]
                
                if not eligible_atoms: continue
                
                attach_idx = int(np.random.choice(eligible_atoms))
                
                new_atom = Chem.Atom('C')
                new_atom_idx = rw_mol.AddAtom(new_atom)
                rw_mol.AddBond(attach_idx, new_atom_idx, Chem.BondType.SINGLE)
                
                Chem.SanitizeMol(rw_mol) # Attempt to sanitize
                new_smiles = Chem.MolToSmiles(rw_mol)
                if new_smiles and new_smiles not in unique_permuted_smiles:
                    unique_permuted_smiles.add(new_smiles)
                    logger.info(f"Generated permuted SMILES (methyl addition): {new_smiles}")
            except Exception as add_e:
                logger.warning(f"Error adding methyl group for {original_smiles}: {add_e}")
                
        # Strategy 2: Replace a hydrogen with a hydroxyl group at various positions
        for _ in range(5): # Attempt up to 5 hydroxyl substitutions
            if len(unique_permuted_smiles) >= 10: break
            try:
                temp_mol = Chem.Mol(mol)
                temp_mol_h = Chem.AddHs(temp_mol) # Add hydrogens to find replaceable ones
                
                # Find a hydrogen attached to a carbon to replace
                replaceable_h_info = [] # (H_atom_idx, C_atom_idx)
                for h_atom in temp_mol_h.GetAtoms():
                    if h_atom.GetAtomicNum() == 1: # Is a hydrogen
                        for neighbor in h_atom.GetNeighbors():
                            if neighbor.GetAtomicNum() == 6: # Attached to a carbon
                                replaceable_h_info.append((h_atom.GetIdx(), neighbor.GetIdx()))
                                break
                
                if not replaceable_h_info: continue
                
                h_to_replace_idx, c_attached_to_h_idx = replaceable_h_info[np.random.choice(len(replaceable_h_info))]
                
                # Create a reaction to replace H with OH
                # This uses a generic reaction SMARTS for H -> OH substitution
                rxn_smarts = f"[C:1]-[H:2]>>[C:1]-[OH:2]"
                rxn = AllChem.ReactionFromSmarts(rxn_smarts)
                
                # Need to apply the reaction to the original molecule, not the one with added Hs
                # This is a bit tricky with generic SMARTS. A more direct way is to use EditableMol.
                
                # Using EditableMol for more control
                rw_mol = Chem.RWMol(temp_mol)
                
                # Find the carbon atom in the original molecule that had the H
                original_c_atom = rw_mol.GetAtomWithIdx(c_attached_to_h_idx)
                
                # Check if replacing H with OH would violate valency
                # This is a heuristic: if carbon has 4 bonds, it can't take another without losing something
                if original_c_atom.GetExplicitValence() >= original_c_atom.GetTotalValence():
                    continue # Cannot replace H with OH without breaking valency

                # Remove a hydrogen (conceptually) and add an oxygen
                # This is a simplified representation. RDKit's ReplaceSubstructs is better for this.
                # For now, let's try a direct bond replacement if possible.
                
                # A more robust way to replace H with OH is to use a specific reaction
                # or to build the molecule step-by-step with EditableMol.
                
                # Let's try a simpler approach: find a C-H bond and replace H with OH
                # This requires finding the bond and then modifying it.
                
                # Simpler: just add OH to a carbon with available valence
                eligible_oh_atoms = [atom.GetIdx() for atom in rw_mol.GetAtoms() if atom.GetAtomicNum() == 6 and atom.GetImplicitValence() > 0]
                if not eligible_oh_atoms: continue
                
                attach_oh_idx = int(np.random.choice(eligible_oh_atoms))
                
                new_o_atom = Chem.Atom('O')
                new_h_atom = Chem.Atom('H')
                
                new_o_idx = rw_mol.AddAtom(new_o_atom)
                rw_mol.AddBond(attach_oh_idx, new_o_idx, Chem.BondType.SINGLE)
                # The H will be added implicitly by SanitizeMol if needed
                
                Chem.SanitizeMol(rw_mol) # Attempt to sanitize
                new_smiles = Chem.MolToSmiles(rw_mol)
                if new_smiles and new_smiles not in unique_permuted_smiles:
                    unique_permuted_smiles.add(new_smiles)
                    logger.info(f"Generated permuted SMILES (hydroxyl addition): {new_smiles}")
            except Exception as sub_e:
                logger.warning(f"Error substituting hydroxyl group for {original_smiles}: {sub_e}")

        # Convert set to list and ensure we have exactly 10 permutations
        permuted_smiles_list = list(unique_permuted_smiles)
        
        # If we have less than 10 unique valid permutations, fill with the original SMILES
        while len(permuted_smiles_list) < 10:
            permuted_smiles_list.append(original_smiles)
        
        # Trim to exactly 10 (if more were generated, e.g., original + 9 permutations)
        permuted_smiles_list = permuted_smiles_list[:10]

        # Add the original SMILES to the list of results as the first entry
        all_smiles_to_process = permuted_smiles_list
        
        for smiles in all_smiles_to_process:
            processed_features_for_result = None
            try:
                logger.info(f"Generating features for permuted SMILES: {smiles}")
                all_features_df = generate_all_features(smiles, mock_alvadesc=settings.MOCK_ALVADESC)
                feature_vector = all_features_df.values.astype(np.float32)

                if feature_vector.ndim == 1:
                    feature_vector = feature_vector.reshape(1, -1)
                elif feature_vector.ndim > 2:
                    feature_vector = feature_vector[0].reshape(1, -1)

                dmatrix_thread = DMatrixCreationThread(feature_vector, FULL_FEATURE_NAMES)
                dmatrix_thread.start()
                dmatrix_thread.join()

                if dmatrix_thread.error:
                    raise dmatrix_thread.error

                dmatrix_feature_vector = dmatrix_thread.dmatrix

                processed_features_for_result = {
                    "morgan_fingerprint": feature_vector[0, :MORGAN_FINGERPRINT_COUNT].tolist(),
                    "descriptors": {"alvadesc_features": feature_vector[0, MORGAN_FINGERPRINT_COUNT:].tolist()}
                }

                features_summary = {}
                summary_feature_names = ["MW", "ALOGP", "nHDon", "TPSA"]
                for sf_name in summary_feature_names:
                    try:
                        feature_index = FULL_FEATURE_NAMES.index(sf_name)
                        features_summary[sf_name] = float(feature_vector[0, feature_index])
                    except (ValueError, IndexError):
                        features_summary[sf_name] = 0.0

                if classifier_model is not None:
                    raw_scores = classifier_model.predict(dmatrix_feature_vector, output_margin=True)
                    classifier_prob_positive = sigmoid(raw_scores)
                    classifier_prob_negative = 1 - classifier_prob_positive
                    classifier_prob = np.array([[classifier_prob_negative[0], classifier_prob_positive[0]]])
                    classifier_pred = (classifier_prob_positive > 0.5).astype(int)[0]
                    confidence_stats = calculate_classification_confidence(classifier_prob[0])
                else:
                    classifier_pred = 0
                    confidence_stats = {"confidence": 0.0, "uncertainty": 1.0, "class_probabilities": [0.5, 0.5]}

                features_summary_list = []
                for name, value in features_summary.items():
                    features_summary_list.append({"name": name, "value": value})

                prediction = 1 if classifier_pred == 1 else 0

                results.append({
                    "smiles": smiles,
                    "prediction": prediction,
                    "confidence": confidence_stats["confidence"],
                    "uncertainty": confidence_stats["uncertainty"],
                    "class_probabilities": confidence_stats["class_probabilities"],
                    "classifier_prediction": int(classifier_pred),
                    "features": processed_features_for_result,
                    "features_summary": features_summary_list,
                    "error": None,
                })

            except Exception as e:
                logger.exception(f"Error processing permuted SMILES {smiles}: {e}")
                results.append({
                    "smiles": smiles,
                    "prediction": 0,
                    "confidence": 0.0,
                    "uncertainty": 1.0,
                    "class_probabilities": [0.5, 0.5],
                    "classifier_prediction": 0,
                    "features": processed_features_for_result,
                    "features_summary": [],
                    "error": str(e),
                })

        return {
            "status": "completed",
            "results": results,
            "total_processed": len(results),
            "successful": len([r for r in results if r["error"] is None]),
            "failed": len([r for r in results if r["error"] is not None]),
        }

    except Exception as e:
        logger.error(f"Exploration task failed for {original_smiles}: {e}")
        traceback.print_exc()
        self.retry(countdown=60, max_retries=3)
        return {"status": "failed", "error": str(e), "results": []}

    logger.info(f"predict_permeability received smiles_list: type={type(smiles_list)}, content={smiles_list}")
    logger.info(f"XGBoost version: {xgb.__version__}")
    xgb.set_config(verbosity=3) # Enable verbose XGBoost logging
    if not isinstance(smiles_list, list):
        logger.error(f"smiles_list is not a list! Type: {type(smiles_list)}")
        raise TypeError("smiles_list must be a list of strings")
    if not all(isinstance(s, str) for s in smiles_list):
        logger.error(f"smiles_list contains non-string elements: {smiles_list}")
        raise TypeError("All elements in smiles_list must be strings")

    try:
        # Load classifier model inside the task for process isolation
        classifier_model = None
        classifier_path = settings.MODEL_CLASSIFIER_PATH
        logger.info(f"Attempting to load classifier model from: {classifier_path} inside task.")
        if os.path.exists(classifier_path):
            with open(classifier_path, "rb") as f:
                classifier_model = xgb.Booster()
                classifier_model.load_model(classifier_path)
            logger.info("Classifier model loaded successfully inside task.")
            if classifier_model.feature_names:
                logger.info(f"Classifier model feature names (first 10): {classifier_model.feature_names[:10]}")
                logger.info(f"Classifier model feature names (last 10): {classifier_model.feature_names[-10:]}")
                logger.info(f"Classifier model feature count: {len(classifier_model.feature_names)}")
            else:
                logger.info("Classifier model has no feature names.")
        else:
            logger.error("Classifier model could not be loaded inside task - check model file path")
            raise FileNotFoundError(f"Model file not found: {classifier_path}")

        results = []

        for smiles in smiles_list:
            processed_features_for_result = None  # Initialize here
            try:
                # Extract features using alvaDesc CLI wrapper and Morgan fingerprints
                logger.info(f"Generating features for SMILES: {smiles}")
                logger.info(f"DEBUG: settings.MOCK_ALVADESC = {settings.MOCK_ALVADESC}")  # New line
                all_features_df = generate_all_features(smiles, mock_alvadesc=settings.MOCK_ALVADESC)
                logger.info(f"Features generated for SMILES: {smiles}")
                # logger.info(f"Head of generated features: {all_features_df.head().to_dict(orient='records')}")

                # Convert DataFrame to numpy array for model input
                logger.info("Converting DataFrame to numpy array...")
                feature_vector = all_features_df.values.astype(np.float32)
                logger.info(f"Numpy array created. Shape: {feature_vector.shape}, Dtype: {feature_vector.dtype}")

                # Ensure feature_vector has the correct shape (1, 10160) for a single sample
                logger.info("Checking feature_vector shape...")
                if feature_vector.ndim == 1:
                    feature_vector = feature_vector.reshape(1, -1)
                    logger.info(f"Reshaped 1D feature_vector to {feature_vector.shape}")
                elif feature_vector.ndim > 2:
                    # For now, we'll take the first row if multiple are returned for a single SMILES
                    feature_vector = feature_vector[0].reshape(1, -1)
                    logger.info(f"Reshaped multi-dimensional feature_vector to {feature_vector.shape}")
                logger.info(f"Final feature_vector shape: {feature_vector.shape}")

                # Convert numpy array to DMatrix for XGBoost Booster model
                try:
                    logger.info(f"Preparing to create actual DMatrix in a separate thread. Feature vector shape: {feature_vector.shape}, FULL_FEATURE_NAMES length: {len(FULL_FEATURE_NAMES)}")
                    dmatrix_thread = DMatrixCreationThread(feature_vector, FULL_FEATURE_NAMES)
                    dmatrix_thread.start()
                    dmatrix_thread.join() # Wait for the thread to complete

                    if dmatrix_thread.error:
                        raise dmatrix_thread.error

                    dmatrix_feature_vector = dmatrix_thread.dmatrix
                    logger.info("Actual DMatrix created successfully in main thread after thread join.")
                    logger.info(f"DMatrix has {dmatrix_feature_vector.num_row()} rows and {dmatrix_feature_vector.num_col()} columns.")

                except Exception as dmatrix_e:
                    logger.exception(f"Error creating DMatrix for SMILES {smiles}: {dmatrix_e}")
                    result = {
                        "smiles": smiles,
                        "prediction": 0,
                        "confidence": 0.0,
                        "uncertainty": 1.0,
                        "class_probabilities": [0.5, 0.5],
                        "classifier_prediction": 0,
                        "features": processed_features_for_result,
                        "error": f"Error creating DMatrix: {str(dmatrix_e)}",
                    }
                    results.append(result)
                    continue  # Skip to the next SMILES in the list

                # Store the generated features for the result
                processed_features_for_result = {
                    "morgan_fingerprint": feature_vector[0, :MORGAN_FINGERPRINT_COUNT].tolist(),
                    "descriptors": {"alvadesc_features": feature_vector[0, MORGAN_FINGERPRINT_COUNT:].tolist()}
                }

                # Extract key features for summary
                features_summary = {}
                summary_feature_names = ["MW", "ALOGP", "nHDon", "TPSA"] # Common features
                
                for sf_name in summary_feature_names:
                    try:
                        # Find the index of the feature name in FULL_FEATURE_NAMES
                        # This assumes FULL_FEATURE_NAMES is a list of strings
                        feature_index = FULL_FEATURE_NAMES.index(sf_name)
                        features_summary[sf_name] = float(feature_vector[0, feature_index])
                    except ValueError:
                        logger.warning(f"Summary feature '{sf_name}' not found in FULL_FEATURE_NAMES.")
                        features_summary[sf_name] = 0.0 # Default or handle as appropriate
                    except IndexError:
                        logger.warning(f"Index error for summary feature '{sf_name}'. Feature vector might be too small.")
                        features_summary[sf_name] = 0.0 # Default or handle as appropriate

                logger.info(f"Feature vector shape: {feature_vector.shape}")
                logger.info("Attempting classification prediction...")
                # Classification prediction
                if classifier_model is not None:
                    logger.info(f"Classifier model type: {type(classifier_model)}")
                    logger.info(f"Feature vector dtype: {feature_vector.dtype}, shape: {feature_vector.shape}")
                    logger.info(f"Feature vector min: {np.min(feature_vector)}")
                    logger.info(f"Feature vector max: {np.max(feature_vector)}")
                    logger.info(f"Feature vector mean: {np.mean(feature_vector)}")
                    logger.info(f"Feature vector std: {np.std(feature_vector)}")
                    logger.info(f"Feature vector sample (first 10): {feature_vector.flatten()[:10]}")
                    np.save("/app/feature_vector_debug.npy", feature_vector)
                    logger.info(f"Feature vector saved to /app/feature_vector_debug.npy")
                    logger.info(f"FULL_FEATURE_NAMES length: {len(FULL_FEATURE_NAMES)}")
                    logger.info(f"FULL_FEATURE_NAMES sample (first 10): {FULL_FEATURE_NAMES[:10]}")
                    logger.info(f"FULL_FEATURE_NAMES sample (last 10): {FULL_FEATURE_NAMES[-10:]}")

                    try:
                        # Get raw scores from predict
                        raw_scores = classifier_model.predict(dmatrix_feature_vector, output_margin=True)
                        logger.info(f"Raw scores from classifier.predict: {raw_scores}")

                        # Apply sigmoid to get probabilities for binary classification
                        classifier_prob_positive = sigmoid(raw_scores)
                        classifier_prob_negative = 1 - classifier_prob_positive
                        classifier_prob = np.array([[classifier_prob_negative[0], classifier_prob_positive[0]]])
                        logger.info(f"Classifier probabilities (after sigmoid): {classifier_prob}")

                        # Determine prediction based on probability (e.g., > 0.5 for positive class)
                        classifier_pred = (classifier_prob_positive > 0.5).astype(int)[0]
                        logger.info(f"Classifier prediction: {classifier_pred}")

                        confidence_stats = calculate_classification_confidence(classifier_prob[0])
                    except Exception as xgb_e:
                        logger.exception(f"XGBoost prediction failed for SMILES {smiles}: {xgb_e}")
                        # Set error for this specific SMILES and continue processing other SMILES
                        result = {
                            "smiles": smiles,
                            "prediction": 0,
                            "confidence": 0.0,
                            "uncertainty": 1.0,
                            "class_probabilities": [0.5, 0.5],
                            "classifier_prediction": 0,
                            "features": processed_features_for_result,
                            "features_summary": [], # Add empty features_summary on error
                            "error": f"XGBoost prediction failed: {str(xgb_e)}",
                        }
                        results.append(result)
                        continue  # Skip to the next SMILES in the list
                    logger.info(f"Confidence stats: {confidence_stats}")
                else:
                    # Fallback if no classifier
                    classifier_pred = 0
                    confidence_stats = {"confidence": 0.0, "uncertainty": 1.0, "class_probabilities": [0.5, 0.5]}
                    logger.warning("No classifier model available - defaulting to non-permeant")

                # Convert features_summary dictionary to a list of objects for frontend consumption
                features_summary_list = []
                for name, value in features_summary.items():
                    features_summary_list.append({"name": name, "value": value})

                logger.info(f"Features summary list being sent to frontend: {features_summary_list}")

                logger.info("Classification prediction completed.")

                # Convert classification to binary prediction
                prediction = 1 if classifier_pred == 1 else 0  # 1 = permeant, 0 = non-permeant

                result = {
                    "smiles": smiles,
                    "prediction": prediction,
                    "confidence": confidence_stats["confidence"],
                    "uncertainty": confidence_stats["uncertainty"],
                    "class_probabilities": confidence_stats["class_probabilities"],
                    "classifier_prediction": int(classifier_pred),
                    "features": processed_features_for_result,
                    "features_summary": features_summary_list,
                    "error": None,
                }

            except Exception as e:
                logger.exception(f"Error processing SMILES {smiles}: {e}")
                result = {
                    "smiles": smiles,
                    "prediction": 0,
                    "confidence": 0.0,
                    "uncertainty": 1.0,
                    "class_probabilities": [0.5, 0.5],
                    "classifier_prediction": 0,
                    "features": processed_features_for_result,
                    "features_summary": [], # Add empty features_summary on error
                    "error": str(e),
                }
            results.append(result)

        return {
            "status": "completed",
            "results": results,
            "total_processed": len(results),
            "successful": len([r for r in results if r["error"] is None]),
            "failed": len([r for r in results if r["error"] is not None]),
        }

    except Exception as e:
        logger.error(f"Task failed: {e}")
        traceback.print_exc()  # Print full traceback
        self.retry(countdown=60, max_retries=3)
        return {"status": "failed", "error": str(e), "results": []}