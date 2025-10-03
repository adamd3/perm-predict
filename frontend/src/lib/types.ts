// GraphQL Types from Backend Schema
export interface MolecularDescriptors {
  molWt: number;
  logP: number;
  tpsa: number;
  numHDonors: number;
  numHAcceptors: number;
  numRotatableBonds: number;
  numAromaticRings: number;
}

export interface PredictionFeatures {
  morganFingerprint: number[];
  descriptors: MolecularDescriptors;
  // For now, we'll also include the raw alvadesc_feature_vector
  alvadescFeatureVector?: number[];
}

export interface PredictionResult {
  smiles: string;
  prediction: number; // Predicted permeability value (float from backend, but we'll display binary)
  confidence: number; // Model confidence score
  permeantProbability: number; // Probability of being permeant (prob_1)
  uncertainty?: number; // Prediction uncertainty from ensemble variance
  ensembleStd?: number; // Standard deviation of ensemble predictions
  classifierPrediction: number; // Binary classifier prediction (0 or 1)
  ensemblePredictions?: number[]; // Individual ensemble model predictions
  classProbabilities: number[]; // Probabilities for each class [prob_0, prob_1]
  features?: PredictionFeatures; // Extracted molecular features
  error?: string; // Error message if prediction failed
}

export interface JobResult {
  results: PredictionResult[];
  totalProcessed: number;
  successful: number;
  failed: number;
  jobId: string;
  createdAt: string;
  completedAt?: string;
  error?: string;
}

export interface JobStatus {
  jobId: string;
  status: 'pending' | 'processing' | 'completed' | 'failed' | 'retrying' | 'cancelled' | 'error' | 'submitted';
  createdAt: string;
  progress?: string;
  error?: string;
}

// Input type for submitting prediction jobs
export interface PredictionJobInput {
  smilesList: string[];
  jobName?: string;
}

// Component Props Types
export interface PredictionResultsProps {
  results: PredictionResult[];
}