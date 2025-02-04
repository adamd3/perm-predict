// API Response Types
export interface PredictionResult {
  smiles: string;
  prediction: number;
  probability: number;
  error?: string;
}

// Form Types
export interface SinglePredictionRequest {
  smiles: string;
}

// Component Props Types
export interface PredictionResultsProps {
  results: PredictionResult[];
}