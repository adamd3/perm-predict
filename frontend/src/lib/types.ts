// GraphQL Types
export interface JobResponse {
  jobId: string;
  status: 'pending' | 'processing' | 'completed' | 'failed';
}

export interface PredictionResult {
  smiles: string;
  prediction: number; // 0 = impermeant, 1 = permeant
  probability: number; // Probability of the predicted class
  confidence: number; // Max probability (confidence score)
  processingTime?: number; // Processing time in seconds
}

export interface JobResult {
  jobId: string;
  status: 'pending' | 'processing' | 'completed' | 'failed';
  result?: PredictionResult | PredictionResult[]; // Single or batch results
  error?: string;
}

// Form Types
export interface SinglePredictionRequest {
  smiles: string;
}

export interface BatchPredictionRequest {
  smilesStrings: string[];
}

// Component Props Types
export interface PredictionResultsProps {
  results: PredictionResult[];
}