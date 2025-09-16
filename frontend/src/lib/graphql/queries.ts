import { gql } from '@apollo/client';

export const SUBMIT_PREDICTION_JOB = gql`
  mutation SubmitPredictionJob($jobInput: PredictionJobInput!) {
    submitPredictionJob(jobInput: $jobInput) {
      jobId
      status
      createdAt
      progress
      error
    }
  }
`;

export const GET_PREDICTION_RESULT = gql`
  query GetPredictionResult($jobId: String!) {
    getPredictionResult(jobId: $jobId) {
      status
      jobId
      createdAt
      completedAt
      totalProcessed
      successful
      failed
      results {
        smiles
        prediction
        confidence
        uncertainty
        classifierPrediction
        classProbabilities
        features {
          # For now, we'll just fetch the alvadesc_feature_vector
          # This will need to be refined once the real alvaDesc output is mapped
          alvadescFeatureVector
        }
        error
      }
    }
  }
`;

export const GET_JOB_STATUS = gql`
  query GetJobStatus($jobId: String!) {
    getJobStatus(jobId: $jobId) {
      jobId
      status
      createdAt
      progress
      error
    }
  }
`;