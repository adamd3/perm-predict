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
        features {
          __typename
        }
        featuresSummary {
          name
          value
        }
        error
        __typename
      }
      __typename
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