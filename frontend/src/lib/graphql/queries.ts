import { gql } from '@apollo/client';

export const SUBMIT_PREDICTION_JOB = gql`
  mutation SubmitPredictionJob($smiles: String!) {
    submitPredictionJob(smiles: $smiles) {
      jobId
      status
    }
  }
`;

export const SUBMIT_BATCH_PREDICTION_JOB = gql`
  mutation SubmitBatchPredictionJob($smilesStrings: [String!]!) {
    submitBatchPredictionJob(smilesStrings: $smilesStrings) {
      jobId
      status
    }
  }
`;

export const GET_PREDICTION_RESULT = gql`
  query GetPredictionResult($jobId: String!) {
    getPredictionResult(jobId: $jobId) {
      jobId
      status
      result {
        smiles
        prediction
        probability
        confidence
        processingTime
      }
      error
    }
  }
`;