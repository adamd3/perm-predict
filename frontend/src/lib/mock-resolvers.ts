import { JobStatus, JobResult, PredictionResult } from './types';

// Mock job storage (in production this would be Redis/database)
const mockJobs: Record<string, JobResult> = {};

// Generate realistic mock predictions based on SMILES complexity
function generateMockPrediction(smiles: string): PredictionResult {
  // Simple heuristic: longer SMILES = more complex = less likely to be permeant
  const complexity = smiles.length;
  const hasRings = smiles.includes('1') || smiles.includes('2');
  const hasAromatics = smiles.includes('c') || smiles.includes('n');

  // Simulate realistic probability distribution
  let baseProb = 0.5;
  if (complexity > 20) baseProb -= 0.2;
  if (hasRings) baseProb += 0.1;
  if (hasAromatics) baseProb += 0.15;

  // Add some randomness
  const randomFactor = (Math.random() - 0.5) * 0.3;
  const probPermeant = Math.max(0.05, Math.min(0.95, baseProb + randomFactor));
  const probImpermeant = 1 - probPermeant;

  const prediction = probPermeant > 0.5 ? 1 : 0;
  const confidence = Math.max(probPermeant, probImpermeant);

  return {
    smiles,
    prediction,
    confidence: Math.round(confidence * 1000) / 1000,
    permeantProbability: probPermeant,
    classifierPrediction: prediction,
    classProbabilities: [probImpermeant, probPermeant],
  };
}

export const mockResolvers = {
  Mutation: {
    submitPredictionJob: (_: any, { smiles }: { smiles: string }): JobStatus => {
      const jobId = `job_${Date.now()}_${Math.random().toString(36).substr(2, 9)}`;
      const createdAt = new Date().toISOString();

      // Create job with pending status
      mockJobs[jobId] = {
        jobId,
        createdAt,
        results: [],
        totalProcessed: 0,
        successful: 0,
        failed: 0,
      };

      // Simulate async processing
      setTimeout(() => {
        // No status update here, as mockJobs stores JobResult, not JobStatus
      }, 500);

      setTimeout(() => {
        if (mockJobs[jobId]) {
          try {
            const result = generateMockPrediction(smiles);
            mockJobs[jobId] = {
              ...mockJobs[jobId],
              results: [result],
              totalProcessed: 1,
              successful: 1,
              completedAt: new Date().toISOString(),
            };
          } catch (error) {
            mockJobs[jobId] = {
              ...mockJobs[jobId],
              error: 'Invalid SMILES string',
              completedAt: new Date().toISOString(),
            };
          }
        }
      }, Math.random() * 3000 + 2000); // 2-5 seconds processing time

      return { jobId, status: 'pending', createdAt };
    },

    submitBatchPredictionJob: (_: any, { smilesStrings }: { smilesStrings: string[] }): JobStatus => {
      const jobId = `batch_${Date.now()}_${Math.random().toString(36).substr(2, 9)}`;
      const createdAt = new Date().toISOString();

      mockJobs[jobId] = {
        jobId,
        createdAt,
        results: [],
        totalProcessed: 0,
        successful: 0,
        failed: 0,
      };

      setTimeout(() => {
        // No status update here, as mockJobs stores JobResult, not JobStatus
      }, 1000);

      setTimeout(() => {
        if (mockJobs[jobId]) {
          try {
            const results = smilesStrings.map(generateMockPrediction);
            mockJobs[jobId] = {
              ...mockJobs[jobId],
              results,
              totalProcessed: results.length,
              successful: results.length,
              completedAt: new Date().toISOString(),
            };
          } catch (error) {
            mockJobs[jobId] = {
              ...mockJobs[jobId],
              error: 'Failed to process batch',
              completedAt: new Date().toISOString(),
            };
          }
        }
      }, Math.random() * 5000 + 3000); // 3-8 seconds for batch

      return { jobId, status: 'pending', createdAt };
    },
  },

  Query: {
    getPredictionResult: (_: any, { jobId }: { jobId: string }): JobResult | null => {
      const job = mockJobs[jobId];
      if (!job) {
        return null;
      }
      return job;
    },
  },
};

// Helper function to clear old jobs (optional cleanup)
export function clearOldMockJobs() {
  const oneHourAgo = Date.now() - 60 * 60 * 1000;
  Object.keys(mockJobs).forEach(jobId => {
    const timestamp = parseInt(jobId.split('_')[1]);
    if (timestamp < oneHourAgo) {
      delete mockJobs[jobId];
    }
  });
}
