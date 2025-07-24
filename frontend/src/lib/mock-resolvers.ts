import { JobResponse, JobResult, PredictionResult } from './types';

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
  const probability = prediction === 1 ? probPermeant : probImpermeant;
  const confidence = Math.max(probPermeant, probImpermeant);
  
  return {
    smiles,
    prediction,
    probability: Math.round(probability * 1000) / 1000,
    confidence: Math.round(confidence * 1000) / 1000,
    processingTime: Math.random() * 5 + 2, // 2-7 seconds
  };
}

export const mockResolvers = {
  Mutation: {
    submitPredictionJob: (_: any, { smiles }: { smiles: string }): JobResponse => {
      const jobId = `job_${Date.now()}_${Math.random().toString(36).substr(2, 9)}`;
      
      // Create job with pending status
      mockJobs[jobId] = {
        jobId,
        status: 'pending',
      };
      
      // Simulate async processing
      setTimeout(() => {
        if (mockJobs[jobId]) {
          mockJobs[jobId].status = 'processing';
        }
      }, 500);
      
      setTimeout(() => {
        if (mockJobs[jobId]) {
          try {
            mockJobs[jobId] = {
              jobId,
              status: 'completed',
              result: generateMockPrediction(smiles),
            };
          } catch (error) {
            mockJobs[jobId] = {
              jobId,
              status: 'failed',
              error: 'Invalid SMILES string',
            };
          }
        }
      }, Math.random() * 3000 + 2000); // 2-5 seconds processing time
      
      return { jobId, status: 'pending' };
    },
    
    submitBatchPredictionJob: (_: any, { smilesStrings }: { smilesStrings: string[] }): JobResponse => {
      const jobId = `batch_${Date.now()}_${Math.random().toString(36).substr(2, 9)}`;
      
      mockJobs[jobId] = {
        jobId,
        status: 'pending',
      };
      
      setTimeout(() => {
        if (mockJobs[jobId]) {
          mockJobs[jobId].status = 'processing';
        }
      }, 1000);
      
      setTimeout(() => {
        if (mockJobs[jobId]) {
          try {
            const results = smilesStrings.map(generateMockPrediction);
            mockJobs[jobId] = {
              jobId,
              status: 'completed',
              result: results,
            };
          } catch (error) {
            mockJobs[jobId] = {
              jobId,
              status: 'failed',
              error: 'Failed to process batch',
            };
          }
        }
      }, Math.random() * 5000 + 3000); // 3-8 seconds for batch
      
      return { jobId, status: 'pending' };
    },
  },
  
  Query: {
    getPredictionResult: (_: any, { jobId }: { jobId: string }): JobResult => {
      const job = mockJobs[jobId];
      if (!job) {
        return {
          jobId,
          status: 'failed',
          error: 'Job not found',
        };
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