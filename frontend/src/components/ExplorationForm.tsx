'use client';

import React, { useState, useEffect } from 'react';
import { useMutation, useLazyQuery } from '@apollo/client';
import { Alert, AlertDescription } from '@/components/ui/alert';
import { Button } from '@/components/ui/button';
import { Input } from '@/components/ui/input';
import { Progress } from '@/components/ui/progress';
import { Loader2, CheckCircle, AlertCircle } from 'lucide-react';
import PredictionResults from './PredictionResults';
import { SUBMIT_EXPLORATION_JOB, GET_EXPLORATION_RESULT, GET_JOB_STATUS } from '@/lib/graphql/queries'; // Will define SUBMIT_EXPLORATION_JOB and GET_EXPLORATION_RESULT

import type { PredictionResult, JobStatus, JobResult } from '@/lib/types'

interface ExplorationFormProps {
  initialSmiles?: string;
  onResultsLoaded?: (hasResults: boolean) => void;
}

const ExplorationForm = ({ initialSmiles = '', onResultsLoaded }: ExplorationFormProps) => {
  const [smilesInput, setSmilesInput] = useState(initialSmiles);
  const [results, setResults] = useState<PredictionResult[]>([]);
  const [error, setError] = useState('');
  const [currentJobId, setCurrentJobId] = useState<string | null>(null);
  const [jobStatus, setJobStatus] = useState<'idle' | 'pending' | 'processing' | 'completed' | 'failed' | 'retrying' | 'cancelled' | 'error' | 'submitted'>('idle');
  const [progress, setProgress] = useState(0);

  // GraphQL hooks
  const [submitExplorationJobMutation] = useMutation(SUBMIT_EXPLORATION_JOB);
  const [getJobStatus, { data: jobStatusData, startPolling: startStatusPolling, stopPolling: stopStatusPolling }] = useLazyQuery(GET_JOB_STATUS, {
    pollInterval: 1000, // Poll every 1 second for job status
    errorPolicy: 'all',
  });
  const [getExplorationResult, { data: explorationResultData }] = useLazyQuery(GET_EXPLORATION_RESULT);

  useEffect(() => {
    if (currentJobId) {
      getJobStatus({ variables: { jobId: currentJobId } });
      startStatusPolling(1000);
    }
  }, [currentJobId, getJobStatus, startStatusPolling]);

  // Effect to handle job polling results
  useEffect(() => {
    if (jobStatusData?.getJobStatus) {
      const currentStatus = jobStatusData.getJobStatus as JobStatus;
      setJobStatus(currentStatus.status);
      
      // Update progress based on backend progress or simulate if not available
      if (currentStatus.progress) {
        const match = currentStatus.progress.match(/(\d+)%/);
        if (match) {
          setProgress(parseInt(match[1], 10));
        } else if (currentStatus.status === 'processing') {
          setProgress(prev => Math.min(prev + 5, 95)); // Simulate progress if no percentage
        } else if (currentStatus.status === 'completed') {
          setProgress(100);
        }
      }

      if (currentStatus.status === 'completed') {
        stopStatusPolling();
        getExplorationResult({ variables: { jobId: currentStatus.jobId } }); // Fetch final results
      } else if (currentStatus.status === 'failed' || currentStatus.status === 'cancelled') {
        setError(currentStatus.error || 'Exploration failed');
        stopStatusPolling();
        setJobStatus('idle');
        setCurrentJobId(null);
        setProgress(0);
      }
    }
  }, [jobStatusData, stopStatusPolling, getExplorationResult]);

  useEffect(() => {
    if (explorationResultData?.getExplorationResult) {
      const result = explorationResultData.getExplorationResult as JobResult;
      console.log("Received explorationResultData:", result);
      if (result.results) {
        const processedResults = result.results.map(pr => ({
          ...pr,
          permeantProbability: pr.classProbabilities && pr.classProbabilities.length > 1 ? pr.classProbabilities[1] : 0,
        }));
        setResults(processedResults);
        console.log("Updated results state with:", processedResults);
        if (onResultsLoaded) {
          onResultsLoaded(processedResults.length > 0);
        }
      }
      // Reset after a delay
      setTimeout(() => {
        setJobStatus('idle');
        setCurrentJobId(null);
        setProgress(0);
      }, 2000);
    }
  }, [explorationResultData, onResultsLoaded]);

  const handleExplorationSubmission = async (e: React.FormEvent) => {
    e.preventDefault();
    setError('');
    setResults([]);
    setProgress(0);
    if (onResultsLoaded) {
      onResultsLoaded(false);
    }
    
    try {
      const { data } = await submitExplorationJobMutation({
        variables: { jobInput: { smilesList: [smilesInput.trim()], jobName: 'Compound Exploration' } }
      });
      
      if (data?.submitExplorationJob) {
        const jobResponse = data.submitExplorationJob as JobStatus;
        console.log("Job Response:", jobResponse);
        setCurrentJobId(jobResponse.jobId);
        setJobStatus('pending');
        setProgress(10);
      }
    } catch (err) {
      setError('Failed to submit exploration. Please try again.');
      setJobStatus('idle');
    }
  };

  const getStatusIcon = () => {
    switch (jobStatus) {
      case 'pending':
      case 'processing':
        return <Loader2 className="h-4 w-4 animate-spin" />;
      case 'completed':
        return <CheckCircle className="h-4 w-4 text-green-600" />;
      case 'failed':
        return <AlertCircle className="h-4 w-4 text-red-600" />;
      default:
        return null;
    }
  };

  const getStatusText = () => {
    switch (jobStatus) {
      case 'pending':
        return 'Queued for exploration...';
      case 'processing':
        return 'Generating and predicting related compounds...';
      case 'completed':
        return 'Exploration completed!';
      case 'failed':
        return 'Exploration failed';
      default:
        return '';
    }
  };

  const isProcessing = jobStatus === 'pending' || jobStatus === 'processing';

  return (
    <>
      <form onSubmit={handleExplorationSubmission} className="space-y-4 w-full">
        <div>
          <Input
            placeholder="Enter SMILES string (e.g., CCO for ethanol)..."
            value={smilesInput}
            onChange={(e) => setSmilesInput(e.target.value)}
            className="w-full font-mono"
            disabled={isProcessing}
          />
        </div>
        <Button type="submit" disabled={isProcessing || !smilesInput} className="w-full">
          {isProcessing ? (
            <>
              <Loader2 className="mr-2 h-4 w-4 animate-spin" />
              Exploring...
            </>
          ) : (
            'Explore Compounds'
          )}
        </Button>
      </form>

      {/* Progress indicator */}
      {isProcessing && (
        <div className="mt-6 space-y-2">
          <div className="flex items-center justify-between text-sm">
            <span className="flex items-center gap-2 text-white">
              {getStatusIcon()}
              {getStatusText()}
            </span>
            <span>{Math.round(progress)}%</span>
          </div>
          <Progress value={progress} className="w-full" />
        </div>
      )}

      {error && (
        <Alert variant="destructive" className="mt-4">
          <AlertCircle className="h-4 w-4" />
          <AlertDescription>{error}</AlertDescription>
        </Alert>
      )}

      {results.length > 0 && <PredictionResults results={results} />}
    </>
  );
};

export default ExplorationForm;
