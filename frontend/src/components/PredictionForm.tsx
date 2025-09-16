import React, { useState, useEffect } from 'react';
import { useMutation, useLazyQuery } from '@apollo/client';
import { Alert, AlertDescription } from '@/components/ui/alert';
import { Button } from '@/components/ui/button';
import { Input } from '@/components/ui/input';
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from '@/components/ui/card';
import { Tabs, TabsContent, TabsList, TabsTrigger } from '@/components/ui/tabs';
import { Progress } from '@/components/ui/progress';
import { Loader2, CheckCircle, AlertCircle } from 'lucide-react';
import PredictionResults from './PredictionResults';
import { SUBMIT_PREDICTION_JOB, GET_PREDICTION_RESULT, GET_JOB_STATUS } from '@/lib/graphql/queries';

import type { PredictionResult, JobStatus, JobResult } from '@/lib/types'

interface PredictionFormProps {
  initialSmiles?: string;
}

const PredictionForm = ({ initialSmiles = '' }: PredictionFormProps) => {
  const [smilesInput, setSmilesInput] = useState(initialSmiles);
  const [batchInput, setBatchInput] = useState('');
  const [results, setResults] = useState<PredictionResult[]>([]);
  const [error, setError] = useState('');
  const [currentJobId, setCurrentJobId] = useState<string | null>(null);
  const [jobStatus, setJobStatus] = useState<'idle' | 'pending' | 'processing' | 'completed' | 'failed'>('idle');
  const [progress, setProgress] = useState(0);

  // GraphQL hooks
  const [submitPredictionJobMutation] = useMutation(SUBMIT_PREDICTION_JOB);
  const [getPredictionResult, { data: jobResult, stopPolling }] = useLazyQuery(GET_PREDICTION_RESULT, {
    pollInterval: 2000,
    errorPolicy: 'all',
  });

  // Effect to handle job polling results
  useEffect(() => {
    if (jobResult?.getPredictionResult) {
      const result = jobResult.getPredictionResult as JobResult;
      setJobStatus(result.status);
      
      if (result.status === 'processing') {
        setProgress(prev => Math.min(prev + 10, 90)); // Simulate progress
      } else if (result.status === 'completed') {
        setProgress(100);
        stopPolling();
        
        if (result.results) {
          setResults(result.results);
        }
        
        // Reset after a delay
        setTimeout(() => {
          setJobStatus('idle');
          setCurrentJobId(null);
          setProgress(0);
        }, 2000);
      } else if (result.status === 'failed') {
        setError(result.error || 'Prediction failed');
        stopPolling();
        setJobStatus('idle');
        setCurrentJobId(null);
        setProgress(0);
      }
    }
  }, [jobResult, stopPolling]);

  // Update input when initialSmiles changes
  useEffect(() => {
    if (initialSmiles) {
      setSmilesInput(initialSmiles);
    }
  }, [initialSmiles]);

  const handleSinglePrediction = async (e: React.FormEvent) => {
    e.preventDefault();
    setError('');
    setResults([]);
    setProgress(0);
    
    try {
      const { data } = await submitPredictionJobMutation({
        variables: { jobInput: { smilesList: [smilesInput], jobName: 'Single Prediction' } }
      });
      
      if (data?.submitPredictionJob) {
        const jobResponse = data.submitPredictionJob as JobStatus;
        setCurrentJobId(jobResponse.jobId);
        setJobStatus('pending');
        setProgress(10);
        
        // Start polling for results
        getPredictionResult({ variables: { jobId: jobResponse.jobId } });
      }
    } catch (err) {
      setError('Failed to submit prediction. Please try again.');
      setJobStatus('idle');
    }
  };

  const handleBatchPrediction = async (e: React.FormEvent) => {
    e.preventDefault();
    if (!batchInput.trim()) {
      setError('Please enter SMILES strings');
      return;
    }

    setError('');
    setResults([]);
    setProgress(0);
    
    // Parse SMILES strings (one per line or comma-separated)
    const smilesStrings = batchInput
      .split(/[\n,]/)
      .map(s => s.trim())
      .filter(s => s.length > 0);
      
    if (smilesStrings.length === 0) {
      setError('No valid SMILES strings found');
      return;
    }

    try {
      const { data } = await submitPredictionJobMutation({
        variables: { jobInput: { smilesList: smilesStrings, jobName: 'Batch Prediction' } }
      });
      
      if (data?.submitBatchPredictionJob) {
        const jobResponse = data.submitBatchPredictionJob as JobStatus;
        setCurrentJobId(jobResponse.jobId);
        setJobStatus('pending');
        setProgress(10);
        
        // Start polling for results
        getPredictionResult({ variables: { jobId: jobResponse.jobId } });
      }
    } catch (err) {
      setError('Failed to submit batch prediction. Please try again.');
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
        return 'Queued for processing...';
      case 'processing':
        return 'Running prediction model...';
      case 'completed':
        return 'Prediction completed!';
      case 'failed':
        return 'Prediction failed';
      default:
        return '';
    }
  };

  const isProcessing = jobStatus === 'pending' || jobStatus === 'processing';

  return (
    <div className="container mx-auto p-4 max-w-4xl">
      <Card>
        <CardHeader>
          <CardTitle className="flex items-center gap-2">
            Chemical Permeability Prediction
            {getStatusIcon()}
          </CardTitle>
          <CardDescription>
            Enter SMILES notation to predict compound permeability using machine learning
          </CardDescription>
        </CardHeader>
        <CardContent>
          <Tabs defaultValue="single" className="w-full">
            <TabsList className="grid w-full grid-cols-2">
              <TabsTrigger value="single">Single Prediction</TabsTrigger>
              <TabsTrigger value="batch">Batch Prediction</TabsTrigger>
            </TabsList>

            <TabsContent value="single" className="space-y-4">
              <form onSubmit={handleSinglePrediction} className="space-y-4">
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
                      Processing...
                    </>
                  ) : (
                    'Predict Permeability'
                  )}
                </Button>
              </form>
            </TabsContent>

            <TabsContent value="batch" className="space-y-4">
              <form onSubmit={handleBatchPrediction} className="space-y-4">
                <div>
                  <textarea
                    placeholder="Enter multiple SMILES strings (one per line or comma-separated)..."
                    value={batchInput}
                    onChange={(e) => setBatchInput(e.target.value)}
                    className="w-full min-h-[120px] p-3 border border-gray-300 rounded-md font-mono text-sm resize-vertical focus:ring-2 focus:ring-blue-500 focus:border-transparent"
                    disabled={isProcessing}
                  />
                </div>
                <Button type="submit" disabled={isProcessing || !batchInput.trim()} className="w-full">
                  {isProcessing ? (
                    <>
                      <Loader2 className="mr-2 h-4 w-4 animate-spin" />
                      Processing Batch...
                    </>
                  ) : (
                    'Predict Batch'
                  )}
                </Button>
              </form>
            </TabsContent>
          </Tabs>

          {/* Progress indicator */}
          {isProcessing && (
            <div className="mt-6 space-y-2">
              <div className="flex items-center justify-between text-sm">
                <span className="flex items-center gap-2">
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
        </CardContent>
      </Card>
    </div>
  );
};

export default PredictionForm;
export { PredictionForm };