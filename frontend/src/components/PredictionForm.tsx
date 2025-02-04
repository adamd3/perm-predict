import React, { useState } from 'react';
import { Alert, AlertDescription } from '@/components/ui/alert';
import { Button } from '@/components/ui/button';
import { Input } from '@/components/ui/input';
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from '@/components/ui/card';
import { Tabs, TabsContent, TabsList, TabsTrigger } from '@/components/ui/tabs';
import PredictionResults from './PredictionResults';

const PredictionForm = () => {
  const [smilesInput, setSmilesInput] = useState('');
  const [file, setFile] = useState(null);
  const [results, setResults] = useState([]);
  const [error, setError] = useState('');
  const [loading, setLoading] = useState(false);

  const handleSinglePrediction = async (e) => {
    e.preventDefault();
    setError('');
    setLoading(true);

    try {
      const response = await fetch('/api/predict/single', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ smiles: smilesInput })
      });

      const data = await response.json();
      if (data.error) {
        setError(data.error);
      } else {
        setResults([data]);
      }
    } catch (err) {
      setError('Failed to get prediction. Please try again.');
    } finally {
      setLoading(false);
    }
  };

  const handleFileUpload = async (e) => {
    e.preventDefault();
    if (!file) {
      setError('Please select a file');
      return;
    }

    setError('');
    setLoading(true);

    const formData = new FormData();
    formData.append('file', file);

    try {
      const response = await fetch('/api/predict/batch', {
        method: 'POST',
        body: formData
      });

      const data = await response.json();
      if (response.ok) {
        setResults(data);
      } else {
        setError(data.detail || 'Failed to process file');
      }
    } catch (err) {
      setError('Failed to upload file. Please try again.');
    } finally {
      setLoading(false);
    }
  };

  return (
    <div className="container mx-auto p-4 max-w-4xl">
      <Card>
        <CardHeader>
          <CardTitle>Chemical Permeability Prediction</CardTitle>
          <CardDescription>
            Enter SMILES notation or upload a CSV file to predict compound permeability
          </CardDescription>
        </CardHeader>
        <CardContent>
          <Tabs defaultValue="single" className="w-full">
            <TabsList>
              <TabsTrigger value="single">Single Prediction</TabsTrigger>
              <TabsTrigger value="batch">Batch Prediction</TabsTrigger>
            </TabsList>

            <TabsContent value="single">
              <form onSubmit={handleSinglePrediction} className="space-y-4">
                <div>
                  <Input
                    placeholder="Enter SMILES string..."
                    value={smilesInput}
                    onChange={(e) => setSmilesInput(e.target.value)}
                    className="w-full"
                  />
                </div>
                <Button type="submit" disabled={loading || !smilesInput}>
                  {loading ? 'Processing...' : 'Predict'}
                </Button>
              </form>
            </TabsContent>

            <TabsContent value="batch">
              <form onSubmit={handleFileUpload} className="space-y-4">
                <div>
                  <Input
                    type="file"
                    accept=".csv"
                    onChange={(e) => setFile(e.target.files[0])}
                    className="w-full"
                  />
                </div>
                <Button type="submit" disabled={loading || !file}>
                  {loading ? 'Processing...' : 'Upload and Predict'}
                </Button>
              </form>
            </TabsContent>
          </Tabs>

          {error && (
            <Alert variant="destructive" className="mt-4">
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