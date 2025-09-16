import React from 'react';
import {
  Table,
  TableBody,
  TableCell,
  TableHead,
  TableHeader,
  TableRow,
} from '@/components/ui/table';
import { Badge } from '@/components/ui/badge';
import { Progress } from '@/components/ui/progress';
import { CheckCircle, XCircle, AlertCircle } from 'lucide-react';

import type { PredictionResultsProps } from '@/lib/types'

const PredictionResults = ({ results }: PredictionResultsProps) => {
  if (!results.length) return null;

  const getPredictionBadge = (classifierPrediction: number) => {
    const isPermeant = classifierPrediction === 1;
    const variant = isPermeant ? 'default' : 'secondary';
    const icon = isPermeant ? <CheckCircle className="w-3 h-3" /> : <XCircle className="w-3 h-3" />;
    
    return (
      <Badge variant={variant} className={`flex items-center gap-1 ${isPermeant ? 'bg-green-100 text-green-800 border-green-300' : 'bg-red-100 text-red-800 border-red-300'}`}>
        {icon}
        {isPermeant ? 'Permeant' : 'Impermeant'}
      </Badge>
    );
  };

  return (
    <div className="mt-8 space-y-4">
      <div className="flex items-center justify-between">
        <h3 className="text-lg font-semibold">Prediction Results</h3>
        <div className="text-sm text-gray-600">
          {results.length} compound{results.length > 1 ? 's' : ''} analyzed
        </div>
      </div>
      
      <div className="overflow-x-auto border rounded-lg">
        <Table>
          <TableHeader>
            <TableRow className="bg-gray-50">
              <TableHead className="w-[300px]">SMILES</TableHead>
              <TableHead className="w-[150px]">Prediction</TableHead>
              <TableHead className="w-[120px]">Permeant Probability</TableHead>
              <TableHead className="w-[200px]">Confidence</TableHead>
              <TableHead className="w-[120px]">Error</TableHead>
            </TableRow>
          </TableHeader>
          <TableBody>
            {results.map((result, index) => (
              <TableRow key={index} className="hover:bg-gray-50">
                <TableCell className="font-mono text-sm break-all max-w-[300px]">
                  {result.smiles}
                </TableCell>
                <TableCell>
                  {result.error ? (
                    <Badge variant="destructive" className="flex items-center gap-1">
                      <AlertCircle className="w-3 h-3" /> Error
                    </Badge>
                  ) : (
                    getPredictionBadge(result.classifierPrediction)
                  )}
                </TableCell>
                <TableCell className="font-semibold">
                  {result.error ? '—' : `${(result.classProbabilities[1] * 100).toFixed(1)}%`}
                </TableCell>
                <TableCell>
                  {result.error ? '—' : (
                    <div className="space-y-1">
                      <div className="flex justify-between text-sm">
                        <span>Confidence</span>
                        <span className="font-medium">{(result.confidence * 100).toFixed(1)}%</span>
                      </div>
                      <Progress 
                        value={result.confidence * 100} 
                        className="h-2"
                      />
                    </div>
                  )}
                </TableCell>
                <TableCell className="text-sm text-red-600">
                  {result.error || '—'}
                </TableCell>
              </TableRow>
            ))}
          </TableBody>
        </Table>
      </div>
      
      {/* Summary Statistics */}
      {results.length > 1 && (
        <div className="grid grid-cols-2 md:grid-cols-4 gap-4 p-4 bg-gray-50 rounded-lg">
          <div className="text-center">
            <div className="text-2xl font-bold text-green-600">
              {results.filter(r => r.classifierPrediction === 1 && !r.error).length}
            </div>
            <div className="text-sm text-gray-600">Permeant</div>
          </div>
          <div className="text-center">
            <div className="text-2xl font-bold text-red-600">
              {results.filter(r => r.classifierPrediction === 0 && !r.error).length}
            </div>
            <div className="text-sm text-gray-600">Impermeant</div>
          </div>
          <div className="text-center">
            <div className="text-2xl font-bold text-blue-600">
              {(results.filter(r => !r.error).reduce((sum, r) => sum + r.confidence, 0) / results.filter(r => !r.error).length * 100).toFixed(1)}%
            </div>
            <div className="text-sm text-gray-600">Avg Confidence</div>
          </div>
          <div className="text-center">
            <div className="text-2xl font-bold text-gray-600">
              {results.filter(r => r.error).length}
            </div>
            <div className="text-sm text-gray-600">Failed</div>
          </div>
        </div>
      )}
    </div>
  );
};

export default PredictionResults;