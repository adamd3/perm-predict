import React from 'react';
import {
  Table,
  TableBody,
  TableCell,
  TableHead,
  TableHeader,
  TableRow,
} from '@/components/ui/table';

interface PredictionResult {
  smiles: string;
  prediction: number;
  probability: number;
  error?: string;
}

interface PredictionResultsProps {
  results: PredictionResult[];
}

const PredictionResults = ({ results }: PredictionResultsProps) => {
  if (!results.length) return null;

  return (
    <div className="mt-8 overflow-x-auto">
      <Table>
        <TableHeader>
          <TableRow>
            <TableHead>SMILES</TableHead>
            <TableHead>Prediction</TableHead>
            <TableHead>Probability</TableHead>
          </TableRow>
        </TableHeader>
        <TableBody>
          {results.map((result, index) => (
            <TableRow key={index}>
              <TableCell className="font-mono">{result.smiles}</TableCell>
              <TableCell>{result.prediction}</TableCell>
              <TableCell>{(result.probability * 100).toFixed(2)}%</TableCell>
            </TableRow>
          ))}
        </TableBody>
      </Table>
    </div>
  );
};

export default PredictionResults;