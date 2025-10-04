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
import dynamic from 'next/dynamic';

const MoleculeViewer = dynamic(() => import('./MoleculeViewer'), {
  ssr: false,
});

import { Bar } from 'react-chartjs-2';
import { Chart as ChartJS, CategoryScale, LinearScale, BarElement, Title, Tooltip, Legend } from 'chart.js';

ChartJS.register(CategoryScale, LinearScale, BarElement, Title, Tooltip, Legend);

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

  const getChartData = (featuresSummary: { [key: string]: number }) => {
    const labels = Object.keys(featuresSummary);
    const data = Object.values(featuresSummary);

    return {
      labels,
      datasets: [
        {
          label: 'Feature Value',
          data,
          backgroundColor: 'rgba(75, 192, 192, 0.6)',
          borderColor: 'rgba(75, 192, 192, 1)',
          borderWidth: 1,
        },
      ],
    };
  };

  const chartOptions = {
    responsive: true,
    maintainAspectRatio: false,
    plugins: {
      legend: { display: false },
      title: { display: false },
    },
    scales: {
      x: { ticks: { font: { size: 8 } } },
      y: { ticks: { font: { size: 8 } } },
    },
  };

  return (
    <div className="mt-8 space-y-4">
      <div className="flex items-center justify-between">
        <h3 className="text-lg font-semibold text-white">Prediction Results</h3>
        <div className="text-sm text-white">
          {results.length} compound{results.length > 1 ? 's' : ''} analyzed
        </div>
      </div>
      
      <div className="overflow-x-auto border rounded-lg">
        <Table>
          <TableHeader>
            <TableRow className="bg-gray-50">
              <TableHead className="w-[300px]">SMILES / Structure</TableHead>
              <TableHead className="w-[150px]">Prediction</TableHead>
              <TableHead className="w-[120px]">Permeant Score</TableHead>
              <TableHead className="w-[200px]">Confidence</TableHead>
              <TableHead className="w-[150px]">Features Summary</TableHead>
              <TableHead className="w-[120px]">Status/Error</TableHead>
            </TableRow>
          </TableHeader>
          <TableBody>
            {results.map((result, index) => (
              <TableRow key={index} className="bg-white hover:bg-gray-50">
                <TableCell className="font-mono text-sm break-all max-w-[300px]">
                  {result.smiles}
                  {!result.error && result.smiles && (
                    <div className="mt-2">
                      <MoleculeViewer smiles={result.smiles} width={150} height={100} />
                    </div>
                  )}
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
                  {result.error ? '—' : `${(result.permeantProbability * 100).toFixed(1)}%`}
                </TableCell>
                <TableCell>
                  {result.error ? '—' : (
                    <div className="space-y-1">
                      <div className="flex justify-between text-sm">
                        <span>Overall Confidence</span>
                        <span className="font-medium">{(result.confidence * 100).toFixed(1)}%</span>
                      </div>
                      <Progress 
                        value={result.confidence * 100} 
                        className="h-2"
                      />
                    </div>
                  )}
                </TableCell>
                <TableCell className="max-w-[150px]">
                  {!result.error && result.featuresSummary && Object.keys(result.featuresSummary).length > 0 ? (
                    <div style={{ width: '100%', height: '100px' }}>
                      <Bar data={getChartData(result.featuresSummary)} options={chartOptions} />
                    </div>
                  ) : (
                    <span className="text-gray-500 text-xs">No summary features</span>
                  )}
                </TableCell>
                <TableCell className="text-sm text-red-600">
                  {result.error || 'Success'}
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