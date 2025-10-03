'use client';

import { useSearchParams } from 'next/navigation';
import { useEffect, useState, Suspense } from 'react';
import { Button } from '@/components/ui/button';
import {
  Card,
  CardContent,
  CardDescription,
  CardHeader,
  CardTitle,
} from '@/components/ui/card';
import { Palette } from 'lucide-react';
import Link from 'next/link';
import PredictionForm from '@/components/PredictionForm';

function PredictionPageComponent() {
  const searchParams = useSearchParams();
  const [initialSmiles, setInitialSmiles] = useState<string>('');

  useEffect(() => {
    const smilesParam = searchParams.get('smiles');
    if (smilesParam) {
      setInitialSmiles(smilesParam);
    }
  }, [searchParams]);

  return (
    <main className="min-h-screen bg-gray-950 py-16">
      <div className="container mx-auto p-4 max-w-4xl">
        <Card className="mb-8 shadow-xl hover:shadow-2xl transition-shadow duration-300 border border-gray-700 bg-gray-800">
          <CardHeader className="text-center">
            <CardTitle className="text-6xl font-extrabold text-blue-300 mb-4">
              Perm-Predict
            </CardTitle>
            <CardDescription className="text-2xl text-gray-400 mb-8">
              Advanced machine learning-based prediction of chemical accumulation in bacteria
            </CardDescription>
          </CardHeader>
        </Card>
        
        <Card className="shadow-xl hover:shadow-2xl transition-shadow duration-300 border border-gray-700 bg-gray-800">
          <CardHeader>
            <CardTitle className="text-4xl font-bold text-blue-300">
              Predict Permeability
            </CardTitle>
            <CardDescription className="text-lg text-gray-400">
              Enter a SMILES string below to predict its permeance in bacterial cells.
            </CardDescription>
          </CardHeader>
          <CardContent>
            <PredictionForm initialSmiles={initialSmiles} />
          </CardContent>
        </Card>
      </div>
    </main>
  );
}

export default function PredictionPage() {
  return (
    <Suspense fallback={<div>Loading...</div>}>
      <PredictionPageComponent />
    </Suspense>
  );
}
