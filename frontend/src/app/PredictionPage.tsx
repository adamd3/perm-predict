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
    <main className="min-h-screen py-16">
      <div className="container mx-auto p-4 max-w-6xl">
        <Card className="mb-8 shadow-xl hover:shadow-2xl transition-shadow duration-300 max-w-2xl mx-auto">
          <CardHeader className="text-center">
            <CardTitle className="text-6xl font-extrabold mb-4 text-foreground dark:text-foreground">
              Perm-Predict
            </CardTitle>
            <CardDescription className="text-2xl mb-8 text-muted-foreground dark:text-muted-foreground">
              AI-based prediction of chemical accumulation in bacteria
            </CardDescription>
          </CardHeader>
        </Card>
        
        <Card className="shadow-xl hover:shadow-2xl transition-shadow duration-300 max-w-2xl mx-auto">
          <CardHeader>
            <CardTitle className="text-4xl font-bold text-foreground dark:text-foreground">
              Predict Permeability
            </CardTitle>
            <CardDescription className="text-lg text-muted-foreground dark:text-muted-foreground">
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
