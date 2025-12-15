'use client';

import Image from 'next/image';
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
  const [hasResults, setHasResults] = useState<boolean>(false);

  const resetKey = searchParams.get('_t') || 'default'; // Use _t as key to force remount

  useEffect(() => {
    const smilesParam = searchParams.get('smiles');
    setInitialSmiles(smilesParam || '');
    // Reset hasResults when the page effectively reloads (due to _t or smiles change)
    setHasResults(false);
  }, [searchParams, resetKey]); // Depend on searchParams AND resetKey to ensure useEffect fires on all relevant changes

  return (
    <main className="min-h-screen py-16">
      <div className="container mx-auto p-4 max-w-full">
        <Card className="mb-8 shadow-xl hover:shadow-2xl transition-shadow duration-300 max-w-2xl mx-auto">
          <CardHeader className="text-center">

            <Image
              src="/images/perm_predict_logo.png"
              alt="Perm-Predict Logo"
              width={300}
              height={300}
              className="mx-auto block mb-4"
            />
            <CardDescription className="text-2xl mb-8 text-muted-foreground dark:text-muted-foreground">
              AI-based prediction of chemical accumulation in bacteria
            </CardDescription>
          </CardHeader>
        </Card>
        
        <Card className={`shadow-xl hover:shadow-2xl transition-shadow duration-300 ${hasResults ? 'w-full' : 'max-w-2xl'} mx-auto`}>
          <CardHeader>
            <CardTitle className="text-4xl font-bold text-foreground dark:text-foreground">
              Predict Permeability
            </CardTitle>
            <CardDescription className="text-lg text-muted-foreground dark:text-muted-foreground">
              Enter a SMILES string below to predict its permeance in bacterial cells.
            </CardDescription>
          </CardHeader>
          <CardContent>
            {/* Use resetKey as the key to force PredictionForm to remount and reset its state */}
            <PredictionForm key={resetKey} initialSmiles={initialSmiles} onResultsLoaded={setHasResults} />
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
