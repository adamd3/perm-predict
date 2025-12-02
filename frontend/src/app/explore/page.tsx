'use client';

import React, { useState, useEffect, Suspense } from 'react'; // Import Suspense
import { useSearchParams } from 'next/navigation';
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from '@/components/ui/card';
import ExplorationForm from '@/components/ExplorationForm';

function ExplorePageComponent() { // Renamed
  const searchParams = useSearchParams();
  const [hasResults, setHasResults] = useState<boolean>(false);

  const resetKey = searchParams.get('_t') || 'default';

  useEffect(() => {
    setHasResults(false);
  }, [resetKey, searchParams]);

  return (
    <div className="flex flex-col items-center justify-center min-h-screen py-2">
      <main className="flex flex-col items-center justify-center w-full flex-1 px-4 sm:px-20 text-center">
        <Card className={`w-full mx-auto ${hasResults ? 'max-w-full' : 'max-w-2xl'}`}>
          <CardHeader>
            <CardTitle className="text-3xl font-bold">Explore compound permeability</CardTitle>
            <CardDescription className="mt-2 text-lg">
              Enter a SMILES string to predict permeability of related compounds
            </CardDescription>
          </CardHeader>
          <CardContent>
            <ExplorationForm key={resetKey} onResultsLoaded={setHasResults} />
          </CardContent>
        </Card>
      </main>
    </div>
  );
}

export default function ExplorePage() { // New wrapper function
  return (
    <Suspense fallback={<div>Loading Explore Page...</div>}>
      <ExplorePageComponent />
    </Suspense>
  );
}

