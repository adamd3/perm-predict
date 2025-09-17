'use client';

import { useSearchParams } from 'next/navigation';
import { useEffect, useState, Suspense } from 'react';
import { Button } from '@/components/ui/button';
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
    <main className="min-h-screen bg-gray-50 py-8">
      <div className="container mx-auto p-4 max-w-4xl">
        {/* Header with navigation */}
        <div className="text-center mb-8">
          <h1 className="text-3xl font-bold text-gray-900 mb-2">
            Perm-Predict
          </h1>
          <p className="text-gray-600 mb-6">
            Advanced machine learning-based prediction of chemical accumulation in bacteria
          </p>
          
          <Link href="/editor">
            <Button variant="outline" className="flex items-center gap-2 mx-auto">
              <Palette className="w-4 h-4" />
              Open Chemical Structure Editor
            </Button>
          </Link>
        </div>
        
        <PredictionForm initialSmiles={initialSmiles} />
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
