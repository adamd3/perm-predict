'use client';

import { useState } from 'react';
import { useRouter } from 'next/navigation';
import ChemicalEditor from '@/components/ChemicalEditor';
import { Button } from '@/components/ui/button';
import { ArrowLeft, ArrowRight } from 'lucide-react';

export default function EditorPage() {
  const router = useRouter();
  const [selectedSmiles, setSelectedSmiles] = useState<string>('');

  const handleSmilesGenerated = (smiles: string) => {
    setSelectedSmiles(smiles);
  };

  const handleProceedToPrediction = () => {
    // Navigate to main page with the SMILES string
    const params = new URLSearchParams({ smiles: selectedSmiles });
    router.push(`/?${params}`);
  };

  return (
    <main className="min-h-screen bg-gray-50 py-8">
      <div className="container mx-auto p-4 max-w-6xl">
        {/* Navigation */}
        <div className="flex items-center justify-between mb-6">
          <Button
            variant="outline"
            onClick={() => router.push('/')}
            className="flex items-center gap-2"
          >
            <ArrowLeft className="w-4 h-4" />
            Back to Predictions
          </Button>
          
          {selectedSmiles && (
            <Button
              onClick={handleProceedToPrediction}
              className="flex items-center gap-2"
            >
              Predict This Compound
              <ArrowRight className="w-4 h-4" />
            </Button>
          )}
        </div>

        {/* Chemical Editor */}
        <ChemicalEditor onSmilesGenerated={handleSmilesGenerated} />

        {/* Selected SMILES Display */}
        {selectedSmiles && (
          <div className="mt-6 p-4 bg-white border rounded-lg">
            <h3 className="font-medium mb-2">Selected Structure:</h3>
            <div className="font-mono text-sm bg-gray-100 p-3 rounded">
              {selectedSmiles}
            </div>
            <div className="mt-3 flex gap-2">
              <Button onClick={handleProceedToPrediction}>
                Predict Permeability
              </Button>
              <Button
                variant="outline"
                onClick={() => navigator.clipboard.writeText(selectedSmiles)}
              >
                Copy SMILES
              </Button>
            </div>
          </div>
        )}
      </div>
    </main>
  );
}