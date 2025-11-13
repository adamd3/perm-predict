'use client';

import React, { useState, useCallback } from 'react';
import dynamic from 'next/dynamic';
import { useRouter } from 'next/navigation'; // Import useRouter
import { Button } from '../../components/ui/button';
import { Card, CardContent, CardHeader, CardTitle } from '../../components/ui/card';

const DynamicKetcherIframeEditor = dynamic(() => import('../../components/KetcherIframeEditor'), { ssr: false });

export default function CreatePage() {
  const [smiles, setSmiles] = useState<string>('');
  const router = useRouter(); // Initialize useRouter

  const handleSmilesChange = useCallback((newSmiles: string) => {
    setSmiles(newSmiles);
  }, []);

  const handleSubmit = () => {
    if (smiles) {
      router.push(`/?smiles=${encodeURIComponent(smiles)}`); // Navigate to predict page with SMILES
    } else {
      alert('Please draw a molecule first.');
    }
  };

  return (
    <div className="container mx-auto p-4">
      <h1 className="text-2xl font-bold mb-4">Molecule Designer</h1>
      <Card className="mb-6">
        <CardHeader>
          <CardTitle>Draw Your Molecule</CardTitle>
        </CardHeader>
        <CardContent>
          <DynamicKetcherIframeEditor onSmilesChange={handleSmilesChange} />
        </CardContent>
      </Card>

      <Card className="mb-6">
        <CardHeader>
          <CardTitle>Generated SMILES</CardTitle>
        </CardHeader>
        <CardContent>
          <p className="break-all font-mono bg-gray-100 p-2 rounded text-gray-900">{smiles || 'Draw a molecule to see its SMILES string here.'}</p>
          <Button onClick={handleSubmit} className="mt-4" disabled={!smiles}>
            Use This Molecule
          </Button>
        </CardContent>
      </Card>
    </div>
  );
}
