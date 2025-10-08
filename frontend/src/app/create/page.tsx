'use client';

import React, { useState, useCallback } from 'react';
import dynamic from 'next/dynamic';
import { Button } from '../../components/ui/button';
import { Card, CardContent, CardHeader, CardTitle } from '../../components/ui/card';

const DynamicJsmeEditor = dynamic(() => import('../../components/JsmeEditor'), {
  ssr: false,
  loading: () => <p>Loading editor...</p>,
});

export default function CreatePage() {
  const [smiles, setSmiles] = useState<string>('');

  const handleSmilesChange = useCallback((newSmiles: string) => {
    setSmiles(newSmiles);
  }, []);

  const handleSubmit = () => {
    if (smiles) {
      alert(`Submitting SMILES: ${smiles}`);
      // In a real application, you would send this SMILES to your backend
    } else {
      alert('Please draw a molecule first.');
    }
  };

  return (
    <div className="container mx-auto p-4">
      <h1 className="text-2xl font-bold mb-4">Create New Molecule</h1>
      <Card className="mb-6">
        <CardHeader>
          <CardTitle>Draw Your Molecule</CardTitle>
        </CardHeader>
        <CardContent>
          <DynamicJsmeEditor onSmilesChange={handleSmilesChange} />
        </CardContent>
      </Card>

      <Card className="mb-6">
        <CardHeader>
          <CardTitle>Generated SMILES</CardTitle>
        </CardHeader>
        <CardContent>
          <p className="break-all font-mono bg-gray-100 p-2 rounded">{smiles || 'Draw a molecule to see its SMILES string here.'}</p>
          <Button onClick={handleSubmit} className="mt-4" disabled={!smiles}>
            Use This Molecule
          </Button>
        </CardContent>
      </Card>
    </div>
  );
}