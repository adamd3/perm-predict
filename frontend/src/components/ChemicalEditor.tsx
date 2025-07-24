import React, { useState } from 'react';
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from '@/components/ui/card';
import { Button } from '@/components/ui/button';
import { Input } from '@/components/ui/input';
import { Tabs, TabsContent, TabsList, TabsTrigger } from '@/components/ui/tabs';
import { Alert, AlertDescription } from '@/components/ui/alert';
import { Copy, Download, Upload, Palette, Info } from 'lucide-react';

interface ChemicalEditorProps {
  onSmilesGenerated: (smiles: string) => void;
}

const ChemicalEditor = ({ onSmilesGenerated }: ChemicalEditorProps) => {
  const [smilesInput, setSmilesInput] = useState('');
  const [showEditor, setShowEditor] = useState(false);

  // Sample common chemical structures for quick testing
  const sampleStructures = [
    { name: 'Ethanol', smiles: 'CCO', description: 'Simple alcohol' },
    { name: 'Caffeine', smiles: 'CN1C=NC2=C1C(=O)N(C(=O)N2C)C', description: 'Stimulant alkaloid' },
    { name: 'Aspirin', smiles: 'CC(=O)OC1=CC=CC=C1C(=O)O', description: 'Pain reliever' },
    { name: 'Glucose', smiles: 'C([C@@H]1[C@H]([C@@H]([C@H]([C@H](O1)O)O)O)O)O', description: 'Simple sugar' },
    { name: 'Benzene', smiles: 'c1ccccc1', description: 'Aromatic hydrocarbon' },
    { name: 'Penicillin G', smiles: 'CC1([C@@H](N2[C@H](S1)[C@@H](C2=O)NC(=O)CC3=CC=CC=C3)C(=O)O)C', description: 'Antibiotic' },
  ];

  const copyToClipboard = (smiles: string) => {
    navigator.clipboard.writeText(smiles);
  };

  const handleUseSample = (smiles: string) => {
    setSmilesInput(smiles);
    onSmilesGenerated(smiles);
  };

  const handleManualInput = () => {
    if (smilesInput.trim()) {
      onSmilesGenerated(smilesInput.trim());
    }
  };

  return (
    <Card className="w-full">
      <CardHeader>
        <CardTitle className="flex items-center gap-2">
          <Palette className="w-5 h-5" />
          Chemical Structure Editor
        </CardTitle>
        <CardDescription>
          Draw or input chemical structures to generate SMILES notation
        </CardDescription>
      </CardHeader>
      <CardContent>
        <Tabs defaultValue="samples" className="w-full">
          <TabsList className="grid w-full grid-cols-3">
            <TabsTrigger value="samples">Sample Structures</TabsTrigger>
            <TabsTrigger value="manual">Manual SMILES</TabsTrigger>
            <TabsTrigger value="editor" disabled>
              Draw Structure
              <span className="ml-1 text-xs text-gray-500">(Coming Soon)</span>
            </TabsTrigger>
          </TabsList>

          <TabsContent value="samples" className="space-y-4">
            <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
              {sampleStructures.map((structure, index) => (
                <div
                  key={index}
                  className="p-4 border rounded-lg hover:bg-gray-50 transition-colors"
                >
                  <div className="flex justify-between items-start mb-2">
                    <div>
                      <h4 className="font-medium">{structure.name}</h4>
                      <p className="text-sm text-gray-600">{structure.description}</p>
                    </div>
                    <Button
                      size="sm"
                      variant="outline"
                      onClick={() => handleUseSample(structure.smiles)}
                    >
                      Use This
                    </Button>
                  </div>
                  <div className="font-mono text-xs bg-gray-100 p-2 rounded flex justify-between items-center">
                    <span className="truncate flex-1 mr-2">{structure.smiles}</span>
                    <Button
                      size="sm"
                      variant="ghost"
                      onClick={() => copyToClipboard(structure.smiles)}
                      className="p-1 h-auto"
                    >
                      <Copy className="w-3 h-3" />
                    </Button>
                  </div>
                </div>
              ))}
            </div>
          </TabsContent>

          <TabsContent value="manual" className="space-y-4">
            <div className="space-y-4">
              <div>
                <Input
                  placeholder="Enter SMILES string (e.g., CCO for ethanol)..."
                  value={smilesInput}
                  onChange={(e) => setSmilesInput(e.target.value)}
                  className="font-mono"
                  onKeyPress={(e) => e.key === 'Enter' && handleManualInput()}
                />
              </div>
              <div className="flex gap-2">
                <Button onClick={handleManualInput} disabled={!smilesInput.trim()}>
                  Use This SMILES
                </Button>
                <Button
                  variant="outline"
                  onClick={() => copyToClipboard(smilesInput)}
                  disabled={!smilesInput.trim()}
                >
                  <Copy className="w-4 h-4 mr-2" />
                  Copy
                </Button>
              </div>
            </div>
          </TabsContent>

          <TabsContent value="editor" className="space-y-4">
            <Alert>
              <Info className="h-4 w-4" />
              <AlertDescription>
                <div className="space-y-2">
                  <p className="font-medium">Interactive Chemical Editor Coming Soon!</p>
                  <p className="text-sm">
                    We're working on integrating Ketcher, a professional chemical structure editor.
                    For now, you can use the sample structures or input SMILES notation manually.
                  </p>
                  <div className="text-xs text-gray-600 space-y-1">
                    <p>• Draw chemical structures with mouse/touch</p>
                    <p>• Automatic SMILES generation</p>
                    <p>• Structure validation and optimization</p>
                    <p>• Import/export various chemical formats</p>
                  </div>
                </div>
              </AlertDescription>
            </Alert>
            
            {/* Placeholder for future Ketcher integration */}
            <div className="h-64 border-2 border-dashed border-gray-300 rounded-lg flex items-center justify-center bg-gray-50">
              <div className="text-center text-gray-500">
                <Palette className="w-12 h-12 mx-auto mb-2 opacity-50" />
                <p className="font-medium">Chemical Structure Editor</p>
                <p className="text-sm">Will be integrated here</p>
              </div>
            </div>
          </TabsContent>
        </Tabs>
      </CardContent>
    </Card>
  );
};

export default ChemicalEditor;