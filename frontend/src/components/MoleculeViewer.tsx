import React, { useEffect, useRef, useState } from 'react';
import initRDKit from '@rdkit/rdkit';
import { Dialog, DialogContent, DialogTrigger } from '@/components/ui/dialog';

interface MoleculeViewerProps {
  smiles: string;
  width?: number;
  height?: number;
}

const MoleculeViewer: React.FC<MoleculeViewerProps> = ({ smiles, width = 300, height = 200 }) => {
  const molRef = useRef<HTMLCanvasElement>(null); // Change to HTMLCanvasElement
  const [RDKit, setRDKit] = useState<any>(null);

  useEffect(() => {
    console.log("MoleculeViewer: Attempting to initialize RDKit...");
    (initRDKit as any)().then((rdkit: any) => {
      console.log("MoleculeViewer: RDKit initialized successfully.");
      console.log("MoleculeViewer: RDKit object:", rdkit);
      setRDKit(rdkit);
    }).catch((error: any) => {
      console.error("MoleculeViewer: Error initializing RDKit:", error);
    });
  }, []);

  useEffect(() => {
    console.log("MoleculeViewer: useEffect for drawing triggered.");
    console.log("MoleculeViewer: RDKit state:", RDKit ? "initialized" : "not initialized");
    console.log("MoleculeViewer: molRef.current:", molRef.current ? "available" : "not available");
    console.log("MoleculeViewer: smiles prop:", smiles);

    if (RDKit && molRef.current && smiles) {
      try {
        console.log("MoleculeViewer: Attempting to get molecule from SMILES:", smiles);
        const mol = RDKit.get_mol(smiles);
        console.log("MoleculeViewer: mol object after get_mol:", mol);
        if (mol) {
          console.log("MoleculeViewer: Molecule obtained, attempting to draw.");
          console.log("MoleculeViewer: typeof mol.draw_to_canvas_with_highlights:", typeof mol.draw_to_canvas_with_highlights);
          
          // Clear the canvas before drawing a new molecule
          const ctx = molRef.current.getContext('2d');
          if (ctx) {
            ctx.clearRect(0, 0, molRef.current.width, molRef.current.height);
          }

          mol.draw_to_canvas_with_highlights(
            molRef.current, // Pass the canvas element
            `{"width": ${width}, "height": ${height}}`
          );
          mol.delete(); // Important for memory management
          console.log("MoleculeViewer: Molecule drawn successfully.");
        } else {
          console.log("MoleculeViewer: Could not generate molecule from SMILES.");
          molRef.current.innerHTML = '<p>Could not generate molecule from SMILES.</p>';
        }
      } catch (e) {
        console.error("MoleculeViewer: Error rendering molecule:", e);
        molRef.current.innerHTML = '<p>Invalid SMILES string.</p>';
      }
    } else if (!smiles) {
      console.log("MoleculeViewer: No SMILES string provided, skipping drawing.");
    }
  }, [RDKit, smiles, width, height]);

  return (
    <Dialog>
      <DialogTrigger asChild>
        <div className="cursor-pointer">
          <canvas ref={molRef} width={width} height={height} style={{ overflow: 'hidden' }}>
            {!RDKit && <p>Loading molecule viewer...</p>}
          </canvas>
        </div>
      </DialogTrigger>
      <DialogContent className="max-w-fit p-0">
        <div className="p-4">
          <h4 className="text-lg font-semibold mb-2 text-black">Molecule Structure</h4>
          <MoleculeViewer smiles={smiles} width={600} height={400} /> {/* Larger version in modal */}
          <p className="text-sm text-gray-500 mt-2 break-all">SMILES: {smiles}</p>
        </div>
      </DialogContent>
    </Dialog>
  );
};

export default MoleculeViewer;