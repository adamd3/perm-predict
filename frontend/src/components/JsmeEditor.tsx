"use client";

import React, { useEffect, useRef, useState, useCallback } from "react";

interface JsmeEditorProps {
  initialSmiles?: string;
  onSmilesChange?: (smiles: string) => void;
}

declare global {
  interface Window {
    JSApplet: {
      JSME: new (
        id: string,
        width: string,
        height: string
      ) => JSMEApplet;
    };
    jsmeOnLoad?: () => void;
  }
}

interface JSMEApplet {
  g: {
    readGenericMolecularInput: (smiles: string) => void;
    smiles: () => string;
    options: (options: string) => void;
    setSmiles: (smiles: string) => void;
    setCallBack: (eventName: string, callback: (event: any) => void) => void;
  };
}

const JsmeEditor: React.FC<JsmeEditorProps> = ({
  initialSmiles = "",
  onSmilesChange,
}) => {
  const jsmeRef = useRef<HTMLDivElement>(null);
  const jsmeApplet = useRef<JSMEApplet | null>(null);
  const [isJsmeLoaded, setIsJsmeLoaded] = useState(false);

  // Function to load the JSME script dynamically
  const loadJsmeScript = useCallback(() => {
    if (document.getElementById("jsme-script")) {
      setIsJsmeLoaded(true);
      return;
    }

    const script = document.createElement("script");
    script.id = "jsme-script";
    script.src = "/jsme/jsme.nocache.js"; // Path to JSME script in public folder
    script.async = true;
    script.onload = () => {
      setIsJsmeLoaded(true);
      console.log("JSME script loaded.");
    };
    script.onerror = () => {
      console.error("Failed to load JSME script.");
    };
    document.head.appendChild(script);
  }, []);

  useEffect(() => {
    loadJsmeScript();

    // Define the global jsmeOnLoad function
    window.jsmeOnLoad = () => {
      if (jsmeRef.current && window.JSApplet && !jsmeApplet.current) {
        console.log("Initializing JSME applet...");
        console.log("window.JSApplet:", window.JSApplet); // Added for debugging
                console.log("jsmeRef.current before JSME init:", jsmeRef.current); // Debugging line
                // Clear the contents of the container before initializing JSME
                if (jsmeRef.current) {
                  jsmeRef.current.innerHTML = '';
                }
                jsmeApplet.current = new window.JSApplet.JSME(
                  jsmeRef.current.id,
                  "100%",
                  "500px"
                );
                console.log("JSME Applet instance:", jsmeApplet.current); // Debugging line
        
                jsmeApplet.current.g.options("query,depict,highlight,structure,autoez,nocanonize");                        
                                
                        
                                        if (initialSmiles) {
                        
                                          jsmeApplet.current.g.readGenericMolecularInput(initialSmiles);
                        
                                        }
                        
                                
                        
                                        // Define the callback function for structural changes
                        
                                        const handleJsmeStructuralChange = (event: any) => {
                        
                                          if (onSmilesChange) {
                        
                                            const currentSmiles = event.src.smiles();
                        
                                            onSmilesChange(currentSmiles);
                        
                                          }
                        
                                        };
                        
                                
                        
                                        // Register the callback with JSME
                        
                                        jsmeApplet.current.g.setCallBack("AfterStructureModified", handleJsmeStructuralChange);

        

              }

            };

        

            return () => {
      // Cleanup: remove script and global function if component unmounts
      const script = document.getElementById("jsme-script");
      if (script) {
        document.head.removeChild(script);
      }
      delete window.jsmeOnLoad;
    };
  }, [loadJsmeScript, initialSmiles, onSmilesChange]);

  return (
    <div
      id="jsme_container"
      ref={jsmeRef}
      className="w-full h-[500px] border border-gray-300 rounded-md"
    >
      {!isJsmeLoaded && (
        <p className="text-center text-gray-500 mt-4">Loading chemical editor...</p>
      )}
    </div>
  );
};

export default JsmeEditor;

