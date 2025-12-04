import React, { useRef, useEffect, useState, useCallback } from 'react';

interface KetcherIframeEditorProps {
  onSmilesChange: (smiles: string) => void;
}

const KETCHER_IFRAME_URL = '/ketcher-iframe/standalone/index.html';
const KETCHER_ORIGIN = window.location.origin; // Assuming same origin for simplicity, adjust if deployed cross-origin

const KetcherIframeEditor: React.FC<KetcherIframeEditorProps> = ({ onSmilesChange }) => {
  const iframeRef = useRef<HTMLIFrameElement>(null);
  const [ketcherReady, setKetcherReady] = useState(false);

  const getSmiles = useCallback(() => {
    if (iframeRef.current && ketcherReady) {
      iframeRef.current.contentWindow?.postMessage({ type: 'GET_SMILES' }, KETCHER_ORIGIN);
    }
  }, [ketcherReady]);

  useEffect(() => {
    const handleMessage = (event: MessageEvent) => {
      // Ensure the message is from the expected origin and is a trusted message
      if (event.origin !== KETCHER_ORIGIN) {
        console.warn('Received message from unknown origin:', event.origin);
        return;
      }

      if (event.data && typeof event.data === 'object') {
        switch (event.data.type) {
          case 'KETCHER_READY':
            setKetcherReady(true);
            console.log('Ketcher iframe is ready.');
            break;
          case 'SMILES_RESULT':
            onSmilesChange(event.data.smiles);
            break;
          case 'ERROR':
            console.error('Error from Ketcher iframe:', event.data.message);
            break;
          default:
            // console.log('Unhandled message type from Ketcher iframe:', event.data.type);
            break;
        }
      }
    };

    window.addEventListener('message', handleMessage);

    return () => {
      window.removeEventListener('message', handleMessage);
    };
  }, [onSmilesChange]);

  return (
    <div>
      <iframe
        ref={iframeRef}
        src={KETCHER_IFRAME_URL}
        width="100%"
        height="600px"
        style={{ border: '1px solid #ccc' }}
        title="Ketcher Editor"
      ></iframe>
      <button onClick={getSmiles} disabled={!ketcherReady} className="mt-2 px-4 py-2 bg-blue-500 text-white rounded disabled:opacity-50">
        Get SMILES
      </button>
    </div>
  );
};

export default KetcherIframeEditor;
