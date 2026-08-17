export default function AboutPage() {
  return (
    <main className="flex min-h-screen flex-col items-center justify-center p-24">
      <h1 className="text-4xl font-bold mb-6">About Perm-Predict</h1>
      <div className="max-w-3xl text-lg text-justify space-y-4">
        <p>
          Perm-Predict is a web application designed for the machine learning-based prediction of chemical accumulation in bacteria. Its goal is to provide a robust tool for researchers and scientists to understand and predict the permeance of chemical compounds.
        </p>
        <p>
          Users can input a single chemical compound (as a SMILES string) or a list of compounds. The application then utilizes a sophisticated machine learning pipeline to predict how likely these compounds are to accumulate in bacterial cells.
        </p>
        <p>
          Beyond just prediction, Perm-Predict is being developed to suggest potential modifications or alternative compounds. By testing similar chemical structures and querying its classification model, the app aims to identify analogs that are more likely to be permeant, thereby aiding in chemical design and optimization.
        </p>
        <h2 className="text-3xl font-semibold mt-8 mb-4">How It Works</h2>
        <p>
          Perm-Predict employs a decoupled, asynchronous architecture to efficiently handle computationally intensive machine learning predictions.
        </p>
        <ul className="list-disc list-inside text-left mx-auto w-fit">
          <li><strong>Prediction Submission:</strong> A user submits a SMILES string via the application's interface.</li>
          <li><strong>Asynchronous Processing:</strong> The application sends a request to its backend, which immediately returns a unique Job ID. A dedicated worker then processes the machine learning prediction in the background.</li>
          <li><strong>Result Retrieval:</strong> The application periodically checks the status of the job using the Job ID and displays the results once available.</li>
        </ul>
      </div>
    </main>
  );
}