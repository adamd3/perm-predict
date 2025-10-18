export default function HelpPage() {
  return (
    <main className="flex min-h-screen flex-col items-center justify-center p-24">
      <h1 className="text-4xl font-bold mb-6">Perm-Predict Help Center</h1>
      <div className="max-w-3xl text-lg text-justify space-y-6">
        <p>
          This guide will walk you through the main features of the application and how to use them effectively.
        </p>

        <h2 className="text-3xl font-semibold mt-8 mb-4">Understanding the Application Pages</h2>

        <h3 className="text-2xl font-medium mt-6 mb-2">Prediction Page</h3>
        <p className="text-left">
          On the Prediction page, users can submit chemical compounds for permeability prediction. Input your chemical compound(s) as SMILES strings. The application will process these inputs using a machine learning model and provide a prediction on their accumulation in bacterial cells. Results are displayed asynchronously, meaning you'll get a job ID and can check back for the completed prediction.
        </p>

        <h3 className="text-2xl font-medium mt-6 mb-2">Explore Page</h3>
        <p className="text-left">
          The Explore page is designed to help users discover similar chemical compounds that might exhibit higher accumulation or permeability. Here, you can input a compound, and the system will suggest potential modifications or alternative compounds that are predicted to be more permeant, aiding in the design and optimization of chemical structures.
        </p>

        <h3 className="text-2xl font-medium mt-6 mb-2">Create Page</h3>
        <p className="text-left">
          The Create page provides a graphical interface for designing chemical compounds. Users can draw or build molecular structures directly within the application and then export them, typically in SMILES format, for use in predictions or exploration. This feature allows for intuitive chemical design without needing external tools.
        </p>

        <h3 className="text-2xl font-medium mt-6 mb-2">About Page</h3>
        <p className="text-left">
          The About page offers a comprehensive overview of Perm-Predict, detailing its purpose, the underlying machine learning methodology, and a high-level explanation of how the asynchronous prediction system works. It's a great place to learn more about the science and technology behind the application.
        </p>
      </div>
    </main>
  );
}