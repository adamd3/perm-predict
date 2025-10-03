import type { Metadata } from "next";
import { Geist, Geist_Mono } from "next/font/google";
import "./globals.css";
import { GraphQLProvider } from "@/lib/apollo-provider";
import Link from "next/link"; // Import Link

const geistSans = Geist({
  variable: "--font-geist-sans",
  subsets: ["latin"],
});

const geistMono = Geist_Mono({
  variable: "--font-geist-mono",
  subsets: ["latin"],
});

export const metadata: Metadata = {
  title: "Perm-Predict - Chemical Permeability Prediction",
  description: "Advanced machine learning-based prediction of chemical accumulation in bacteria",
};

export default function RootLayout({
  children,
}: Readonly<{
  children: React.ReactNode;
}>) {
  return (
    <html lang="en">
      <body
        className={`${geistSans.variable} ${geistMono.variable} antialiased flex flex-col min-h-screen bg-gray-950 text-gray-200`}
      >
        <GraphQLProvider>
          <header className="bg-gray-900 text-gray-200 p-4 shadow-md">
            <nav className="container mx-auto flex justify-between items-center px-8">
              <Link href="/" className="text-2xl font-bold">
                Perm-Predict
              </Link>
              <ul className="flex space-x-8">
                <li>
                  <Link href="/" className="hover:text-gray-50">
                    Predict
                  </Link>
                </li>
                <li>
                  <Link href="/create" className="hover:text-gray-50">
                    Create
                  </Link>
                </li>
                <li>
                  <Link href="/explore" className="hover:text-gray-50">
                    Explore
                  </Link>
                </li>
                <li>
                  <Link href="/help" className="hover:text-gray-50">
                    Help
                  </Link>
                </li>
                <li>
                  <Link href="/about" className="hover:text-gray-50">
                    About
                  </Link>
                </li>
              </ul>
            </nav>
          </header>

          <main className="flex-grow">
            {children}
          </main>

          <footer className="bg-gray-900 text-gray-200 p-6 mt-8">
            <div className="container mx-auto flex flex-col md:flex-row justify-between items-center">
              <div className="text-sm">
                &copy; {new Date().getFullYear()} Perm-Predict. All rights reserved.
              </div>
              <div className="flex space-x-4 mt-4 md:mt-0">
                <Link href="/contact" className="hover:text-gray-50">
                  Contact Us
                </Link>
                {/* Placeholder images */}
                <img src="/placeholder-1.svg" alt="Placeholder 1" className="h-6 w-6" />
                <img src="/placeholder-2.svg" alt="Placeholder 2" className="h-6 w-6" />
              </div>
            </div>
          </footer>
        </GraphQLProvider>
      </body>
    </html>
  );
}
