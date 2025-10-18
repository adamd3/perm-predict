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
  description: "AI-based prediction of chemical accumulation in bacteria",
};

import { ThemeProvider } from '@/components/theme-provider';
import { ThemeToggle } from '@/components/ThemeToggle';

export default function RootLayout({
  children,
}: Readonly<{
  children: React.ReactNode;
}>) {
  return (
    <html lang="en" suppressHydrationWarning>
      <body
        className={`${geistSans.variable} ${geistMono.variable} antialiased flex flex-col min-h-screen bg-gray-50 text-gray-900 dark:bg-gray-950 dark:text-gray-200`}
      >
        <ThemeProvider
          attribute="class"
          defaultTheme="system"
          enableSystem
          disableTransitionOnChange
        >
          <GraphQLProvider>
            <header className="bg-gray-100 text-gray-900 p-4 shadow-md dark:bg-gray-900 dark:text-gray-200">
              <nav className="container mx-auto flex justify-between items-center px-8">
                <Link href="/" className="text-2xl font-bold">
                  Perm-Predict
                </Link>
                <ul className="flex space-x-8 items-center">
                  <li>
                    <Link href="/" className="hover:text-gray-700 dark:hover:text-gray-50">
                      Predict
                    </Link>
                  </li>
                  <li>
                    <Link href="/create" className="hover:text-gray-700 dark:hover:text-gray-50">
                      Create
                    </Link>
                  </li>
                  <li>
                    <Link href="/explore" className="hover:text-gray-700 dark:hover:text-gray-50">
                      Explore
                    </Link>
                  </li>
                  <li>
                    <Link href="/help" className="hover:text-gray-700 dark:hover:text-gray-50">
                      Help
                    </Link>
                  </li>
                  <li>
                    <Link href="/about" className="hover:text-gray-700 dark:hover:text-gray-50">
                      About
                    </Link>
                  </li>
                  <li>
                    <ThemeToggle />
                  </li>
                </ul>
              </nav>
            </header>

            <main className="flex-grow">
              {children}
            </main>

            <footer className="bg-gray-100 text-gray-900 p-6 mt-8 dark:bg-gray-900 dark:text-gray-200">
              <div className="container mx-auto flex flex-col md:flex-row justify-between items-center">
                <div className="text-sm">
                  &copy; {new Date().getFullYear()} Perm-Predict. All rights reserved.
                </div>
                <div className="flex space-x-4 mt-4 md:mt-0">
                  <Link href="/contact" className="hover:text-gray-700 dark:hover:text-gray-50">
                    Contact Us
                  </Link>
                  <Link href="https://github.com/adamd3/perm-predict" target="_blank" rel="noopener noreferrer" className="hover:text-gray-700 dark:hover:text-gray-50">
                    Source Code
                  </Link>
                  {/* Placeholder images */}
                  <img src="/placeholder-1.svg" alt="Placeholder 1" className="h-6 w-6" />
                  <img src="/placeholder-2.svg" alt="Placeholder 2" className="h-6 w-6" />
                </div>
              </div>
            </footer>
          </GraphQLProvider>
        </ThemeProvider>
      </body>
    </html>
  );
}
