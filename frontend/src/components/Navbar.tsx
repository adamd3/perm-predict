'use client';

import Link from "next/link";
import { useRouter, usePathname } from "next/navigation";
import { ThemeToggle } from "./ThemeToggle";

export function Navbar() {
  const router = useRouter();
  const pathname = usePathname();

  const handlePredictClick = (e: React.MouseEvent) => {
    e.preventDefault();
    router.push(`/predict?_t=${new Date().getTime()}`);
  };

  const handleExploreClick = (e: React.MouseEvent) => {
    e.preventDefault();
    router.push(`/explore?_t=${new Date().getTime()}`);
  };

  return (
    <header className="bg-gray-100 text-gray-900 p-4 shadow-md dark:bg-gray-900 dark:text-gray-200">
      <nav className="container mx-auto flex justify-between items-center px-8">
        <Link href="/" className="text-2xl font-bold">
          Perm-Predict
        </Link>
        <ul className="flex space-x-8 items-center">
          <li>
            <a onClick={handlePredictClick} className="hover:text-gray-700 dark:hover:text-gray-50 cursor-pointer">
              Predict
            </a>
          </li>
          <li>
            <Link href="/create" className="hover:text-gray-700 dark:hover:text-gray-50">
              Create
            </Link>
          </li>
          <li>
            <a onClick={handleExploreClick} className="hover:text-gray-700 dark:hover:text-gray-50 cursor-pointer">
              Explore
            </a>
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
  );
}
