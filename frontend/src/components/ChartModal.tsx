import React, { useRef, useMemo } from 'react';
import {
  Dialog,
  DialogContent,
  DialogHeader,
  DialogTitle,
  DialogDescription,
  DialogFooter,
  DialogClose,
} from '@/components/ui/dialog';
import { Button } from '@/components/ui/button';
import { Bar } from 'react-chartjs-2';
import { ChartData, ChartOptions } from 'chart.js';

interface ChartModalProps {
  isOpen: boolean;
  onClose: () => void;
  chartData: ChartData<'bar'> | null;
  chartOptions: ChartOptions<'bar'>;
  title: string;
}

const ChartModal: React.FC<ChartModalProps> = ({ isOpen, onClose, chartData, chartOptions, title }) => {
  const chartRef = useRef<any>(null); // Ref to access the Chart.js instance

  // Create a modified chartOptions for the modal with larger fonts
  const modalChartOptions = useMemo(() => {
    const newOptions = JSON.parse(JSON.stringify(chartOptions)); // Deep copy
    if (newOptions.scales && newOptions.scales.x && newOptions.scales.x.ticks && newOptions.scales.x.ticks.font) {
      newOptions.scales.x.ticks.font.size = 12; // Larger font for x-axis ticks
    }
    if (newOptions.scales && newOptions.scales.y && newOptions.scales.y.ticks && newOptions.scales.y.ticks.font) {
      newOptions.scales.y.ticks.font.size = 12; // Larger font for y-axis ticks
    }
    if (newOptions.scales && newOptions.scales.y && newOptions.scales.y.title && newOptions.scales.y.title.font) {
      newOptions.scales.y.title.font.size = 14; // Larger font for y-axis title
    } else if (newOptions.scales && newOptions.scales.y && newOptions.scales.y.title) {
        // If font object doesn't exist, create it
        newOptions.scales.y.title.font = { size: 14 };
    }
    return newOptions;
  }, [chartOptions]);

  const handleSaveChart = () => {
    if (chartRef.current) {
      // Chart.js 3.x and later: chartInstance is directly on the ref.current
      const chartInstance = chartRef.current.chart; 
      if (chartInstance) {
        const link = document.createElement('a');
        link.href = chartInstance.toBase64Image('image/png', 1); // Get data URL of the chart
        link.download = `${title.replace(/[^a-z0-9]/gi, '_').toLowerCase()}_chart.png`; // Sanitize title for filename
        document.body.appendChild(link);
        link.click();
        document.body.removeChild(link);
      }
    }
  };

  if (!chartData) return null; // Don't render if no data

  return (
    <Dialog open={isOpen} onOpenChange={onClose}>
      <DialogContent className="max-w-3xl">
        <DialogHeader>
          <DialogTitle className="text-black">{title}</DialogTitle>
          <DialogDescription>
            Detailed view of the feature summary.
          </DialogDescription>
        </DialogHeader>
        <div className="relative h-[400px] w-full">
          <Bar ref={chartRef} data={chartData} options={modalChartOptions} /> {/* Use modalChartOptions here */}
        </div>
        <DialogFooter>
          <DialogClose asChild>
            <Button variant="outline" onClick={onClose}>Close</Button>
          </DialogClose>
          <Button onClick={handleSaveChart}>Save Chart</Button>
        </DialogFooter>
      </DialogContent>
    </Dialog>
  );
};

export default ChartModal;
