import { NextRequest, NextResponse } from 'next/server';
import { mockResolvers } from '@/lib/mock-resolvers';

export async function POST(request: NextRequest) {
  try {
    const body = await request.json();
    const { query, variables } = body;
    
    // Simple GraphQL query parsing and routing
    if (query.includes('submitPredictionJob') && !query.includes('Batch')) {
      const result = mockResolvers.Mutation.submitPredictionJob(null, variables);
      return NextResponse.json({
        data: { submitPredictionJob: result }
      });
    }
    
    if (query.includes('submitBatchPredictionJob')) {
      const result = mockResolvers.Mutation.submitBatchPredictionJob(null, variables);
      return NextResponse.json({
        data: { submitBatchPredictionJob: result }
      });
    }
    
    if (query.includes('getPredictionResult')) {
      const result = mockResolvers.Query.getPredictionResult(null, variables);
      return NextResponse.json({
        data: { getPredictionResult: result }
      });
    }
    
    return NextResponse.json({
      errors: [{ message: 'Unknown query' }]
    }, { status: 400 });
    
  } catch (error) {
    return NextResponse.json({
      errors: [{ message: 'Invalid request' }]
    }, { status: 400 });
  }
}

export async function GET() {
  return NextResponse.json({
    message: 'GraphQL endpoint - use POST requests'
  });
}