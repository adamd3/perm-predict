import { NextResponse } from 'next/server'

export async function POST(request: Request) {
  const body = await request.json()
  
  const response = await fetch(`${process.env.BACKEND_URL}/predict/single`, {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify(body)
  })
  
  const data = await response.json()
  return NextResponse.json(data)
}

// export async function POST(request: Request) {
//   try {
//     const formData = await request.formData()
    
//     // Forward the request to FastAPI backend
//     const response = await fetch('http://your-fastapi-url/predict', {
//       method: 'POST',
//       body: formData,
//     })
    
//     if (!response.ok) {
//       throw new Error(`Prediction failed: ${response.statusText}`)
//     }
    
//     const data = await response.json()
//     return NextResponse.json(data)
//   } catch (error) {
//     return NextResponse.json(
//       { error: error.message },
//       { status: 500 }
//     )
//   }
// }