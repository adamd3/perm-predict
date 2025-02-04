import { NextResponse } from 'next/server'

export async function POST(request: Request) {
  const formData = await request.formData()
  
  const response = await fetch(`${process.env.BACKEND_URL}/predict/batch`, {
    method: 'POST',
    body: formData
  })
  
  const data = await response.json()
  return NextResponse.json(data)
}