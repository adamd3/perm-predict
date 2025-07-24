"""
Simple launcher for the FastAPI app.
This allows running uvicorn main:app from the backend directory.
"""

from app.main import app

if __name__ == "__main__":
    import uvicorn
    uvicorn.run("main:app", host="0.0.0.0", port=8000, reload=True)