from fastapi import FastAPI, Request
from starlette.middleware.base import BaseHTTPMiddleware
from strawberry.fastapi import GraphQLRouter

import time
import logging
import datetime

from app.utils.logger import setup_logging
from app.schema import schema


class RequestLoggingMiddleware(BaseHTTPMiddleware):
    async def dispatch(self, request: Request, call_next):
        start_time = time.time()
        if request.method == "POST" and request.url.path == "/graphql":
            body = await request.body()
            logging.info(f"GraphQL Request Body: {body.decode()}")
        response = await call_next(request)
        process_time = time.time() - start_time
        logging.info(
            f"Method: {request.method} Path: {request.url.path} "
            f"Status: {response.status_code} Duration: {process_time:.3f}s"
        )
        return response


# Initialize logging
setup_logging()

# Create FastAPI app
app = FastAPI(
    title="Perm-Predict GraphQL API",
    description="Machine learning-based prediction of chemical accumulation in bacteria",
    version="2.0.0"
)

# Add middleware
app.add_middleware(RequestLoggingMiddleware)

# Create GraphQL router
graphql_app = GraphQLRouter(schema)

# Include GraphQL router
app.include_router(graphql_app, prefix="/graphql")

@app.get("/health")
async def health_check():
    """Health check endpoint."""
    return {
        "status": "healthy",
        "service": "perm-predict-api",
        "version": "2.0.0",
        "timestamp": datetime.datetime.now().isoformat(),
    }
