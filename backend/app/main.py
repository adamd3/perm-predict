from fastapi import FastAPI, Request, status, Depends # Import status and Depends for HTTP status codes
from starlette.middleware.base import BaseHTTPMiddleware
from strawberry.fastapi import GraphQLRouter
from fastapi.responses import JSONResponse # Import JSONResponse
from slowapi.errors import RateLimitExceeded # Import RateLimitExceeded

import time
import logging
import datetime

# Import for fastapi-limiter
from fastapi_limiter import FastAPILimiter
from redis.asyncio import Redis as Aioredis # Use Aioredis for async support
from fastapi_limiter.depends import RateLimiter # Import RateLimiter for dependencies

from app.utils.logger import logger
from app.schema import schema
from app.config import settings # Import settings


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


# Create FastAPI app
app = FastAPI(
    title="Perm-Predict GraphQL API",
    description="Machine learning-based prediction of chemical accumulation in bacteria",
    version="2.0.0"
)

# Add middleware
app.add_middleware(RequestLoggingMiddleware)

@app.on_event("startup")
async def startup():
    redis_connection = Aioredis.from_url(settings.CELERY_BROKER_URL, encoding="utf-8", decode_responses=True)
    await FastAPILimiter.init(redis_connection)
    logger.info("FastAPILimiter initialized.")

@app.on_event("shutdown")
async def shutdown():
    if FastAPILimiter.redis:
        await FastAPILimiter.redis.close()
    logger.info("FastAPILimiter shut down.")

@app.exception_handler(RateLimitExceeded)
async def rate_limit_exception_handler(request: Request, exc: RateLimitExceeded):
    """
    Handles RateLimitExceeded exceptions, returning a GraphQL-formatted error.
    """
    logger.warning(f"Rate limit exceeded for client {request.client.host}: {exc.detail}")
    return JSONResponse(
        status_code=status.HTTP_200_OK, # GraphQL typically returns 200 OK with errors in payload
        content={
            "data": None,
            "errors": [
                {
                    "message": f"Too many requests: {exc.detail}",
                    "extensions": {
                        "code": "RATE_LIMIT_EXCEEDED",
                        "status": status.HTTP_429_TOO_MANY_REQUESTS,
                        "retryAfter": exc.retry_after,
                    }
                }
            ]
        },
        headers={"Retry-After": str(exc.retry_after)}
    )

# Create GraphQL router
graphql_app = GraphQLRouter(schema, dependencies=[Depends(RateLimiter(times=30, seconds=60))])

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
