"""
core/auth.py - X-API-Key middleware, timing-safe comparison.

Auth is opt-in: if ATTOFAB_API_KEY is unset (the default for local/dev
use), the middleware is a no-op. Set ATTOFAB_API_KEY to require every
request to carry a matching X-API-Key header - useful once the backend is
exposed beyond localhost.
"""
from __future__ import annotations

import hmac
import os

from fastapi import Request
from fastapi.responses import JSONResponse
from starlette.middleware.base import BaseHTTPMiddleware

# Paths that never require an API key, even when one is configured.
_PUBLIC_PATHS = {"/api/health", "/docs", "/openapi.json", "/redoc"}


class ApiKeyMiddleware(BaseHTTPMiddleware):
    async def dispatch(self, request: Request, call_next):
        expected = os.environ.get("ATTOFAB_API_KEY")
        if not expected or request.url.path in _PUBLIC_PATHS:
            return await call_next(request)

        provided = request.headers.get("X-API-Key", "")
        # hmac.compare_digest is constant-time - avoids leaking key length/
        # prefix information through response-time side channels.
        if not hmac.compare_digest(provided, expected):
            return JSONResponse(status_code=401, content={"detail": "Invalid or missing X-API-Key"})

        return await call_next(request)
