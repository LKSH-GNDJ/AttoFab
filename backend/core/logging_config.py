"""
core/logging_config.py - structured JSON logging.

Every log line is a single JSON object (timestamp, level, logger, message,
plus any extra fields passed via `logger.info(..., extra={...})`), so logs
are grep/jq-able and ready to ship to a log aggregator without a separate
parsing step.
"""
from __future__ import annotations

import json
import logging
import sys
import time
from pathlib import Path

LOG_DIR = Path(__file__).resolve().parent.parent.parent / "logs"

_RESERVED = set(logging.LogRecord(
    "", 0, "", 0, "", (), None
).__dict__.keys()) | {"message", "asctime"}


class JsonFormatter(logging.Formatter):
    def format(self, record: logging.LogRecord) -> str:
        payload = {
            "ts": time.strftime("%Y-%m-%dT%H:%M:%S", time.gmtime(record.created)),
            "level": record.levelname,
            "logger": record.name,
            "message": record.getMessage(),
        }
        if record.exc_info:
            payload["exc_info"] = self.formatException(record.exc_info)
        extras = {k: v for k, v in record.__dict__.items() if k not in _RESERVED}
        if extras:
            payload.update(extras)
        return json.dumps(payload, default=str)


def configure_logging(level: str = "INFO") -> None:
    LOG_DIR.mkdir(parents=True, exist_ok=True)

    root = logging.getLogger()
    root.setLevel(level)
    root.handlers.clear()

    formatter = JsonFormatter()

    stream_handler = logging.StreamHandler(sys.stdout)
    stream_handler.setFormatter(formatter)
    root.addHandler(stream_handler)

    file_handler = logging.FileHandler(LOG_DIR / "attofab-backend.log")
    file_handler.setFormatter(formatter)
    root.addHandler(file_handler)

    # Quiet down noisy third-party loggers unless something's actually wrong.
    logging.getLogger("uvicorn.access").setLevel("WARNING")
