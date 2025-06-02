"""
Enhanced logging configuration for quantmsio.
Provides structured logging, rotation, and environment-based configuration.
"""

import os
import sys
import json
import logging
import logging.handlers
from typing import Dict, Any, Optional
from pathlib import Path
import uuid
from functools import wraps
from datetime import datetime

# Default configuration
DEFAULT_CONFIG = {
    "level": "INFO",
    "format": "%(asctime)s - %(name)s - %(levelname)s - %(request_id)s - %(message)s",
    "date_format": "%Y-%m-%d %H:%M:%S",
    "file": None,
    "max_bytes": 10 * 1024 * 1024,  # 10MB
    "backup_count": 5,
    "log_dir": "logs",
}

class StructuredLogger(logging.Logger):
    """
    Enhanced logger that supports structured logging and request tracking.
    """
    def __init__(self, name: str):
        super().__init__(name)
        self.request_id = None

    def _log_structured(self, level: int, msg: str, *args, **kwargs) -> None:
        """
        Log a message with structured data.
        """
        extra = kwargs.get('extra', {})
        if not extra:
            kwargs['extra'] = extra
        
        if self.request_id:
            extra['request_id'] = self.request_id
        else:
            extra['request_id'] = 'NO_REQUEST'

        if isinstance(msg, dict):
            msg = json.dumps(msg)
        
        super()._log(level, msg, args, **kwargs)

    def set_request_id(self, request_id: Optional[str] = None) -> None:
        """Set a request ID for the current context."""
        self.request_id = request_id or str(uuid.uuid4())

    def clear_request_id(self) -> None:
        """Clear the current request ID."""
        self.request_id = None

class JsonFormatter(logging.Formatter):
    """
    Formatter that outputs JSON strings.
    """
    def format(self, record: logging.LogRecord) -> str:
        """Format the log record as JSON."""
        data = {
            'timestamp': datetime.fromtimestamp(record.created).isoformat(),
            'name': record.name,
            'level': record.levelname,
            'message': record.getMessage(),
            'request_id': getattr(record, 'request_id', 'NO_REQUEST'),
        }
        
        if record.exc_info:
            data['exception'] = self.formatException(record.exc_info)
            
        return json.dumps(data)

def with_request_tracking(func):
    """
    Decorator to automatically track requests with a unique ID.
    """
    @wraps(func)
    def wrapper(*args, **kwargs):
        logger = get_logger(func.__module__)
        request_id = str(uuid.uuid4())
        logger.set_request_id(request_id)
        try:
            return func(*args, **kwargs)
        finally:
            logger.clear_request_id()
    return wrapper

def configure_from_env() -> Dict[str, Any]:
    """
    Get logging configuration from environment variables.
    
    Environment variables:
    - QUANTMSIO_LOG_LEVEL: Logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL)
    - QUANTMSIO_LOG_FILE: Log file path
    - QUANTMSIO_LOG_FORMAT: Log format string
    - QUANTMSIO_LOG_DATE_FORMAT: Date format string
    - QUANTMSIO_LOG_JSON: Use JSON formatting if set to "true"
    """
    config = DEFAULT_CONFIG.copy()
    
    if os.getenv("QUANTMSIO_LOG_LEVEL"):
        config["level"] = os.getenv("QUANTMSIO_LOG_LEVEL")
        
    if os.getenv("QUANTMSIO_LOG_FILE"):
        config["file"] = os.getenv("QUANTMSIO_LOG_FILE")
        
    if os.getenv("QUANTMSIO_LOG_FORMAT"):
        config["format"] = os.getenv("QUANTMSIO_LOG_FORMAT")
        
    if os.getenv("QUANTMSIO_LOG_DATE_FORMAT"):
        config["date_format"] = os.getenv("QUANTMSIO_LOG_DATE_FORMAT")
        
    if os.getenv("QUANTMSIO_LOG_JSON", "").lower() == "true":
        config["json"] = True
    
    return config

def setup_logging(
    level: str = None,
    log_file: str = None,
    max_bytes: int = None,
    backup_count: int = None,
    use_json: bool = False
) -> None:
    """
    Configure logging for the application.
    
    Args:
        level: Logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL)
        log_file: Path to log file
        max_bytes: Maximum size of log file before rotation
        backup_count: Number of backup files to keep
        use_json: Whether to use JSON formatting
    """
    # Register custom logger class
    logging.setLoggerClass(StructuredLogger)
    
    # Get config from environment or use defaults
    config = configure_from_env()
    
    # Override with provided parameters
    if level:
        config["level"] = level
    if log_file:
        config["file"] = log_file
    if max_bytes:
        config["max_bytes"] = max_bytes
    if backup_count:
        config["backup_count"] = backup_count
    
    # Set up root logger
    root_logger = logging.getLogger()
    root_logger.setLevel(config["level"])
    
    # Clear existing handlers
    root_logger.handlers = []
    
    # Create formatter
    if use_json or config.get("json", False):
        formatter = JsonFormatter()
    else:
        formatter = logging.Formatter(
            fmt=config["format"],
            datefmt=config["date_format"]
        )
    
    # Console handler
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setFormatter(formatter)
    root_logger.addHandler(console_handler)
    
    # File handler (if configured)
    if config["file"]:
        log_path = Path(config["file"])
        
        # Create log directory if it doesn't exist
        log_path.parent.mkdir(parents=True, exist_ok=True)
        
        file_handler = logging.handlers.RotatingFileHandler(
            filename=str(log_path),
            maxBytes=config["max_bytes"],
            backupCount=config["backup_count"]
        )
        file_handler.setFormatter(formatter)
        root_logger.addHandler(file_handler)

def get_logger(name: str) -> StructuredLogger:
    """
    Get a logger instance with the given name.
    
    Args:
        name: Name of the logger (typically __name__)
        
    Returns:
        A StructuredLogger instance
    """
    return logging.getLogger(name)