"""
Logger configuration for quantmsio.
"""

import logging


def get_logger(name: str) -> logging.Logger:
    """Get a logger with the given name."""
    return logging.getLogger(name)
