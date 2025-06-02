"""
Base command class providing common functionality for all quantmsio commands.
"""

import click
from typing import Any, Dict, Optional, List
from pathlib import Path
import time
from functools import wraps
import logging

from quantmsio.utils.logger import get_logger, with_request_tracking
from quantmsio.core.common import QUANTMSIO_VERSION


class BaseCommand:
    """
    Base class for all quantmsio commands providing common functionality.

    Features:
    - Logging with request tracking
    - Progress reporting
    - Input validation
    - Error handling
    - Performance metrics
    """

    def __init__(self, name: str):
        """
        Initialize the command.

        Args:
            name: Name of the command
        """
        self.logger = get_logger(f"quantmsio.commands.{name}")
        self.start_time = None
        self.name = name

    def validate_input_file(
        self, file_path: str, required_extensions: Optional[List[str]] = None
    ) -> Path:
        """
        Validate that an input file exists and has the correct extension.

        Args:
            file_path: Path to the input file
            required_extensions: List of allowed file extensions (e.g., ['.txt', '.csv'])

        Returns:
            Path object for the validated file

        Raises:
            click.BadParameter: If the file doesn't exist or has wrong extension
        """
        path = Path(file_path)
        if not path.exists():
            raise click.BadParameter(f"Input file does not exist: {file_path}")

        if required_extensions and path.suffix.lower() not in required_extensions:
            raise click.BadParameter(
                f"Input file must have one of these extensions: {', '.join(required_extensions)}"
            )

        return path

    def validate_output_file(
        self, file_path: str, required_extension: Optional[str] = None
    ) -> Path:
        """
        Validate that an output file path is valid and create parent directories if needed.

        Args:
            file_path: Path where the output file will be written
            required_extension: Required file extension (e.g., '.parquet')

        Returns:
            Path object for the validated file

        Raises:
            click.BadParameter: If the path is invalid or has wrong extension
        """
        path = Path(file_path)

        if required_extension and path.suffix.lower() != required_extension.lower():
            raise click.BadParameter(
                f"Output file must have extension: {required_extension}"
            )

        # Create parent directories if they don't exist
        path.parent.mkdir(parents=True, exist_ok=True)

        return path

    def start_progress(self, message: str) -> None:
        """Start a progress operation with the given message."""
        self.start_time = time.time()
        self.logger.info(f"Starting {message}")

    def end_progress(self, message: str) -> None:
        """End a progress operation with the given message."""
        if self.start_time:
            elapsed = time.time() - self.start_time
            self.logger.info(f"Completed {message} in {elapsed:.2f} seconds")
            self.start_time = None

    @staticmethod
    def with_error_handling(func):
        """
        Decorator to add error handling to command methods.
        Catches exceptions, logs them, and shows user-friendly error messages.
        """

        @wraps(func)
        def wrapper(self, *args, **kwargs):
            try:
                return func(self, *args, **kwargs)
            except Exception as e:
                self.logger.exception(f"Error in {self.name} command: {str(e)}")
                raise click.ClickException(
                    f"Error in {self.name} command: {str(e)}\n"
                    f"Check the logs for more details."
                )

        return wrapper

    def get_metadata(self) -> Dict[str, Any]:
        """Get metadata for the command execution."""
        return {
            "command": self.name,
            "version": QUANTMSIO_VERSION,
            "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
        }

    @staticmethod
    def common_options(f):
        """Common options for all commands."""
        f = click.option(
            "--verbose",
            is_flag=True,
            help="Enable verbose logging",
        )(f)
        f = click.option(
            "--quiet",
            is_flag=True,
            help="Suppress all output except errors",
        )(f)
        f = click.option(
            "--log-file",
            help="Log file path",
        )(f)
        
        @wraps(f)
        def wrapper(*args, **kwargs):
            """Wrapper for all commands."""
            try:
                # Extract and remove common options from kwargs
                verbose = kwargs.pop("verbose", False)
                quiet = kwargs.pop("quiet", False)
                log_file = kwargs.pop("log_file", None)
                
                # Configure logging based on options
                logger = get_logger("quantmsio.cli")  # Use a default logger name for CLI operations
                if verbose:
                    logger.setLevel(logging.DEBUG)
                elif quiet:
                    logger.setLevel(logging.WARNING)
                
                if log_file:
                    file_handler = logging.FileHandler(log_file)
                    logger.addHandler(file_handler)
                    
                return f(*args, **kwargs)
            except Exception as e:
                logger = get_logger("quantmsio.cli")  # Use same logger name for consistency
                logger.exception(str(e))
                raise click.ClickException(
                    f"Error: {str(e)}\nCheck the logs for more details."
                )
        
        return wrapper
