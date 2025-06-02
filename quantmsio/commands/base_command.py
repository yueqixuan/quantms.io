"""
Base command class providing common functionality for all quantmsio commands.
"""

import click
from typing import Any, Dict, Optional, List
from pathlib import Path
import time
from functools import wraps

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
    
    def validate_input_file(self, file_path: str, required_extensions: Optional[List[str]] = None) -> Path:
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
    
    def validate_output_file(self, file_path: str, required_extension: Optional[str] = None) -> Path:
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
            raise click.BadParameter(f"Output file must have extension: {required_extension}")
            
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
            "timestamp": time.strftime("%Y-%m-%d %H:%M:%S")
        }
    
    @staticmethod
    def common_options(func):
        """
        Decorator to add common command options.
        """
        @click.option(
            "--verbose", "-v",
            is_flag=True,
            help="Enable verbose logging"
        )
        @click.option(
            "--quiet", "-q",
            is_flag=True,
            help="Suppress all output except errors"
        )
        @click.option(
            "--log-file",
            type=click.Path(),
            help="Log file path"
        )
        @wraps(func)
        def wrapper(*args, **kwargs):
            return func(*args, **kwargs)
        return wrapper 