# Logging System for quantmsio

This document describes the logging system for the quantmsio project and provides guidelines on how to use it effectively.

## Overview

The quantmsio project uses a centralized logging system built on top of Python's built-in `logging` module. The system provides:

- Consistent logging format across the entire project
- Support for different logging levels (DEBUG, INFO, WARNING, ERROR, CRITICAL)
- Colored console output for better readability
- File logging with rotation to prevent log files from growing too large
- Configuration via environment variables
- Structured logging with additional context

## Basic Usage

### Getting a Logger

To use logging in your module, import the `get_logger` function from `quantmsio.utils.logger` and create a logger instance:

```python
from quantmsio.utils.logger import get_logger

# Create a logger for this module
logger = get_logger(__name__)

# Use the logger
logger.info("This is an informational message")
logger.error("An error occurred")
```

### Logging Levels

The logging system supports the following levels (in order of increasing severity):

1. **DEBUG**: Detailed information, typically useful only for diagnosing problems
2. **INFO**: Confirmation that things are working as expected
3. **WARNING**: An indication that something unexpected happened, but the program can continue
4. **ERROR**: Due to a more serious problem, the program may not be able to perform some function
5. **CRITICAL**: A serious error, indicating that the program itself may be unable to continue running

Example:

```python
logger.debug("This is a debug message")
logger.info("This is an info message")
logger.warning("This is a warning message")
logger.error("This is an error message")
logger.critical("This is a critical message")
```

### Logging Exceptions

To log an exception with its traceback:

```python
try:
    # Some code that might raise an exception
    result = 1 / 0
except Exception as e:
    # Log the exception with traceback
    logger.exception(f"An error occurred: {e}")
    
    # Or log without traceback
    logger.error(f"An error occurred: {e}")
```

## Advanced Usage

### Custom Configuration

The logging system is automatically configured when the application starts, but you can customize it if needed:

```python
from quantmsio.utils.logger import setup_logging

# Configure logging with custom settings
setup_logging(
    level="DEBUG",  # Set logging level to DEBUG
    log_file="app.log",  # Log to a file
    max_file_size=5 * 1024 * 1024,  # 5 MB
    backup_count=3,  # Keep 3 backup files
    log_format="[%(asctime)s] %(levelname)s | %(name)s:%(lineno)d | %(message)s",
    date_format="%Y-%m-%d %H:%M:%S",
)
```

### Environment Variables

The logging system can be configured using environment variables:

- `QUANTMSIO_LOG_LEVEL`: Set the logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL)
- `QUANTMSIO_LOG_FILE`: Path to the log file

Example:

```bash
# Set logging level to DEBUG
export QUANTMSIO_LOG_LEVEL=DEBUG

# Log to a file
export QUANTMSIO_LOG_FILE=/path/to/quantmsio.log

# Run your application
python -m quantmsio
```

## Best Practices

1. **Use the appropriate logging level**: Use DEBUG for detailed information, INFO for general information, WARNING for unexpected events, ERROR for errors, and CRITICAL for serious errors.

2. **Include context in log messages**: Make log messages informative by including relevant context.

3. **Be consistent**: Use the same logging style throughout the project.

4. **Don't log sensitive information**: Avoid logging sensitive information like passwords, API keys, etc.

5. **Use structured logging when appropriate**: For complex data, consider using structured logging.

6. **Log at entry and exit points**: Log at the beginning and end of important functions to help with debugging.

7. **Include error details**: When logging errors, include details about what went wrong and how to fix it.

## Example

See `utils/logging_example.py` for a complete example of how to use the logging system.