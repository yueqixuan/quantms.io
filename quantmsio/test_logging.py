#!/usr/bin/env python
"""
Test script for the new logging system.
Run this script to verify that the logging system works as expected.
"""

import os
import sys
import time

# Add the project root to the Python path
sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)))

from quantmsio.utils.logger import get_logger, setup_logging

# Get a logger for this module
logger = get_logger("test_logging")


def test_basic_logging():
    """Test basic logging functionality"""
    logger.info("Testing basic logging")
    logger.debug("This is a DEBUG message")
    logger.info("This is an INFO message")
    logger.warning("This is a WARNING message")
    logger.error("This is an ERROR message")
    logger.critical("This is a CRITICAL message")


def test_exception_logging():
    """Test exception logging"""
    logger.info("Testing exception logging")
    try:
        # Deliberately cause an exception
        result = 1 / 0
    except Exception as e:
        # Log the exception with traceback
        logger.exception(f"An exception occurred: {e}")


def test_file_logging():
    """Test logging to a file"""
    log_file = "test_logging.log"
    logger.info(f"Testing file logging to {log_file}")
    
    # Configure logging to write to a file
    setup_logging(
        level="DEBUG",
        log_file=log_file,
        max_file_size=1024 * 1024,  # 1 MB
        backup_count=3,
    )
    
    logger.debug("This message should be written to the log file")
    logger.info("Check the log file to verify that this message was written")
    
    return log_file


def test_environment_variables():
    """Test logging configuration via environment variables"""
    logger.info("Testing environment variable configuration")
    
    # Save original environment variables
    original_level = os.environ.get("QUANTMSIO_LOG_LEVEL")
    original_file = os.environ.get("QUANTMSIO_LOG_FILE")
    
    try:
        # Set environment variables
        os.environ["QUANTMSIO_LOG_LEVEL"] = "DEBUG"
        os.environ["QUANTMSIO_LOG_FILE"] = "env_test.log"
        
        # Import the module again to trigger configuration from environment variables
        from quantmsio.utils.logger import configure_from_env
        
        # Get configuration from environment variables
        config = configure_from_env()
        logger.info(f"Environment configuration: {config}")
        
        # Apply configuration
        setup_logging(**config)
        
        # Log a message
        logger.debug("This message should be logged at DEBUG level to env_test.log")
        
    finally:
        # Restore original environment variables
        if original_level is not None:
            os.environ["QUANTMSIO_LOG_LEVEL"] = original_level
        else:
            os.environ.pop("QUANTMSIO_LOG_LEVEL", None)
            
        if original_file is not None:
            os.environ["QUANTMSIO_LOG_FILE"] = original_file
        else:
            os.environ.pop("QUANTMSIO_LOG_FILE", None)


def main():
    """Main function"""
    logger.info("Starting logging system test")
    
    # Test basic logging
    test_basic_logging()
    
    # Test exception logging
    test_exception_logging()
    
    # Test file logging
    log_file = test_file_logging()
    
    # Test environment variables
    test_environment_variables()
    
    logger.info("Logging system test completed")
    
    # Print instructions for checking the log file
    print(f"\nCheck the log files to verify that logging is working correctly:")
    print(f"- {log_file}")
    print(f"- env_test.log")


if __name__ == "__main__":
    main()