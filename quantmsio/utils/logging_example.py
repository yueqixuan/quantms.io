"""
Example script demonstrating how to use the logging module in quantmsio.
This serves as a reference for developers to understand how to use the logging module consistently.
"""

from quantmsio.utils.logger import get_logger, setup_logging

# Get a logger for this module
logger = get_logger(__name__)


def demonstrate_logging_levels():
    """
    Demonstrate different logging levels
    """
    logger.debug("This is a DEBUG message - detailed information for debugging")
    logger.info("This is an INFO message - confirmation that things are working as expected")
    logger.warning("This is a WARNING message - something unexpected happened but the program can continue")
    logger.error("This is an ERROR message - something failed and the program may not be able to continue")
    logger.critical("This is a CRITICAL message - a serious error that may prevent the program from continuing")


def demonstrate_exception_logging():
    """
    Demonstrate logging exceptions
    """
    try:
        # Deliberately cause an exception
        result = 1 / 0
    except Exception as e:
        # Log the exception with traceback
        logger.exception(f"An error occurred: {e}")
        # Or log without traceback
        logger.error(f"An error occurred: {e}")


def demonstrate_custom_configuration():
    """
    Demonstrate custom logging configuration
    """
    # Configure logging with custom settings
    setup_logging(
        level="DEBUG",  # Set logging level to DEBUG
        log_file="example.log",  # Log to a file
        log_format="[%(asctime)s] %(levelname)s | %(name)s:%(lineno)d | %(message)s",  # Custom format
        date_format="%Y-%m-%d %H:%M:%S",  # Custom date format
    )
    
    logger.debug("This message will now be logged to the file as well")


def demonstrate_structured_logging():
    """
    Demonstrate structured logging with additional context
    """
    # Log with additional context
    logger.info(
        "Processing file",
        extra={
            "file_name": "example.txt",
            "file_size": 1024,
            "processing_time": 0.5,
        },
    )


if __name__ == "__main__":
    # Default logging configuration is applied when the logger module is imported
    # You can override it here if needed
    
    # Demonstrate different logging levels
    demonstrate_logging_levels()
    
    # Demonstrate exception logging
    demonstrate_exception_logging()
    
    # Note: The following would change the logging configuration for the entire application
    # Uncomment to see the effect
    # demonstrate_custom_configuration()
    
    # Demonstrate structured logging
    demonstrate_structured_logging()
    
    logger.info("Example script completed")