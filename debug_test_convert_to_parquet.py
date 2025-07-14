#!/usr/bin/env python3
"""
Debug script for test_convert_to_parquet function.
Run this script to manually debug the parquet conversion with verbose logging.
"""

import logging
import sys
import os
from pathlib import Path

# Add the project root to the Python path
project_root = Path(__file__).parent
sys.path.insert(0, str(project_root))

# Set up verbose logging
logging.basicConfig(
    level=logging.DEBUG,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler("debug_test_convert_to_parquet.log"),
    ],
)

# Import the test function
from tests.core.idxml.test_idxml_psm import test_convert_to_parquet


def main():
    """Run the test_convert_to_parquet function with debug logging."""
    print("=== Starting debug session for test_convert_to_parquet ===")
    print("Logs will be written to both console and debug_test_convert_to_parquet.log")
    print()

    try:
        # Run the test function
        test_convert_to_parquet()
        print("\n=== Debug session completed successfully! ===")
    except Exception as e:
        print(f"\n=== Debug session failed with error: {e} ===")
        import traceback

        print(f"=== Full traceback: ===")
        traceback.print_exc()
        return 1

    return 0


if __name__ == "__main__":
    exit_code = main()
    sys.exit(exit_code)
