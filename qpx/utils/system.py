import logging
import os
from pathlib import Path
import psutil
import shutil


def log_memory_usage(logger, step: str):
    """Log memory usage for a given step."""
    process = psutil.Process(os.getpid())
    mem_info = process.memory_info()
    logger.debug(
        f"Memory usage after {step}: RSS={mem_info.rss / 1024 / 1024:.1f}MB, VMS={mem_info.vms / 1024 / 1024:.1f}MB"
    )


def check_disk_space(output_dir: Path, required_gb: float = 1.0) -> bool:
    """
    Check if there's sufficient disk space for processing.

    Args:
        output_dir: Directory to check space for
        required_gb: Required space in GB (default: 1GB)

    Returns:
        True if sufficient space is available, False otherwise
    """
    try:
        total, used, free = shutil.disk_usage(output_dir)
        free_gb = free / (1024**3)  # Convert to GB
        return free_gb >= required_gb
    except Exception as e:
        logging.warning(f"Could not check disk space: {e}")
        return True  # Assume sufficient space if we can't check
