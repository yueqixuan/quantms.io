"""
Caching utilities for quantms.io.
Provides memory and disk caching for frequently accessed data.
"""

import os
import hashlib
import pickle
import time
from typing import Any, Optional, Dict, Callable, TypeVar, cast
from pathlib import Path
from functools import wraps
import pyarrow as pa
import pyarrow.parquet as pq
import pandas as pd
from datetime import datetime, timedelta

from quantmsio.utils.logger import get_logger

logger = get_logger(__name__)

# Type variable for generic function return type
T = TypeVar("T")


class Cache:
    """
    Cache manager for storing frequently accessed data.
    Supports both memory and disk caching with TTL.
    """

    def __init__(
        self,
        cache_dir: Optional[str] = None,
        memory_size: int = 1000,
        ttl: Optional[int] = None,
    ):
        """
        Initialize the cache manager.

        Args:
            cache_dir: Directory for disk cache (None for memory-only)
            memory_size: Maximum number of items in memory cache
            ttl: Time to live in seconds (None for no expiration)
        """
        self.memory_cache: Dict[str, Any] = {}
        self.memory_size = memory_size
        self.ttl = ttl

        # Set up disk cache if directory provided
        if cache_dir:
            self.cache_dir = Path(cache_dir)
            self.cache_dir.mkdir(parents=True, exist_ok=True)
        else:
            self.cache_dir = None

        # Track access times for LRU eviction
        self.access_times: Dict[str, float] = {}

    def _get_cache_key(self, *args: Any, **kwargs: Any) -> str:
        """Generate a unique cache key from arguments."""
        # Convert args and kwargs to strings and sort for consistency
        key_parts = [str(arg) for arg in args]
        key_parts.extend(f"{k}={v}" for k, v in sorted(kwargs.items()))

        # Create hash of the key parts
        key = hashlib.sha256("|".join(key_parts).encode()).hexdigest()
        return key

    def _get_cache_path(self, key: str) -> Optional[Path]:
        """Get the disk cache path for a key."""
        if not self.cache_dir:
            return None
        return self.cache_dir / f"{key}.cache"

    def _is_expired(self, timestamp: float) -> bool:
        """Check if a cache entry has expired."""
        if self.ttl is None:
            return False
        return time.time() - timestamp > self.ttl

    def _evict_old_entries(self) -> None:
        """Remove old entries when cache is full."""
        if len(self.memory_cache) >= self.memory_size:
            # Sort by access time and remove oldest
            oldest_key = min(self.access_times.items(), key=lambda x: x[1])[0]
            del self.memory_cache[oldest_key]
            del self.access_times[oldest_key]

    def get(self, key: str) -> Optional[Any]:
        """
        Get a value from the cache.

        Args:
            key: Cache key

        Returns:
            Cached value or None if not found/expired
        """
        # Check memory cache first
        if key in self.memory_cache:
            timestamp = self.access_times[key]
            if not self._is_expired(timestamp):
                # Update access time
                self.access_times[key] = time.time()
                return self.memory_cache[key]
            else:
                # Remove expired entry
                del self.memory_cache[key]
                del self.access_times[key]

        # Check disk cache if enabled
        cache_path = self._get_cache_path(key)
        if cache_path and cache_path.exists():
            try:
                with open(cache_path, "rb") as f:
                    timestamp, value = pickle.load(f)
                    if not self._is_expired(timestamp):
                        # Store in memory cache for faster access
                        self._evict_old_entries()
                        self.memory_cache[key] = value
                        self.access_times[key] = time.time()
                        return value
                    else:
                        # Remove expired file
                        cache_path.unlink()
            except Exception as e:
                logger.warning(f"Failed to read cache file {cache_path}: {e}")

        return None

    def set(self, key: str, value: Any) -> None:
        """
        Store a value in the cache.

        Args:
            key: Cache key
            value: Value to store
        """
        timestamp = time.time()

        # Store in memory cache
        self._evict_old_entries()
        self.memory_cache[key] = value
        self.access_times[key] = timestamp

        # Store in disk cache if enabled
        cache_path = self._get_cache_path(key)
        if cache_path:
            try:
                with open(cache_path, "wb") as f:
                    pickle.dump((timestamp, value), f)
            except Exception as e:
                logger.warning(f"Failed to write cache file {cache_path}: {e}")

    def clear(self) -> None:
        """Clear all cache entries."""
        self.memory_cache.clear()
        self.access_times.clear()

        if self.cache_dir and self.cache_dir.exists():
            for cache_file in self.cache_dir.glob("*.cache"):
                try:
                    cache_file.unlink()
                except Exception as e:
                    logger.warning(f"Failed to delete cache file {cache_file}: {e}")


# Global cache instance
_cache = Cache(
    cache_dir=os.environ.get("QUANTMSIO_CACHE_DIR"),
    memory_size=int(os.environ.get("QUANTMSIO_CACHE_SIZE", "1000")),
    ttl=int(os.environ.get("QUANTMSIO_CACHE_TTL", "3600")),  # 1 hour default
)


def cached(ttl: Optional[int] = None) -> Callable:
    """
    Decorator for caching function results.

    Args:
        ttl: Time to live in seconds (None for default)

    Returns:
        Decorated function
    """

    def decorator(func: Callable[..., T]) -> Callable[..., T]:
        @wraps(func)
        def wrapper(*args: Any, **kwargs: Any) -> T:
            # Generate cache key
            cache_key = _cache._get_cache_key(
                func.__module__, func.__name__, *args, **kwargs
            )

            # Try to get from cache
            cached_value = _cache.get(cache_key)
            if cached_value is not None:
                return cast(T, cached_value)

            # Call function and cache result
            result = func(*args, **kwargs)
            _cache.set(cache_key, result)
            return result

        return wrapper

    return decorator


def cached_property(ttl: Optional[int] = None) -> Callable:
    """
    Decorator for caching class property results.
    Similar to @property but with caching.

    Args:
        ttl: Time to live in seconds (None for default)

    Returns:
        Decorated property
    """

    def decorator(func: Callable[..., T]) -> property:
        cache_key_prefix = f"{func.__module__}.{func.__name__}"

        @property
        @wraps(func)
        def wrapper(self: Any) -> T:
            # Generate cache key including instance id
            cache_key = _cache._get_cache_key(cache_key_prefix, id(self))

            # Try to get from cache
            cached_value = _cache.get(cache_key)
            if cached_value is not None:
                return cast(T, cached_value)

            # Call function and cache result
            result = func(self)
            _cache.set(cache_key, result)
            return result

        return wrapper

    return decorator


def clear_cache() -> None:
    """Clear all cache entries."""
    _cache.clear()


def get_cache_stats() -> Dict[str, Any]:
    """Get cache statistics."""
    return {
        "memory_entries": len(_cache.memory_cache),
        "memory_size": _cache.memory_size,
        "ttl": _cache.ttl,
        "cache_dir": str(_cache.cache_dir) if _cache.cache_dir else None,
    }
