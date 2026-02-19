"""Tests for resolve_columns() utility."""
import pytest


def test_resolve_columns_basic():
    from qpx.converters.base import resolve_columns
    mappings = {
        "intensity": ["Precursor.Quantity"],
        "rt": ["RT"],
    }
    available = {"Precursor.Quantity", "RT", "Run", "Sequence"}
    result = resolve_columns(mappings, available)
    assert result == {"intensity": "Precursor.Quantity", "rt": "RT"}


def test_resolve_columns_ordered_fallback():
    from qpx.converters.base import resolve_columns
    mappings = {
        "rt": ["RT", "RetentionTime"],
    }
    # Only the fallback name exists
    available = {"RetentionTime", "Sequence"}
    result = resolve_columns(mappings, available)
    assert result == {"rt": "RetentionTime"}


def test_resolve_columns_first_match_wins():
    from qpx.converters.base import resolve_columns
    mappings = {
        "rt": ["RT", "RetentionTime"],
    }
    # Both exist — first candidate wins
    available = {"RT", "RetentionTime"}
    result = resolve_columns(mappings, available)
    assert result == {"rt": "RT"}


def test_resolve_columns_missing_skipped():
    from qpx.converters.base import resolve_columns
    mappings = {
        "intensity": ["Precursor.Quantity"],
        "predicted_rt": ["Predicted.RT"],
    }
    available = {"Precursor.Quantity"}
    result = resolve_columns(mappings, available)
    assert result == {"intensity": "Precursor.Quantity"}
    assert "predicted_rt" not in result


def test_resolve_columns_empty():
    from qpx.converters.base import resolve_columns
    result = resolve_columns({}, {"col1", "col2"})
    assert result == {}
