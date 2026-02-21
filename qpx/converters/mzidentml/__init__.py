"""mzIdentML converter — converts mzIdentML (mzid) files to QPX Parquet.

Requires lxml (optional dependency). Install with: pip install qpx[mzidentml]
"""


def __getattr__(name: str):
    if name == "MzIdentMLPsmAdapter":
        from qpx.converters.mzidentml.psm_adapter import MzIdentMLPsmAdapter

        return MzIdentMLPsmAdapter
    if name == "MzIdentMLConverter":
        from qpx.converters.mzidentml.converter import MzIdentMLConverter

        return MzIdentMLConverter
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


__all__ = ["MzIdentMLPsmAdapter", "MzIdentMLConverter"]
