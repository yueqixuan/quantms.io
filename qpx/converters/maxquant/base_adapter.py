"""Shared base adapter for MaxQuant converters.

Provides common data-loading helpers used by
:class:`MaxQuantFeatureAdapter` and :class:`MaxQuantPgAdapter`.
"""

from __future__ import annotations

from typing import Optional

from qpx.converters.base import BaseConverter


class MaxQuantBaseAdapter(BaseConverter):
    """Base adapter with MaxQuant-specific loading utilities.

    Subclasses inherit ``_load_sdrf()`` so the SDRF loading logic is
    defined once rather than duplicated across the feature and PG adapters.
    """

    def _load_sdrf(self, sdrf_path: Optional[str]) -> tuple[dict, str, list]:
        """Load SDRF and return ``(sample_map, experiment_type, tmt_channels)``.

        When *sdrf_path* is ``None`` the method returns sensible LFQ defaults.

        Args:
            sdrf_path: Filesystem path to the SDRF file, or ``None``.

        Returns:
            A 3-tuple of:
              - *sample_map*: ``{run_file -> sample_accession}``
              - *experiment_type*: e.g. ``"LFQ"``, ``"TMT11"``
              - *tmt_channels*: sorted list of TMT channel labels (empty for LFQ)
        """
        if not sdrf_path:
            return {}, "LFQ", []

        from qpx.core.sdrf import SDRFHandler

        handler = SDRFHandler(sdrf_path)
        sample_map = handler.get_sample_map_run()
        experiment_type = handler.get_experiment_type_from_sdrf()

        tmt_channels: list[str] = []
        if experiment_type and "TMT" in experiment_type.upper():
            labels = handler.sdrf_table.get("comment[label]")
            if labels is not None:
                tmt_labels = [
                    l for l in labels.unique() if l and "TMT" in str(l).upper()
                ]
                tmt_channels = sorted(tmt_labels)

        return sample_map, experiment_type or "LFQ", tmt_channels
