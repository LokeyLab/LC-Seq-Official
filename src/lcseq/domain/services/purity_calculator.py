"""
Purity calculation service for synthesis validation.

This service provides a single implementation of purity calculation used by
all validation services (BayesianValidator, AdaptiveValidator, etc.).

Implementation based on THEORY.md Section 6.3.
"""

from typing import List, Dict, Any
import numpy as np
from ..entities.compound import Compound
from ..entities.peak import PeakType


class PurityCalculator:
    """
    Domain service for calculating compound purity.

    Purity is defined as the fraction of total chromatographic signal
    attributed to the putative product peak vs. all other peaks (impurities).

    This is a stateless service with a single static method to ensure
    consistent purity calculation across all validation workflows.

    Notes
    -----
    Purity = Area(PUTATIVE_PRODUCT) / [Area(PUTATIVE_PRODUCT) + Area(all others)]

    Only ACCEPTED peaks (peaks that passed detection thresholds) are included.
    Rejected peaks (failed significance, prominence, baseline, or SNR filters)
    are excluded from purity calculation.

    All non-product ACCEPTED peaks count as impurities:
    - TRUNCATION peaks (incomplete synthesis)
    - NULL peaks (unknown origin)
    - UNKNOWN peaks (unclassified)

    Purity is always in [0, 1]:
    - purity = 1.0: 100% pure product
    - purity = 0.5: 50% product, 50% impurities
    - purity = 0.0: No product detected

    References
    ----------
    THEORY.md Section 6.3: Purity Definition and Calculation

    Examples
    --------
    >>> from lcseq.domain.services import PurityCalculator
    >>> purity = PurityCalculator.calculate(compound)
    >>> purity
    0.85
    """

    @staticmethod
    def calculate(compound: Compound) -> float:
        """
        Calculate purity as product area / total area.

        Parameters
        ----------
        compound : Compound
            Compound with classified peaks (must have detected_peaks populated)

        Returns
        -------
        float
            Purity value in [0, 1]

        Notes
        -----
        Calculation steps:
        1. Sum areas of all peaks with peak_type == PUTATIVE_PRODUCT
        2. Sum areas of ALL peaks (product + impurities)
        3. Return product_area / total_area
        4. Clip result to [0, 1] for numerical safety

        Edge cases:
        - No detected peaks → purity = 0.0
        - Total area = 0 → purity = 0.0
        - All peaks are product → purity = 1.0

        References
        ----------
        THEORY.md Section 6.3: Purity Definition and Calculation

        Examples
        --------
        >>> compound = Compound(...)
        >>> compound.detected_peaks = [
        ...     Peak(area=100, peak_type=PeakType.PUTATIVE_PRODUCT),
        ...     Peak(area=15, peak_type=PeakType.TRUNCATION),
        ...     Peak(area=5, peak_type=PeakType.NULL)
        ... ]
        >>> PurityCalculator.calculate(compound)
        0.833  # 100 / (100 + 15 + 5)
        """
        if not compound.detected_peaks:
            return 0.0

        product_area = 0.0
        total_area = 0.0

        for peak in compound.detected_peaks:
            # Only include accepted peaks in purity calculation
            if not peak.is_accepted:
                continue
            total_area += peak.area
            if peak.peak_type == PeakType.PUTATIVE_PRODUCT:
                product_area += peak.area

        if total_area == 0:
            return 0.0

        purity = product_area / total_area
        return float(np.clip(purity, 0.0, 1.0))

    @staticmethod
    def calculate_from_peaks(peaks: List[Dict[str, Any]]) -> float:
        """
        Calculate purity from serialized peak data (JSONL records).

        This method enables purity calculation from peak dictionaries
        without requiring full Compound objects, useful for analysis
        scripts that read JSONL files directly.

        Parameters
        ----------
        peaks : List[Dict]
            Peak dictionaries with 'area', 'classification', and optionally 'is_accepted'

        Returns
        -------
        float
            Purity value in [0, 1]

        Notes
        -----
        Peak dictionaries must contain:
        - 'area': float (peak area)
        - 'classification': str (e.g., "PUTATIVE_PRODUCT", "TRUNCATION", "NULL")
        - 'is_accepted': bool (required)

        Only accepted peaks with area > 0 are included in the calculation.

        Edge cases:
        - Empty peaks list → purity = 0.0
        - No accepted peaks → purity = 0.0
        - Total area = 0 → purity = 0.0

        Examples
        --------
        >>> peaks = [
        ...     {"area": 100, "classification": "PUTATIVE_PRODUCT", "is_accepted": True},
        ...     {"area": 15, "classification": "TRUNCATION", "is_accepted": True},
        ...     {"area": 5, "classification": "NULL", "is_accepted": True}
        ... ]
        >>> PurityCalculator.calculate_from_peaks(peaks)
        0.833  # 100 / (100 + 15 + 5)
        """
        if not peaks:
            return 0.0

        total_area = 0.0
        product_area = 0.0

        for p in peaks:
            area = p.get("area", 0)
            if area <= 0:
                continue
            is_accepted = p["is_accepted"]  # Required field
            if not is_accepted:
                continue
            total_area += area
            if p.get("classification") == "PUTATIVE_PRODUCT":
                product_area += area

        if total_area == 0:
            return 0.0

        purity = product_area / total_area
        return float(np.clip(purity, 0.0, 1.0))
