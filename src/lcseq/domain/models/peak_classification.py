"""
PeakClassification - result of classifying a peak for a compound.

Implementation based on THEORY.md Section 5.3.
"""

from dataclasses import dataclass
from typing import Optional
from ..entities.compound import Compound
from ..entities.peak import Peak, PeakType


@dataclass(frozen=True)
class PeakClassification:
    """
    Result of classifying a detected peak for a compound.

    Associates a peak with its compound and classification type
    (NULL/TRUNCATION/PUTATIVE_PRODUCT/UNKNOWN) with optional confidence.

    Attributes
    ----------
    compound : Compound
        The compound being analyzed
    peak : Peak
        The detected peak being classified
    classification : PeakType
        Classification type (NULL, TRUNCATION, PUTATIVE_PRODUCT, UNKNOWN)
    confidence : Optional[float]
        Classification confidence score [0.0, 1.0] (if available)

    Notes
    -----
    - Immutable value object (frozen dataclass)
    - Classification is positional hypothesis (NOT chemical validation)
    - Confidence represents classification certainty, not synthesis success
    - PUTATIVE_PRODUCT means "positionally consistent", not "chemically confirmed"
    - Validation is separate step (THEORY.md Section 6)

    Examples
    --------
    >>> from lcseq.domain.entities.building_block import BuildingBlock
    >>> from lcseq.domain.entities.chromatogram import Chromatogram
    >>> import numpy as np
    >>>
    >>> # Create compound and peak
    >>> chromatogram = Chromatogram(
    ...     time_points=np.array([1.0, 2.0, 3.0]),
    ...     counts=np.array([100.0, 200.0, 150.0])
    ... )
    >>> bb0 = BuildingBlock.from_code(0, "Pro")
    >>> bb1 = BuildingBlock.from_code(1, "Leu")
    >>> compound = Compound([bb0, bb1], chromatogram)
    >>>
    >>> peak = Peak(
    ...     position=2.0,
    ...     left_base=1.5,
    ...     right_base=2.5,
    ...     height=200.0,
    ...     area=150.0,
    ...     peak_type=PeakType.PUTATIVE_PRODUCT
    ... )
    >>>
    >>> # Classify peak
    >>> classification = PeakClassification(
    ...     compound=compound,
    ...     peak=peak,
    ...     classification=PeakType.PUTATIVE_PRODUCT,
    ...     confidence=0.95
    ... )

    References
    ----------
    THEORY.md Section 5.3: Peak Type Classification
    THEORY.md Section 5.6: Classification Scope and Limitations
    THEORY.md Section 6.13: Classification vs Validation
    """

    compound: Compound
    peak: Peak
    classification: PeakType
    confidence: Optional[float] = None

    def __post_init__(self) -> None:
        """Validate classification properties."""
        # Validate confidence range
        if self.confidence is not None:
            if not 0.0 <= self.confidence <= 1.0:
                raise ValueError(
                    f"Confidence must be in [0.0, 1.0], got {self.confidence}"
                )

    @property
    def is_product(self) -> bool:
        """
        Check if classified as putative product.

        Returns
        -------
        bool
            True if classification is PUTATIVE_PRODUCT

        Notes
        -----
        This is a positional hypothesis, NOT chemical confirmation.
        """
        return self.classification == PeakType.PUTATIVE_PRODUCT

    @property
    def is_truncation(self) -> bool:
        """
        Check if classified as truncation.

        Returns
        -------
        bool
            True if classification is TRUNCATION
        """
        return self.classification == PeakType.TRUNCATION

    @property
    def is_null(self) -> bool:
        """
        Check if classified as NULL (L₀).

        Returns
        -------
        bool
            True if classification is NULL
        """
        return self.classification == PeakType.NULL

    @property
    def is_unknown(self) -> bool:
        """
        Check if classification is unknown.

        Returns
        -------
        bool
            True if classification is UNKNOWN
        """
        return self.classification == PeakType.UNKNOWN

    def __repr__(self) -> str:
        """Detailed representation for debugging."""
        conf_str = f", confidence={self.confidence:.3f}" if self.confidence else ""
        return (
            f"PeakClassification("
            f"compound='{self.compound.positional_sequence}', "
            f"peak_pos={self.peak.position:.2f}, "
            f"type={self.classification.value}"
            f"{conf_str})"
        )

    def __str__(self) -> str:
        """String representation for display."""
        conf_str = f" (conf={self.confidence:.2f})" if self.confidence else ""
        return (
            f"{self.compound.positional_sequence}: "
            f"{self.classification.value} at {self.peak.position:.2f}{conf_str}"
        )
