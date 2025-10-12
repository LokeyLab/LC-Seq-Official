"""
Peak entity - represents detected feature in chromatogram.

Implementation based on THEORY.md Section 2.1, 5.3, 6.10.
"""

from dataclasses import dataclass
from enum import Enum
from typing import Optional


class PeakType(Enum):
    """
    Peak type classification based on retention time and DAG constraints.

    Classification types from THEORY.md Section 5.3.

    Attributes
    ----------
    NULL : str
        Peak at L₀ (full-null) retention time - DNA tag only
    TRUNCATION : str
        Peak matching ancestor retention time - incomplete synthesis
    PUTATIVE_PRODUCT : str
        Peak positionally consistent with expected product elution
        (NOT chemically validated - positional hypothesis only)
    UNKNOWN : str
        Peak that cannot be classified (late-eluting, unmatched)
    """
    NULL = "NULL"
    TRUNCATION = "TRUNCATION"
    PUTATIVE_PRODUCT = "PUTATIVE_PRODUCT"
    UNKNOWN = "UNKNOWN"


class ValidationStatus(Enum):
    """
    Synthesis validation status from adaptive validation framework.

    Validation categories from THEORY.md Section 6.10.

    Attributes
    ----------
    VALIDATED : str
        Synthesis succeeded with high confidence
    LIKELY_SUCCESS : str
        Synthesis probably succeeded (moderate confidence)
    UNCERTAIN : str
        Ambiguous result - cannot determine success/failure
    LIKELY_FAILURE : str
        Synthesis probably failed (moderate confidence)
    FAILED : str
        Synthesis failed with high confidence
    NOT_VALIDATED : str
        Validation not yet performed
    """
    VALIDATED = "VALIDATED"
    LIKELY_SUCCESS = "LIKELY_SUCCESS"
    UNCERTAIN = "UNCERTAIN"
    LIKELY_FAILURE = "LIKELY_FAILURE"
    FAILED = "FAILED"
    NOT_VALIDATED = "NOT_VALIDATED"


@dataclass
class Peak:
    """
    Represents a detected peak in chromatogram.

    A peak is defined by its retention time (position), integration boundaries,
    classification (NULL/TRUNCATION/PUTATIVE_PRODUCT/UNKNOWN), and optional
    validation status.

    Attributes
    ----------
    position : float
        Retention time in seconds/minutes (scalar absolute time)
    left_base : float
        Left integration boundary
    right_base : float
        Right integration boundary
    height : float
        Peak height (maximum intensity)
    area : float
        Integrated peak area
    peak_type : PeakType
        Classification: NULL, TRUNCATION, PUTATIVE_PRODUCT, or UNKNOWN
    validation_status : ValidationStatus
        Synthesis validation result (default: NOT_VALIDATED)
    left_valley : Optional[float]
        Left valley position (if detected)
    right_valley : Optional[float]
        Right valley position (if detected)
    prominence : Optional[float]
        Chromatographic prominence (height above surrounding valleys)

    Notes
    -----
    - Position is absolute time (not array index)
    - Peak boundaries defined by valley detection or 5% threshold
    - Classification is positional hypothesis (NOT chemical validation)
    - Validation is separate from classification (THEORY.md Section 6.13)

    References
    ----------
    THEORY.md Section 2.1: Core Entities
    THEORY.md Section 5.3: Peak Type Classification
    THEORY.md Section 6.10: Validation Classification
    THEORY.md Section 2.3.1: Absolute Time Representation
    """

    position: float
    left_base: float
    right_base: float
    height: float
    area: float
    peak_type: PeakType = PeakType.UNKNOWN
    validation_status: ValidationStatus = ValidationStatus.NOT_VALIDATED
    left_valley: Optional[float] = None
    right_valley: Optional[float] = None
    prominence: Optional[float] = None

    def __post_init__(self) -> None:
        """Validate peak properties."""
        if self.height < 0:
            raise ValueError(f"Peak height must be non-negative, got {self.height}")
        
        if self.area < 0:
            raise ValueError(f"Peak area must be non-negative, got {self.area}")
        
        if self.left_base >= self.right_base:
            raise ValueError(
                f"Left base ({self.left_base}) must be < right base ({self.right_base})"
            )
        
        if not (self.left_base <= self.position <= self.right_base):
            raise ValueError(
                f"Peak position ({self.position}) must be within boundaries "
                f"[{self.left_base}, {self.right_base}]"
            )

        if self.prominence is not None and self.prominence < 0:
            raise ValueError(f"Prominence must be non-negative, got {self.prominence}")

    @property
    def width(self) -> float:
        """
        Peak width (base-to-base).

        Returns
        -------
        float
            Right base - left base
        """
        return self.right_base - self.left_base

    @property
    def is_product_peak(self) -> bool:
        """
        Check if this is a putative product peak.

        Returns
        -------
        bool
            True if peak_type is PUTATIVE_PRODUCT

        Notes
        -----
        This is a positional hypothesis, NOT chemical confirmation.
        See THEORY.md Section 5.6 for scope and limitations.
        """
        return self.peak_type == PeakType.PUTATIVE_PRODUCT

    @property
    def is_validated(self) -> bool:
        """
        Check if synthesis is validated as successful.

        Returns
        -------
        bool
            True if validation_status is VALIDATED

        Notes
        -----
        Validation is separate from classification.
        A peak can be PUTATIVE_PRODUCT but not VALIDATED.
        """
        return self.validation_status == ValidationStatus.VALIDATED

    def __repr__(self) -> str:
        """Detailed representation for debugging."""
        return (
            f"Peak(position={self.position:.2f}, "
            f"height={self.height:.1f}, "
            f"area={self.area:.1f}, "
            f"type={self.peak_type.value}, "
            f"validation={self.validation_status.value})"
        )
