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
        Peak matching descendant's product retention time - incomplete synthesis
    TRUNCATION_UNKNOWN : str
        Peak matching descendant's non-product peak - propagated species
        (e.g., oligomer, side product that originated from a truncation)
    PUTATIVE_PRODUCT : str
        Peak positionally consistent with expected product elution
        (NOT chemically validated - positional hypothesis only)
    UNKNOWN : str
        Peak that cannot be classified (unmatched to any descendant)
    """
    NULL = "NULL"
    TRUNCATION = "TRUNCATION"
    TRUNCATION_UNKNOWN = "TRUNCATION_UNKNOWN"
    PUTATIVE_PRODUCT = "PUTATIVE_PRODUCT"
    UNKNOWN = "UNKNOWN"


class RejectionReason(Enum):
    """
    Reason why a peak candidate was rejected during detection.

    Used to track why local maxima were not accepted as peaks,
    enabling diagnostic visualization of rejected candidates.

    Attributes
    ----------
    NONE : str
        Peak was accepted (not rejected)
    SIGNIFICANCE : str
        Rejected due to failing significance test (p-value >= alpha)
    PROMINENCE : str
        Rejected due to prominence below threshold
    BASELINE : str
        Rejected due to height below baseline threshold
    SNR : str
        Rejected due to local SNR below threshold
    NOT_MAXIMUM : str
        Rejected because not a verified local maximum
    """
    NONE = "none"
    SIGNIFICANCE = "significance"
    PROMINENCE = "prominence"
    BASELINE = "baseline"
    SNR = "snr"
    NOT_MAXIMUM = "not_maximum"


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
        Classification: NULL, TRUNCATION, TRUNCATION_UNKNOWN, PUTATIVE_PRODUCT, or UNKNOWN
    validation_status : ValidationStatus
        Synthesis validation result (default: NOT_VALIDATED)
    rejection_reason : RejectionReason
        Why this peak was rejected during detection (default: NONE = accepted).
        Rejected peaks are kept for diagnostic visualization.
    left_valley : Optional[float]
        Left valley position (if detected)
    right_valley : Optional[float]
        Right valley position (if detected)
    prominence : Optional[float]
        Chromatographic prominence (height above surrounding valleys)
    matched_compound_sequence : Optional[str]
        Block support sequence of the descendant compound whose peak matched this one.
        Used for tracing peak origin through the hierarchy.
    matched_peak_position : Optional[float]
        Retention time of the matched peak in the descendant compound.
        Together with matched_compound_sequence, enables chain tracing.
    matched_peak_type : Optional[PeakType]
        Peak type of the matched descendant peak (TRUNCATION, UNKNOWN, etc.).
        Indicates whether this peak matched a product or non-product peak.

    Notes
    -----
    - Position is absolute time (not array index)
    - Peak boundaries defined by valley detection or 5% threshold
    - Classification is positional hypothesis (NOT chemical validation)
    - Validation is separate from classification (THEORY.md Section 6.13)
    - Match tracking fields enable purity breakdown by peak origin

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
    rejection_reason: RejectionReason = RejectionReason.NONE
    left_valley: Optional[float] = None
    right_valley: Optional[float] = None
    prominence: Optional[float] = None
    # Match tracking - traces peak origin through hierarchy
    matched_compound_sequence: Optional[str] = None
    matched_peak_position: Optional[float] = None
    matched_peak_type: Optional["PeakType"] = None
    # Significance testing - p-value from detection
    # Used for filtering product candidates with stricter alpha_product
    p_value: Optional[float] = None

    # cLPE (chromatographic Linear Peptide Equation) validation fields
    # See clpe_validator.py for details
    clpe_residual: Optional[float] = None  # LogK(observed) - LogK(predicted)
    clpe_z_score: Optional[float] = None  # residual / model_std
    clpe_is_outlier: Optional[bool] = None  # abs(z_score) > threshold
    clpe_reselected: bool = False  # True if this peak was selected by cLPE re-selection

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

    @property
    def is_rejected(self) -> bool:
        """
        Check if this peak was rejected during detection.

        Returns
        -------
        bool
            True if rejection_reason is not NONE

        Notes
        -----
        Rejected peaks are local maxima that failed one of the
        detection filters (significance, prominence, baseline, SNR).
        They are kept for diagnostic purposes to show why
        expected peaks might not appear in the final results.
        """
        return self.rejection_reason != RejectionReason.NONE

    @property
    def is_accepted(self) -> bool:
        """
        Check if this peak was accepted during detection.

        Returns
        -------
        bool
            True if rejection_reason is NONE
        """
        return self.rejection_reason == RejectionReason.NONE

    def __repr__(self) -> str:
        """Detailed representation for debugging."""
        rejection_str = f", rejected={self.rejection_reason.value}" if self.is_rejected else ""
        return (
            f"Peak(position={self.position:.2f}, "
            f"height={self.height:.1f}, "
            f"area={self.area:.1f}, "
            f"type={self.peak_type.value}, "
            f"validation={self.validation_status.value}{rejection_str})"
        )
