"""
LC-Seq Configuration - Centralized parameter defaults.

All default parameters for peak detection, filtering, and analysis are defined here.
Modify values in this file to change behavior across the entire codebase.

USAGE
-----
To change default parameters, simply edit the values in the classes below.
Changes will automatically propagate to:
  - PeakDetector.detect_peaks()
  - ProcessChromatogramsUseCase.execute()
  - ProcessChromatogramsWithIntegrationUseCase.execute()
  - LineageOffsetPlotter.plot()

EXAMPLE
-------
To disable local SNR filtering globally:
    MIN_SNR = 0.0

To make baseline filter more strict:
    MIN_BASELINE_SDS = 2.0  # Require 2 SDs above minimum

To retain only top 90% most prominent peaks:
    PROMINENCE_PERCENTILE = 0.1  # Lower = more strict

CURRENT CONFIGURATION
---------------------
See PeakDetectionConfig and VisualizationConfig classes below.
"""

# =============================================================================
# Peak Detection Parameters
# =============================================================================


class PeakDetectionConfig:
    """Peak detection and filtering parameters."""

    # Poisson Z-score threshold for statistical significance
    # Z > 3.0 corresponds to p < 0.001 (highly significant)
    Z_THRESHOLD = 3.0

    # Prominence percentile threshold
    # 0.2 = retain top 80% most prominent peaks
    # Lower value = more strict (fewer peaks retained)
    # Set to 0.0 to disable (retain all peaks)
    PROMINENCE_PERCENTILE = 0.5

    # Local SNR filter (adaptive threshold)
    # Threshold = min(SNRs) + MIN_SNR * std(SNRs)
    # Set to 0.0 to disable local SNR filtering
    # NOTE: This threshold is NOT visualized (it would be misleading)
    # DISABLED for now since it can't be accurately visualized
    MIN_SNR = 0.001

    # Global baseline filter (in standard deviations)
    # Peak height must exceed: min(signal) + MIN_BASELINE_SDS * std(signal)
    # Set to 0.0 to disable global baseline filtering
    MIN_BASELINE_SDS = 1.0

    # Signal variant to use for detection
    # "raw" = no baseline correction, no smoothing (per THEORY.md)
    SIGNAL_VARIANT = "raw"

    # Truncation boundary margin (in seconds)
    # Product peaks must be this far beyond max(truncation_positions)
    # Accounts for retention time variability in peak matching
    TRUNCATION_MARGIN = 60.0

    # Peak matching tolerance (relative fraction)
    # Used for Hungarian algorithm and position matching
    # Relative tolerance: |measured - expected| / expected
    PEAK_MATCHING_TOLERANCE = 0.01

    # Minimum threshold for Hungarian algorithm (fraction of signal length)
    # Accounts for discrete fractionation in LC-Seq
    HUNGARIAN_MIN_THRESHOLD = 0.02


# =============================================================================
# Visualization Parameters
# =============================================================================


class VisualizationConfig:
    """Visualization parameters for lineage plots."""

    # Layout parameters
    SECONDS_PER_MINUTE = 60.0  # Time conversion factor
    OFFSET_SPACING = 0.5  # Normal vertical spacing between compounds
    GROUP_SPACING_EXTRA = 1.0  # Extra spacing between different canonical sequence groups

    # Color parameters
    USE_COLORMAP = True  # True = rainbow colormap by canonical sequence, False = all black
    # When True, chemically identical compounds (same canonical sequence) get the same color
    # Colors assigned using rainbow colormap, organized by level then sorted by sequence
    # When False, all signals and markers are black (publication-ready)

    # Line styling parameters
    LINEWIDTH_DEFAULT = 1.5  # Default trace line width
    LINEWIDTH_REFERENCE = 3.0  # Reference compound line width (thicker)
    ALPHA_TRACE = 1.0  # Trace transparency (1.0 = opaque)

    # Marker parameters
    MARKER_SIZE = 6  # Peak marker size

    # Text parameters
    LABEL_FONTSIZE = 10  # Sequence label font size
    LABEL_MAX_LENGTH = 40  # Maximum label length before truncation
    LABEL_TRUNCATE_LENGTH = 4  # Number of characters to show when truncating

    # Figure dimensions
    FIG_WIDTH = 16  # Figure width in inches
    FIG_HEIGHT_BASE = 4  # Base figure height
    FIG_HEIGHT_PER_TRACE = 0.4  # Additional height per trace
    FIG_HEIGHT_MIN = 10  # Minimum figure height


# =============================================================================
# Consensus Mode Parameters
# =============================================================================


class ConsensusConfig:
    """Consensus mode parameters for positional variant aggregation."""

    # Correlation threshold for validity
    # Variants must have min(pairwise correlation) > threshold
    # Recommended: 0.8 (strong similarity required)
    CORRELATION_THRESHOLD = 0.0

    # Aggregation method for consensus signal
    # "mean" = faster, smoother (recommended)
    # "median" = more robust to outliers
    AGGREGATION_METHOD = "mean"


# =============================================================================
# Easy access to default parameters
# =============================================================================

# Peak detection defaults (for convenience)
DEFAULT_Z_THRESHOLD = PeakDetectionConfig.Z_THRESHOLD
DEFAULT_PROMINENCE_PERCENTILE = PeakDetectionConfig.PROMINENCE_PERCENTILE
DEFAULT_MIN_SNR = PeakDetectionConfig.MIN_SNR
DEFAULT_MIN_BASELINE_SDS = PeakDetectionConfig.MIN_BASELINE_SDS
DEFAULT_SIGNAL_VARIANT = PeakDetectionConfig.SIGNAL_VARIANT
DEFAULT_TRUNCATION_MARGIN = PeakDetectionConfig.TRUNCATION_MARGIN
DEFAULT_PEAK_MATCHING_TOLERANCE = PeakDetectionConfig.PEAK_MATCHING_TOLERANCE
DEFAULT_HUNGARIAN_MIN_THRESHOLD = PeakDetectionConfig.HUNGARIAN_MIN_THRESHOLD

# Consensus mode defaults (for convenience)
DEFAULT_CORRELATION_THRESHOLD = ConsensusConfig.CORRELATION_THRESHOLD
DEFAULT_AGGREGATION_METHOD = ConsensusConfig.AGGREGATION_METHOD
