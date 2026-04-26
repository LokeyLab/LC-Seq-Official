"""
LC-Seq Configuration - Single Source of Truth from configs/default.yaml

This module loads configuration from configs/default.yaml and provides
convenient access to parameters throughout the codebase.

ALL configuration values come from the YAML file - no hardcoded defaults.
If a key is missing from the YAML, a ConfigurationError will be raised.
"""

# Direct import to avoid circular dependency through infrastructure/__init__.py
# The circular import chain was:
#   config -> infrastructure -> exporters -> application -> config
from lcseq.infrastructure.configuration.yaml_loader import ConfigurationLoader

# Load configuration from Single Source of Truth (configs/default.yaml)
# This happens once on module import
_DEFAULT_CONFIG = ConfigurationLoader.get_default_config()

# Peak Detection Parameters (from configs/default.yaml - no fallbacks)
DEFAULT_ALPHA = _DEFAULT_CONFIG.peak_detection_params["alpha"]
DEFAULT_ALPHA_PRODUCT = _DEFAULT_CONFIG.peak_detection_params["alpha_product"]
DEFAULT_PROMINENCE_PERCENTILE = _DEFAULT_CONFIG.peak_detection_params["prominence_percentile"]
DEFAULT_MIN_SNR = _DEFAULT_CONFIG.peak_detection_params["min_snr"]
DEFAULT_MIN_BASELINE_SDS = _DEFAULT_CONFIG.peak_detection_params["min_baseline_sds"]
DEFAULT_SIGNAL_VARIANT = _DEFAULT_CONFIG.peak_detection_params["signal_variant"]
DEFAULT_MIN_DISPERSION_R = _DEFAULT_CONFIG.peak_detection_params["min_dispersion_r"]
DEFAULT_BOUNDARY_THRESHOLD_FRACTION = _DEFAULT_CONFIG.peak_detection_params["boundary_threshold_fraction"]

# Peak Classification Parameters (from configs/default.yaml - no fallbacks)
DEFAULT_TRUNCATION_MARGIN = _DEFAULT_CONFIG.classification_params["truncation_margin"]
DEFAULT_PEAK_MATCHING_TOLERANCE = _DEFAULT_CONFIG.classification_params["peak_matching_tolerance"]
DEFAULT_HUNGARIAN_MIN_THRESHOLD = _DEFAULT_CONFIG.classification_params["hungarian_min_threshold"]

# Pooling Parameters (from configs/default.yaml - no fallbacks)
DEFAULT_CORRELATION_THRESHOLD = _DEFAULT_CONFIG.validation_params["correlation_threshold"]
DEFAULT_AGGREGATION_METHOD = _DEFAULT_CONFIG.validation_params["aggregation_method"]


# Configuration Classes for convenient access
class PeakDetectionConfig:
    """Peak detection parameters from configs/default.yaml"""
    ALPHA = DEFAULT_ALPHA
    ALPHA_PRODUCT = DEFAULT_ALPHA_PRODUCT
    PROMINENCE_PERCENTILE = DEFAULT_PROMINENCE_PERCENTILE
    MIN_SNR = DEFAULT_MIN_SNR
    MIN_BASELINE_SDS = DEFAULT_MIN_BASELINE_SDS
    SIGNAL_VARIANT = DEFAULT_SIGNAL_VARIANT
    TRUNCATION_MARGIN = DEFAULT_TRUNCATION_MARGIN
    PEAK_MATCHING_TOLERANCE = DEFAULT_PEAK_MATCHING_TOLERANCE
    HUNGARIAN_MIN_THRESHOLD = DEFAULT_HUNGARIAN_MIN_THRESHOLD


class VisualizationConfig:
    """Visualization parameters from configs/default.yaml"""
    # Layout parameters
    SECONDS_PER_MINUTE = _DEFAULT_CONFIG.peak_detection_params["viz_seconds_per_minute"]
    OFFSET_SPACING = _DEFAULT_CONFIG.peak_detection_params["viz_offset_spacing"]
    GROUP_SPACING_EXTRA = _DEFAULT_CONFIG.peak_detection_params["viz_group_spacing_extra"]

    # Color parameters
    USE_COLORMAP = _DEFAULT_CONFIG.peak_detection_params["viz_use_colormap"]

    # Line styling
    LINEWIDTH_DEFAULT = _DEFAULT_CONFIG.peak_detection_params["viz_linewidth_default"]
    LINEWIDTH_REFERENCE = _DEFAULT_CONFIG.peak_detection_params["viz_linewidth_reference"]
    ALPHA_TRACE = _DEFAULT_CONFIG.peak_detection_params["viz_alpha_trace"]

    # Marker parameters
    MARKER_SIZE = _DEFAULT_CONFIG.peak_detection_params["viz_marker_size"]

    # Text parameters
    LABEL_FONTSIZE = _DEFAULT_CONFIG.peak_detection_params["viz_label_fontsize"]
    LABEL_MAX_LENGTH = _DEFAULT_CONFIG.peak_detection_params["viz_label_max_length"]
    LABEL_TRUNCATE_LENGTH = _DEFAULT_CONFIG.peak_detection_params["viz_label_truncate_length"]

    # Figure dimensions
    FIG_WIDTH = _DEFAULT_CONFIG.peak_detection_params["viz_fig_width"]
    FIG_HEIGHT_BASE = _DEFAULT_CONFIG.peak_detection_params["viz_fig_height_base"]
    FIG_HEIGHT_PER_TRACE = _DEFAULT_CONFIG.peak_detection_params["viz_fig_height_per_trace"]
    FIG_HEIGHT_MIN = _DEFAULT_CONFIG.peak_detection_params["viz_fig_height_min"]


class PeakAppearanceConfig:
    """
    Peak appearance configuration for consistent visualization across plotters.

    Centralizes peak colors, markers, and legend creation to eliminate DRY violations.
    All plotters should use this config instead of defining their own PEAK_COLORS dicts.
    """
    from lcseq.domain.entities.peak import PeakType
    from matplotlib.lines import Line2D
    from typing import Dict, Tuple, List

    # Peak type colors
    COLORS: Dict[PeakType, str] = {
        PeakType.NULL: "gray",
        PeakType.TRUNCATION: "orange",
        PeakType.TRUNCATION_UNKNOWN: "darkorange",
        PeakType.PUTATIVE_PRODUCT: "green",
        PeakType.UNKNOWN: "red",
    }

    # Peak type markers: (marker_char, size_multiplier)
    MARKERS: Dict[PeakType, Tuple[str, float]] = {
        PeakType.NULL: ('s', 1.5),           # square
        PeakType.TRUNCATION: ('D', 1.5),     # diamond
        PeakType.TRUNCATION_UNKNOWN: ('d', 1.5),  # small diamond
        PeakType.PUTATIVE_PRODUCT: ('o', 2.0),    # circle (larger)
        PeakType.UNKNOWN: ('^', 1.5),        # triangle
    }

    @classmethod
    def get_color(cls, peak_type: 'PeakType') -> str:
        """Get color for a peak type."""
        return cls.COLORS.get(peak_type, "gray")

    @classmethod
    def get_marker(cls, peak_type: 'PeakType') -> Tuple[str, float]:
        """
        Get marker shape and size for a peak type.

        Returns
        -------
        Tuple[str, float]
            (marker_char, size) where size = base_size * size_multiplier
        """
        marker_char, multiplier = cls.MARKERS.get(peak_type, ('o', 1.0))
        size = VisualizationConfig.MARKER_SIZE * multiplier
        return marker_char, size

    @classmethod
    def create_legend_handles(cls, include_rejected: bool = True) -> 'List[Line2D]':
        """
        Create legend handles for peak markers.

        Parameters
        ----------
        include_rejected : bool
            Whether to include the "Rejected" marker in legend

        Returns
        -------
        List[Line2D]
            List of Line2D objects for legend
        """
        handles = [
            cls.Line2D([0], [0], marker='s', color='gray', linestyle='None',
                   markersize=7, markeredgecolor='white', markeredgewidth=1,
                   label='NULL'),
            cls.Line2D([0], [0], marker='D', color='orange', linestyle='None',
                   markersize=7, markeredgecolor='white', markeredgewidth=1,
                   label='Truncation'),
            cls.Line2D([0], [0], marker='d', color='darkorange', linestyle='None',
                   markersize=7, markeredgecolor='white', markeredgewidth=1,
                   label='Trunc. (unknown)'),
            cls.Line2D([0], [0], marker='o', color='green', linestyle='None',
                   markersize=9, markeredgecolor='white', markeredgewidth=1,
                   label='Product'),
            cls.Line2D([0], [0], marker='^', color='red', linestyle='None',
                   markersize=7, markeredgecolor='white', markeredgewidth=1,
                   label='Unknown'),
        ]

        if include_rejected:
            handles.append(
                cls.Line2D([0], [0], marker='o', color='gray', linestyle='None',
                       markersize=7, markerfacecolor='none', markeredgewidth=1.5,
                       label='Rejected')
            )

        return handles

    @classmethod
    def create_linestyle_handles(cls, include_truncation: bool = True) -> 'List[Line2D]':
        """
        Create legend handles for signal line styles.

        Parameters
        ----------
        include_truncation : bool
            Whether to include truncation-specific line styles

        Returns
        -------
        List[Line2D]
            List of Line2D objects for legend
        """
        handles = [
            cls.Line2D([0], [0], color='black', linestyle='-', linewidth=2,
                   label='Valid (above threshold)'),
        ]

        if include_truncation:
            handles.extend([
                cls.Line2D([0], [0], color='black', linestyle='--', linewidth=2,
                       label='Truncation region'),
                cls.Line2D([0], [0], color='black', linestyle='-.', linewidth=2,
                       label='Truncation + below'),
            ])

        handles.append(
            cls.Line2D([0], [0], color='black', linestyle=':', linewidth=2,
                   label='Below threshold')
        )

        return handles
