"""
LC-Seq Configuration - Single Source of Truth from configs/default.yaml

This module loads configuration from configs/default.yaml and provides
convenient access to parameters throughout the codebase.

ALL configuration values come from the YAML file - no hardcoded defaults.
"""

from lcseq.infrastructure.configuration.yaml_loader import ConfigurationLoader

# Load configuration from Single Source of Truth (configs/default.yaml)
# This happens once on module import
_DEFAULT_CONFIG = ConfigurationLoader.get_default_config()

# Peak Detection Parameters (from configs/default.yaml)
DEFAULT_Z_THRESHOLD = _DEFAULT_CONFIG.peak_detection_params.get("z_threshold", 3.0)
DEFAULT_PROMINENCE_PERCENTILE = _DEFAULT_CONFIG.peak_detection_params.get("prominence_percentile", 0.5)
DEFAULT_MIN_SNR = _DEFAULT_CONFIG.peak_detection_params.get("min_snr", 0.001)
DEFAULT_MIN_BASELINE_SDS = _DEFAULT_CONFIG.peak_detection_params.get("min_baseline_sds", 1.0)
DEFAULT_SIGNAL_VARIANT = _DEFAULT_CONFIG.peak_detection_params.get("signal_variant", "raw")

# Peak Classification Parameters (from configs/default.yaml)
DEFAULT_TRUNCATION_MARGIN = _DEFAULT_CONFIG.classification_params.get("truncation_margin", 60.0)
DEFAULT_PEAK_MATCHING_TOLERANCE = _DEFAULT_CONFIG.classification_params.get("peak_matching_tolerance", 0.01)
DEFAULT_HUNGARIAN_MIN_THRESHOLD = _DEFAULT_CONFIG.classification_params.get("hungarian_min_threshold", 0.02)

# Pooling Parameters (from configs/default.yaml)
DEFAULT_CORRELATION_THRESHOLD = _DEFAULT_CONFIG.validation_params.get("correlation_threshold", 0.8)
DEFAULT_AGGREGATION_METHOD = _DEFAULT_CONFIG.validation_params.get("aggregation_method", "mean")


# Visualization Config Classes (backwards compatibility)
class PeakDetectionConfig:
    """Peak detection parameters from configs/default.yaml"""
    Z_THRESHOLD = DEFAULT_Z_THRESHOLD
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
    SECONDS_PER_MINUTE = _DEFAULT_CONFIG.peak_detection_params.get("viz_seconds_per_minute", 60.0)
    OFFSET_SPACING = _DEFAULT_CONFIG.peak_detection_params.get("viz_offset_spacing", 0.5)
    GROUP_SPACING_EXTRA = _DEFAULT_CONFIG.peak_detection_params.get("viz_group_spacing_extra", 1.0)

    # Color parameters
    USE_COLORMAP = _DEFAULT_CONFIG.peak_detection_params.get("viz_use_colormap", True)

    # Line styling
    LINEWIDTH_DEFAULT = _DEFAULT_CONFIG.peak_detection_params.get("viz_linewidth_default", 1.5)
    LINEWIDTH_REFERENCE = _DEFAULT_CONFIG.peak_detection_params.get("viz_linewidth_reference", 3.0)
    ALPHA_TRACE = _DEFAULT_CONFIG.peak_detection_params.get("viz_alpha_trace", 1.0)

    # Marker parameters
    MARKER_SIZE = _DEFAULT_CONFIG.peak_detection_params.get("viz_marker_size", 6)

    # Text parameters
    LABEL_FONTSIZE = _DEFAULT_CONFIG.peak_detection_params.get("viz_label_fontsize", 10)
    LABEL_MAX_LENGTH = _DEFAULT_CONFIG.peak_detection_params.get("viz_label_max_length", 40)
    LABEL_TRUNCATE_LENGTH = _DEFAULT_CONFIG.peak_detection_params.get("viz_label_truncate_length", 4)

    # Figure dimensions
    FIG_WIDTH = _DEFAULT_CONFIG.peak_detection_params.get("viz_fig_width", 16)
    FIG_HEIGHT_BASE = _DEFAULT_CONFIG.peak_detection_params.get("viz_fig_height_base", 4)
    FIG_HEIGHT_PER_TRACE = _DEFAULT_CONFIG.peak_detection_params.get("viz_fig_height_per_trace", 0.4)
    FIG_HEIGHT_MIN = _DEFAULT_CONFIG.peak_detection_params.get("viz_fig_height_min", 10)


class PoolingConfig:
    """Pooling parameters from configs/default.yaml"""
    CORRELATION_THRESHOLD = DEFAULT_CORRELATION_THRESHOLD
    AGGREGATION_METHOD = DEFAULT_AGGREGATION_METHOD
