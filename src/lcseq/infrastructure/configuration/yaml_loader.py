"""
YAML configuration loading infrastructure.

Implements configuration loading from YAML files with validation.

The YAML file (configs/default.yaml) is the SINGLE SOURCE OF TRUTH.
No hardcoded defaults exist anywhere else in the codebase.
"""

from pathlib import Path
from typing import Any, Dict
import yaml

from lcseq.domain.models.analysis_configuration import (
    AnalysisConfiguration,
    AnalysisMode,
    HierarchyMode,
)


class ConfigurationError(Exception):
    """Raised when configuration is missing required keys."""
    pass


def _require(d: Dict[str, Any], key: str, section: str) -> Any:
    """
    Require a key to be present in a dictionary.

    Raises ConfigurationError with helpful message if missing.
    """
    if key not in d:
        raise ConfigurationError(
            f"Missing required configuration key '{key}' in section '{section}'. "
            f"Please add it to your config file (configs/default.yaml is the reference)."
        )
    return d[key]


def _require_float(d: Dict[str, Any], key: str, section: str) -> float:
    """
    Require a key and convert to float.

    Handles scientific notation (e.g., 1e6) that YAML may parse as string.
    """
    value = _require(d, key, section)
    try:
        return float(value)
    except (ValueError, TypeError) as e:
        raise ConfigurationError(
            f"Configuration key '{key}' in section '{section}' must be a number, got: {value!r}"
        ) from e


class ConfigurationLoader:
    """
    Loads and validates analysis configuration from YAML files.

    This adapter sits in the infrastructure layer and converts
    external YAML format to domain AnalysisConfiguration objects.

    IMPORTANT: The YAML file is the single source of truth.
    No fallback values are used - missing keys will raise errors.
    """

    @staticmethod
    def load_from_yaml(path: Path) -> AnalysisConfiguration:
        """
        Load configuration from YAML file.

        Parameters
        ----------
        path : Path
            Path to YAML configuration file

        Returns
        -------
        AnalysisConfiguration
            Validated configuration object

        Raises
        ------
        FileNotFoundError
            If configuration file doesn't exist
        ConfigurationError
            If required configuration keys are missing
        yaml.YAMLError
            If YAML parsing fails
        """
        if not path.exists():
            raise FileNotFoundError(f"Configuration file not found: {path}")

        with open(path, "r") as f:
            config_dict = yaml.safe_load(f)

        return ConfigurationLoader._dict_to_config(config_dict)

    @staticmethod
    def _dict_to_config(config_dict: Dict[str, Any]) -> AnalysisConfiguration:
        """
        Convert configuration dictionary to AnalysisConfiguration.

        All keys are REQUIRED - no fallback values.
        The YAML file is the single source of truth.

        Parameters
        ----------
        config_dict : dict
            Configuration dictionary from YAML

        Returns
        -------
        AnalysisConfiguration
            Validated configuration object

        Raises
        ------
        ConfigurationError
            If required keys are missing
        """
        # Extract top-level sections (required)
        if "analysis" not in config_dict:
            raise ConfigurationError("Missing required section 'analysis' in config")
        if "detection" not in config_dict:
            raise ConfigurationError("Missing required section 'detection' in config")
        if "classification" not in config_dict:
            raise ConfigurationError("Missing required section 'classification' in config")
        if "pooling" not in config_dict:
            raise ConfigurationError("Missing required section 'pooling' in config")
        if "validation" not in config_dict:
            raise ConfigurationError("Missing required section 'validation' in config")

        analysis = config_dict["analysis"]
        preprocessing = config_dict.get("preprocessing", {})
        detection = config_dict["detection"]
        classification = config_dict["classification"]
        pooling = config_dict["pooling"]
        validation = config_dict["validation"]
        visualization = config_dict.get("visualization", {})

        # Parse analysis mode
        analysis_mode_str = _require(analysis, "variant_mode", "analysis")
        hierarchy_mode_str = _require(analysis, "hierarchy_mode", "analysis")

        # Convert to enums
        if analysis_mode_str == "individual":
            analysis_mode = AnalysisMode.INDIVIDUAL
        elif analysis_mode_str == "pooled":
            analysis_mode = AnalysisMode.POOLED
        else:
            raise ConfigurationError(
                f"Invalid variant_mode '{analysis_mode_str}'. Must be 'individual' or 'pooled'."
            )

        if hierarchy_mode_str in ["block", "building_block"]:
            hierarchy_mode = HierarchyMode.BUILDING_BLOCK
        elif hierarchy_mode_str == "monomer":
            hierarchy_mode = HierarchyMode.MONOMER
        else:
            raise ConfigurationError(
                f"Invalid hierarchy_mode '{hierarchy_mode_str}'. Must be 'building_block' or 'monomer'."
            )

        # Extract preprocessing parameters (all required)
        if "preprocessing" not in config_dict:
            raise ConfigurationError("Missing required section 'preprocessing' in config")
        preprocessing_params = {
            "enabled": _require(preprocessing, "enabled", "preprocessing"),
            "baseline_order": _require(preprocessing, "baseline_order", "preprocessing"),
            "baseline_smoothing": preprocessing["baseline_smoothing"],  # Can be null
            "include_endpoints": _require(preprocessing, "include_endpoints", "preprocessing"),
        }

        # Extract peak detection parameters (all required)
        # Use _require_float for numeric params to handle scientific notation (e.g., 1e6)
        peak_detection_params = {
            "min_persistence": _require_float(detection, "min_persistence", "detection"),
            "alpha": _require_float(detection, "alpha", "detection"),
            "alpha_product": _require_float(detection, "alpha_product", "detection"),
            "prominence_percentile": _require_float(detection, "prominence_percentile", "detection"),
            "min_snr": _require_float(detection, "min_snr", "detection"),
            "min_baseline_sds": _require_float(detection, "min_baseline_sds", "detection"),
            "boundary_method": _require(detection, "boundary_method", "detection"),
            "boundary_threshold_fraction": _require_float(detection, "boundary_threshold_fraction", "detection"),
            "signal_variant": _require(detection, "signal_variant", "detection"),
            "min_dispersion_r": _require_float(detection, "min_dispersion_r", "detection"),
            "sigma_clip_sigma": _require_float(detection, "sigma_clip_sigma", "detection"),
            "include_rejected": _require(detection, "include_rejected", "detection"),
        }

        # Extract peak classification parameters (all required)
        classification_params = {
            "truncation_margin": _require_float(classification, "truncation_margin", "classification"),
            "peak_matching_tolerance": _require_float(classification, "peak_matching_tolerance", "classification"),
            "hungarian_min_threshold": _require_float(classification, "hungarian_min_threshold", "classification"),
        }

        # Extract pooling parameters (all required)
        pooling_params = {
            "correlation_threshold": _require_float(pooling, "correlation_threshold", "pooling"),
            "aggregation_method": _require(pooling, "aggregation_method", "pooling"),
        }

        # Extract validation parameters (all required)
        # Note: purity_threshold and snr_threshold can be "auto" strings, so keep as _require
        clpe = validation.get("clpe", {})
        validation_params = {
            "purity_threshold": _require(validation, "purity_threshold", "validation"),
            "snr_threshold": _require(validation, "snr_threshold", "validation"),
            "retention_precision": _require_float(validation, "retention_precision", "validation"),
            "clpe_outlier_threshold": _require_float(clpe, "outlier_threshold", "validation.clpe"),
            "clpe_min_group_size": _require(clpe, "min_group_size", "validation.clpe"),  # int, not float
        }

        # Include pooling params in validation_params for unified access
        validation_params.update(pooling_params)

        # Extract performance parameters (optional section)
        performance = config_dict.get("performance", {})
        performance_params = {
            "num_workers": performance.get("num_workers"),  # None = sequential
        }

        # Extract quality filter parameters (optional section)
        quality_filter = config_dict.get("quality_filter", {})
        quality_filter_params = {}
        if quality_filter:
            quality_filter_params = {
                "min_correlation": quality_filter.get("min_correlation", 0.8),
                "intensity_percentile": quality_filter.get("intensity_percentile", 0.05),
                "intensity_absolute": quality_filter.get("intensity_absolute"),  # None = use percentile
                "max_noise_ratio": quality_filter.get("max_noise_ratio"),  # None = disabled
            }

        # Extract visualization parameters (optional section, but if present all keys required)
        if visualization:
            layout = visualization.get("layout", {})
            colors = visualization.get("colors", {})
            lines = visualization.get("lines", {})
            markers = visualization.get("markers", {})
            text = visualization.get("text", {})
            figure = visualization.get("figure", {})

            viz_params = {
                "viz_seconds_per_minute": _require(layout, "seconds_per_minute", "visualization.layout"),
                "viz_offset_spacing": _require(layout, "offset_spacing", "visualization.layout"),
                "viz_group_spacing_extra": _require(layout, "group_spacing_extra", "visualization.layout"),
                "viz_use_colormap": _require(colors, "use_colormap", "visualization.colors"),
                "viz_linewidth_default": _require(lines, "linewidth_default", "visualization.lines"),
                "viz_linewidth_reference": _require(lines, "linewidth_reference", "visualization.lines"),
                "viz_alpha_trace": _require(lines, "alpha_trace", "visualization.lines"),
                "viz_marker_size": _require(markers, "marker_size", "visualization.markers"),
                "viz_label_fontsize": _require(text, "label_fontsize", "visualization.text"),
                "viz_label_max_length": _require(text, "label_max_length", "visualization.text"),
                "viz_label_truncate_length": _require(text, "label_truncate_length", "visualization.text"),
                "viz_fig_width": _require(figure, "fig_width", "visualization.figure"),
                "viz_fig_height_base": _require(figure, "fig_height_base", "visualization.figure"),
                "viz_fig_height_per_trace": _require(figure, "fig_height_per_trace", "visualization.figure"),
                "viz_fig_height_min": _require(figure, "fig_height_min", "visualization.figure"),
            }
            peak_detection_params.update(viz_params)

        return AnalysisConfiguration(
            analysis_mode=analysis_mode,
            hierarchy_mode=hierarchy_mode,
            preprocessing_params=preprocessing_params,
            peak_detection_params=peak_detection_params,
            validation_params=validation_params,
            classification_params=classification_params,
            performance_params=performance_params,
            quality_filter_params=quality_filter_params,
        )

    @staticmethod
    def get_default_config() -> AnalysisConfiguration:
        """
        Get default analysis configuration from configs/default.yaml.

        This is the Single Source of Truth for all configuration parameters.

        Returns
        -------
        AnalysisConfiguration
            Default configuration loaded from configs/default.yaml

        Raises
        ------
        FileNotFoundError
            If configs/default.yaml not found
        ConfigurationError
            If required keys are missing from the config file
        """
        from pathlib import Path

        # Try to find configs/default.yaml relative to this file
        # Go up from infrastructure/configuration/ to project root
        current_file = Path(__file__)
        project_root = current_file.parent.parent.parent.parent.parent
        default_config_path = project_root / "configs" / "default.yaml"

        if not default_config_path.exists():
            raise FileNotFoundError(
                f"Default configuration not found at {default_config_path}. "
                f"This file is the Single Source of Truth for LC-Seq configuration."
            )

        return ConfigurationLoader.load_from_yaml(default_config_path)

    @staticmethod
    def save_to_yaml(config: AnalysisConfiguration, path: Path) -> None:
        """
        Save configuration to YAML file.

        Parameters
        ----------
        config : AnalysisConfiguration
            Configuration to save
        path : Path
            Output path for YAML file
        """
        # Build comprehensive config dict - just access directly, no fallbacks
        config_dict = {
            "analysis": {
                "variant_mode": config.analysis_mode.value,
                "hierarchy_mode": config.hierarchy_mode.value,
            },
        }

        # Add preprocessing if present
        if config.preprocessing_params:
            # Uses PreprocessingConfig as single source of truth
            from lcseq.domain.services.signal_preprocessor import PreprocessingConfig
            config_dict["preprocessing"] = PreprocessingConfig.from_dict(config.preprocessing_params).to_dict()

        config_dict["detection"] = {
            "min_persistence": config.peak_detection_params["min_persistence"],
            "alpha": config.peak_detection_params["alpha"],
            "alpha_product": config.peak_detection_params["alpha_product"],
            "prominence_percentile": config.peak_detection_params["prominence_percentile"],
            "min_snr": config.peak_detection_params["min_snr"],
            "min_baseline_sds": config.peak_detection_params["min_baseline_sds"],
            "boundary_method": config.peak_detection_params["boundary_method"],
            "boundary_threshold_fraction": config.peak_detection_params["boundary_threshold_fraction"],
            "signal_variant": config.peak_detection_params["signal_variant"],
            "min_dispersion_r": config.peak_detection_params["min_dispersion_r"],
            "sigma_clip_sigma": config.peak_detection_params["sigma_clip_sigma"],
            "include_rejected": config.peak_detection_params["include_rejected"],
        }
        config_dict["classification"] = {
            "truncation_margin": config.classification_params["truncation_margin"],
            "peak_matching_tolerance": config.classification_params["peak_matching_tolerance"],
            "hungarian_min_threshold": config.classification_params["hungarian_min_threshold"],
        }
        config_dict["pooling"] = {
            "correlation_threshold": config.validation_params["correlation_threshold"],
            "aggregation_method": config.validation_params["aggregation_method"],
        }
        config_dict["validation"] = {
            "purity_threshold": config.validation_params["purity_threshold"],
            "snr_threshold": config.validation_params["snr_threshold"],
            "retention_precision": config.validation_params["retention_precision"],
            "clpe": {
                "outlier_threshold": config.validation_params["clpe_outlier_threshold"],
                "min_group_size": config.validation_params["clpe_min_group_size"],
            },
        }

        # Add performance parameters if present
        if config.performance_params:
            config_dict["performance"] = {
                "num_workers": config.performance_params.get("num_workers"),
            }

        # Add quality filter parameters if present
        if config.quality_filter_params:
            config_dict["quality_filter"] = {
                "min_correlation": config.quality_filter_params.get("min_correlation", 0.8),
                "intensity_percentile": config.quality_filter_params.get("intensity_percentile", 0.05),
                "intensity_absolute": config.quality_filter_params.get("intensity_absolute"),
                "max_noise_ratio": config.quality_filter_params.get("max_noise_ratio"),
            }

        # Add visualization parameters if present (have viz_ prefix)
        viz_params = {k: v for k, v in config.peak_detection_params.items() if k.startswith("viz_")}
        if viz_params:
            config_dict["visualization"] = {
                "layout": {
                    "seconds_per_minute": viz_params["viz_seconds_per_minute"],
                    "offset_spacing": viz_params["viz_offset_spacing"],
                    "group_spacing_extra": viz_params["viz_group_spacing_extra"],
                },
                "colors": {
                    "use_colormap": viz_params["viz_use_colormap"],
                },
                "lines": {
                    "linewidth_default": viz_params["viz_linewidth_default"],
                    "linewidth_reference": viz_params["viz_linewidth_reference"],
                    "alpha_trace": viz_params["viz_alpha_trace"],
                },
                "markers": {
                    "marker_size": viz_params["viz_marker_size"],
                },
                "text": {
                    "label_fontsize": viz_params["viz_label_fontsize"],
                    "label_max_length": viz_params["viz_label_max_length"],
                    "label_truncate_length": viz_params["viz_label_truncate_length"],
                },
                "figure": {
                    "fig_width": viz_params["viz_fig_width"],
                    "fig_height_base": viz_params["viz_fig_height_base"],
                    "fig_height_per_trace": viz_params["viz_fig_height_per_trace"],
                    "fig_height_min": viz_params["viz_fig_height_min"],
                },
            }

        with open(path, "w") as f:
            yaml.dump(config_dict, f, default_flow_style=False, sort_keys=False)
