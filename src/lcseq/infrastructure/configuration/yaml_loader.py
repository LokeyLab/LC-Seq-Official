"""
YAML configuration loading infrastructure.

Implements configuration loading from YAML files with validation.
"""

from pathlib import Path
from typing import Any, Dict
import yaml

from lcseq.domain.models.analysis_configuration import (
    AnalysisConfiguration,
    AnalysisMode,
    HierarchyMode,
)


class ConfigurationLoader:
    """
    Loads and validates analysis configuration from YAML files.

    This adapter sits in the infrastructure layer and converts
    external YAML format to domain AnalysisConfiguration objects.
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
        ValueError
            If configuration is invalid
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

        Parameters
        ----------
        config_dict : dict
            Configuration dictionary from YAML

        Returns
        -------
        AnalysisConfiguration
            Validated configuration object
        """
        # Extract top-level settings
        analysis = config_dict.get("analysis", {})
        analysis_mode_str = analysis.get("variant_mode", "individual")
        hierarchy_mode_str = analysis.get("hierarchy_mode", "building_block")

        # Convert to enums
        analysis_mode = AnalysisMode.INDIVIDUAL if analysis_mode_str == "individual" else AnalysisMode.POOLED
        hierarchy_mode = HierarchyMode.BUILDING_BLOCK if hierarchy_mode_str in ["block", "building_block"] else HierarchyMode.MONOMER

        # Extract peak detection parameters (comprehensive)
        detection = config_dict.get("detection", {})
        peak_detection_params = {
            "min_persistence": detection.get("min_persistence", 0.05),
            "z_threshold": detection.get("z_threshold", 3.0),
            "prominence_percentile": detection.get("prominence_percentile", 0.5),
            "min_snr": detection.get("min_snr", 0.001),
            "min_baseline_sds": detection.get("min_baseline_sds", 1.0),
            "boundary_method": detection.get("boundary_method", "valley_or_5pct"),
            "signal_variant": detection.get("signal_variant", "raw"),
        }

        # Extract peak classification parameters
        classification = config_dict.get("classification", {})
        classification_params = {
            "truncation_margin": classification.get("truncation_margin", 60.0),
            "peak_matching_tolerance": classification.get("peak_matching_tolerance", 0.01),
            "hungarian_min_threshold": classification.get("hungarian_min_threshold", 0.02),
        }

        # Extract pooling parameters
        pooling = config_dict.get("pooling", {})
        pooling_params = {
            "correlation_threshold": pooling.get("correlation_threshold", 0.8),
            "aggregation_method": pooling.get("aggregation_method", "mean"),
        }

        # Extract validation parameters
        validation = config_dict.get("validation", {})
        validation_params = {
            "purity_threshold": validation.get("purity_threshold", "auto"),
            "snr_threshold": validation.get("snr_threshold", "auto"),
        }

        # Merge pooling params into validation_params for backwards compatibility
        # (correlation_threshold was historically in validation_params)
        validation_params.update(pooling_params)

        # Extract visualization parameters (stored in peak_detection_params for now)
        visualization = config_dict.get("visualization", {})
        if visualization:
            # Store visualization config in peak_detection_params with "viz_" prefix
            # This keeps it accessible without changing AnalysisConfiguration structure
            layout = visualization.get("layout", {})
            colors = visualization.get("colors", {})
            lines = visualization.get("lines", {})
            markers = visualization.get("markers", {})
            text = visualization.get("text", {})
            figure = visualization.get("figure", {})

            viz_params = {
                "viz_seconds_per_minute": layout.get("seconds_per_minute", 60.0),
                "viz_offset_spacing": layout.get("offset_spacing", 0.5),
                "viz_group_spacing_extra": layout.get("group_spacing_extra", 1.0),
                "viz_use_colormap": colors.get("use_colormap", True),
                "viz_linewidth_default": lines.get("linewidth_default", 1.5),
                "viz_linewidth_reference": lines.get("linewidth_reference", 3.0),
                "viz_alpha_trace": lines.get("alpha_trace", 1.0),
                "viz_marker_size": markers.get("marker_size", 6),
                "viz_label_fontsize": text.get("label_fontsize", 10),
                "viz_label_max_length": text.get("label_max_length", 40),
                "viz_label_truncate_length": text.get("label_truncate_length", 4),
                "viz_fig_width": figure.get("fig_width", 16),
                "viz_fig_height_base": figure.get("fig_height_base", 4),
                "viz_fig_height_per_trace": figure.get("fig_height_per_trace", 0.4),
                "viz_fig_height_min": figure.get("fig_height_min", 10),
            }
            peak_detection_params.update(viz_params)

        return AnalysisConfiguration(
            analysis_mode=analysis_mode,
            hierarchy_mode=hierarchy_mode,
            peak_detection_params=peak_detection_params,
            validation_params=validation_params,
            classification_params=classification_params,
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
        # Build comprehensive config dict
        config_dict = {
            "analysis": {
                "variant_mode": config.analysis_mode.value,
                "hierarchy_mode": config.hierarchy_mode.value,
            },
            "detection": {
                "min_persistence": config.peak_detection_params.get("min_persistence", 0.05),
                "z_threshold": config.peak_detection_params.get("z_threshold", 3.0),
                "prominence_percentile": config.peak_detection_params.get("prominence_percentile", 0.5),
                "min_snr": config.peak_detection_params.get("min_snr", 0.001),
                "min_baseline_sds": config.peak_detection_params.get("min_baseline_sds", 1.0),
                "boundary_method": config.peak_detection_params.get("boundary_method", "valley_or_5pct"),
                "signal_variant": config.peak_detection_params.get("signal_variant", "raw"),
            },
            "classification": {
                "truncation_margin": config.classification_params.get("truncation_margin", 60.0),
                "peak_matching_tolerance": config.classification_params.get("peak_matching_tolerance", 0.01),
                "hungarian_min_threshold": config.classification_params.get("hungarian_min_threshold", 0.02),
            },
            "pooling": {
                "correlation_threshold": config.validation_params.get("correlation_threshold", 0.8),
                "aggregation_method": config.validation_params.get("aggregation_method", "mean"),
            },
            "validation": {
                "purity_threshold": config.validation_params.get("purity_threshold", "auto"),
                "snr_threshold": config.validation_params.get("snr_threshold", "auto"),
            },
        }

        # Add visualization parameters if present (have viz_ prefix)
        viz_params = {k: v for k, v in config.peak_detection_params.items() if k.startswith("viz_")}
        if viz_params:
            config_dict["visualization"] = {
                "layout": {
                    "seconds_per_minute": viz_params.get("viz_seconds_per_minute", 60.0),
                    "offset_spacing": viz_params.get("viz_offset_spacing", 0.5),
                    "group_spacing_extra": viz_params.get("viz_group_spacing_extra", 1.0),
                },
                "colors": {
                    "use_colormap": viz_params.get("viz_use_colormap", True),
                },
                "lines": {
                    "linewidth_default": viz_params.get("viz_linewidth_default", 1.5),
                    "linewidth_reference": viz_params.get("viz_linewidth_reference", 3.0),
                    "alpha_trace": viz_params.get("viz_alpha_trace", 1.0),
                },
                "markers": {
                    "marker_size": viz_params.get("viz_marker_size", 6),
                },
                "text": {
                    "label_fontsize": viz_params.get("viz_label_fontsize", 10),
                    "label_max_length": viz_params.get("viz_label_max_length", 40),
                    "label_truncate_length": viz_params.get("viz_label_truncate_length", 4),
                },
                "figure": {
                    "fig_width": viz_params.get("viz_fig_width", 16),
                    "fig_height_base": viz_params.get("viz_fig_height_base", 4),
                    "fig_height_per_trace": viz_params.get("viz_fig_height_per_trace", 0.4),
                    "fig_height_min": viz_params.get("viz_fig_height_min", 10),
                },
            }

        with open(path, "w") as f:
            yaml.dump(config_dict, f, default_flow_style=False, sort_keys=False)
