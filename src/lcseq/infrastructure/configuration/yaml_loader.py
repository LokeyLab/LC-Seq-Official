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
        hierarchy_mode_str = analysis.get("hierarchy_mode", "block")

        # Convert to enums
        analysis_mode = AnalysisMode.INDIVIDUAL if analysis_mode_str == "individual" else AnalysisMode.CONSENSUS
        hierarchy_mode = HierarchyMode.BUILDING_BLOCK if hierarchy_mode_str == "block" else HierarchyMode.MONOMER

        # Extract peak detection parameters
        detection = config_dict.get("detection", {})
        peak_detection_params = {
            "min_persistence": detection.get("min_persistence", 0.05),
            "boundary_method": detection.get("boundary_method", "valley_or_5pct"),
        }

        # Extract validation parameters
        validation = config_dict.get("validation", {})
        validation_params = {
            "purity_threshold": validation.get("purity_threshold", "auto"),
            "snr_threshold": validation.get("snr_threshold", "auto"),
        }

        return AnalysisConfiguration(
            analysis_mode=analysis_mode,
            hierarchy_mode=hierarchy_mode,
            peak_detection_params=peak_detection_params,
            validation_params=validation_params,
        )

    @staticmethod
    def get_default_config() -> AnalysisConfiguration:
        """
        Get default analysis configuration.

        Returns
        -------
        AnalysisConfiguration
            Default configuration matching THEORY.md specifications
        """
        return AnalysisConfiguration.default()

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
        config_dict = {
            "analysis": {
                "variant_mode": config.analysis_mode.value,
                "hierarchy_mode": config.hierarchy_mode.value,
            },
            "detection": {
                "min_persistence": config.peak_detection_params.get("min_persistence", 0.05),
                "boundary_method": config.peak_detection_params.get("boundary_method", "valley_or_5pct"),
            },
            "validation": {
                "purity_threshold": config.validation_params.get("purity_threshold", "auto"),
                "snr_threshold": config.validation_params.get("snr_threshold", "auto"),
            },
        }

        with open(path, "w") as f:
            yaml.dump(config_dict, f, default_flow_style=False, sort_keys=False)
