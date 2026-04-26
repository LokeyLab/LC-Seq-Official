"""
AnalysisConfiguration - configuration for analysis run.

Implementation based on THEORY.md Part 5, 6.
"""

from dataclasses import dataclass, field
from typing import Dict, Any, Optional
from enum import Enum
from .compound_hierarchy import HierarchyMode


class AnalysisMode(Enum):
    """
    Analysis mode for compound processing.

    Attributes
    ----------
    INDIVIDUAL : str
        Process each positional variant independently (default)
    POOLED : str
        Aggregate positional variants by equivalence class (optional optimization)
    """

    INDIVIDUAL = "individual"
    POOLED = "pooled"


@dataclass
class AnalysisConfiguration:
    """
    Configuration parameters for LC-Seq analysis run.

    Stores all parameters for preprocessing, peak detection, classification,
    and validation. Ensures consistent parameter application across entire dataset.

    Attributes
    ----------
    analysis_mode : AnalysisMode
        Individual or pooled mode
    hierarchy_mode : HierarchyMode
        Building-block or monomer hierarchy
    peak_detection_params : Dict[str, Any]
        Morse theory + persistence parameters
    validation_params : Dict[str, Any]
        Bayesian validation framework parameters
    classification_params : Dict[str, Any]
        Peak classification parameters (optional)

    Notes
    -----
    - All parameters applied uniformly across dataset
    - Peak detection operates directly on raw signals
    - Validation thresholds: adaptive (computed from dataset)
    - Mode selection affects graph structure and processing

    Examples
    --------
    >>> # Standard configuration
    >>> config = AnalysisConfiguration.default()
    >>> config.analysis_mode
    <AnalysisMode.INDIVIDUAL: 'individual'>

    >>> # Pooled mode configuration
    >>> config_pooled = AnalysisConfiguration.default()
    >>> config_pooled.analysis_mode = AnalysisMode.POOLED
    >>> config_pooled.validation_params['correlation_threshold'] = 0.8

    References
    ----------
    THEORY.md Section 5.2: Persistent Homology Parameters
    THEORY.md Section 6.4: Adaptive Thresholds
    THEORY.md Section 4.2.2: Analysis Modes
    """

    analysis_mode: AnalysisMode = AnalysisMode.INDIVIDUAL
    hierarchy_mode: HierarchyMode = HierarchyMode.BUILDING_BLOCK
    preprocessing_params: Dict[str, Any] = field(default_factory=dict)
    peak_detection_params: Dict[str, Any] = field(default_factory=dict)
    validation_params: Dict[str, Any] = field(default_factory=dict)
    classification_params: Dict[str, Any] = field(default_factory=dict)
    performance_params: Dict[str, Any] = field(default_factory=dict)
    quality_filter_params: Dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        """Validate and set default parameters."""
        # Set peak detection defaults if not provided (merge with existing)
        peak_defaults = {
            "min_persistence": 0.05,
            "boundary_method": "valley_or_5pct",
        }
        for key, value in peak_defaults.items():
            if key not in self.peak_detection_params:
                self.peak_detection_params[key] = value

        # Validate persistence threshold
        if "min_persistence" in self.peak_detection_params:
            min_pers = self.peak_detection_params["min_persistence"]
            if not 0.0 <= min_pers <= 1.0:
                raise ValueError(
                    f"min_persistence must be in [0, 1], got {min_pers}"
                )

        # Set validation defaults if not provided (merge with existing)
        validation_defaults = {
            "purity_threshold": "auto",  # Use P₂₅ from dataset
            "snr_threshold": "auto",  # Use P₅₀ from dataset
        }
        for key, value in validation_defaults.items():
            if key not in self.validation_params:
                self.validation_params[key] = value

        # Add correlation threshold for pooled mode
        if self.analysis_mode == AnalysisMode.POOLED:
            if "correlation_threshold" not in self.validation_params:
                self.validation_params["correlation_threshold"] = 0.8

            # Validate correlation threshold
            corr_thresh = self.validation_params["correlation_threshold"]
            if not 0.0 <= corr_thresh <= 1.0:
                raise ValueError(
                    f"correlation_threshold must be in [0, 1], got {corr_thresh}"
                )

    @classmethod
    def default(cls) -> "AnalysisConfiguration":
        """
        Create configuration with standard default parameters.

        Returns
        -------
        AnalysisConfiguration
            Configuration with standard defaults

        Notes
        -----
        Standard parameters from THEORY.md:
        - Persistence: min=0.05
        - Validation: adaptive thresholds
        - Mode: individual, building-block

        Examples
        --------
        >>> config = AnalysisConfiguration.default()
        >>> config.peak_detection_params['min_persistence']
        0.05
        """
        return cls()

    @classmethod
    def for_pooled_mode(
        cls, correlation_threshold: float = 0.8
    ) -> "AnalysisConfiguration":
        """
        Create configuration for pooled mode analysis.

        Parameters
        ----------
        correlation_threshold : float
            Minimum correlation for pooling validity (default 0.8)

        Returns
        -------
        AnalysisConfiguration
            Configuration with pooled mode enabled

        Notes
        -----
        Pooled mode requires high correlation between positional variants
        (THEORY.md Section 4.2.4).

        Examples
        --------
        >>> config = AnalysisConfiguration.for_pooled_mode(correlation_threshold=0.85)
        >>> config.analysis_mode
        <AnalysisMode.POOLED: 'pooled'>
        >>> config.validation_params['correlation_threshold']
        0.85
        """
        config = cls(analysis_mode=AnalysisMode.POOLED)
        config.validation_params["correlation_threshold"] = correlation_threshold
        return config

    @classmethod
    def for_monomer_mode(cls) -> "AnalysisConfiguration":
        """
        Create configuration for monomer-level hierarchy.

        Returns
        -------
        AnalysisConfiguration
            Configuration with monomer hierarchy mode

        Notes
        -----
        Monomer mode creates DAG with convergence patterns
        (THEORY.md Section 1.5.4).

        Examples
        --------
        >>> config = AnalysisConfiguration.for_monomer_mode()
        >>> config.hierarchy_mode
        <HierarchyMode.MONOMER: 'monomer'>
        """
        return cls(hierarchy_mode=HierarchyMode.MONOMER)

    def copy(self) -> "AnalysisConfiguration":
        """
        Create a copy of this configuration.

        Returns
        -------
        AnalysisConfiguration
            Deep copy of configuration
        """
        return AnalysisConfiguration(
            analysis_mode=self.analysis_mode,
            hierarchy_mode=self.hierarchy_mode,
            preprocessing_params=self.preprocessing_params.copy(),
            peak_detection_params=self.peak_detection_params.copy(),
            validation_params=self.validation_params.copy(),
            classification_params=self.classification_params.copy(),
            performance_params=self.performance_params.copy(),
            quality_filter_params=self.quality_filter_params.copy(),
        )

    def __repr__(self) -> str:
        """Detailed representation for debugging."""
        return (
            f"AnalysisConfiguration("
            f"analysis_mode={self.analysis_mode.value}, "
            f"hierarchy_mode={self.hierarchy_mode.value}, "
            f"min_persistence={self.peak_detection_params.get('min_persistence', 0.05)})"
        )
