"""
Signal preprocessing service for chromatogram data.

Applies baseline correction to prepare signals for peak detection.
"""

from dataclasses import dataclass
from typing import Optional

import numpy as np
from numpy.typing import NDArray

from lcseq.domain.entities.chromatogram import Chromatogram
from lcseq.domain.services.baseline_estimator import (
    MinimaSplineParams,
    minima_spline_baseline,
)


@dataclass
class PreprocessingConfig:
    """
    Configuration for signal preprocessing.

    Attributes
    ----------
    enabled : bool
        Whether preprocessing is enabled.
    baseline_order : int
        Order for local minima detection (higher = fewer minima).
    baseline_smoothing : float or None
        Spline smoothing factor (None = auto).
    include_start : bool
        Include start point as virtual minimum for baseline fitting.
    include_end : bool
        Include end point as virtual minimum for baseline fitting.
    """

    enabled: bool = True
    baseline_order: int = 3
    baseline_smoothing: Optional[float] = None
    include_start: bool = True
    include_end: bool = True

    @classmethod
    def from_dict(cls, params: dict) -> "PreprocessingConfig":
        """
        Create PreprocessingConfig from a dictionary.

        This is the single source of truth for mapping dict keys to config fields.
        Handles the include_endpoints -> include_start/include_end mapping.

        Parameters
        ----------
        params : dict
            Dictionary with preprocessing parameters (e.g., from YAML config).
            If empty/None, returns default config.
            If provided, must contain all required keys.

        Returns
        -------
        PreprocessingConfig
            Configuration object with values from dict

        Raises
        ------
        KeyError
            If params is provided but missing required keys
        """
        if not params:
            # Explicit request for defaults
            return cls()

        # All keys required when config is provided
        include_endpoints = params["include_endpoints"]
        return cls(
            enabled=params["enabled"],
            baseline_order=params["baseline_order"],
            baseline_smoothing=params["baseline_smoothing"],
            include_start=params.get("include_start", include_endpoints),
            include_end=params.get("include_end", include_endpoints),
        )

    def to_dict(self) -> dict:
        """
        Convert config to dictionary for serialization.

        Returns
        -------
        dict
            Dictionary representation of this config
        """
        return {
            "enabled": self.enabled,
            "baseline_order": self.baseline_order,
            "baseline_smoothing": self.baseline_smoothing,
            "include_endpoints": self.include_start and self.include_end,
        }


class SignalPreprocessor:
    """
    Service for preprocessing chromatogram signals.

    Applies baseline correction to improve signal quality before peak detection.

    The preprocessing pipeline:
    1. Baseline correction using local minima spline fitting
    2. Store corrected signal as a variant on the chromatogram

    Examples
    --------
    >>> preprocessor = SignalPreprocessor()
    >>> preprocessor.preprocess(chromatogram)
    >>> corrected = chromatogram.get_signal("corrected")
    """

    def __init__(self, config: Optional[PreprocessingConfig] = None):
        """
        Initialize preprocessor with configuration.

        Parameters
        ----------
        config : PreprocessingConfig, optional
            Preprocessing parameters. Uses defaults if not provided.
        """
        self.config = config or PreprocessingConfig()

    def preprocess(
        self,
        chromatogram: Chromatogram,
        store_variant: str = "corrected",
    ) -> NDArray[np.float64]:
        """
        Apply preprocessing to a chromatogram's signal.

        Applies baseline correction, stores result as a signal variant,
        and returns the corrected signal.

        Parameters
        ----------
        chromatogram : Chromatogram
            Chromatogram to preprocess (modified in-place with new variant)
        store_variant : str
            Name for the corrected signal variant (default: "corrected")

        Returns
        -------
        NDArray[np.float64]
            Baseline-corrected signal
        """
        # Build baseline correction parameters
        params = MinimaSplineParams(
            order=self.config.baseline_order,
            smoothing=self.config.baseline_smoothing,
            include_start=self.config.include_start,
            include_end=self.config.include_end,
        )

        # Apply baseline correction
        result = minima_spline_baseline(
            signal=chromatogram.counts,
            x=chromatogram.time_points,
            params=params,
        )

        # Store corrected signal as variant
        chromatogram.add_signal_variant(store_variant, result.corrected)

        # Also store baseline for reference
        chromatogram.add_signal_variant("baseline", result.baseline)

        return result.corrected

    def preprocess_batch(
        self,
        chromatograms: list[Chromatogram],
        store_variant: str = "corrected",
    ) -> None:
        """
        Apply preprocessing to multiple chromatograms.

        Parameters
        ----------
        chromatograms : list[Chromatogram]
            Chromatograms to preprocess (modified in-place)
        store_variant : str
            Name for the corrected signal variant
        """
        for chrom in chromatograms:
            self.preprocess(chrom, store_variant)
