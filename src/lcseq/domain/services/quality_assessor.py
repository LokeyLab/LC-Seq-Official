"""
QualityAssessor - Assesses signal quality for library filtering.

Computes quality metrics for chromatograms and equivalence classes
to enable pre-filtering of low-quality data before analysis.
"""

import numpy as np
from dataclasses import dataclass
from typing import List, Optional, Tuple
from scipy.interpolate import interp1d

from ..entities.compound import Compound
from ..entities.chromatogram import Chromatogram


@dataclass
class SignalQualityMetrics:
    """Quality metrics for a single chromatogram."""

    total_signal: float  # Sum of all counts
    max_signal: float  # Maximum signal value
    median_signal: float  # Median signal value
    noise_std: float  # Standard deviation of baseline region
    noise_ratio: float  # noise_std / median_signal (lower is better)
    baseline: float  # Estimated baseline level


@dataclass
class EquivalenceClassQuality:
    """Quality metrics for an equivalence class (group of positional variants)."""

    block_support_sequence: str
    n_variants: int
    min_correlation: float  # Minimum pairwise correlation (1.0 for single variant)
    mean_total_signal: float
    mean_noise_ratio: float

    # Filter results
    passes_correlation: bool
    passes_intensity: bool
    passes_noise: bool
    passes_all: bool

    # Individual variant metrics
    variant_metrics: List[SignalQualityMetrics]


class QualityAssessor:
    """
    Assesses signal quality for library filtering.

    Computes metrics at both individual compound and equivalence class levels
    to enable principled pre-filtering of low-quality data.

    Quality Metrics
    ---------------
    1. Signal Intensity: Total/max signal (low = failed synthesis)
    2. Noise Level: Baseline std / median signal (high = unreliable)
    3. Replicate Correlation: Min pairwise correlation in equivalence class

    Notes
    -----
    - Stateless service
    - Uses sigma-clipping for robust baseline/noise estimation
    - Correlation computed on aligned signals (handles different time grids)

    References
    ----------
    THEORY.md Section 4.2.8: Validity Requirements
    """

    def __init__(self, sigma_clip_sigma: float = 2.0):
        """
        Initialize quality assessor.

        Parameters
        ----------
        sigma_clip_sigma : float
            Sigma for sigma-clipping baseline estimation (default 2.0 = 95% CI)
        """
        self._sigma = sigma_clip_sigma

    def assess_signal(self, chromatogram: Chromatogram) -> SignalQualityMetrics:
        """
        Compute quality metrics for a single chromatogram.

        Parameters
        ----------
        chromatogram : Chromatogram
            Chromatogram to assess

        Returns
        -------
        SignalQualityMetrics
            Quality metrics for the signal
        """
        signal = chromatogram.get_signal()

        # Basic signal statistics
        total_signal = float(np.sum(signal))
        max_signal = float(np.max(signal))
        median_signal = float(np.median(signal))

        # Estimate baseline and noise using sigma-clipping
        baseline, noise_std = self._estimate_baseline_and_noise(signal)

        # Noise ratio (lower is better)
        noise_ratio = noise_std / median_signal if median_signal > 0 else float('inf')

        return SignalQualityMetrics(
            total_signal=total_signal,
            max_signal=max_signal,
            median_signal=median_signal,
            noise_std=noise_std,
            noise_ratio=noise_ratio,
            baseline=baseline,
        )

    def assess_equivalence_class(
        self,
        variants: List[Compound],
        block_support_sequence: str,
        correlation_threshold: float = 0.8,
        intensity_percentile_threshold: Optional[float] = None,
        intensity_absolute_threshold: Optional[float] = None,
        max_noise_ratio: Optional[float] = None,
        library_intensity_percentile: Optional[float] = None,
    ) -> EquivalenceClassQuality:
        """
        Assess quality of an equivalence class (group of positional variants).

        Parameters
        ----------
        variants : List[Compound]
            Positional variants in the equivalence class
        block_support_sequence : str
            Block support sequence identifying this class
        correlation_threshold : float
            Minimum acceptable correlation between variants
        intensity_percentile_threshold : float, optional
            Intensity threshold from library percentile
        intensity_absolute_threshold : float, optional
            Absolute minimum intensity threshold
        max_noise_ratio : float, optional
            Maximum acceptable noise ratio
        library_intensity_percentile : float, optional
            Pre-computed library intensity percentile for comparison

        Returns
        -------
        EquivalenceClassQuality
            Quality assessment for the equivalence class
        """
        # Compute individual variant metrics
        variant_metrics = []
        for variant in variants:
            if variant.chromatogram is not None:
                metrics = self.assess_signal(variant.chromatogram)
                variant_metrics.append(metrics)

        if not variant_metrics:
            # No valid chromatograms
            return EquivalenceClassQuality(
                block_support_sequence=block_support_sequence,
                n_variants=len(variants),
                min_correlation=0.0,
                mean_total_signal=0.0,
                mean_noise_ratio=float('inf'),
                passes_correlation=False,
                passes_intensity=False,
                passes_noise=False,
                passes_all=False,
                variant_metrics=[],
            )

        # Aggregate metrics
        mean_total_signal = np.mean([m.total_signal for m in variant_metrics])
        mean_noise_ratio = np.mean([m.noise_ratio for m in variant_metrics])

        # Compute replicate correlation
        if len(variants) > 1:
            min_correlation = self._compute_min_correlation(variants)
        else:
            min_correlation = 1.0  # Single variant, perfect "correlation"

        # Apply filters
        passes_correlation = min_correlation >= correlation_threshold

        # Intensity filter
        if intensity_percentile_threshold is not None:
            passes_intensity = mean_total_signal >= intensity_percentile_threshold
        elif intensity_absolute_threshold is not None:
            passes_intensity = mean_total_signal >= intensity_absolute_threshold
        else:
            passes_intensity = True  # No intensity filter applied

        # Noise filter
        if max_noise_ratio is not None:
            passes_noise = mean_noise_ratio <= max_noise_ratio
        else:
            passes_noise = True  # No noise filter applied

        passes_all = passes_correlation and passes_intensity and passes_noise

        return EquivalenceClassQuality(
            block_support_sequence=block_support_sequence,
            n_variants=len(variants),
            min_correlation=min_correlation,
            mean_total_signal=mean_total_signal,
            mean_noise_ratio=mean_noise_ratio,
            passes_correlation=passes_correlation,
            passes_intensity=passes_intensity,
            passes_noise=passes_noise,
            passes_all=passes_all,
            variant_metrics=variant_metrics,
        )

    def compute_library_intensity_percentile(
        self,
        compounds: List[Compound],
        percentile: float,
    ) -> float:
        """
        Compute intensity threshold from library percentile.

        Parameters
        ----------
        compounds : List[Compound]
            All compounds in library
        percentile : float
            Percentile to compute (e.g., 0.05 for bottom 5%)

        Returns
        -------
        float
            Intensity value at the given percentile
        """
        intensities = []
        for compound in compounds:
            if compound.chromatogram is not None:
                signal = compound.chromatogram.get_signal()
                intensities.append(float(np.sum(signal)))

        if not intensities:
            return 0.0

        return float(np.percentile(intensities, percentile * 100))

    def _estimate_baseline_and_noise(
        self, signal: np.ndarray, max_iter: int = 10
    ) -> Tuple[float, float]:
        """
        Estimate baseline and noise using sigma-clipping.

        Parameters
        ----------
        signal : np.ndarray
            Signal array
        max_iter : int
            Maximum iterations for sigma-clipping

        Returns
        -------
        Tuple[float, float]
            (baseline, noise_std)
        """
        mask = np.ones(len(signal), dtype=bool)

        for _ in range(max_iter):
            masked_signal = signal[mask]
            if len(masked_signal) == 0:
                break

            mean = np.mean(masked_signal)
            std = np.std(masked_signal)

            new_mask = signal <= (mean + self._sigma * std)

            if np.array_equal(mask, new_mask):
                break

            mask = new_mask

        remaining = signal[mask]
        if len(remaining) == 0:
            return float(np.median(signal)), float(np.std(signal))

        baseline = float(np.median(remaining))
        noise_std = float(np.std(remaining))

        return baseline, noise_std

    def _compute_min_correlation(self, variants: List[Compound]) -> float:
        """
        Compute minimum pairwise Pearson correlation between variant signals.

        Parameters
        ----------
        variants : List[Compound]
            Variants to compute correlation for

        Returns
        -------
        float
            Minimum pairwise correlation
        """
        if len(variants) < 2:
            return 1.0

        # Get chromatograms with valid signals
        valid_variants = [v for v in variants if v.chromatogram is not None]
        if len(valid_variants) < 2:
            return 1.0

        # Compute common time grid
        all_times = set()
        for variant in valid_variants:
            all_times.update(variant.chromatogram.time_points)
        time_grid = np.array(sorted(all_times))

        # Align signals to common grid
        aligned_signals = []
        for variant in valid_variants:
            aligned = self._align_signal(variant.chromatogram, time_grid)
            aligned_signals.append(aligned)

        # Compute correlation matrix
        signal_array = np.array(aligned_signals)
        corr_matrix = np.corrcoef(signal_array)

        # Extract upper triangle (pairwise correlations)
        upper_indices = np.triu_indices(len(aligned_signals), k=1)
        pairwise_correlations = corr_matrix[upper_indices]

        # Handle NaN (constant signals)
        valid_correlations = pairwise_correlations[~np.isnan(pairwise_correlations)]
        if len(valid_correlations) == 0:
            return 1.0

        return float(np.min(valid_correlations))

    def _align_signal(
        self, chromatogram: Chromatogram, time_grid: np.ndarray
    ) -> np.ndarray:
        """
        Align chromatogram to common time grid via interpolation.

        Parameters
        ----------
        chromatogram : Chromatogram
            Chromatogram to align
        time_grid : np.ndarray
            Common time grid

        Returns
        -------
        np.ndarray
            Aligned signal
        """
        signal = chromatogram.get_signal()

        interpolator = interp1d(
            chromatogram.time_points,
            signal,
            kind="linear",
            bounds_error=False,
            fill_value=(signal[0], signal[-1]),
        )

        return interpolator(time_grid)
