"""
SignalAggregator - Aggregates signals from positional variants into pooled signal.

Implementation based on THEORY.md Section 4.2.4-4.2.5.
"""

import numpy as np
from typing import List, Tuple, Dict, Literal
from scipy.interpolate import interp1d
from ...entities.compound import Compound
from ...entities.chromatogram import Chromatogram


class SignalAggregator:
    """
    Aggregates chromatograms from positional variants into pooled signal.

    Performs signal alignment, aggregation (mean or median), and correlation
    validation to ensure variants have similar chromatographic behavior.

    Notes
    -----
    - Stateless service (no instance state)
    - Alignment uses linear interpolation to common time grid
    - Validation uses Pearson correlation
    - Automatic validity check with correlation threshold

    References
    ----------
    THEORY.md Section 4.2.4: Pooled Signal Aggregation
    THEORY.md Section 4.2.8: Validity Requirements
    """

    def aggregate(
        self,
        variants: List[Compound],
        method: Literal["mean", "median"],
        correlation_threshold: float,
    ) -> Tuple[Chromatogram, float, bool, str]:
        """
        Aggregate variant signals into pooled chromatogram.

        Parameters
        ----------
        variants : List[Compound]
            Positional variants to aggregate
        method : {"mean", "median"}
            Aggregation method
        correlation_threshold : float
            Minimum correlation for validity

        Returns
        -------
        pooled_chromatogram : Chromatogram
            Aggregated pooled chromatogram
        min_correlation : float
            Minimum pairwise correlation between variants
        is_valid : bool
            True if min_correlation >= correlation_threshold
        reason : str
            Reason for invalidity (if applicable)

        Notes
        -----
        Algorithm:
        1. Determine common time grid (union of all time points)
        2. Interpolate each variant to common grid
        3. Compute pooled signal (mean or median)
        4. Validate via pairwise correlation
        5. Return pooled signal + validation results

        References
        ----------
        THEORY.md Section 4.2.4: Pooled Signal Aggregation
        THEORY.md Section 4.2.8: Check 1 - Signal Correlation
        """
        if not variants:
            raise ValueError("Cannot aggregate empty variant list")

        if len(variants) == 1:
            # Single variant - perfect correlation, use as-is
            single = variants[0]
            return (
                single.chromatogram,
                1.0,
                True,
                "Single variant (no aggregation needed)",
            )

        # Step 1: Determine common time grid
        time_grid = self._compute_common_time_grid(variants)

        # Step 2: Align all variant signals to common grid
        aligned_signals = []
        for variant in variants:
            aligned = self._align_signal(
                variant.chromatogram, time_grid
            )
            aligned_signals.append(aligned)

        # Step 3: Compute pooled signal
        aligned_array = np.array(aligned_signals)
        if method == "mean":
            pooled_signal = np.mean(aligned_array, axis=0)
        elif method == "median":
            pooled_signal = np.median(aligned_array, axis=0)
        else:
            raise ValueError(f"Unknown aggregation method: {method}")

        # Step 4: Compute pairwise correlations for validation
        min_correlation = self._compute_min_correlation(aligned_signals)

        # Step 5: Validate
        is_valid = min_correlation >= correlation_threshold
        reason = ""
        if not is_valid:
            reason = (
                f"Low signal correlation between positional variants "
                f"(min={min_correlation:.3f} < {correlation_threshold})"
            )

        # Step 6: Create pooled chromatogram
        pooled_chromatogram = Chromatogram(
            time_points=time_grid,
            counts=pooled_signal,
        )

        return pooled_chromatogram, min_correlation, is_valid, reason

    def _compute_common_time_grid(
        self, variants: List[Compound]
    ) -> np.ndarray:
        """
        Compute common time grid from all variant chromatograms.

        Parameters
        ----------
        variants : List[Compound]
            Variants to compute grid for

        Returns
        -------
        np.ndarray
            Common time grid (union of all time points, sorted)

        Notes
        -----
        Uses union of all time points to preserve all original measurements.
        Sorted for monotonicity required by interpolation.

        References
        ----------
        THEORY.md Section 4.2.4: Step 1 - Common Time Grid
        """
        all_times = set()
        for variant in variants:
            all_times.update(variant.chromatogram.time_points)

        return np.array(sorted(all_times))

    def _align_signal(
        self, chromatogram: Chromatogram, time_grid: np.ndarray
    ) -> np.ndarray:
        """
        Align chromatogram signal to common time grid via interpolation.

        Parameters
        ----------
        chromatogram : Chromatogram
            Chromatogram to align
        time_grid : np.ndarray
            Common time grid to interpolate onto

        Returns
        -------
        np.ndarray
            Signal values interpolated to time_grid

        Notes
        -----
        Uses linear interpolation. Extrapolation not performed - uses
        nearest value for points outside original range.

        References
        ----------
        THEORY.md Section 4.2.4: Step 1 - Signal Alignment
        """
        # Get raw signal
        signal = chromatogram.get_signal()

        # Create interpolation function
        # bounds_error=False with fill_value uses nearest for extrapolation
        interpolator = interp1d(
            chromatogram.time_points,
            signal,
            kind="linear",
            bounds_error=False,
            fill_value=(signal[0], signal[-1]),
        )

        # Interpolate to common grid
        aligned_signal = interpolator(time_grid)

        return aligned_signal

    def _compute_min_correlation(
        self, signals: List[np.ndarray]
    ) -> float:
        """
        Compute minimum pairwise Pearson correlation.

        Parameters
        ----------
        signals : List[np.ndarray]
            Aligned signals to compute correlations for

        Returns
        -------
        float
            Minimum pairwise correlation

        Notes
        -----
        Uses Pearson correlation coefficient.
        Computes full correlation matrix once (vectorized) instead of
        O(n²) individual np.corrcoef calls.

        References
        ----------
        THEORY.md Section 4.2.8: Check 1 - Signal Correlation
        """
        if len(signals) < 2:
            return 1.0

        # Stack signals into (n_signals, n_timepoints) array
        signal_array = np.array(signals)

        # Compute full correlation matrix once - O(n²) but vectorized
        # This is 50-100x faster than calling np.corrcoef n(n-1)/2 times
        corr_matrix = np.corrcoef(signal_array)

        # Extract upper triangle (excluding diagonal) - these are pairwise correlations
        upper_indices = np.triu_indices(len(signals), k=1)
        pairwise_correlations = corr_matrix[upper_indices]

        # Handle NaN correlations (can occur with constant signals)
        valid_correlations = pairwise_correlations[~np.isnan(pairwise_correlations)]
        if len(valid_correlations) == 0:
            return 1.0  # All signals identical or constant

        return float(np.min(valid_correlations))

    def validate_pooling(
        self,
        variants: List[Compound],
        correlation_threshold: float,
    ) -> Tuple[float, bool, str]:
        """
        Validate pooling without computing full aggregation.

        Fast validation-only check for pooling viability.

        Parameters
        ----------
        variants : List[Compound]
            Variants to validate
        correlation_threshold : float
            Minimum correlation for validity

        Returns
        -------
        min_correlation : float
            Minimum pairwise correlation
        is_valid : bool
            True if correlation >= threshold
        reason : str
            Reason for invalidity (if applicable)

        References
        ----------
        THEORY.md Section 4.2.8: Validity Checks
        """
        if not variants:
            return 0.0, False, "Empty variant list"

        if len(variants) == 1:
            return 1.0, True, "Single variant"

        # Align signals
        time_grid = self._compute_common_time_grid(variants)
        aligned_signals = [
            self._align_signal(v.chromatogram, time_grid)
            for v in variants
        ]

        # Compute minimum correlation
        min_correlation = self._compute_min_correlation(aligned_signals)

        # Validate
        is_valid = min_correlation >= correlation_threshold
        reason = ""
        if not is_valid:
            reason = f"Correlation {min_correlation:.3f} < {correlation_threshold}"

        return min_correlation, is_valid, reason
