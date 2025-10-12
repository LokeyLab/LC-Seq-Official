"""
Peak detection service using Discrete Morse theory and Poisson statistics.

Implementation based on THEORY.md Section 5.1-5.2.

This module implements peak detection for LC-Seq discrete fraction counts
using Discrete Morse Theory for local maxima detection and Poisson statistics
for significance testing, as described in THEORY.md Section 5.2.
"""

import numpy as np
from numpy.typing import NDArray
from typing import List, Tuple
from ..entities.chromatogram import Chromatogram
from ..entities.peak import Peak
from ...config import (
    DEFAULT_Z_THRESHOLD,
    DEFAULT_PROMINENCE_PERCENTILE,
    DEFAULT_MIN_SNR,
    DEFAULT_MIN_BASELINE_SDS,
    DEFAULT_SIGNAL_VARIANT,
)


class PeakDetector:
    """
    Detects significant peaks using Discrete Morse theory and Poisson statistics.

    Peak detection is based on finding local maxima (Discrete Morse theory) and
    filtering by statistical significance (Poisson Z-score) and chromatographic
    prominence.

    This is a stateless service - all methods are operations on input data with
    no instance state.

    Notes
    -----
    The detection pipeline:
    1. Find local maxima using Discrete Morse theory (no smoothing)
    2. Find shoulder peaks (inflection points for co-eluting compounds)
    3. Compute peak boundaries using valley detection
    4. Compute prominence (height above surrounding valleys)
    5. Compute Poisson Z-score for statistical significance
    6. Filter by Z-score threshold (Z > 3)
    7. Filter by prominence percentile threshold

    Statistical significance uses Poisson count statistics:
    Z = (c[i] - μ_bg) / √(μ_bg + ε)
    where μ_bg is the background level (10th percentile of all counts).

    Prominence measures chromatographic significance:
    prominence = height - max(valley_left, valley_right)

    References
    ----------
    THEORY.md Section 5.1: Discrete Morse Theory Framework
    THEORY.md Section 5.2: Statistical Significance Testing (Poisson + Prominence)
    THEORY.md Section 5.2.4: Peak Boundary Determination

    Examples
    --------
    >>> detector = PeakDetector()
    >>> chromatogram = Chromatogram(time_points=[...], counts=[...])
    >>> peaks = detector.detect_peaks(
    ...     chromatogram,
    ...     z_threshold=3.0,
    ...     prominence_percentile=0.2
    ... )
    """

    def detect_peaks(
        self,
        chromatogram: Chromatogram,
        z_threshold: float = DEFAULT_Z_THRESHOLD,
        prominence_percentile: float = DEFAULT_PROMINENCE_PERCENTILE,
        min_snr: float = DEFAULT_MIN_SNR,
        min_baseline_sds: float = DEFAULT_MIN_BASELINE_SDS,
        signal_variant: str = DEFAULT_SIGNAL_VARIANT,
    ) -> List[Peak]:
        """
        Detect significant peaks in chromatogram.

        Parameters
        ----------
        chromatogram : Chromatogram
            Chromatogram with signal data
        z_threshold : float, optional
            Poisson Z-score threshold for detection (default: 3.0, p < 0.001)
        prominence_percentile : float, optional
            Prominence percentile threshold - retain top (1-percentile) peaks
            e.g., 0.2 retains top 80% most prominent peaks (default: 0.2)
        min_snr : float, optional
            Standard deviations above minimum SNR for adaptive threshold (default: 0.0, disabled)
            Threshold = min(SNRs) + min_snr * std(SNRs)
            Data-adaptive filtering that removes peaks in bottom tail of SNR distribution
            Set to 0.0 to disable local SNR filtering
        min_baseline_sds : float, optional
            Global baseline threshold in standard deviations (default: 1.0)
            Peak height must exceed: min(signal) + min_baseline_sds * std(signal)
            Ensures peaks are above noise floor
        signal_variant : str, optional
            Signal variant to use for detection. Default is "raw" (per THEORY.md).

        Returns
        -------
        List[Peak]
            Detected peaks with positions, boundaries, prominence, and Z-scores

        Notes
        -----
        Detection pipeline:
        1. Find local maxima using Discrete Morse Theory
        2. Find shoulder peaks (for co-eluting compounds)
        3. Compute properties: prominence, Z-score, local SNR
        4. Filter by global baseline threshold (peak height > min + N*std)
        5. Filter by Z-score ≥ z_threshold (global significance)
        6. Filter by adaptive local SNR (data-driven quality threshold)
        7. Filter by prominence percentile (chromatographic significance)
        8. Verify local maximum property (height ≥ both neighbors)

        Local SNR computed as:
            SNR = (peak_height - median(surrounding)) / std(surrounding)
            where surrounding = regions AROUND peak (not including peak itself)

        Adaptive SNR threshold computed as:
            threshold = min(SNRs) + min_snr * std(SNRs)

        This data-driven approach removes peaks in the bottom tail of the
        SNR distribution, automatically adapting to each compound's signal quality.
        Noise is estimated from windows adjacent to peaks, following literature
        recommendations (PMC4324452, MDPI Processes 2022)

        References
        ----------
        THEORY.md Section 5.2: Statistical Significance Testing
        THEORY.md Section 5.0.1: LC-Seq Signal Characteristics
        PMC10612323: Picky with peakpicking (SNR via residuals)
        PMC4324452: Wavelet-Based Peak Detection (SNR threshold = 1.0)

        Examples
        --------
        >>> detector = PeakDetector()
        >>> peaks = detector.detect_peaks(
        ...     chromatogram,
        ...     z_threshold=3.0,
        ...     prominence_percentile=0.2,
        ...     min_snr=1.0
        ... )
        >>> len(peaks)
        5
        >>> peaks[0].prominence  # Height above local baseline
        150.3
        """
        if not chromatogram.has_signal_variant(signal_variant):
            raise ValueError(
                f"Signal variant '{signal_variant}' not found in chromatogram. "
                f"Available: {list(chromatogram.signal_variants.keys())}"
            )

        # Get signal (using raw signal per THEORY.md - no baseline correction, no smoothing)
        signal = chromatogram.get_signal(signal_variant)
        time_points = chromatogram.time_points

        # Estimate background level (10th percentile)
        background = self._estimate_background(signal)

        # Find all local maxima using Discrete Morse Theory
        peak_indices = self._find_local_maxima(signal)

        # Also find shoulder peaks (statistically significant but not strict maxima)
        # This handles co-eluting compounds where peaks overlap
        shoulder_indices = self._find_shoulder_peaks(signal, background, z_threshold, peak_indices)

        # Combine local maxima and shoulder peaks
        all_peak_indices = sorted(set(peak_indices + shoulder_indices))

        if not all_peak_indices:
            return []  # No peaks found

        # Compute properties for each peak
        peak_candidates = []
        for idx in all_peak_indices:
            # Get peak properties
            peak_time = time_points[idx]
            peak_height = signal[idx]

            # Skip peaks with non-positive height
            if peak_height <= 0:
                continue

            # Find boundaries (valley detection or 5% threshold)
            left_idx, right_idx = self._find_peak_boundaries(signal, idx)
            left_base = time_points[left_idx]
            right_base = time_points[right_idx]

            # Find valleys
            left_valley_idx = self._find_valley_left(signal, idx, left_idx)
            right_valley_idx = self._find_valley_right(signal, idx, right_idx)

            # Compute prominence (height above surrounding valleys)
            prominence = self._compute_prominence(
                signal, idx, left_valley_idx, right_valley_idx
            )

            # Compute Poisson Z-score
            z_score = self._compute_poisson_z_score(peak_height, background)

            # Compute local SNR
            local_snr = self._compute_local_snr(
                signal, idx, left_idx, right_idx, background
            )

            # Integrate area
            area = self._integrate_peak_area(signal, left_idx, right_idx)
            area = max(0.0, area)

            # Store candidate with all properties
            peak_candidates.append({
                'idx': idx,
                'time': peak_time,
                'height': peak_height,
                'prominence': prominence,
                'z_score': z_score,
                'local_snr': local_snr,
                'left_idx': left_idx,
                'right_idx': right_idx,
                'left_base': left_base,
                'right_base': right_base,
                'left_valley_idx': left_valley_idx,
                'right_valley_idx': right_valley_idx,
                'area': area
            })

        if not peak_candidates:
            return []

        # DIAGNOSTIC: Track filtering stages
        n_initial = len(peak_candidates)

        # Filter by global baseline threshold (peak height > min + N*std)
        signal_min = np.min(signal)
        signal_std = np.std(signal)
        baseline_threshold = signal_min + (min_baseline_sds * signal_std)

        baseline_filtered = [
            p for p in peak_candidates if p['height'] > baseline_threshold
        ]
        n_baseline = len(baseline_filtered)

        if not baseline_filtered:
            return []

        # Filter by Z-score threshold (statistical significance)
        significant_peaks = [
            p for p in baseline_filtered if p['z_score'] >= z_threshold
        ]
        n_zscore = len(significant_peaks)

        if not significant_peaks:
            return []

        # Filter by adaptive local SNR threshold (data-driven quality)
        # Skip if min_snr=0.0 (disabled)
        if min_snr > 0.0 and significant_peaks:
            # Compute SNR statistics from significant peaks
            snr_values = [p['local_snr'] for p in significant_peaks]
            snr_min = np.min(snr_values)
            snr_std = np.std(snr_values)

            # Adaptive threshold: min + N*std (N controlled by min_snr parameter)
            # This removes peaks in the bottom tail of the SNR distribution
            snr_threshold = snr_min + (min_snr * snr_std)

            snr_filtered = [p for p in significant_peaks if p['local_snr'] >= snr_threshold]
        else:
            # SNR filtering disabled
            snr_filtered = significant_peaks
        n_snr = len(snr_filtered)

        if not snr_filtered:
            return []

        # Filter by prominence percentile threshold
        prominences = [p['prominence'] for p in snr_filtered]
        prominence_cutoff = np.percentile(prominences, prominence_percentile * 100)
        prominence_filtered = [
            p for p in snr_filtered if p['prominence'] >= prominence_cutoff
        ]
        n_prominence = len(prominence_filtered)

        # Verify each peak is an actual local maximum (higher than or equal to both neighbors)
        verified_peaks = []
        for p in prominence_filtered:
            idx = p['idx']

            # Check left neighbor
            if idx > 0:
                left_neighbor = signal[idx - 1]
                if signal[idx] < left_neighbor:
                    continue  # Not a peak - left neighbor is higher

            # Check right neighbor
            if idx < len(signal) - 1:
                right_neighbor = signal[idx + 1]
                if signal[idx] < right_neighbor:
                    continue  # Not a peak - right neighbor is higher

            # Peak verified - higher than or equal to both neighbors
            verified_peaks.append(p)
        n_verified = len(verified_peaks)

        # DIAGNOSTIC: Print filtering stages for debugging
        # Uncomment to see filtering progression:
        # print(f"  Filtering: {n_initial} → baseline:{n_baseline} → zscore:{n_zscore} → snr:{n_snr} → prom:{n_prominence} → verified:{n_verified}")

        # Build Peak objects
        peaks = []
        for p in verified_peaks:
            left_valley = time_points[p['left_valley_idx']] if p['left_valley_idx'] is not None else None
            right_valley = time_points[p['right_valley_idx']] if p['right_valley_idx'] is not None else None

            peak = Peak(
                position=p['time'],
                left_base=p['left_base'],
                right_base=p['right_base'],
                height=p['height'],
                area=p['area'],
                left_valley=left_valley,
                right_valley=right_valley,
                prominence=p['prominence']
            )
            peaks.append(peak)

        return peaks

    def _estimate_background(self, signal: NDArray[np.float64]) -> float:
        """
        Estimate background level using 10th percentile.

        Parameters
        ----------
        signal : NDArray[np.float64]
            Signal array

        Returns
        -------
        float
            Background level (μ_bg)

        Notes
        -----
        Uses 10th percentile to capture low-count baseline.
        Per THEORY.md Section 5.2: μ_bg = percentile(all_counts, 10)
        """
        return float(np.percentile(signal, 10))

    def _compute_poisson_z_score(
        self,
        peak_height: float,
        background: float,
        epsilon: float = 1.0
    ) -> float:
        """
        Compute Poisson Z-score for statistical significance.

        Parameters
        ----------
        peak_height : float
            Peak height (counts)
        background : float
            Background level (μ_bg)
        epsilon : float, optional
            Regularization for low counts (default: 1.0)

        Returns
        -------
        float
            Z-score = (c[i] - μ_bg) / √(μ_bg + ε)

        Notes
        -----
        For Poisson counts, variance ≈ mean, so σ ≈ √mean.
        Z-score measures how many standard deviations above background.
        Z > 3 corresponds to p < 0.001 (highly significant).

        References
        ----------
        THEORY.md Section 5.2: Statistical Hypothesis Testing
        """
        return (peak_height - background) / np.sqrt(background + epsilon)

    def _compute_prominence(
        self,
        signal: NDArray[np.float64],
        peak_idx: int,
        left_valley_idx: int | None,
        right_valley_idx: int | None
    ) -> float:
        """
        Compute prominence (height above surrounding valleys).

        Parameters
        ----------
        signal : NDArray[np.float64]
            Signal array
        peak_idx : int
            Peak maximum index
        left_valley_idx : int or None
            Left valley index
        right_valley_idx : int or None
            Right valley index

        Returns
        -------
        float
            Prominence = height - max(valley_left, valley_right)

        Notes
        -----
        Prominence measures how much a peak rises above its local baseline.
        If no valleys found, uses signal start/end as baseline.

        References
        ----------
        THEORY.md Section 5.2: Prominence (Local Significance)
        """
        peak_height = signal[peak_idx]

        # Get valley heights (use signal edges if valleys not found)
        if left_valley_idx is not None:
            left_valley_height = signal[left_valley_idx]
        else:
            # No left valley - use signal start as baseline
            left_valley_height = signal[0]

        if right_valley_idx is not None:
            right_valley_height = signal[right_valley_idx]
        else:
            # No right valley - use signal end as baseline
            right_valley_height = signal[-1]

        # Prominence = height above highest surrounding valley
        local_baseline = max(left_valley_height, right_valley_height)
        prominence = peak_height - local_baseline

        # Ensure non-negative (can be negative for edge peaks if signal edges are high)
        # This shouldn't happen with proper local maxima, but handle edge cases gracefully
        return float(max(0.0, prominence))

    def _compute_local_snr(
        self,
        signal: NDArray[np.float64],
        peak_idx: int,
        left_idx: int,
        right_idx: int,
        background: float,
        window_size: int = 5
    ) -> float:
        """
        Compute local signal-to-noise ratio using surrounding regions.

        Parameters
        ----------
        signal : NDArray[np.float64]
            Signal array
        peak_idx : int
            Peak maximum index
        left_idx : int
            Left boundary index
        right_idx : int
            Right boundary index
        background : float
            Global background level (fallback for noise estimate)
        window_size : int, optional
            Size of surrounding windows for noise estimation (default: 5)

        Returns
        -------
        float
            Local SNR = signal_amplitude / local_noise

        Notes
        -----
        Signal amplitude = peak height - median(surrounding regions)
        Local noise = std(surrounding regions)

        Noise is estimated from windows AROUND the peak (not including peak):
        - Left window: [left_idx-window_size, left_idx)
        - Right window: (right_idx, right_idx+window_size]

        This follows literature recommendations:
        - PMC4324452: Noise from local bins excluding peaks
        - MDPI Processes 2022: Separate windows for signal/noise
        - Avoids inflating SNR by including peak in noise calculation

        If surrounding regions have < 3 points, falls back to global background.

        References
        ----------
        PMC4324452: Wavelet-Based Peak Detection
        MDPI Processes 10(6):1098: Multi-Sliding Window Method
        """
        peak_height = signal[peak_idx]

        # Define noise estimation windows (regions AROUND peak, not including it)
        # Left noise window: [left_idx-window_size, left_idx)
        left_start = max(0, left_idx - window_size)
        left_noise = signal[left_start:left_idx] if left_idx > 0 else np.array([])

        # Right noise window: (right_idx, right_idx+window_size]
        right_end = min(len(signal), right_idx + 1 + window_size)
        right_noise = signal[right_idx+1:right_end] if right_idx < len(signal)-1 else np.array([])

        # Combine surrounding regions
        surrounding = np.concatenate([left_noise, right_noise])

        if len(surrounding) < 3:
            # Not enough surrounding data, fall back to global background
            local_baseline = background
            local_noise = np.sqrt(background + 1.0)
        else:
            # Use median of surrounding as baseline (robust to outliers)
            local_baseline = np.median(surrounding)
            # Use standard deviation as noise
            local_noise = np.std(surrounding)

            # Ensure minimum noise level
            if local_noise < 1e-6:
                local_noise = np.sqrt(background + 1.0)

        # Signal amplitude above local baseline
        signal_amplitude = peak_height - local_baseline

        return float(signal_amplitude / local_noise)

    def _find_shoulder_peaks(
        self,
        signal: NDArray[np.float64],
        background: float,
        z_threshold: float,
        existing_maxima: List[int],
        min_separation: int = 2
    ) -> List[int]:
        """
        Find shoulder peaks - statistically significant points that aren't strict local maxima.

        Shoulder peaks occur when compounds co-elute, creating overlapping signals.
        These appear as inflection points or plateaus within broader peaks.

        Parameters
        ----------
        signal : NDArray[np.float64]
            Signal array
        background : float
            Background level for Z-score calculation
        z_threshold : float
            Z-score threshold for significance
        existing_maxima : List[int]
            Indices of already-detected local maxima
        min_separation : int, optional
            Minimum index separation from existing maxima (default: 2)

        Returns
        -------
        List[int]
            Indices of shoulder peaks

        Notes
        -----
        Detects shoulder peaks by finding:
        1. Points with Z-score >= threshold (statistically significant)
        2. NOT within min_separation of existing local maxima
        3. Local maxima in the signal derivative (inflection points)

        This handles overlapping peaks from co-eluting compounds, which is
        common in LC-MS where multiple peptides may have similar retention times.
        """
        if len(signal) < 5:
            return []

        shoulders = []

        # Compute first derivative (discrete approximation)
        derivative = np.zeros(len(signal))
        for i in range(1, len(signal) - 1):
            derivative[i] = (signal[i+1] - signal[i-1]) / 2.0

        # Find points with high Z-score that are derivative maxima
        for i in range(2, len(signal) - 2):
            # Check Z-score
            z_score = (signal[i] - background) / np.sqrt(background + 1.0)
            if z_score < z_threshold:
                continue

            # Skip if too close to existing local maximum
            too_close = any(abs(i - m) <= min_separation for m in existing_maxima)
            if too_close:
                continue

            # Check if this is a local maximum in the derivative
            # (indicates inflection point in original signal)
            is_derivative_max = (
                derivative[i] > derivative[i-1] and
                derivative[i] >= derivative[i+1]
            )

            # Or check if signal is in a plateau region (derivative near zero)
            is_plateau = abs(derivative[i]) < 0.5 and signal[i] > signal[i-2] and signal[i] > signal[i+2]

            if is_derivative_max or is_plateau:
                shoulders.append(i)

        return shoulders

    def _find_local_maxima(self, signal: NDArray[np.float64]) -> List[int]:
        """
        Find local maxima using Discrete Morse theory.

        A local maximum in a discrete sequence is defined as:
        c[i] > c[i-1]  AND  c[i] >= c[i+1]

        Parameters
        ----------
        signal : NDArray[np.float64]
            Signal array

        Returns
        -------
        List[int]
            Indices of local maxima

        Notes
        -----
        This is the discrete analogue of continuous Morse theory.
        Uses >= for right comparison to handle flat peaks consistently.

        References
        ----------
        THEORY.md Section 5.1: Discrete Morse Theory Framework
        """
        if len(signal) < 3:
            return []

        maxima = []

        # Check each interior point
        for i in range(1, len(signal) - 1):
            # Simple local maximum check: higher than both neighbors
            # Use >= for right comparison to handle flat peaks
            if signal[i] > signal[i-1] and signal[i] >= signal[i+1]:
                # Ensure we're at a true maximum, not on a plateau
                # If equal to right neighbor, only add if next value is lower
                if signal[i] == signal[i+1]:
                    # Check if this is start of flat region that then descends
                    if i+2 < len(signal) and signal[i+1] > signal[i+2]:
                        maxima.append(i)
                else:
                    maxima.append(i)

        return maxima

    def _find_peak_boundaries(
        self,
        signal: NDArray[np.float64],
        peak_idx: int,
        threshold_fraction: float = 0.05
    ) -> Tuple[int, int]:
        """
        Find peak boundaries using valley detection or 5% threshold.

        Scans left and right from peak until:
        1. Valley (local minimum) is found, OR
        2. Signal drops below 5% of peak height, OR
        3. Signal edge is reached

        Parameters
        ----------
        signal : NDArray[np.float64]
            Signal array
        peak_idx : int
            Index of peak maximum
        threshold_fraction : float, optional
            Fraction of peak height for threshold. Default is 0.05 (5%).

        Returns
        -------
        Tuple[int, int]
            (left_boundary_idx, right_boundary_idx)

        Notes
        -----
        The 5% threshold is conservative - captures full peak including tails.
        Valley detection provides natural separation between adjacent peaks.

        References
        ----------
        THEORY.md Section 5.2.4: Peak Boundary Determination
        """
        peak_height = signal[peak_idx]
        threshold = threshold_fraction * peak_height

        # Find left boundary
        left_idx = peak_idx
        for i in range(peak_idx - 1, -1, -1):
            if signal[i] < threshold:
                left_idx = i
                break
            # Check for valley (local minimum)
            if i > 0 and signal[i] < signal[i-1] and signal[i] < signal[i+1]:
                left_idx = i
                break
            left_idx = i  # Update to current position

        # Find right boundary
        right_idx = peak_idx
        for i in range(peak_idx + 1, len(signal)):
            if signal[i] < threshold:
                right_idx = i
                break
            # Check for valley (local minimum)
            if i < len(signal) - 1 and signal[i] < signal[i-1] and signal[i] < signal[i+1]:
                right_idx = i
                break
            right_idx = i  # Update to current position

        return left_idx, right_idx

    def _find_valley_left(
        self,
        signal: NDArray[np.float64],
        peak_idx: int,
        left_boundary: int
    ) -> int | None:
        """
        Find valley (local minimum) to the left of peak.

        Parameters
        ----------
        signal : NDArray[np.float64]
            Signal array
        peak_idx : int
            Peak maximum index
        left_boundary : int
            Left boundary index

        Returns
        -------
        int or None
            Valley index if found, None otherwise
        """
        if left_boundary >= peak_idx - 1:
            return None

        for i in range(peak_idx - 1, left_boundary, -1):
            if i > 0 and signal[i] < signal[i-1] and signal[i] < signal[i+1]:
                return i

        return None

    def _find_valley_right(
        self,
        signal: NDArray[np.float64],
        peak_idx: int,
        right_boundary: int
    ) -> int | None:
        """
        Find valley (local minimum) to the right of peak.

        Parameters
        ----------
        signal : NDArray[np.float64]
            Signal array
        peak_idx : int
            Peak maximum index
        right_boundary : int
            Right boundary index

        Returns
        -------
        int or None
            Valley index if found, None otherwise
        """
        if right_boundary <= peak_idx + 1:
            return None

        for i in range(peak_idx + 1, right_boundary):
            if i < len(signal) - 1 and signal[i] < signal[i-1] and signal[i] < signal[i+1]:
                return i

        return None

    def _integrate_peak_area(
        self,
        signal: NDArray[np.float64],
        left_idx: int,
        right_idx: int
    ) -> float:
        """
        Integrate peak area using simple summation.

        Parameters
        ----------
        signal : NDArray[np.float64]
            Signal array
        left_idx : int
            Left boundary index
        right_idx : int
            Right boundary index

        Returns
        -------
        float
            Integrated area

        Notes
        -----
        Uses simple summation of signal values in the peak region.
        Works on raw signals (per THEORY.md - no baseline correction).

        References
        ----------
        THEORY.md Section 5.0.7: Peak Area Integration
        """
        if right_idx <= left_idx:
            return 0.0

        return float(np.sum(signal[left_idx:right_idx+1]))
