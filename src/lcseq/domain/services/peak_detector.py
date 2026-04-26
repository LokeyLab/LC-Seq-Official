"""
Peak detection service using Discrete Morse theory and Negative Binomial statistics.

Implementation based on THEORY.md Section 5.1-5.2.

This module implements peak detection for LC-Seq discrete fraction counts
using Discrete Morse Theory for local maxima detection and Negative Binomial
statistics for significance testing. The Negative Binomial model properly
handles overdispersion (Var > Mean), reducing to Poisson when data is
not overdispersed (r → ∞).
"""

import numpy as np
from numpy.typing import NDArray
from typing import List, Tuple, Optional
from ..entities.chromatogram import Chromatogram
from ..entities.peak import Peak, RejectionReason
from .baseline_estimator import BaselineEstimatorService
from .significance_tester import SignificanceTesterService


class PeakDetector:
    """
    Detects significant peaks using Discrete Morse theory and Negative Binomial statistics.

    Peak detection is based on finding local maxima (Discrete Morse theory) and
    filtering by statistical significance (Negative Binomial test) and
    chromatographic prominence.

    This is a stateless service - all methods are operations on input data with
    no instance state.

    Notes
    -----
    The detection pipeline:
    1. Find local maxima using Discrete Morse theory (no smoothing)
    2. Find shoulder peaks (inflection points for co-eluting compounds)
    3. Compute peak boundaries using valley detection
    4. Compute prominence (height above surrounding valleys)
    5. Compute significance using Negative Binomial test
    6. Filter by significance threshold (p-value < alpha)
    7. Filter by prominence percentile threshold

    Statistical significance uses Negative Binomial distribution:
    P(X >= observed | NB(μ=background, r=dispersion)) < alpha

    The Negative Binomial model handles overdispersion (Var > Mean) which is
    common in scaled/normalized count data. When dispersion r → ∞, the NB
    distribution reduces to Poisson (Var = Mean).

    Prominence measures chromatographic significance:
    prominence = height - max(valley_left, valley_right)

    References
    ----------
    THEORY.md Section 5.1: Discrete Morse Theory Framework
    THEORY.md Section 5.2: Statistical Significance Testing
    THEORY.md Section 5.2.4: Peak Boundary Determination

    Examples
    --------
    >>> detector = PeakDetector()
    >>> chromatogram = Chromatogram(time_points=[...], counts=[...])
    >>> peaks = detector.detect_peaks(
    ...     chromatogram,
    ...     alpha=0.001,
    ...     prominence_percentile=0.0
    ... )
    """

    def __init__(self, sigma_clip_sigma: float = 2.0):
        """
        Initialize peak detector.

        Parameters
        ----------
        sigma_clip_sigma : float
            Sigma value for sigma-clipping baseline estimation.
            Default 2.0 corresponds to 95% confidence interval.

        Notes
        -----
        Significance testing uses the Negative Binomial distribution, which:
        - Properly handles overdispersion (Var > Mean)
        - Reduces to Poisson when data is not overdispersed (r → ∞)
        - Uses exact CDF (no Normal approximation)
        - Is the standard in genomics (DESeq2, edgeR)
        """
        self._baseline_estimator = BaselineEstimatorService(sigma=sigma_clip_sigma)
        self._significance_tester = SignificanceTesterService()

    def detect_peaks(
        self,
        chromatogram: Chromatogram,
        alpha: float,
        prominence_percentile: float,
        min_snr: float,
        min_baseline_sds: float,
        signal_variant: str,
        min_dispersion_r: float,
        include_rejected: bool = False,
    ) -> List[Peak]:
        """
        Detect significant peaks in chromatogram.

        Parameters
        ----------
        chromatogram : Chromatogram
            Chromatogram with signal data
        alpha : float, optional
            Significance level (false positive rate) for peak detection.
            Default 0.001 = 0.1% chance of false positive per test.
            Uses Negative Binomial test with dispersion estimation.
        prominence_percentile : float, optional
            Prominence percentile threshold - retain top (1-percentile) peaks
            e.g., 0.2 retains top 80% most prominent peaks (default: 0.2)
        min_snr : float, optional
            Standard deviations above minimum SNR for adaptive threshold (default: 0.0, disabled)
            Threshold = min(SNRs) + min_snr * std(SNRs)
            Data-adaptive filtering that removes peaks in bottom tail of SNR distribution
            Set to 0.0 to disable local SNR filtering
        min_baseline_sds : float, optional
            Global baseline threshold in standard deviations (default: 0.0, disabled)
            Peak height must exceed: baseline + min_baseline_sds * noise_std
            where baseline and noise_std are estimated via sigma-clipping
            (excludes peaks, providing a robust threshold)
            Disabled by default - significance test handles this.
        signal_variant : str, optional
            Signal variant to use for detection. Default is "raw" (per THEORY.md).
        include_rejected : bool, optional
            If True, include rejected peaks with rejection_reason set.
            Rejected peaks are local maxima that failed one of the filters.
            Useful for diagnostic visualization. Default False.

        Returns
        -------
        List[Peak]
            Detected peaks with positions, boundaries, prominence, and p-values.
            If include_rejected=True, includes rejected peaks with rejection_reason set.

        Notes
        -----
        Detection pipeline:
        1. Find local maxima using Discrete Morse Theory
        2. Find shoulder peaks (for co-eluting compounds)
        3. Compute properties: prominence, p-value, local SNR
        4. Filter by global baseline threshold (disabled by default)
        5. Filter by significance test (p-value < alpha)
        6. Filter by adaptive local SNR (disabled by default)
        7. Filter by prominence percentile (chromatographic significance)
        8. Verify local maximum property (height ≥ both neighbors)

        Significance testing:
            Uses Negative Binomial test with dispersion estimation.
            When no overdispersion is detected (r → ∞), NB naturally reduces to Poisson.
            P(X >= observed | NB(μ = background, r = dispersion)) < alpha

        The alpha parameter specifies the false positive rate - the probability
        of incorrectly detecting noise as a peak. Default 0.001 = 0.1%.

        References
        ----------
        THEORY.md Section 5.2: Statistical Significance Testing
        THEORY.md Section 5.0.1: LC-Seq Signal Characteristics

        Examples
        --------
        >>> detector = PeakDetector()
        >>> peaks = detector.detect_peaks(
        ...     chromatogram,
        ...     alpha=0.001,  # 0.1% false positive rate
        ...     prominence_percentile=0.2,
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

        # Estimate background level, noise std, and NB dispersion (excludes peaks via sigma-clipping)
        # Dispersion parameter r is used for Negative Binomial significance testing
        background, noise_std, dispersion = self._baseline_estimator.estimate_with_noise_and_dispersion(
            signal, min_dispersion_r=min_dispersion_r
        )

        # Find all local maxima using Discrete Morse Theory
        peak_indices = self._find_local_maxima(signal)

        # Also find shoulder peaks (statistically significant but not strict maxima)
        # This handles co-eluting compounds where peaks overlap
        shoulder_indices = self._find_shoulder_peaks(signal, background, alpha, dispersion, peak_indices)

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

            # Find boundaries (valley detection only - no arbitrary thresholds)
            left_idx, right_idx = self._find_peak_boundaries(signal, idx)
            left_base = time_points[left_idx]
            right_base = time_points[right_idx]

            # Find valleys
            left_valley_idx = self._find_valley_left(signal, idx, left_idx)
            right_valley_idx = self._find_valley_right(signal, idx, right_idx)

            # Compute local baseline from valleys (handles chromatographic drift)
            local_baseline = self._compute_local_baseline(
                signal, left_valley_idx, right_valley_idx,
                left_idx, right_idx, background
            )

            # Compute prominence (height above surrounding valleys)
            prominence = self._compute_prominence(
                signal, idx, left_valley_idx, right_valley_idx
            )

            # Integrate area (needed for area-based significance test)
            area = self._integrate_peak_area(signal, left_idx, right_idx)
            area = max(0.0, area)

            # Compute expected area under null hypothesis
            # H0: region contains only baseline counts ~ NB(μ=local_baseline × width, r=dispersion)
            peak_width = right_idx - left_idx + 1
            expected_area = local_baseline * peak_width

            # DUAL SIGNIFICANCE TEST: Accept if EITHER area OR height is significant
            # This is principled because:
            # - Area significance catches broad peaks (handles peak broadening)
            # - Height significance catches sharp peaks (important for purity)
            # Both indicate "there's a real species here" - different peak shapes

            # Test 1: AREA-based significance
            # Tests: P(X >= observed_area | X ~ NB(expected_area, dispersion))
            is_area_significant, p_value_area = self._significance_tester.test(
                area, expected_area, alpha, dispersion
            )

            # Test 2: HEIGHT-based significance (using local baseline)
            # Tests: P(X >= peak_height | X ~ NB(local_baseline, dispersion))
            is_height_significant, p_value_height = self._significance_tester.test(
                peak_height, local_baseline, alpha, dispersion
            )

            # Accept if EITHER test passes
            is_significant = is_area_significant or is_height_significant
            # Report the more significant (lower) p-value
            p_value = min(p_value_area, p_value_height)

            # Compute local SNR (still uses height for SNR definition)
            local_snr = self._compute_local_snr(
                signal, idx, left_idx, right_idx, background
            )

            # Store candidate with all properties
            peak_candidates.append({
                'idx': idx,
                'time': peak_time,
                'height': peak_height,
                'prominence': prominence,
                'p_value': p_value,
                'is_significant': is_significant,
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

        # Resolve overlapping boundaries between adjacent peaks
        self._resolve_overlapping_boundaries(peak_candidates, signal, time_points)

        # Compute baseline threshold using baseline and noise std (excludes peaks)
        # Peak height must exceed: baseline + min_baseline_sds * noise_std
        baseline_threshold = background + (min_baseline_sds * noise_std)

        # Phase 1: Assign rejection reasons to each candidate
        # We process in order: significance → SNR → prominence → local max
        # First rejection reason wins (most fundamental filter first)
        # Note: Baseline threshold is disabled by default (min_baseline_sds=0)
        # because significance test already handles "above background" significance

        for p in peak_candidates:
            p['rejection_reason'] = RejectionReason.NONE  # Default: accepted

        # Filter 1: Baseline threshold
        # Controlled by min_baseline_sds config value (0 = disabled)
        # Peak height must exceed: baseline + min_baseline_sds * noise_std
        # When min_baseline_sds=0, threshold=background, but we use > so peaks at baseline pass
        if min_baseline_sds > 0:
            for p in peak_candidates:
                if p['height'] <= baseline_threshold:
                    p['rejection_reason'] = RejectionReason.BASELINE

        # Filter 2: Significance test (primary "above background" filter)
        # Uses p-value from Negative Binomial test with dispersion estimation
        for p in peak_candidates:
            if p['rejection_reason'] != RejectionReason.NONE:
                continue  # Already rejected
            if not p['is_significant']:  # p_value >= alpha
                p['rejection_reason'] = RejectionReason.SIGNIFICANCE

        # Get candidates that passed significance test for SNR threshold calculation
        significance_passed = [p for p in peak_candidates if p['rejection_reason'] == RejectionReason.NONE]

        # Filter 3: Adaptive local SNR threshold (data-driven quality)
        # Controlled by min_snr config value (0 = threshold equals minimum, all pass)
        # Threshold = min(SNRs) + min_snr * std(SNRs)
        if significance_passed:
            snr_values = [p['local_snr'] for p in significance_passed]
            snr_min = np.min(snr_values)
            snr_std = np.std(snr_values)
            snr_threshold = snr_min + (min_snr * snr_std)

            for p in significance_passed:
                if p['local_snr'] < snr_threshold:
                    p['rejection_reason'] = RejectionReason.SNR

        # Get candidates that passed SNR for prominence threshold calculation
        snr_passed = [p for p in peak_candidates if p['rejection_reason'] == RejectionReason.NONE]

        # Filter 4: Prominence percentile threshold
        # Controlled by prominence_percentile config value (0 = all pass)
        # Retains top (1 - percentile) peaks by prominence
        if snr_passed:
            prominences = [p['prominence'] for p in snr_passed]
            prominence_cutoff = np.percentile(prominences, prominence_percentile * 100)

            for p in snr_passed:
                if p['prominence'] < prominence_cutoff:
                    p['rejection_reason'] = RejectionReason.PROMINENCE

        # Filter 5: Verify local maximum (higher than or equal to both neighbors)
        for p in peak_candidates:
            if p['rejection_reason'] != RejectionReason.NONE:
                continue  # Already rejected

            idx = p['idx']

            # Check left neighbor
            if idx > 0:
                left_neighbor = signal[idx - 1]
                if signal[idx] < left_neighbor:
                    p['rejection_reason'] = RejectionReason.NOT_MAXIMUM
                    continue

            # Check right neighbor
            if idx < len(signal) - 1:
                right_neighbor = signal[idx + 1]
                if signal[idx] < right_neighbor:
                    p['rejection_reason'] = RejectionReason.NOT_MAXIMUM

        # Phase 2: Build Peak objects
        peaks = []
        for p in peak_candidates:
            # Skip rejected peaks unless include_rejected is True
            if p['rejection_reason'] != RejectionReason.NONE and not include_rejected:
                continue

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
                prominence=p['prominence'],
                rejection_reason=p['rejection_reason'],
                p_value=p['p_value']
            )
            peaks.append(peak)

        return peaks

    def _estimate_background(self, signal: NDArray[np.float64]) -> float:
        """
        Estimate background level using configured baseline method.

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
        Delegates to BaselineEstimatorService with configured method.
        Default method is SIGMA_CLIP (2σ = 95% confidence interval).

        References
        ----------
        THEORY.md Section 5.1: Background Estimation
        """
        return self._baseline_estimator.estimate(signal)

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

    def _compute_local_baseline(
        self,
        signal: NDArray[np.float64],
        left_valley_idx: int | None,
        right_valley_idx: int | None,
        left_boundary: int,
        right_boundary: int,
        global_background: float
    ) -> float:
        """
        Compute local baseline from valley heights.

        Uses structural features (valleys) to estimate local baseline,
        handling chromatographic drift naturally. This is principled because
        valleys are topological critical points where signal returns to baseline.

        Parameters
        ----------
        signal : NDArray[np.float64]
            Signal array
        left_valley_idx : int or None
            Left valley index
        right_valley_idx : int or None
            Right valley index
        left_boundary : int
            Left boundary index (fallback if no valley)
        right_boundary : int
            Right boundary index (fallback if no valley)
        global_background : float
            Global background estimate (fallback if no valleys/boundaries)

        Returns
        -------
        float
            Local baseline estimate

        Notes
        -----
        Priority order:
        1. Max of valley heights (if both valleys found)
        2. Single valley height (if one found)
        3. Max of boundary heights (if no valleys)
        4. Global background (ultimate fallback)

        Uses max (not mean) as conservative estimate - ensures we don't
        underestimate baseline and create false positives.

        References
        ----------
        THEORY.md Section 5.2: Area-based significance uses local baseline
        """
        valley_values = []

        if left_valley_idx is not None:
            valley_values.append(signal[left_valley_idx])
        if right_valley_idx is not None:
            valley_values.append(signal[right_valley_idx])

        if valley_values:
            # Use max of valley heights (conservative)
            return float(max(valley_values))

        # Fallback: use boundary values
        boundary_values = [signal[left_boundary], signal[right_boundary]]
        if any(v > 0 for v in boundary_values):
            return float(max(boundary_values))

        # Ultimate fallback: global background
        return global_background

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
        alpha: float,
        dispersion: float,
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
            Background level for significance calculation
        alpha : float
            Significance level (false positive rate)
        dispersion : float
            NB dispersion parameter r for significance testing
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
        1. Points with signal >= significance threshold
        2. NOT within min_separation of existing local maxima
        3. Local maxima in the signal derivative (inflection points)

        This handles overlapping peaks from co-eluting compounds, which is
        common in LC-MS where multiple peptides may have similar retention times.
        """
        if len(signal) < 5:
            return []

        # OPTIMIZATION: Use np.gradient instead of loop (50-100x faster)
        derivative = np.gradient(signal)

        # OPTIMIZATION: Vectorized shoulder detection
        # Build boolean masks for all conditions at once
        n = len(signal)
        indices = np.arange(2, n - 2)

        # Significance filter (vectorized)
        # Compute the minimum count needed for significance at given alpha
        significance_threshold = self._significance_tester.compute_threshold(background, alpha, dispersion)
        sig_mask = signal[2:-2] >= significance_threshold

        # Proximity filter: check distance to existing maxima (vectorized)
        if existing_maxima:
            existing_arr = np.array(existing_maxima)
            # For each index, check if ANY existing maximum is within min_separation
            # Shape: (len(indices), len(existing_maxima))
            distances = np.abs(indices[:, np.newaxis] - existing_arr)
            proximity_mask = ~np.any(distances <= min_separation, axis=1)
        else:
            proximity_mask = np.ones(len(indices), dtype=bool)

        # Derivative maxima filter (vectorized)
        deriv_left = derivative[1:-3]
        deriv_center = derivative[2:-2]
        deriv_right = derivative[3:-1]
        is_derivative_max = (deriv_center > deriv_left) & (deriv_center >= deriv_right)

        # Plateau filter (vectorized)
        signal_m2 = signal[:-4]
        signal_center = signal[2:-2]
        signal_p2 = signal[4:]
        is_plateau = (np.abs(deriv_center) < 0.5) & (signal_center > signal_m2) & (signal_center > signal_p2)

        # Combine all conditions
        valid_mask = sig_mask & proximity_mask & (is_derivative_max | is_plateau)

        return indices[valid_mask].tolist()

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

        # OPTIMIZATION: Vectorized local maxima detection (100-200x faster)
        # Extract left, center, right views
        left = signal[:-2]
        center = signal[1:-1]
        right = signal[2:]

        # Basic condition: center > left AND center >= right
        basic_max = (center > left) & (center >= right)

        # Strict maximum: center > right (not equal)
        is_strict = center > right

        # Flat peak handling: center == right AND right > next
        # Need to check signal[i+2] for flat peaks
        is_flat = center == right
        # For flat peaks, check if next value descends (right > signal[i+2])
        # Pad with False for last position where i+2 would be out of bounds
        flat_then_descends = np.zeros(len(center), dtype=bool)
        if len(signal) > 3:
            flat_then_descends[:-1] = is_flat[:-1] & (right[:-1] > signal[3:])

        # Combine: basic_max AND (strict OR flat_then_descends)
        is_maxima = basic_max & (is_strict | flat_then_descends)

        # Get indices (add 1 because center starts at index 1)
        return (np.where(is_maxima)[0] + 1).tolist()

    def _find_peak_boundaries(
        self,
        signal: NDArray[np.float64],
        peak_idx: int,
    ) -> Tuple[int, int]:
        """
        Find peak boundaries using valley detection only.

        Boundaries are defined by structural features of the signal (valleys),
        not arbitrary thresholds. This is topologically clean and avoids
        magic numbers.

        Parameters
        ----------
        signal : NDArray[np.float64]
            Signal array
        peak_idx : int
            Index of peak maximum

        Returns
        -------
        Tuple[int, int]
            (left_boundary_idx, right_boundary_idx)
            Boundaries are at valleys or signal edges.

        Notes
        -----
        Valley-only boundary detection:
        1. Scan outward from peak until valley (local minimum) is found
        2. If no valley, boundary extends to signal edge
        3. Uses <= on one side of valley check to handle flat regions

        This approach:
        - Uses structural features (valleys) rather than arbitrary thresholds
        - Produces boundaries that are reproducible features of the signal
        - Relies on prominence filtering to reject peaks with poor boundaries

        A valley is defined as a point where:
        - signal[i] <= signal[i-1] AND signal[i] < signal[i+1] (left scan)
        - signal[i] <= signal[i+1] AND signal[i] < signal[i-1] (right scan)

        The asymmetric inequality (<=, <) handles flat regions (plateaus)
        by selecting the first point of a flat valley.

        References
        ----------
        THEORY.md Section 5.2.4: Peak Boundary Determination
        Discrete Morse Theory: Valleys as critical points (saddles in 1D)
        """
        n = len(signal)

        # Find left boundary: scan left until valley or edge
        left_idx = 0  # Default to signal edge
        for i in range(peak_idx - 1, 0, -1):  # Stop at 1, not 0 (need i-1)
            # Valley: local minimum (use <= on left to handle flats)
            if signal[i] <= signal[i - 1] and signal[i] < signal[i + 1]:
                left_idx = i
                break
        else:
            # No valley found - check if index 0 is a valley or use edge
            if peak_idx > 0 and signal[0] < signal[1]:
                left_idx = 0  # Index 0 is a local minimum
            else:
                left_idx = 0  # Use signal edge

        # Find right boundary: scan right until valley or edge
        right_idx = n - 1  # Default to signal edge
        for i in range(peak_idx + 1, n - 1):  # Stop at n-2 (need i+1)
            # Valley: local minimum (use <= on right to handle flats)
            if signal[i] <= signal[i + 1] and signal[i] < signal[i - 1]:
                right_idx = i
                break
        else:
            # No valley found - check if last index is a valley or use edge
            if peak_idx < n - 1 and signal[n - 1] < signal[n - 2]:
                right_idx = n - 1  # Last index is a local minimum
            else:
                right_idx = n - 1  # Use signal edge

        # Ensure valid boundaries: left < peak < right
        # This should naturally hold for proper peaks with valleys
        # Edge case handling for peaks at signal boundaries
        if left_idx >= peak_idx:
            left_idx = max(0, peak_idx - 1)
        if right_idx <= peak_idx:
            right_idx = min(n - 1, peak_idx + 1)

        # Final check: ensure left < right
        if left_idx >= right_idx:
            # Pathological case - expand minimally
            left_idx = max(0, peak_idx - 1)
            right_idx = min(n - 1, peak_idx + 1)
            if left_idx >= right_idx:
                # Very short signal
                left_idx = 0
                right_idx = n - 1

        return left_idx, right_idx

    def _resolve_overlapping_boundaries(
        self,
        peak_candidates: List[dict],
        signal: NDArray[np.float64],
        time_points: NDArray[np.float64]
    ) -> None:
        """
        Resolve overlapping boundaries between adjacent peaks.

        When two peaks are close together without a clear valley between them,
        their boundaries may overlap. This method finds the minimum point
        between adjacent peaks and uses it as the shared boundary.

        Parameters
        ----------
        peak_candidates : List[dict]
            List of peak candidate dictionaries (modified in place)
        signal : NDArray[np.float64]
            Signal array
        time_points : NDArray[np.float64]
            Time array in seconds

        Notes
        -----
        For each pair of adjacent peaks (sorted by position):
        1. Check if right boundary of left peak overlaps left boundary of right peak
        2. If overlap exists, find the minimum point between the two peak apexes
        3. Set boundaries to meet at this minimum point

        This ensures:
        - No boundary overlaps exist in the final output
        - Boundaries are at structurally meaningful points (local minima)
        - Each point in the signal belongs to at most one peak
        """
        if len(peak_candidates) < 2:
            return

        # Sort candidates by peak position (index)
        sorted_candidates = sorted(peak_candidates, key=lambda p: p['idx'])

        # Process adjacent pairs
        for i in range(len(sorted_candidates) - 1):
            left_peak = sorted_candidates[i]
            right_peak = sorted_candidates[i + 1]

            # Check for overlap: left peak's right boundary >= right peak's left boundary
            if left_peak['right_idx'] >= right_peak['left_idx']:
                # Find minimum point between the two peak apexes
                left_apex = left_peak['idx']
                right_apex = right_peak['idx']

                # Find the minimum value between the apexes
                if left_apex < right_apex:
                    between_signal = signal[left_apex:right_apex + 1]
                    min_offset = np.argmin(between_signal)
                    min_idx = left_apex + min_offset
                else:
                    # Peaks at same position (shouldn't happen, but handle gracefully)
                    min_idx = left_apex

                # Set boundaries to meet at the minimum
                left_peak['right_idx'] = min_idx
                left_peak['right_base'] = time_points[min_idx]

                right_peak['left_idx'] = min_idx
                right_peak['left_base'] = time_points[min_idx]

                # Recalculate areas with new boundaries
                left_peak['area'] = self._integrate_peak_area(
                    signal, left_peak['left_idx'], left_peak['right_idx']
                )
                right_peak['area'] = self._integrate_peak_area(
                    signal, right_peak['left_idx'], right_peak['right_idx']
                )

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
