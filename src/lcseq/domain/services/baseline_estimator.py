"""
Baseline estimation service with swappable methods.

Provides configurable baseline estimation for peak detection.
Each method has a clear statistical or physical justification.

Includes:
- Scalar baseline estimation (constant baseline assumption)
- Curve-based baseline fitting via pybaselines (handles drift)

References
----------
THEORY.md Section 5.1: Background Estimation

pybaselines Reference:
    Erb, D. (2022). pybaselines: A Python library of algorithms for
    the baseline correction of experimental data.
    https://github.com/derb12/pybaselines
"""

from dataclasses import dataclass, field
from enum import Enum
from typing import Optional

import numpy as np
from numpy.typing import NDArray
from scipy.signal import argrelmin
from scipy.interpolate import UnivariateSpline
from pybaselines import Baseline


class BaselineEstimatorService:
    """
    Estimate baseline using sigma-clipping.

    Uses iterative σ-clipping outlier removal to identify and exclude peaks,
    then estimates baseline from remaining points.

    The default σ=2.0 is statistically principled - it corresponds to the 95%
    confidence interval, the standard threshold for outlier detection.

    Notes
    -----
    Sigma-clipping is commonly used in astronomy for baseline estimation
    and is well-established in signal processing literature.

    References
    ----------
    THEORY.md Section 5.1: Background Estimation
    """

    def __init__(self, sigma: float = 2.0):
        """
        Initialize baseline estimator.

        Parameters
        ----------
        sigma : float
            Number of standard deviations for sigma-clipping threshold.
            Default 2.0 corresponds to 95% confidence interval.
            - 1σ = 68.27% CI (too permissive, includes many peaks)
            - 2σ = 95.45% CI (standard outlier threshold)
            - 3σ = 99.73% CI (too strict, may mask legitimate baseline)
        """
        self._sigma = sigma

    def estimate(self, signal: NDArray[np.float64]) -> float:
        """
        Estimate baseline level from signal using sigma-clipping.

        Parameters
        ----------
        signal : NDArray[np.float64]
            Signal array (counts)

        Returns
        -------
        float
            Estimated baseline level

        Raises
        ------
        ValueError
            If signal is empty
        """
        if len(signal) == 0:
            raise ValueError("Cannot estimate baseline from empty signal")

        return self._sigma_clip(signal)

    def _sigma_clip(
        self,
        signal: NDArray[np.float64],
        max_iter: int = 10,
    ) -> float:
        """
        Iterative sigma-clipping baseline estimation.

        Removes outliers (peaks) iteratively until convergence.
        Uses self._sigma for clipping threshold.

        Parameters
        ----------
        signal : NDArray[np.float64]
            Signal array
        max_iter : int
            Maximum iterations (default: 10)

        Returns
        -------
        float
            Median of remaining (non-peak) points

        Notes
        -----
        Iteratively masks points above mean + sigma*std until no more
        points are removed or max_iter is reached.

        Uses median of remaining points (more robust than mean) as
        the final baseline estimate.
        """
        mask = np.ones(len(signal), dtype=bool)

        for _ in range(max_iter):
            masked_signal = signal[mask]
            if len(masked_signal) == 0:
                break

            mean = np.mean(masked_signal)
            std = np.std(masked_signal)

            # Threshold: keep points within mean + sigma*std
            # Only clip above (peaks are positive outliers)
            new_mask = signal <= (mean + self._sigma * std)

            # Check for convergence
            if np.array_equal(mask, new_mask):
                break

            mask = new_mask

        # Return median of remaining points (robust to any remaining outliers)
        remaining = signal[mask]
        if len(remaining) == 0:
            # Fallback if all points were clipped (shouldn't happen with σ=2)
            return float(np.median(signal))

        return float(np.median(remaining))

    def estimate_with_noise(
        self, signal: NDArray[np.float64]
    ) -> tuple[float, float]:
        """
        Estimate baseline level and noise standard deviation.

        Parameters
        ----------
        signal : NDArray[np.float64]
            Signal array (counts)

        Returns
        -------
        tuple[float, float]
            (baseline, noise_std) where:
            - baseline: Estimated baseline level
            - noise_std: Standard deviation of baseline region (excluding peaks)

        Notes
        -----
        Uses sigma-clipping to identify baseline points, then computes
        both the median (baseline) and std (noise) of those points.

        This provides a consistent baseline and noise estimate that
        excludes peak contributions.
        """
        if len(signal) == 0:
            raise ValueError("Cannot estimate baseline from empty signal")

        # Use sigma-clipping to identify baseline points (regardless of method)
        # This ensures noise_std is always computed from non-peak regions
        mask = np.ones(len(signal), dtype=bool)
        max_iter = 10

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

        # Compute baseline and noise from remaining (non-peak) points
        remaining = signal[mask]
        if len(remaining) == 0:
            return float(np.median(signal)), float(np.std(signal))

        baseline = float(np.median(remaining))
        noise_std = float(np.std(remaining))

        return baseline, noise_std

    def estimate_with_noise_and_dispersion(
        self, signal: NDArray[np.float64], min_dispersion_r: float
    ) -> tuple[float, float, float]:
        """
        Estimate baseline level, noise std, and NB dispersion parameter.

        This method extends estimate_with_noise() to also compute the
        Negative Binomial dispersion parameter r, enabling more principled
        significance testing that properly handles overdispersion.

        Parameters
        ----------
        signal : NDArray[np.float64]
            Signal array (counts)
        min_dispersion_r : float, optional
            Minimum dispersion parameter floor. Prevents numerical issues
            when r is very small. Default 0.1.

        Returns
        -------
        tuple[float, float, float]
            (baseline, noise_std, dispersion_r) where:
            - baseline: Estimated baseline level
            - noise_std: Standard deviation of baseline region (excluding peaks)
            - dispersion_r: NB dispersion parameter (larger = closer to Poisson)

        Notes
        -----
        The dispersion parameter r is estimated from the baseline region using
        the method-of-moments estimator. For Negative Binomial:

            Var = μ + μ²/r

        Solving for r:

            r = μ² / (Var - μ)

        When Var ≈ μ (Poisson-like), r → ∞ (we return 1e6).
        When Var > μ (overdispersed), r reflects the degree of overdispersion.

        This is the standard approach used in genomics (DESeq2, edgeR) for
        handling overdispersed count data.

        References
        ----------
        Robinson, M.D. and Smyth, G.K. (2007). Moderated statistical tests for
        assessing differences in tag abundance. Bioinformatics 23(21):2881-2887.

        Anders, S. and Huber, W. (2010). Differential expression analysis for
        sequence count data. Genome Biology 11:R106.
        """
        if len(signal) == 0:
            raise ValueError("Cannot estimate baseline from empty signal")

        # Use sigma-clipping to identify baseline points (regardless of method)
        # This ensures estimates are computed from non-peak regions
        mask = np.ones(len(signal), dtype=bool)
        max_iter = 10

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

        # Compute baseline and noise from remaining (non-peak) points
        remaining = signal[mask]
        if len(remaining) == 0:
            return float(np.median(signal)), float(np.std(signal)), 1e6

        baseline = float(np.median(remaining))
        noise_std = float(np.std(remaining))

        # Compute NB dispersion parameter from baseline values
        # From Var = μ + μ²/r, solving for r: r = μ² / (Var - μ)
        mu = np.mean(remaining)
        var = np.var(remaining, ddof=1)  # Sample variance (unbiased)

        if var <= mu or mu <= 0:
            # No overdispersion detected (or invalid data) → Poisson limit
            dispersion_r = 1e6
        else:
            dispersion_r = mu**2 / (var - mu)
            # Apply minimum floor for numerical stability
            dispersion_r = max(dispersion_r, min_dispersion_r)

        return baseline, noise_std, float(dispersion_r)


# =============================================================================
# Baseline Curve Fitting via pybaselines
# =============================================================================


class BaselineCurveMethod(Enum):
    """Available baseline curve fitting methods from pybaselines.

    Each method has different characteristics:
    - ASLS: Standard asymmetric least squares
    - AIRPLS: Adaptive iteratively reweighted PLS (better edge handling)
    - ARPLS: Asymmetrically reweighted PLS (automatic weight adjustment)
    - SNIP: Statistics-sensitive non-linear iterative peak-clipping
    - MOR: Morphological baseline (rolling ball equivalent)
    - IMODPOLY: Improved modified polynomial
    """

    ASLS = "asls"
    AIRPLS = "airpls"
    ARPLS = "arpls"
    SNIP = "snip"
    MOR = "mor"
    IMODPOLY = "imodpoly"


@dataclass
class BaselineCurveParams:
    """Parameters for baseline curve fitting.

    Attributes
    ----------
    method : BaselineCurveMethod
        Algorithm to use. Default AIRPLS handles edges well.
    lam : float
        Smoothness parameter (λ) for ALS-based methods.
        Larger = smoother. Typical range: 1e3 to 1e7.
    p : float
        Asymmetry parameter for ASLS. Range: 0 < p < 0.5.
        Smaller = baseline pushed further below peaks.
    max_half_window : int
        Half-window size for SNIP and MOR methods.
    poly_order : int
        Polynomial order for IMODPOLY.
    max_iter : int
        Maximum iterations.
    tol : float
        Convergence tolerance.
    """

    method: BaselineCurveMethod = BaselineCurveMethod.AIRPLS
    lam: float = 1e5
    p: float = 0.01
    max_half_window: int = 30
    poly_order: int = 2
    max_iter: int = 50
    tol: float = 1e-3


@dataclass
class BaselineCurveResult:
    """Result from baseline curve fitting.

    Attributes
    ----------
    baseline : NDArray[np.float64]
        Fitted baseline curve.
    corrected : NDArray[np.float64]
        Baseline-corrected signal (original - baseline).
    params : dict
        Parameters returned by the fitting algorithm.
    method : str
        Method used for fitting.
    """

    baseline: NDArray[np.float64]
    corrected: NDArray[np.float64]
    params: dict
    method: str


def fit_baseline_curve(
    signal: NDArray[np.float64],
    x: Optional[NDArray[np.float64]] = None,
    params: Optional[BaselineCurveParams] = None,
) -> BaselineCurveResult:
    """
    Fit a baseline curve using pybaselines.

    Parameters
    ----------
    signal : NDArray[np.float64]
        Input signal (1D array)
    x : NDArray[np.float64], optional
        X-axis values (e.g., time points). If None, uses indices.
    params : BaselineCurveParams, optional
        Fitting parameters. Uses defaults if not provided.

    Returns
    -------
    BaselineCurveResult
        Contains baseline, corrected signal, and fit info.

    Notes
    -----
    Recommended methods:
    - AIRPLS: Good default, handles edges well, adaptive
    - ARPLS: Automatic weight selection, robust
    - SNIP: No smoothness parameter, good for well-separated peaks
    - MOR: Fast morphological approach
    """
    if params is None:
        params = BaselineCurveParams()

    # Create baseline fitter with x-axis if provided
    fitter = Baseline(x_data=x)

    # Dispatch to appropriate method
    if params.method == BaselineCurveMethod.ASLS:
        baseline, fit_params = fitter.asls(
            signal,
            lam=params.lam,
            p=params.p,
            max_iter=params.max_iter,
            tol=params.tol,
        )
    elif params.method == BaselineCurveMethod.AIRPLS:
        baseline, fit_params = fitter.airpls(
            signal,
            lam=params.lam,
            max_iter=params.max_iter,
            tol=params.tol,
        )
    elif params.method == BaselineCurveMethod.ARPLS:
        baseline, fit_params = fitter.arpls(
            signal,
            lam=params.lam,
            max_iter=params.max_iter,
            tol=params.tol,
        )
    elif params.method == BaselineCurveMethod.SNIP:
        baseline, fit_params = fitter.snip(
            signal,
            max_half_window=params.max_half_window,
        )
    elif params.method == BaselineCurveMethod.MOR:
        baseline, fit_params = fitter.mor(
            signal,
            half_window=params.max_half_window,
        )
    elif params.method == BaselineCurveMethod.IMODPOLY:
        baseline, fit_params = fitter.imodpoly(
            signal,
            poly_order=params.poly_order,
            max_iter=params.max_iter,
            tol=params.tol,
        )
    else:
        raise ValueError(f"Unknown method: {params.method}")

    return BaselineCurveResult(
        baseline=baseline,
        corrected=signal - baseline,
        params=fit_params,
        method=params.method.value,
    )


def compare_baseline_methods(
    signal: NDArray[np.float64],
    x: Optional[NDArray[np.float64]] = None,
    methods: Optional[list[BaselineCurveMethod]] = None,
    lam: float = 1e5,
) -> dict[str, BaselineCurveResult]:
    """
    Compare multiple baseline fitting methods on the same signal.

    Parameters
    ----------
    signal : NDArray[np.float64]
        Input signal
    x : NDArray[np.float64], optional
        X-axis values
    methods : list[BaselineCurveMethod], optional
        Methods to compare. Defaults to [AIRPLS, ARPLS, SNIP].
    lam : float
        Smoothness parameter for ALS-based methods.

    Returns
    -------
    dict[str, BaselineCurveResult]
        Results keyed by method name.
    """
    if methods is None:
        methods = [
            BaselineCurveMethod.AIRPLS,
            BaselineCurveMethod.ARPLS,
            BaselineCurveMethod.SNIP,
        ]

    results = {}
    for method in methods:
        params = BaselineCurveParams(method=method, lam=lam)
        results[method.value] = fit_baseline_curve(signal, x, params)

    return results


# =============================================================================
# Local-Minima Piecewise ALS (from SERDS paper)
# =============================================================================


@dataclass
class PiecewiseALSParams:
    """Parameters for local-minima piecewise ALS baseline fitting.

    This approach splits the signal at local minima (valleys between peaks)
    rather than at fixed intervals. This is more principled for chromatography
    because local minima are where the baseline is most clearly visible.

    Reference:
        Zhao et al. (2015). A shifted-excitation Raman difference spectroscopy
        (SERDS) evaluation strategy for the efficient isolation of Raman spectra
        from extreme fluorescence interference.

    Attributes
    ----------
    max_window : int
        Maximum number of points per segment. If no local minimum is found
        within this range, the global minimum in the range is used as the
        split point. Default 200 (from the paper).
    min_window : int
        Minimum segment size. Segments smaller than this are merged with
        neighbors. Prevents over-segmentation from noise.
    lam : float
        Smoothness parameter for ALS fitting within each segment.
    p : float
        Asymmetry parameter for ALS. Smaller = baseline pushed lower.
    order : int
        Order of local minima detection (argrelmin). Higher = smoother
        minima detection, ignores small wiggles.
    method : BaselineCurveMethod
        Which pybaselines method to use within each segment.
    """

    max_window: int = 200
    min_window: int = 20
    lam: float = 1e4
    p: float = 0.01
    order: int = 5
    method: BaselineCurveMethod = BaselineCurveMethod.ASLS


@dataclass
class PiecewiseALSResult:
    """Result from local-minima piecewise ALS baseline fitting.

    Attributes
    ----------
    baseline : NDArray[np.float64]
        Final concatenated baseline curve.
    corrected : NDArray[np.float64]
        Baseline-corrected signal (original - baseline).
    split_points : list[int]
        Indices where the signal was split (local minima).
    segment_baselines : list[tuple[int, int, NDArray[np.float64]]]
        Per-segment results as (start_idx, end_idx, baseline) tuples.
    """

    baseline: NDArray[np.float64]
    corrected: NDArray[np.float64]
    split_points: list[int]
    segment_baselines: list[tuple[int, int, NDArray[np.float64]]]


def find_split_points(
    signal: NDArray[np.float64],
    max_window: int = 200,
    min_window: int = 20,
    order: int = 5,
) -> list[int]:
    """
    Find split points at local minima for piecewise baseline fitting.

    The algorithm:
    1. Find all local minima in the signal
    2. Starting from index 0, find the first local minimum within max_window
    3. If no local minimum exists, use the global minimum in that range
    4. Repeat from that point until the end of the signal
    5. Merge segments that are too small (< min_window)

    Parameters
    ----------
    signal : NDArray[np.float64]
        Input signal
    max_window : int
        Maximum distance to search for a local minimum
    min_window : int
        Minimum segment size (smaller segments are merged)
    order : int
        Order parameter for scipy.signal.argrelmin

    Returns
    -------
    list[int]
        Indices of split points (local minima between segments)
    """
    n = len(signal)
    if n <= max_window:
        return []  # Signal too short, no splitting needed

    # Find all local minima
    local_min_indices = argrelmin(signal, order=order)[0]
    local_min_set = set(local_min_indices)

    # Build split points iteratively
    split_points = []
    current_pos = 0

    while current_pos < n - min_window:
        # Define search window
        window_end = min(current_pos + max_window, n)

        # Look for local minima in this window (excluding current_pos itself)
        candidates = [
            idx for idx in local_min_indices
            if current_pos + min_window <= idx < window_end
        ]

        if candidates:
            # Use the first local minimum in the window
            split_idx = candidates[0]
        else:
            # No local minimum found - use global minimum in window
            search_start = current_pos + min_window
            if search_start >= window_end:
                break
            window_signal = signal[search_start:window_end]
            split_idx = search_start + np.argmin(window_signal)

        # Don't add split point too close to the end
        if split_idx >= n - min_window:
            break

        split_points.append(split_idx)
        current_pos = split_idx

    return split_points


def piecewise_als_baseline(
    signal: NDArray[np.float64],
    x: Optional[NDArray[np.float64]] = None,
    params: Optional[PiecewiseALSParams] = None,
) -> PiecewiseALSResult:
    """
    Local-minima piecewise ALS baseline fitting.

    Splits the signal at local minima (valleys between peaks) and fits
    ALS baseline independently to each segment. This is more principled
    than fixed-window approaches because:

    1. Local minima are where baseline is most clearly visible
    2. Each segment naturally contains one or a few peaks
    3. ALS within each segment isn't "pulled" by distant peaks

    Parameters
    ----------
    signal : NDArray[np.float64]
        Input signal (1D array)
    x : NDArray[np.float64], optional
        X-axis values (e.g., time points). Used for proper interpolation.
    params : PiecewiseALSParams, optional
        Fitting parameters. Uses defaults if not provided.

    Returns
    -------
    PiecewiseALSResult
        Contains baseline, corrected signal, split points, and segment info.

    References
    ----------
    Zhao et al. (2015). A shifted-excitation Raman difference spectroscopy
    (SERDS) evaluation strategy...
    """
    if params is None:
        params = PiecewiseALSParams()

    n = len(signal)

    # Find split points at local minima
    split_points = find_split_points(
        signal,
        max_window=params.max_window,
        min_window=params.min_window,
        order=params.order,
    )

    # Create segment boundaries: [0, split1, split2, ..., n]
    boundaries = [0] + split_points + [n]

    # Fit baseline to each segment
    baseline = np.zeros(n)
    segment_baselines = []

    for i in range(len(boundaries) - 1):
        start = boundaries[i]
        end = boundaries[i + 1]

        segment_signal = signal[start:end]

        # Create new fitter for each segment (pybaselines caches x_data length)
        fitter = Baseline(x_data=np.arange(len(segment_signal)))

        # Fit baseline to this segment
        if params.method == BaselineCurveMethod.ASLS:
            seg_baseline, _ = fitter.asls(
                segment_signal,
                lam=params.lam,
                p=params.p,
            )
        elif params.method == BaselineCurveMethod.AIRPLS:
            seg_baseline, _ = fitter.airpls(
                segment_signal,
                lam=params.lam,
            )
        elif params.method == BaselineCurveMethod.ARPLS:
            seg_baseline, _ = fitter.arpls(
                segment_signal,
                lam=params.lam,
            )
        else:
            # Default to ASLS
            seg_baseline, _ = fitter.asls(
                segment_signal,
                lam=params.lam,
                p=params.p,
            )

        baseline[start:end] = seg_baseline
        segment_baselines.append((start, end, seg_baseline))

    return PiecewiseALSResult(
        baseline=baseline,
        corrected=signal - baseline,
        split_points=split_points,
        segment_baselines=segment_baselines,
    )


# =============================================================================
# Local Minima + Smoothed Spline Baseline
# =============================================================================


@dataclass
class MinimaSplineParams:
    """Parameters for local-minima spline baseline fitting.

    This approach finds local minima (valleys between peaks) and fits a
    smoothed spline only through those points. Since peaks are not included
    in the fit, the baseline cannot cut into them.

    Attributes
    ----------
    order : int
        Order parameter for scipy.signal.argrelmin. Higher values detect
        fewer, more significant minima. Default 5.
    smoothing : float or None
        Spline smoothing factor for UnivariateSpline. None = automatic
        based on number of points. 0 = interpolate exactly through all
        minima. Higher values = smoother approximation.
    spline_degree : int
        Degree of the spline. 3 = cubic (smooth, default), 1 = linear.
    include_start : bool
        Whether to include the first signal point as a virtual minimum.
        Useful for capturing true baseline before the null peak.
    include_end : bool
        Whether to include the last signal point as a virtual minimum.
        Usually False since the last minimum better represents baseline.
    min_value : float or None
        Minimum allowed baseline value. None = no constraint.
        0.0 = non-negative baseline.
    max_iterations : int
        Maximum iterations for constrained fitting. The algorithm
        iteratively adds constraint points where spline > signal.
    constraint_tolerance : float
        Stop iterating when max overshoot is below this value.
    """

    order: int = 5
    smoothing: Optional[float] = None
    spline_degree: int = 3
    include_start: bool = True
    include_end: bool = True
    min_value: Optional[float] = 0.0
    max_iterations: int = 10
    constraint_tolerance: float = 1.0


@dataclass
class MinimaSplineResult:
    """Result from local-minima spline baseline fitting.

    Attributes
    ----------
    baseline : NDArray[np.float64]
        Fitted baseline curve.
    corrected : NDArray[np.float64]
        Baseline-corrected signal (signal - baseline).
    original : NDArray[np.float64]
        Original input signal.
    x : NDArray[np.float64]
        X-axis for signals (time points).
    minima_indices : NDArray[np.int64]
        Indices of detected local minima.
    minima_x : NDArray[np.float64]
        X-coordinates of minima points used for spline fitting.
    minima_y : NDArray[np.float64]
        Y-coordinates (signal values) of minima points.
    """

    baseline: NDArray[np.float64]
    corrected: NDArray[np.float64]
    original: NDArray[np.float64]
    x: NDArray[np.float64]
    minima_indices: NDArray[np.int64]
    minima_x: NDArray[np.float64]
    minima_y: NDArray[np.float64]


def minima_spline_baseline(
    signal: NDArray[np.float64],
    x: Optional[NDArray[np.float64]] = None,
    params: Optional[MinimaSplineParams] = None,
) -> MinimaSplineResult:
    """
    Baseline estimation via smoothed spline through local minima.

    This approach:
    1. Find local minima (valleys between peaks) in the signal
    2. Fit a smoothed spline only through those minima points
    3. Subtract baseline from signal

    Since peaks are not included in the fit, the baseline cannot be
    "pulled up" into peaks - it naturally follows the lower envelope.

    Parameters
    ----------
    signal : NDArray[np.float64]
        Input signal (1D array)
    x : NDArray[np.float64], optional
        X-axis values (e.g., time points). If None, uses indices.
    params : MinimaSplineParams, optional
        Fitting parameters. Uses defaults if not provided.

    Returns
    -------
    MinimaSplineResult
        Contains original, baseline, corrected signals and minima info.

    Notes
    -----
    The smoothing parameter controls the trade-off between:
    - s=0: Exact interpolation through all minima (may be wiggly if noisy)
    - s>0: Smooth approximation (may not pass through all minima exactly)

    For noisy signals, use higher smoothing or increase the order parameter
    to detect fewer, more significant minima.
    """
    if params is None:
        params = MinimaSplineParams()

    n = len(signal)
    if x is None:
        x = np.arange(n, dtype=np.float64)

    # Store original signal
    original = signal.copy()

    # Find local minima in the signal
    minima_idx = argrelmin(signal, order=params.order)[0]

    # Include start/end points if requested
    extra_points = []
    if params.include_start and 0 not in minima_idx:
        extra_points.append(0)
    if params.include_end and n - 1 not in minima_idx:
        extra_points.append(n - 1)
    if extra_points:
        minima_idx = np.concatenate([extra_points, minima_idx])
        minima_idx = np.sort(minima_idx)

    # Start with minima as constraint points
    constraint_idx = minima_idx.copy()

    # Iteratively fit spline and add constraint points where it overshoots
    for iteration in range(params.max_iterations):
        # Handle edge case: too few points for cubic spline
        if len(constraint_idx) < params.spline_degree + 1:
            if len(constraint_idx) >= 2:
                degree = len(constraint_idx) - 1
            else:
                # Only one point - use constant baseline
                baseline = np.full(n, signal[constraint_idx[0]] if len(constraint_idx) > 0 else np.min(signal))
                return MinimaSplineResult(
                    baseline=baseline,
                    corrected=signal - baseline,
                    original=original,
                    x=x,
                    minima_indices=minima_idx,
                    minima_x=x[minima_idx] if len(minima_idx) > 0 else np.array([]),
                    minima_y=signal[minima_idx] if len(minima_idx) > 0 else np.array([]),
                )
        else:
            degree = params.spline_degree

        # Extract constraint point coordinates
        x_pts = x[constraint_idx]
        y_pts = signal[constraint_idx]

        # Determine smoothing factor
        if params.smoothing is None:
            smoothing = len(y_pts) * 0.1
        else:
            smoothing = params.smoothing

        # Fit spline through constraint points
        spline = UnivariateSpline(x_pts, y_pts, s=smoothing, k=degree)

        # Evaluate at all x positions
        baseline = spline(x)

        # Handle edge extrapolation: use constant value outside minima range
        # This prevents wild spline extrapolation at edges
        if len(x_pts) >= 2:
            # Constant extrapolation at start (if start not included)
            if not params.include_start:
                first_constraint_x = x_pts[0]
                first_constraint_val = spline(first_constraint_x)
                baseline[x < first_constraint_x] = first_constraint_val

            # Constant extrapolation at end (if end not included)
            if not params.include_end:
                last_constraint_x = x_pts[-1]
                last_constraint_val = spline(last_constraint_x)
                baseline[x > last_constraint_x] = last_constraint_val

        # Apply minimum value constraint if specified
        if params.min_value is not None:
            baseline = np.maximum(baseline, params.min_value)

        # Find where baseline exceeds signal (violations)
        overshoot = baseline - signal
        max_overshoot = np.max(overshoot)

        # Check convergence
        if max_overshoot <= params.constraint_tolerance:
            break

        # Find violation regions and add their local maxima as constraints
        # This pulls the spline down where it overshoots
        violation_mask = overshoot > params.constraint_tolerance

        if not np.any(violation_mask):
            break

        # Find local maxima of overshoot in violation regions
        # These are the worst offenders - add them as constraints
        new_constraint_idx = []

        # Find contiguous violation regions
        violation_indices = np.where(violation_mask)[0]
        if len(violation_indices) > 0:
            # Split into contiguous regions
            splits = np.where(np.diff(violation_indices) > 1)[0] + 1
            regions = np.split(violation_indices, splits)

            for region in regions:
                if len(region) > 0:
                    # Find the index with maximum overshoot in this region
                    region_overshoots = overshoot[region]
                    max_idx_in_region = region[np.argmax(region_overshoots)]
                    new_constraint_idx.append(max_idx_in_region)

        if not new_constraint_idx:
            break

        # Add new constraint points
        constraint_idx = np.unique(np.concatenate([constraint_idx, new_constraint_idx]))
        constraint_idx = np.sort(constraint_idx)

    # Final baseline - should now respect signal constraint (or be very close)
    # Apply a final soft clamp just in case
    baseline = np.minimum(baseline, signal)

    # Compute corrected signal (signal - baseline)
    corrected = signal - baseline

    return MinimaSplineResult(
        baseline=baseline,
        corrected=corrected,
        original=original,
        x=x,
        minima_indices=minima_idx,
        minima_x=x[minima_idx],
        minima_y=signal[minima_idx],
    )
