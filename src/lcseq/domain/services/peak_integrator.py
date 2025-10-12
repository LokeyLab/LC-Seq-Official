"""
Peak integration service for area calculation and boundary refinement.

Implementation based on THEORY.md Section 5.0.7, 5.2.4.
"""

import numpy as np
from numpy.typing import NDArray
from typing import Tuple
from ..entities.chromatogram import Chromatogram


class PeakIntegrator:
    """
    Integrates peak areas and determines integration boundaries.

    This service handles precise peak area integration and boundary determination,
    typically used after initial peak detection or when integrating on individual
    variants using pooled boundaries.

    This is a stateless service - all methods are operations on input data with
    no instance state.

    Notes
    -----
    Peak integration workflow:
    1. Locate peak maximum near given position
    2. Find boundaries using valley detection or threshold method
    3. Integrate area using trapezoidal rule
    4. Return boundaries and area for further analysis

    Used in:
    - Pooled mode: Integrate individual variants using pooled boundaries
    - Peak refinement: Re-integrate peaks with different boundary criteria
    - Area calculation: Compute areas for purity analysis

    References
    ----------
    THEORY.md Section 5.0.7: Peak Area Integration
    THEORY.md Section 5.2.4: Peak Boundary Determination
    THEORY.md Section 4.2.6: Area Integration on Individual Variants

    Examples
    --------
    >>> integrator = PeakIntegrator()
    >>> chromatogram = Chromatogram(time_points=[...], counts=[...])
    >>> left_base, right_base, area = integrator.integrate_peak(
    ...     chromatogram, peak_position=125.5, signal_variant="corrected"
    ... )
    """

    def integrate_peak(
        self,
        chromatogram: Chromatogram,
        peak_position: float,
        signal_variant: str = "raw",
        left_boundary: float | None = None,
        right_boundary: float | None = None
    ) -> Tuple[float, float, float]:
        """
        Integrate peak area at given position.

        Parameters
        ----------
        chromatogram : Chromatogram
            Chromatogram containing signal data
        peak_position : float
            Approximate peak retention time (seconds or minutes)
        signal_variant : str, optional
            Signal variant to use. Default is "raw" (per THEORY.md - no baseline correction).
        left_boundary : float or None, optional
            Fixed left boundary time. If None, auto-detect.
        right_boundary : float or None, optional
            Fixed right boundary time. If None, auto-detect.

        Returns
        -------
        Tuple[float, float, float]
            (left_base, right_base, area)
            - left_base: Left integration boundary (time)
            - right_base: Right integration boundary (time)
            - area: Integrated peak area

        Notes
        -----
        If boundaries are provided (e.g., from pooled detection), they are
        used directly. Otherwise, boundaries are detected using valley method
        or 5% threshold.

        Area integration uses simple summation on raw signal (per THEORY.md - no baseline correction).

        References
        ----------
        THEORY.md Section 4.2.6: Pooled Mode Area Integration
        THEORY.md Section 5.0.7: Peak Area Integration

        Examples
        --------
        >>> integrator = PeakIntegrator()
        >>> # Auto-detect boundaries
        >>> left, right, area = integrator.integrate_peak(chrom, peak_position=120.0)
        >>> # Use fixed boundaries (e.g., from pooled mode)
        >>> left, right, area = integrator.integrate_peak(
        ...     chrom, peak_position=120.0,
        ...     left_boundary=115.0, right_boundary=125.0
        ... )
        """
        if not chromatogram.has_signal_variant(signal_variant):
            raise ValueError(
                f"Signal variant '{signal_variant}' not found. "
                f"Available: {list(chromatogram.signal_variants.keys())}"
            )

        signal = chromatogram.get_signal(signal_variant)
        time_points = chromatogram.time_points

        # Find peak index closest to given position
        peak_idx = self._find_nearest_index(time_points, peak_position)

        # Determine boundaries
        if left_boundary is not None and right_boundary is not None:
            # Use provided boundaries (e.g., from pooled mode)
            left_idx = self._find_nearest_index(time_points, left_boundary)
            right_idx = self._find_nearest_index(time_points, right_boundary)
        else:
            # Auto-detect boundaries
            left_idx, right_idx = self._find_valley_boundaries(signal, peak_idx)

        # Get boundary times
        left_base = float(time_points[left_idx])
        right_base = float(time_points[right_idx])

        # Integrate area
        area = self._integrate_area(signal, left_idx, right_idx)

        return left_base, right_base, area

    def _find_valley_boundaries(
        self,
        signal: NDArray[np.float64],
        peak_idx: int,
        threshold_fraction: float = 0.05
    ) -> Tuple[int, int]:
        """
        Find peak boundaries using valley detection or 5% threshold.

        Scans outward from peak maximum until valley or threshold is reached.

        Parameters
        ----------
        signal : NDArray[np.float64]
            Signal array (raw signal per THEORY.md)
        peak_idx : int
            Index of peak maximum
        threshold_fraction : float, optional
            Fraction of peak height for fallback threshold. Default is 0.05.

        Returns
        -------
        Tuple[int, int]
            (left_boundary_idx, right_boundary_idx)

        Notes
        -----
        Boundary detection priority:
        1. Valley (local minimum) - natural separation between peaks
        2. 5% threshold - conservative fallback for overlapping peaks
        3. Signal edge - if neither valley nor threshold found

        References
        ----------
        THEORY.md Section 5.3.1: Valley Detection
        THEORY.md Section 5.3.2: 5% Threshold Method
        """
        peak_height = signal[peak_idx]
        threshold = threshold_fraction * peak_height

        # Find left boundary
        left_idx = peak_idx
        for i in range(peak_idx - 1, -1, -1):
            # Check for 5% threshold
            if signal[i] < threshold:
                left_idx = i
                break
            # Check for valley (local minimum)
            if i > 0 and signal[i] <= signal[i-1] and signal[i] <= signal[i+1]:
                left_idx = i
                break
            left_idx = i

        # Find right boundary
        right_idx = peak_idx
        for i in range(peak_idx + 1, len(signal)):
            # Check for 5% threshold
            if signal[i] < threshold:
                right_idx = i
                break
            # Check for valley (local minimum)
            if i < len(signal) - 1 and signal[i] <= signal[i-1] and signal[i] <= signal[i+1]:
                right_idx = i
                break
            right_idx = i

        return left_idx, right_idx

    def _integrate_area(
        self,
        signal: NDArray[np.float64],
        left_idx: int,
        right_idx: int
    ) -> float:
        """
        Integrate peak area between boundaries.

        Uses simple summation of signal values (appropriate for raw signals
        per THEORY.md - no baseline correction).

        Parameters
        ----------
        signal : NDArray[np.float64]
            Signal array (raw signal per THEORY.md)
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
        For raw signals (per THEORY.md - no baseline correction), simple summation
        gives accurate results. The area represents total signal intensity in the peak region.

        For purity calculations, areas from multiple peaks are combined:
            Purity = Area(product) / [Area(product) + Area(truncations) + ...]

        References
        ----------
        THEORY.md Section 5.0.7: Peak Area Integration
        """
        if right_idx <= left_idx:
            return 0.0

        # Sum signal values in integration region
        area = float(np.sum(signal[left_idx:right_idx+1]))

        # Ensure non-negative (raw signals can have noise but should be mostly positive)
        return max(0.0, area)

    def _find_nearest_index(
        self,
        time_points: NDArray[np.float64],
        target_time: float
    ) -> int:
        """
        Find index of time point closest to target.

        Parameters
        ----------
        time_points : NDArray[np.float64]
            Array of time points (strictly increasing)
        target_time : float
            Target time value

        Returns
        -------
        int
            Index of closest time point

        Notes
        -----
        Uses absolute difference to find closest point.
        For time points at equal distance, returns the earlier index.
        """
        differences = np.abs(time_points - target_time)
        return int(np.argmin(differences))
