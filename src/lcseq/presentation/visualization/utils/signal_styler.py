"""
Signal Styling Utility for Multi-Style Signal Plotting.

Analyzes signal for multi-style plotting based on truncation region and thresholds,
providing a shared implementation for lineage_plotter and compound_diagnostic_plotter.
"""

from dataclasses import dataclass
from typing import Dict, Optional
import numpy as np
from matplotlib.axes import Axes


@dataclass
class StyledSignalSegments:
    """
    Result of signal styling analysis.

    Attributes
    ----------
    time_interp : np.ndarray
        2x interpolated time array for smooth style transitions
    signal_interp : np.ndarray
        2x interpolated signal array
    segments : Dict[str, np.ndarray]
        Mapping from linestyle to array of indices in interpolated arrays
    truncation_boundary_idx : Optional[int]
        Index of truncation boundary in interpolated arrays (None if no boundary)
    """
    time_interp: np.ndarray
    signal_interp: np.ndarray
    segments: Dict[str, np.ndarray]
    truncation_boundary_idx: Optional[int]


class SignalStyler:
    """
    Analyzes signal for multi-style plotting based on truncation region and thresholds.

    Line styles:
    - Solid (-): After truncation, above threshold (valid peaks)
    - Dashed (--): In truncation region, above threshold
    - Dot-dash (-.): In truncation region, below threshold
    - Dotted (:): After truncation, below threshold

    This provides a shared implementation for both lineage_plotter and
    compound_diagnostic_plotter, replacing their duplicated styling logic.
    """

    def analyze(
        self,
        time: np.ndarray,
        signal: np.ndarray,
        truncation_boundary: Optional[float],
        baseline_threshold: float,
        snr_threshold: Optional[np.ndarray] = None,
    ) -> StyledSignalSegments:
        """
        Interpolate signal to 2x resolution and classify regions by linestyle.

        Parameters
        ----------
        time : np.ndarray
            Time array (same units as truncation_boundary)
        signal : np.ndarray
            Signal array
        truncation_boundary : float or None
            Time position of truncation boundary (None = no truncation region)
        baseline_threshold : float
            Global baseline threshold value
        snr_threshold : np.ndarray, optional
            Local SNR threshold array (same length as signal).
            If provided, point is "above threshold" only if it exceeds BOTH
            baseline_threshold AND snr_threshold (logical AND).

        Returns
        -------
        StyledSignalSegments
            Interpolated signal with segment classification
        """
        # 2x interpolation for smooth transitions between line styles
        # Creates intermediate points between each pair of consecutive time points
        time_interp = np.zeros(2 * len(time) - 1)
        time_interp[::2] = time  # Original points at even indices
        time_interp[1::2] = (time[:-1] + time[1:]) / 2  # Midpoints at odd indices

        # Interpolate signal to match
        signal_interp = np.interp(time_interp, time, signal)

        # Create threshold mask on interpolated data
        # Point is "above threshold" only if it exceeds ALL thresholds
        above_baseline = signal_interp >= baseline_threshold

        if snr_threshold is not None:
            # Interpolate SNR threshold array
            snr_threshold_interp = np.interp(time_interp, time, snr_threshold)
            above_snr = signal_interp >= snr_threshold_interp
            # Above threshold = above BOTH baseline AND SNR
            above_threshold = above_baseline & above_snr
        else:
            # No SNR threshold - only baseline matters
            above_threshold = above_baseline

        # Create truncation region mask on interpolated data
        if truncation_boundary is not None:
            boundary_idx = np.searchsorted(time_interp, truncation_boundary, side='right') - 1
            in_truncation = time_interp <= truncation_boundary
        else:
            boundary_idx = None
            in_truncation = np.zeros(len(time_interp), dtype=bool)

        # Classify into 4 regions based on truncation status and threshold
        segments = {}
        if truncation_boundary is not None:
            # Case 1: In truncation region, below threshold -> dot-dash ("-.")
            segments['-.'] = np.where(in_truncation & ~above_threshold)[0]

            # Case 2: In truncation region, above threshold -> dashed ("--")
            segments['--'] = np.where(in_truncation & above_threshold)[0]

            # Case 3: After truncation, below threshold -> dotted (":")
            segments[':'] = np.where(~in_truncation & ~above_threshold)[0]

            # Case 4: After truncation, above threshold -> solid ("-")
            segments['-'] = np.where(~in_truncation & above_threshold)[0]
        else:
            # No truncation region - simpler logic
            # Below threshold -> dotted (":")
            segments[':'] = np.where(~above_threshold)[0]

            # Above threshold -> solid ("-")
            segments['-'] = np.where(above_threshold)[0]

        return StyledSignalSegments(
            time_interp=time_interp,
            signal_interp=signal_interp,
            segments=segments,
            truncation_boundary_idx=boundary_idx,
        )

    def plot_styled(
        self,
        ax: Axes,
        styled: StyledSignalSegments,
        color: str = 'black',
        linewidth: float = 1.0,
        offset: float = 0.0,
        alpha: float = 0.9,
        zorder: int = 4,
        white_outline: bool = False,
    ) -> None:
        """
        Plot signal with linestyles based on classification.

        Parameters
        ----------
        ax : Axes
            Matplotlib axes to plot on
        styled : StyledSignalSegments
            Styled signal segments from analyze()
        color : str
            Line color
        linewidth : float
            Line width
        offset : float
            Vertical offset to add to signal (for stacked plots)
        alpha : float
            Line transparency
        zorder : int
            Z-order for layering
        white_outline : bool
            If True, draw white outline behind colored line
        """
        # Helper to plot continuous segments with given linestyle
        def plot_segments(indices, linestyle, is_truncation_region):
            """
            Plot continuous segments with boundary-respecting overlap.

            Parameters
            ----------
            indices : np.ndarray
                Indices of points in this segment
            linestyle : str
                Matplotlib linestyle code
            is_truncation_region : bool
                True if this segment is in truncation region (dashed/dot-dash),
                False if post-truncation (solid/dotted)
            """
            if len(indices) == 0:
                return

            # Find continuous segments (split where diff > 1)
            segment_splits = np.where(np.diff(indices) > 1)[0] + 1
            segments = np.split(indices, segment_splits)

            for seg in segments:
                if len(seg) == 0:
                    continue

                # Extend segment by one point on each side for continuity,
                # but respect truncation boundary to avoid visual artifacts
                start = max(0, seg[0] - 1) if seg[0] > 0 else seg[0]
                end = min(len(styled.time_interp), seg[-1] + 2)

                # Cap extension at truncation boundary, but allow overlap at boundary point
                # for visual continuity between segments
                if styled.truncation_boundary_idx is not None:
                    if is_truncation_region:
                        # Truncation segments: don't extend past boundary
                        end = min(end, styled.truncation_boundary_idx + 1)
                    else:
                        # Post-truncation segments: allow extension to boundary for seamless connection
                        start = max(start, styled.truncation_boundary_idx)

                # Draw white outline first (thicker, behind) if requested
                if white_outline:
                    ax.plot(
                        styled.time_interp[start:end],
                        styled.signal_interp[start:end] + offset,
                        color='white',
                        linewidth=linewidth + 1.5,  # Thicker for outline effect
                        linestyle=linestyle,
                        alpha=1.0,
                        zorder=zorder - 0.1,  # Behind the colored line
                        solid_capstyle='round',
                        solid_joinstyle='round',
                    )

                # Draw colored line on top
                ax.plot(
                    styled.time_interp[start:end],
                    styled.signal_interp[start:end] + offset,
                    color=color,
                    linewidth=linewidth,
                    linestyle=linestyle,
                    alpha=alpha,
                    zorder=zorder,
                )

        # Plot each region with appropriate linestyle
        if styled.truncation_boundary_idx is not None:
            # Case 1: In truncation region, below threshold -> dot-dash
            plot_segments(styled.segments.get('-.', np.array([])), "-.", is_truncation_region=True)

            # Case 2: In truncation region, above threshold -> dashed
            plot_segments(styled.segments.get('--', np.array([])), "--", is_truncation_region=True)

            # Case 3: After truncation, below threshold -> dotted
            plot_segments(styled.segments.get(':', np.array([])), ":", is_truncation_region=False)

            # Case 4: After truncation, above threshold -> solid
            plot_segments(styled.segments.get('-', np.array([])), "-", is_truncation_region=False)
        else:
            # No truncation region - simpler logic
            # Below threshold -> dotted
            plot_segments(styled.segments.get(':', np.array([])), ":", is_truncation_region=False)

            # Above threshold -> solid
            plot_segments(styled.segments.get('-', np.array([])), "-", is_truncation_region=False)
