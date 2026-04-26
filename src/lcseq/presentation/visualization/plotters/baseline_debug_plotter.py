"""
Debug visualization for baseline estimation.

Shows original signal, fitted baseline(s), and corrected signal
on the same axes with distinct colors.
"""

from pathlib import Path
from typing import Optional

import numpy as np
from matplotlib.figure import Figure
from numpy.typing import NDArray

from lcseq.domain.services.baseline_estimator import (
    BaselineCurveResult,
    BaselineCurveParams,
    BaselineCurveMethod,
    fit_baseline_curve,
    compare_baseline_methods,
    PiecewiseALSParams,
    PiecewiseALSResult,
    piecewise_als_baseline,
    MinimaSplineParams,
    MinimaSplineResult,
    minima_spline_baseline,
)
from lcseq.presentation.visualization.plotters.base_plotter import BasePlotter


class BaselineDebugPlotter(BasePlotter):
    """
    Debug plotter for baseline curve fitting.

    Creates a figure showing:
    - Original signal
    - Fitted baseline curve(s)
    - Corrected signal (optional)

    All curves on the same axes with distinct colors.
    """

    # Color scheme
    COLORS = {
        "signal": "#2c3e50",  # Dark blue-gray
        "baseline": "#e74c3c",  # Red
        "corrected": "#27ae60",  # Green
    }

    def __init__(
        self,
        figsize: tuple[float, float] = (14, 8),
        dpi: int = 150,
        style: str = "seaborn-v0_8-whitegrid",
    ):
        super().__init__(figsize=figsize, dpi=dpi, style=style)

    def plot(
        self,
        signal: NDArray[np.float64],
        result: BaselineCurveResult,
        time_axis: Optional[NDArray[np.float64]] = None,
        title: str = "Baseline Debug",
        show_corrected: bool = True,
    ) -> Figure:
        """
        Create debug visualization of baseline fitting.

        Parameters
        ----------
        signal : NDArray[np.float64]
            Original signal
        result : BaselineCurveResult
            Result from fit_baseline_curve()
        time_axis : NDArray[np.float64], optional
            Time/retention axis. If None, uses point indices.
        title : str
            Plot title
        show_corrected : bool
            Show corrected signal

        Returns
        -------
        Figure
            Matplotlib figure
        """
        if time_axis is None:
            time_axis = np.arange(len(signal))

        fig, ax = self.create_figure()

        # Plot original signal
        ax.plot(
            time_axis,
            signal,
            color=self.COLORS["signal"],
            linewidth=1.0,
            alpha=0.8,
            label="Original signal",
            zorder=1,
        )

        # Plot baseline
        ax.plot(
            time_axis,
            result.baseline,
            color=self.COLORS["baseline"],
            linewidth=2.0,
            alpha=0.9,
            label=f"Baseline ({result.method})",
            zorder=3,
        )

        # Plot corrected signal
        if show_corrected:
            ax.plot(
                time_axis,
                result.corrected,
                color=self.COLORS["corrected"],
                linewidth=1.0,
                alpha=0.7,
                label="Corrected signal",
                zorder=2,
            )

        # Styling
        self.apply_common_styling(
            ax,
            title=title,
            xlabel="Time" if time_axis is not None else "Point index",
            ylabel="Intensity",
            grid=True,
        )

        ax.legend(loc="upper right", fontsize=10, framealpha=0.9)

        return fig

    def plot_comparison(
        self,
        signal: NDArray[np.float64],
        results: dict[str, BaselineCurveResult],
        time_axis: Optional[NDArray[np.float64]] = None,
        title: str = "Baseline Method Comparison",
    ) -> Figure:
        """
        Compare multiple baseline fitting results.

        Parameters
        ----------
        signal : NDArray[np.float64]
            Original signal
        results : dict[str, BaselineCurveResult]
            Dictionary mapping method name to result
        time_axis : NDArray[np.float64], optional
            Time axis
        title : str
            Plot title

        Returns
        -------
        Figure
            Matplotlib figure with comparison
        """
        if time_axis is None:
            time_axis = np.arange(len(signal))

        fig, ax = self.create_figure()

        # Plot original signal
        ax.plot(
            time_axis,
            signal,
            color=self.COLORS["signal"],
            linewidth=1.5,
            alpha=0.8,
            label="Original signal",
        )

        # Plot each baseline with different color
        import matplotlib.cm as cm

        cmap = cm.get_cmap("Set1")
        for i, (name, result) in enumerate(results.items()):
            color = cmap(i / max(len(results) - 1, 1))
            ax.plot(
                time_axis,
                result.baseline,
                color=color,
                linewidth=2.0,
                alpha=0.8,
                label=f"{name}",
            )

        self.apply_common_styling(
            ax,
            title=title,
            xlabel="Time",
            ylabel="Intensity",
            grid=True,
        )

        ax.legend(loc="upper right", fontsize=10)

        return fig


def plot_baseline_debug(
    signal: NDArray[np.float64],
    params: Optional[BaselineCurveParams] = None,
    time_axis: Optional[NDArray[np.float64]] = None,
    output_path: Optional[Path] = None,
    title: str = "Baseline Debug",
    show_corrected: bool = True,
) -> tuple[Figure, BaselineCurveResult]:
    """
    Convenience function to fit baseline and create debug plot.

    Parameters
    ----------
    signal : NDArray[np.float64]
        Input signal
    params : BaselineCurveParams, optional
        Baseline fitting parameters
    time_axis : NDArray[np.float64], optional
        Time axis for x-axis
    output_path : Path, optional
        If provided, save figure to this path
    title : str
        Plot title
    show_corrected : bool
        Show corrected signal

    Returns
    -------
    tuple[Figure, BaselineCurveResult]
        Figure and baseline fitting result
    """
    # Fit baseline
    result = fit_baseline_curve(signal, time_axis, params)

    # Create plot
    plotter = BaselineDebugPlotter()
    fig = plotter.plot(
        signal=signal,
        result=result,
        time_axis=time_axis,
        title=title,
        show_corrected=show_corrected,
    )

    # Save if path provided
    if output_path:
        plotter.save(fig, output_path)

    return fig, result


def plot_baseline_comparison(
    signal: NDArray[np.float64],
    time_axis: Optional[NDArray[np.float64]] = None,
    methods: Optional[list[BaselineCurveMethod]] = None,
    lam: float = 1e5,
    output_path: Optional[Path] = None,
    title: str = "Baseline Method Comparison",
) -> tuple[Figure, dict[str, BaselineCurveResult]]:
    """
    Compare multiple baseline methods and create comparison plot.

    Parameters
    ----------
    signal : NDArray[np.float64]
        Input signal
    time_axis : NDArray[np.float64], optional
        Time axis
    methods : list[BaselineCurveMethod], optional
        Methods to compare. Defaults to [AIRPLS, ARPLS, SNIP].
    lam : float
        Smoothness parameter for ALS-based methods.
    output_path : Path, optional
        If provided, save figure to this path
    title : str
        Plot title

    Returns
    -------
    tuple[Figure, dict[str, BaselineCurveResult]]
        Figure and results dictionary
    """
    # Compare methods
    results = compare_baseline_methods(signal, time_axis, methods, lam)

    # Create comparison plot
    plotter = BaselineDebugPlotter()
    fig = plotter.plot_comparison(
        signal=signal,
        results=results,
        time_axis=time_axis,
        title=title,
    )

    # Save if path provided
    if output_path:
        plotter.save(fig, output_path)

    return fig, results


def plot_piecewise_baseline(
    signal: NDArray[np.float64],
    params: Optional[PiecewiseALSParams] = None,
    time_axis: Optional[NDArray[np.float64]] = None,
    output_path: Optional[Path] = None,
    title: str = "Piecewise ALS Baseline (Local Minima)",
    show_corrected: bool = True,
    show_splits: bool = True,
    show_segments: bool = True,
) -> tuple[Figure, PiecewiseALSResult]:
    """
    Fit piecewise ALS baseline and create debug plot.

    Parameters
    ----------
    signal : NDArray[np.float64]
        Input signal
    params : PiecewiseALSParams, optional
        Piecewise ALS parameters
    time_axis : NDArray[np.float64], optional
        Time axis for x-axis
    output_path : Path, optional
        If provided, save figure to this path
    title : str
        Plot title
    show_corrected : bool
        Show corrected signal
    show_splits : bool
        Show vertical lines at split points
    show_segments : bool
        Show individual segment baselines in different colors

    Returns
    -------
    tuple[Figure, PiecewiseALSResult]
        Figure and baseline fitting result
    """
    import matplotlib.cm as cm

    # Fit baseline
    result = piecewise_als_baseline(signal, time_axis, params)

    if time_axis is None:
        time_axis = np.arange(len(signal))

    # Create plot
    plotter = BaselineDebugPlotter()
    fig, ax = plotter.create_figure()

    # Plot original signal
    ax.plot(
        time_axis,
        signal,
        color=plotter.COLORS["signal"],
        linewidth=1.0,
        alpha=0.8,
        label="Original signal",
        zorder=1,
    )

    # Plot segment baselines in different colors
    if show_segments and len(result.segment_baselines) > 1:
        cmap = cm.get_cmap("tab10")
        n_segments = len(result.segment_baselines)
        for i, (start, end, seg_baseline) in enumerate(result.segment_baselines):
            color = cmap(i / max(n_segments - 1, 1))
            seg_time = time_axis[start:end]
            label = f"Segment {i + 1}" if i < 3 else None
            ax.plot(
                seg_time,
                seg_baseline,
                color=color,
                linewidth=2.0,
                alpha=0.6,
                linestyle="--",
                label=label,
                zorder=2,
            )

    # Plot final baseline
    ax.plot(
        time_axis,
        result.baseline,
        color=plotter.COLORS["baseline"],
        linewidth=2.5,
        alpha=0.9,
        label="Baseline (piecewise)",
        zorder=3,
    )

    # Plot corrected signal
    if show_corrected:
        ax.plot(
            time_axis,
            result.corrected,
            color=plotter.COLORS["corrected"],
            linewidth=1.0,
            alpha=0.7,
            label="Corrected signal",
            zorder=2,
        )

    # Plot split points
    if show_splits and result.split_points:
        for i, split_idx in enumerate(result.split_points):
            split_time = time_axis[split_idx]
            label = "Split points (local minima)" if i == 0 else None
            ax.axvline(
                split_time,
                color="#95a5a6",
                linewidth=1.5,
                linestyle=":",
                alpha=0.7,
                label=label,
                zorder=0,
            )
            # Mark the actual minimum point
            ax.scatter(
                [split_time],
                [signal[split_idx]],
                color="#e74c3c",
                s=50,
                marker="v",
                zorder=4,
            )

    # Styling
    plotter.apply_common_styling(
        ax,
        title=title,
        xlabel="Time",
        ylabel="Intensity",
        grid=True,
    )

    ax.legend(loc="upper right", fontsize=9, framealpha=0.9)

    # Add segment count to title
    n_segs = len(result.segment_baselines)
    ax.set_title(f"{title}\n({n_segs} segments, {len(result.split_points)} split points)")

    # Save if path provided
    if output_path:
        plotter.save(fig, output_path)

    return fig, result


def plot_minima_spline_baseline(
    signal: NDArray[np.float64],
    params: Optional[MinimaSplineParams] = None,
    time_axis: Optional[NDArray[np.float64]] = None,
    output_path: Optional[Path] = None,
    title: str = "Minima Spline Baseline",
    show_corrected: bool = True,
    show_minima: bool = True,
) -> tuple[Figure, MinimaSplineResult]:
    """
    Fit baseline via smoothed spline through local minima and create debug plot.

    Shows three traces:
    - Original signal (green)
    - Baseline (red)
    - Corrected signal (black)

    Parameters
    ----------
    signal : NDArray[np.float64]
        Input signal
    params : MinimaSplineParams, optional
        Spline fitting parameters
    time_axis : NDArray[np.float64], optional
        Time axis for x-axis
    output_path : Path, optional
        If provided, save figure to this path
    title : str
        Plot title
    show_corrected : bool
        Show corrected signal
    show_minima : bool
        Show detected minima as scatter points

    Returns
    -------
    tuple[Figure, MinimaSplineResult]
        Figure and baseline fitting result
    """
    # Fit baseline
    result = minima_spline_baseline(signal, time_axis, params)

    # Color scheme
    COLORS = {
        "original": "#27ae60",    # Green
        "baseline": "#e74c3c",    # Red
        "corrected": "#2c3e50",   # Black/dark
    }

    # Create plot
    plotter = BaselineDebugPlotter()
    fig, ax = plotter.create_figure()

    # Plot original signal (green)
    ax.plot(
        result.x,
        result.original,
        color=COLORS["original"],
        linewidth=1.0,
        alpha=0.7,
        label="Original",
        zorder=1,
        marker='o' if len(result.original) < 200 else None,
        markersize=3,
    )

    # Plot baseline (red)
    ax.plot(
        result.x,
        result.baseline,
        color=COLORS["baseline"],
        linewidth=2.5,
        alpha=0.9,
        label="Baseline",
        zorder=2,
    )

    # Plot corrected signal (black)
    if show_corrected:
        ax.plot(
            result.x,
            result.corrected,
            color=COLORS["corrected"],
            linewidth=1.5,
            alpha=0.9,
            label="Corrected",
            zorder=3,
        )

    # Plot detected minima on original signal
    if show_minima and len(result.minima_x) > 0:
        ax.scatter(
            result.minima_x,
            result.minima_y,
            color=COLORS["baseline"],
            s=60,
            marker="v",
            label=f"Minima ({len(result.minima_x)})",
            zorder=4,
            edgecolors="white",
            linewidths=1,
        )

    # Styling
    plotter.apply_common_styling(
        ax,
        title=title,
        xlabel="Time",
        ylabel="Intensity",
        grid=True,
    )

    ax.legend(loc="upper right", fontsize=10, framealpha=0.9)

    # Save if path provided
    if output_path:
        plotter.save(fig, output_path)

    return fig, result
