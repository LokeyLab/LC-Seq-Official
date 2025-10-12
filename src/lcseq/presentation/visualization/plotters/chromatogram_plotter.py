"""
Chromatogram plotter for visualizing time-series chromatogram data.

Implementation based on THEORY.md Section 5.0-5.3.

Note: Baseline correction removed per THEORY.md - raw signals perform better.
"""

from pathlib import Path
from typing import List, Optional, Tuple, Dict
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.axes import Axes

from .base_plotter import BasePlotter
from ....domain.entities.chromatogram import Chromatogram
from ....domain.entities.peak import Peak, PeakType


class ChromatogramPlotter(BasePlotter):
    """
    Plot chromatograms with peaks and signal variants.

    Visualizes:
    - Time-series chromatogram data
    - Detected peaks with markers and labels
    - Multiple signal variants (raw, derivatives)
    - Peak classifications (NULL, TRUNCATION, PUTATIVE_PRODUCT, UNKNOWN)
    """

    # Color scheme for peak types
    PEAK_COLORS = {
        PeakType.NULL: "gray",
        PeakType.TRUNCATION: "orange",
        PeakType.PUTATIVE_PRODUCT: "green",
        PeakType.UNKNOWN: "red",
    }

    def plot(
        self,
        chromatogram: Chromatogram,
        peaks: Optional[List[Peak]] = None,
        signal_variant: str = "raw",
        title: Optional[str] = None,
    ) -> Figure:
        """
        Plot chromatogram with optional peaks.

        Parameters
        ----------
        chromatogram : Chromatogram
            Chromatogram to plot
        peaks : List[Peak], optional
            Detected peaks to overlay
        signal_variant : str, optional
            Signal variant to plot (default: "raw")
        title : str, optional
            Plot title

        Returns
        -------
        Figure
            Matplotlib figure
        """
        fig, ax = self.create_figure()

        # Get signal to plot
        signal = chromatogram.get_signal(signal_variant)
        time = chromatogram.time_points

        # Plot main signal
        ax.plot(time, signal, "b-", linewidth=1.5, label=f"{signal_variant.capitalize()} signal")

        # Plot peaks
        if peaks:
            self._plot_peaks(ax, peaks, signal, time)

        # Style
        plot_title = title or "Chromatogram"
        self.apply_common_styling(
            ax,
            title=plot_title,
            xlabel="Time (min)",
            ylabel="Intensity (counts)",
            grid=True,
        )

        ax.legend(loc="best", framealpha=0.9)

        return fig

    def _plot_peaks(self, ax: Axes, peaks: List[Peak], signal: np.ndarray, time: np.ndarray) -> None:
        """Plot detected peaks as markers."""
        plotted_types = set()

        for peak in peaks:
            # Find closest time index to peak position
            idx = np.argmin(np.abs(time - peak.position))
            height = signal[idx]

            # Get color for peak type
            color = self.PEAK_COLORS.get(peak.peak_type, "gray")

            # Add label only for first peak of each type
            label = None
            if peak.peak_type not in plotted_types:
                label = f"{peak.peak_type.value} peak"
                plotted_types.add(peak.peak_type)

            # Plot peak marker
            ax.plot(
                peak.position,
                height,
                "o",
                color=color,
                markersize=8,
                markeredgecolor="black",
                markeredgewidth=0.5,
                label=label,
            )

            # Plot peak boundaries
            ax.axvline(peak.left_base, color=color, linestyle=":", alpha=0.3, linewidth=1)
            ax.axvline(peak.right_base, color=color, linestyle=":", alpha=0.3, linewidth=1)

    def plot_multi_signal(
        self,
        chromatogram: Chromatogram,
        variants: List[str],
        title: Optional[str] = None,
    ) -> Figure:
        """
        Plot multiple signal variants in subplots.

        Parameters
        ----------
        chromatogram : Chromatogram
            Chromatogram with multiple signal variants
        variants : List[str]
            List of signal variant names to plot
        title : str, optional
            Overall plot title

        Returns
        -------
        Figure
            Matplotlib figure with subplots
        """
        n_variants = len(variants)
        fig, axes = plt.subplots(n_variants, 1, figsize=(self.figsize[0], self.figsize[1] * n_variants // 2), sharex=True)

        if n_variants == 1:
            axes = [axes]

        time = chromatogram.time_points

        for i, (ax, variant) in enumerate(zip(axes, variants)):
            signal = chromatogram.get_signal(variant)
            ax.plot(time, signal, "b-", linewidth=1.5)

            ax.set_ylabel(f"{variant.capitalize()}\nIntensity", fontsize=10)
            ax.grid(True, alpha=0.3)

        axes[-1].set_xlabel("Time (min)", fontsize=12)

        if title:
            fig.suptitle(title, fontsize=14, fontweight="bold")

        return fig
