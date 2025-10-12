"""
Validation plotter for dashboard-style validation result visualization.

Implementation based on THEORY.md Section 6.10-6.11.
"""

from typing import List, Optional, Dict
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.figure import Figure

from .base_plotter import BasePlotter
from ....domain.models.analysis_result import AnalysisResult
from ....domain.entities.compound import Compound
from ....domain.entities.peak import ValidationStatus


class ValidationPlotter(BasePlotter):
    """
    Create dashboard visualizations for validation results.

    Shows:
    - Validation status distribution (pie chart)
    - Purity distribution (histogram)
    - SNR distribution (histogram)
    - Top/bottom performers (bar charts)
    """

    def plot(
        self,
        analysis_result: AnalysisResult,
        title: Optional[str] = None,
    ) -> Figure:
        """
        Create validation dashboard with multiple subplots.

        Parameters
        ----------
        analysis_result : AnalysisResult
            Analysis results with compounds and validation data
        title : str, optional
            Overall title

        Returns
        -------
        Figure
            Matplotlib figure with dashboard
        """
        fig = plt.figure(figsize=(14, 10))

        # Create grid layout
        gs = fig.add_gridspec(3, 2, hspace=0.3, wspace=0.3)
        ax1 = fig.add_subplot(gs[0, 0])  # Status distribution
        ax2 = fig.add_subplot(gs[0, 1])  # SNR distribution
        ax3 = fig.add_subplot(gs[1, 0])  # Purity distribution
        ax4 = fig.add_subplot(gs[1, 1])  # Peak count distribution
        ax5 = fig.add_subplot(gs[2, :])  # Top performers

        # Plot each component
        self._plot_status_distribution(ax1, analysis_result)
        self._plot_snr_distribution(ax2, analysis_result)
        self._plot_purity_distribution(ax3, analysis_result)
        self._plot_peak_count_distribution(ax4, analysis_result)
        self._plot_top_performers(ax5, analysis_result)

        # Overall title
        if title:
            fig.suptitle(title, fontsize=16, fontweight="bold")

        return fig

    def _plot_status_distribution(self, ax, analysis_result: AnalysisResult) -> None:
        """Plot pie chart of validation status."""
        if not analysis_result.validation_results:
            ax.text(0.5, 0.5, "No data", ha="center", va="center")
            return

        # Get status summary
        summary = analysis_result.get_validation_summary()

        if not summary:
            ax.text(0.5, 0.5, "No data", ha="center", va="center")
            return

        # Map ValidationStatus to colors
        colors = {
            ValidationStatus.VALIDATED: "#90EE90",
            ValidationStatus.LIKELY_SUCCESS: "#98FB98",
            ValidationStatus.UNCERTAIN: "#FFD700",
            ValidationStatus.LIKELY_FAILURE: "#FFA07A",
            ValidationStatus.FAILED: "#FFB6C1",
        }

        labels = [status.value for status in summary.keys()]
        values = list(summary.values())
        pie_colors = [colors.get(status, "gray") for status in summary.keys()]

        ax.pie(values, labels=labels, autopct="%1.1f%%", colors=pie_colors, startangle=90)
        ax.set_title("Validation Status Distribution", fontsize=12, fontweight="bold")

    def _plot_snr_distribution(self, ax, analysis_result: AnalysisResult) -> None:
        """Plot histogram of SNR values (placeholder)."""
        # TODO: Implement when SNR metrics are added to domain model
        ax.text(0.5, 0.5, "SNR metrics\nnot yet implemented", ha="center", va="center", fontsize=14)
        ax.set_title("SNR Distribution", fontsize=12, fontweight="bold")
        ax.axis("off")

    def _plot_purity_distribution(self, ax, analysis_result: AnalysisResult) -> None:
        """Plot histogram of purity values (placeholder)."""
        # TODO: Implement when purity metrics are added to domain model
        ax.text(0.5, 0.5, "Purity metrics\nnot yet implemented", ha="center", va="center", fontsize=14)
        ax.set_title("Purity Distribution", fontsize=12, fontweight="bold")
        ax.axis("off")

    def _plot_peak_count_distribution(self, ax, analysis_result: AnalysisResult) -> None:
        """Plot histogram of peak counts."""
        # Get peak counts per compound
        peak_counts = [len(peaks) for peaks in analysis_result.peak_classifications.values()]

        if not peak_counts:
            ax.text(0.5, 0.5, "No peak data", ha="center", va="center")
            return

        max_peaks = max(peak_counts)
        ax.hist(peak_counts, bins=range(max_peaks + 2), color="orange", edgecolor="black", alpha=0.7)
        ax.set_xlabel("Number of Detected Peaks")
        ax.set_ylabel("Count")
        ax.set_title("Peak Count Distribution", fontsize=12, fontweight="bold")
        ax.grid(True, alpha=0.3)

    def _plot_top_performers(self, ax, analysis_result: AnalysisResult, top_n: int = 10) -> None:
        """Plot bar chart of validated compounds."""
        # Get validated compounds
        validated = analysis_result.get_validated_compounds()

        if not validated:
            ax.text(0.5, 0.5, "No validated compounds", ha="center", va="center")
            return

        # Sort by compound ID and take top N
        validated_sorted = sorted(validated, key=lambda c: c.compound_id or "")[:top_n]

        ids = [c.compound_id[:15] if c.compound_id else "Unknown" for c in validated_sorted]
        # Use peak count as a simple metric
        peak_counts = [
            len(analysis_result.peak_classifications.get(c, [])) for c in validated_sorted
        ]

        ax.barh(ids, peak_counts, color="seagreen", edgecolor="black")
        ax.set_xlabel("Number of Peaks")
        ax.set_title(f"Top {top_n} Validated Compounds", fontsize=12, fontweight="bold")
        ax.grid(True, alpha=0.3, axis="x")
        ax.invert_yaxis()  # First at top
