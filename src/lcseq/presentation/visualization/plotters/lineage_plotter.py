"""
Lineage offset chromatogram plotter.

Visualizes lineage hierarchies as vertically offset chromatograms with peaks.
All visualization logic (color assignment, normalization, plotting) contained here.

Terminology (THEORY.md Section 3.1):
    - Reference compound: The compound currently being analyzed (not "parent")
    - Lineage: All ancestors + descendants + self (not "family")
"""

from pathlib import Path
from typing import List, Dict, Optional
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.lines import Line2D

from .base_plotter import BasePlotter
from ....domain.entities import Compound, Peak
from ....domain.entities.peak import PeakType
from ....domain.models import HierarchyMode, CompoundHierarchy
from ....domain.services import SequenceAligner
from ....domain.services.baseline_estimator import BaselineEstimatorService
from ....config import (
    VisualizationConfig,
    PeakDetectionConfig,
    PeakAppearanceConfig,
    DEFAULT_TRUNCATION_MARGIN,
)
from ..utils.signal_styler import SignalStyler


class LineageOffsetPlotter(BasePlotter):
    """
    Plot lineage chromatograms as vertically offset traces.

    Handles all visualization concerns:
    - Color assignment by canonical sequence
    - Signal normalization for display
    - Label formatting and truncation
    - Grid layout and styling
    - Peak markers

    Terminology (THEORY.md Section 3.1):
    - Uses "reference" instead of "parent"
    - Uses "lineage" instead of "family"

    THEORY.md Compliance:
    - Uses CompoundOrderingService for similarity-based ordering (Section 8.6)
    - No business logic - only presentation concerns
    """

    def __init__(self, **kwargs):
        """Initialize plotter with custom defaults for lineage plots."""
        # Override default figsize - will be recalculated per plot
        # All parameters come from VisualizationConfig
        super().__init__(
            figsize=(VisualizationConfig.FIG_WIDTH, VisualizationConfig.FIG_HEIGHT_MIN),
            dpi=150,
            **kwargs
        )
        # Reuse SequenceAligner across all plots (single source of truth)
        self._aligner = SequenceAligner()

    def plot(
        self,
        lineage: List[Compound],
        peaks_dict: Dict[Compound, List[Peak]],
        output_path: Optional[Path] = None,
        reference: Optional[Compound] = None,
        hierarchy_mode: HierarchyMode = HierarchyMode.MONOMER,
        hierarchy: Optional[CompoundHierarchy] = None,
        min_baseline_sds: float = PeakDetectionConfig.MIN_BASELINE_SDS,
        min_snr: float = PeakDetectionConfig.MIN_SNR,
        truncation_margin: float = DEFAULT_TRUNCATION_MARGIN,
    ) -> Figure:
        """
        Plot lineage as offset chromatograms.

        Parameters
        ----------
        lineage : List[Compound]
            List of compounds in lineage (reference + descendants)
        peaks_dict : Dict[Compound, List[Peak]]
            Detected peaks for each compound
        output_path : Path, optional
            Path to save plot (if None, returns figure without saving)
        reference : Compound, optional
            Reference compound (for highlighting) - THEORY.md Section 3.1:
            "The compound currently being analyzed"
        hierarchy_mode : HierarchyMode, optional
            Hierarchy mode for level display (default: MONOMER per THEORY.md)
        hierarchy : CompoundHierarchy, optional
            Hierarchy for determining truncation regions (enables dashed lines)
        min_baseline_sds : float, optional
            Global baseline threshold in standard deviations (default: 1.0)
            Signal below min(signal) + min_baseline_sds * std(signal) shown as dotted
        min_snr : float, optional
            Local SNR threshold in standard deviations (default: 0.5)
            Signal below local SNR threshold shown as dotted/dot-dash
        truncation_margin : float, optional
            Margin beyond truncation positions (in seconds).
            Dashed region extended by this margin to show buffer zone

        Returns
        -------
        Figure
            Matplotlib figure object

        References
        ----------
        THEORY.md Section 3.1: Terminology - Ancestry and Lineage
        THEORY.md Section 3.3: Truncation Level
        """
        # Sort by level (maximal → minimal, THEORY.md Section 3.3)
        # Use appropriate level attribute based on mode
        # Secondary sort by block support sequence to group chemically identical compounds
        level_attr = "monomer_level" if hierarchy_mode == HierarchyMode.MONOMER else "level"
        lineage_sorted = sorted(lineage, key=lambda c: (-getattr(c, level_attr), c.block_support_sequence))

        # Assign colors (presentation concern)
        # Colors grouped by level for visual clarity
        color_map = self._assign_lineage_colors(lineage_sorted, reference, hierarchy_mode)

        # Create figure with dynamic height
        n_compounds = len(lineage_sorted)
        fig_height = max(
            VisualizationConfig.FIG_HEIGHT_MIN,
            VisualizationConfig.FIG_HEIGHT_BASE + n_compounds * VisualizationConfig.FIG_HEIGHT_PER_TRACE
        )
        fig, ax = plt.subplots(figsize=(VisualizationConfig.FIG_WIDTH, fig_height))

        # Use reference (maximal) as alignment template (MSA-style, THEORY.md Section 2.3.4)
        alignment_reference = reference if reference else lineage_sorted[0]

        # Calculate common label position based on reference compound
        ref_time = alignment_reference.chromatogram.time_points
        label_x_right = (ref_time[-1] / VisualizationConfig.SECONDS_PER_MINUTE) + 1.0

        # Calculate offsets with extra spacing between block support sequence groups
        offsets = []
        prev_block_support = None
        cumulative_offset = 0.0

        for compound in lineage_sorted:
            block_support = compound.block_support_sequence

            # Add extra space when transitioning to a new block support sequence group
            if prev_block_support is not None and block_support != prev_block_support:
                cumulative_offset += VisualizationConfig.GROUP_SPACING_EXTRA

            offsets.append(cumulative_offset)
            cumulative_offset += VisualizationConfig.OFFSET_SPACING
            prev_block_support = block_support

        # Plot each compound
        # Z-order: maximal (bottom of plot, index 0) should be on top
        # Reverse z-order so later compounds (higher index) are drawn underneath
        for i, compound in enumerate(lineage_sorted):
            self._plot_compound_trace(
                ax, compound, offsets[i], i, peaks_dict, color_map, reference, hierarchy_mode,
                alignment_reference, n_compounds, label_x_right, hierarchy, min_baseline_sds, min_snr,
                truncation_margin
            )

        # Format axes and title
        # Pass max offset for proper y-axis scaling
        max_offset = offsets[-1] if offsets else 0.0
        self._format_plot(fig, ax, lineage_sorted, reference, n_compounds, max_offset)

        plt.tight_layout()

        # Save or return
        if output_path:
            self.save(fig, output_path)
            print(f"✓ Plot saved: {output_path.name}")

        return fig

    def _assign_lineage_colors(
        self,
        compounds: List[Compound],
        reference: Optional[Compound] = None,
        hierarchy_mode: HierarchyMode = HierarchyMode.MONOMER
    ) -> Dict[str, str]:
        """
        Assign colors to chemically identical compounds (same block support sequence).

        All positional isomers (same block support sequence) get the same color.
        Colors are assigned using a rainbow colormap, organized by level groups so that
        compounds at the same level get similar colors from the same region of the
        color spectrum. Reference and null compounds are always black.

        Color mode controlled by VisualizationConfig.USE_COLORMAP:
        - True: Rainbow colormap by block support sequence, organized by level
        - False: All black (publication-ready)

        Parameters
        ----------
        compounds : List[Compound]
            Compounds to assign colors to
        reference : Compound, optional
            Reference compound (assigned black) - THEORY.md Section 3.1
        hierarchy_mode : HierarchyMode, optional
            Hierarchy mode for level grouping (default: MONOMER)

        Returns
        -------
        Dict[str, str]
            Mapping from positional_block_sequence to hex color
        """
        # Check if colormap is enabled
        if not VisualizationConfig.USE_COLORMAP:
            # All black mode
            return {compound.positional_block_sequence: "black" for compound in compounds}

        # Determine level attribute based on hierarchy mode
        level_attr = "monomer_level" if hierarchy_mode == HierarchyMode.MONOMER else "level"

        # Group by level, then by block support sequence within each level
        level_groups: Dict[int, Dict[str, List[Compound]]] = {}

        for compound in compounds:
            level = getattr(compound, level_attr)
            block_support = compound.block_support_sequence

            if level not in level_groups:
                level_groups[level] = {}
            if block_support not in level_groups[level]:
                level_groups[level][block_support] = []
            level_groups[level][block_support].append(compound)

        # Sort levels (descending: maximal to minimal)
        sorted_levels = sorted(level_groups.keys(), reverse=True)

        # Generate colors using rainbow colormap
        # Divide spectrum into chunks, one per level
        n_levels = len(sorted_levels)
        cmap = plt.colormaps.get_cmap("rainbow")

        # Build color mapping
        color_map: Dict[str, str] = {}
        reference_block_support = reference.block_support_sequence if reference else None

        for level_idx, level in enumerate(sorted_levels):
            block_support_seqs = sorted(level_groups[level].keys())
            n_block_support_at_level = len(block_support_seqs)

            # Calculate color range for this level
            # Each level gets an equal chunk of the rainbow spectrum
            level_start = level_idx / n_levels
            level_end = (level_idx + 1) / n_levels

            # Generate colors for this level's block support sequences
            if n_block_support_at_level == 1:
                # Single block support sequence - use midpoint of level's color range
                level_colors = [cmap((level_start + level_end) / 2)]
            else:
                # Multiple sequences - spread across level's color range
                level_colors = [cmap(level_start + (i / (n_block_support_at_level - 1)) * (level_end - level_start))
                               for i in range(n_block_support_at_level)]

            # Assign colors to block support sequences at this level
            for block_support_idx, block_support in enumerate(block_support_seqs):
                is_reference_block_support = block_support == reference_block_support
                is_null = block_support == ""

                if is_reference_block_support or is_null:
                    color = "black"
                else:
                    # Convert RGBA to hex
                    rgba = level_colors[block_support_idx]
                    color = f"#{int(rgba[0]*255):02x}{int(rgba[1]*255):02x}{int(rgba[2]*255):02x}"

                # Assign same color to all positional variants
                for compound in level_groups[level][block_support]:
                    color_map[compound.positional_block_sequence] = color

        return color_map

    def _plot_compound_trace(
        self,
        ax,
        compound: Compound,
        offset: float,
        index: int,
        peaks_dict: Dict[Compound, List[Peak]],
        color_map: Dict[str, str],
        reference: Optional[Compound],
        hierarchy_mode: HierarchyMode,
        alignment_reference: Compound,
        n_compounds: int,
        label_x_right: float,
        hierarchy: Optional[CompoundHierarchy] = None,
        min_baseline_sds: float = 0.0,
        min_snr: float = 0.5,
        truncation_margin: float = 0.02,
    ) -> None:
        """
        Plot a single compound trace with peaks and labels.

        Parameters
        ----------
        ax : Axes
            Matplotlib axes to plot on
        compound : Compound
            Compound to plot
        offset : float
            Vertical offset for this trace (includes group spacing)
        index : int
            Index in sorted list (for z-order calculation)
        peaks_dict : Dict[Compound, List[Peak]]
            Peak dictionary
        color_map : Dict[str, str]
            Color mapping
        reference : Compound, optional
            Reference compound for highlighting - THEORY.md Section 3.1
        hierarchy_mode : HierarchyMode
            Hierarchy mode for level display
        alignment_reference : Compound
            Reference compound to align sequences to (MSA-style)
        n_compounds : int
            Total number of compounds (for z-order calculation)
        label_x_right : float
            X-position for sequence labels in minutes
        hierarchy : CompoundHierarchy, optional
            Hierarchy for determining truncation regions (enables dashed lines)
        min_baseline_sds : float, optional
            Global baseline threshold in standard deviations (default: 1.0)
            Signal below threshold shown as dotted
        min_snr : float, optional
            Local SNR threshold in standard deviations (default: 0.5)
            Signal below local SNR threshold shown as dotted/dot-dash
        truncation_margin : float, optional
            Margin beyond truncation positions (in seconds).
            Extends dashed region to show buffer zone where products cannot be selected
        """

        # Get baseline-corrected signal
        time = compound.chromatogram.time_points
        signal = compound.chromatogram.get_signal("corrected")

        # Calculate global baseline threshold (for dotted/solid split)
        # Uses same sigma-clipping method as peak detector for consistency
        baseline_estimator = BaselineEstimatorService()
        background, noise_std = baseline_estimator.estimate_with_noise(signal)
        baseline_threshold = background + (min_baseline_sds * noise_std)

        # NOTE: Local SNR threshold visualization is disabled because it doesn't match
        # the peak detector's threshold (detector uses SNR distribution of peak candidates,
        # not all points). Visualizing it would be misleading.
        # TODO: Fix by either computing SNR only for peak candidates or using a different approach
        local_snr_threshold_signal = np.zeros_like(signal)  # Disabled for now

        # Normalize signal for display (presentation concern)
        signal_max = np.max(signal)
        if signal_max > 0:
            normalized = signal / signal_max
            normalized_baseline_threshold = baseline_threshold / signal_max
            normalized_snr_threshold = local_snr_threshold_signal / signal_max
        else:
            normalized = signal
            normalized_baseline_threshold = 0.0
            normalized_snr_threshold = np.zeros_like(signal)

        # Get color and styling
        color = color_map.get(compound.positional_block_sequence, "gray")
        is_reference = compound == reference
        linewidth = VisualizationConfig.LINEWIDTH_REFERENCE if is_reference else VisualizationConfig.LINEWIDTH_DEFAULT

        # Reverse z-order: maximal (index 0) on top, minimal (index n-1) on bottom
        trace_zorder = 5 + (n_compounds - index)
        # Peaks should be above their own trace but below the next trace
        peak_zorder = trace_zorder + 0.5

        # Determine truncation time for dashed line region
        # Uses selected_peak.position (same as diagnostic plotter) for consistency
        max_truncation_time = None
        if hierarchy is not None:
            truncation_times = []

            # Get all descendants (truncation products) of this compound
            descendants = hierarchy.get_descendants(compound)
            if descendants:
                # Find the latest (maximum) truncation product retention time
                # Use selected_peak (same as diagnostic plotter)
                for desc in descendants:
                    if desc.selected_peak is not None:
                        truncation_times.append(desc.selected_peak.position)

                # Use the maximum truncation time
                if truncation_times:
                    max_truncation_time_base = max(truncation_times)
                    # Apply truncation margin to extend the dashed region
                    max_truncation_time = max_truncation_time_base + truncation_margin

        # Use SignalStyler for consistent multi-style signal plotting
        # Supports both global baseline AND local SNR thresholds
        # Convert time to minutes (x-axis is in minutes)
        time_minutes_arr = time / VisualizationConfig.SECONDS_PER_MINUTE
        max_truncation_time_minutes = (
            max_truncation_time / VisualizationConfig.SECONDS_PER_MINUTE
            if max_truncation_time is not None else None
        )
        styler = SignalStyler()
        styled = styler.analyze(
            time=time_minutes_arr,
            signal=normalized,
            truncation_boundary=max_truncation_time_minutes,
            baseline_threshold=normalized_baseline_threshold,
            snr_threshold=normalized_snr_threshold,
        )

        # Plot with color and white outline
        styler.plot_styled(
            ax=ax,
            styled=styled,
            color=color,
            linewidth=linewidth,
            offset=offset,
            alpha=VisualizationConfig.ALPHA_TRACE,
            zorder=trace_zorder,
            white_outline=True,
        )

        # Keep interpolated data for peak plotting (now in minutes)
        time_interp = styled.time_interp
        signal_interp = styled.signal_interp
        time_minutes = time / VisualizationConfig.SECONDS_PER_MINUTE

        # Plot peaks with different markers based on peak type
        # Matches diagnostic plotter exactly: hollow gray for rejected, colored for accepted
        peaks = peaks_dict.get(compound, [])

        for peak in peaks:
            # Find peak position in interpolated data (both in minutes)
            peak_time = peak.position / VisualizationConfig.SECONDS_PER_MINUTE
            idx = np.argmin(np.abs(time_interp - peak_time))

            # Get marker shape and size based on peak type
            marker, marker_size = PeakAppearanceConfig.get_marker(peak.peak_type)
            peak_color = PeakAppearanceConfig.get_color(peak.peak_type)

            # Use hollow markers for rejected peaks (same as diagnostic plotter)
            if peak.is_rejected:
                ax.plot(
                    peak_time,
                    signal_interp[idx] + offset,
                    marker,
                    color='gray',
                    markersize=marker_size * 0.8,
                    markerfacecolor='none',
                    markeredgecolor='gray',
                    markeredgewidth=1.5,
                    zorder=peak_zorder - 0.1,
                )
            else:
                ax.plot(
                    peak_time,
                    signal_interp[idx] + offset,
                    marker,
                    color=peak_color,
                    markersize=marker_size,
                    markeredgecolor="white",
                    markeredgewidth=1.0,
                    zorder=peak_zorder,
                )

        # Add validation indicator on left side
        # Green check if expected peak found:
        #   - For L0 (null compound): NULL peak is expected
        #   - For others: PUTATIVE_PRODUCT peak is expected
        is_null_compound = compound.is_null_compound
        has_expected_peak = False

        if compound.selected_peak is not None:
            if is_null_compound:
                # L0 compound: NULL peak is the expected signal (DNA tag only)
                has_expected_peak = compound.selected_peak.peak_type == PeakType.NULL
            else:
                # Non-null compound: PUTATIVE_PRODUCT is expected
                has_expected_peak = compound.selected_peak.peak_type == PeakType.PUTATIVE_PRODUCT

        indicator = "✓" if has_expected_peak else "✗"
        indicator_color = "green" if has_expected_peak else "red"

        # Position indicator to the left of the signal start
        indicator_x = time_minutes[0] - 0.5

        ax.text(
            indicator_x,
            offset,
            indicator,
            verticalalignment="center",
            horizontalalignment="center",
            fontsize=VisualizationConfig.LABEL_FONTSIZE + 2,
            color=indicator_color,
            fontweight="bold",
            alpha=1.0,
            zorder=trace_zorder + 1,
        )

        # Add label (MSA-style alignment with monospace font)
        label = self._format_compound_label(compound, is_reference, hierarchy_mode, alignment_reference)
        fontweight = "bold" if is_reference else "normal"

        ax.text(
            label_x_right,  # Consistent x-position for all sequence labels
            offset,
            label,
            verticalalignment="center",
            horizontalalignment="left",  # Left-align text (padding already applied)
            fontsize=VisualizationConfig.LABEL_FONTSIZE,
            fontfamily="monospace",  # Monospaced font for proper character alignment
            color="black",  # All text black
            fontweight=fontweight,
            alpha=1.0,
        )

        # Horizontal gridline
        ax.axhline(
            y=offset,
            color="lightgrey",
            linestyle="--",
            linewidth=0.8,
            alpha=0.5,
            zorder=0,
        )

    def _format_compound_label(
        self,
        compound: Compound,
        is_reference: bool,
        hierarchy_mode: HierarchyMode,
        alignment_reference: Compound
    ) -> str:
        """
        Format compound label with MSA-style alignment to reference.

        Aligns sequences like Multiple Sequence Alignment, showing gaps (----)
        where AgxNull appears. All sequences align to the reference/maximal compound.

        Parameters
        ----------
        compound : Compound
            Compound to format label for
        is_reference : bool
            Whether this is the reference compound
        hierarchy_mode : HierarchyMode
            Hierarchy mode for level display (THEORY.md Section 3.3)
        alignment_reference : Compound
            Reference compound to align to (defines column positions)

        Returns
        -------
        str
            Formatted label with MSA-style alignment

        References
        ----------
        THEORY.md Section 2.3.4: Sequence alignment for visualization
        THEORY.md Section 3.3: Truncation Level
        """
        # Get level for annotation
        if hierarchy_mode == HierarchyMode.MONOMER:
            level = compound.monomer_level
        else:  # BUILDING_BLOCK
            level = compound.level

        # Perform MSA-style alignment using SequenceAligner (single source of truth)
        aligned_seq = self._aligner.align_to_reference(
            compound, alignment_reference, hierarchy_mode
        )

        label = f"{aligned_seq}  (L{level})"
        return label

    def _format_plot(
        self,
        fig: Figure,
        ax,
        lineage_sorted: List[Compound],
        reference: Optional[Compound],
        n_compounds: int,
        max_offset: float
    ) -> None:
        """
        Format plot axes, labels, and title.

        Parameters
        ----------
        fig : Figure
            Matplotlib figure
        ax : Axes
            Matplotlib axes to format
        lineage_sorted : List[Compound]
            Sorted lineage list
        reference : Compound, optional
            Reference compound - THEORY.md Section 3.1
        n_compounds : int
            Number of compounds
        max_offset : float
            Maximum vertical offset (accounts for group spacing)
        """
        # X-axis label (black text)
        ax.set_xlabel("Time (minutes)", fontsize=12, fontweight="bold", color="black")

        # Generate title (black text)
        reference_label = reference.block_support_sequence if reference else lineage_sorted[-1].block_support_sequence
        ax.set_title(
            f"Lineage Analysis - {reference_label} (n={n_compounds})",
            fontsize=14,
            fontweight="bold",
            color="black"
        )

        # Configure axes
        ax.grid(False)  # Disable grid (we have horizontal lines at each compound)
        ax.set_ylim(-0.2, max_offset + 1.5)
        ax.set_yticks([])

        # Set x-axis limits to end exactly at first and last ticks
        if lineage_sorted:
            # First, set data limits
            time_points = lineage_sorted[0].chromatogram.time_points
            time_min = time_points[0] / VisualizationConfig.SECONDS_PER_MINUTE
            time_max = time_points[-1] / VisualizationConfig.SECONDS_PER_MINUTE
            ax.set_xlim(time_min, time_max)

            # Now let matplotlib generate ticks based on data range
            fig.canvas.draw()
            ticks = ax.get_xticks()

            # Set xlim to span exactly from first to last tick
            if len(ticks) > 0:
                ax.set_xlim(ticks[0], ticks[-1])

        # Remove left, top, and right spines
        ax.spines['left'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        # Make bottom spine black
        ax.spines['bottom'].set_color('black')

        # Make tick labels black
        ax.tick_params(axis='x', colors='black')

        # Add legend at bottom center
        self._add_legend(ax)

    def _add_legend(self, ax) -> None:
        """
        Add legend showing peak markers and line styles at bottom center.

        Creates a horizontal legend positioned below the x-axis, centered on the plot.
        Shows both peak type markers and signal region line styles.
        """
        # Create legend handles using PeakAppearanceConfig
        marker_handles = PeakAppearanceConfig.create_legend_handles(include_rejected=True)
        line_handles = PeakAppearanceConfig.create_linestyle_handles(include_truncation=True)

        all_handles = marker_handles + line_handles

        # Legend at bottom center, below x-axis label
        ax.legend(
            handles=all_handles,
            loc='upper center',
            bbox_to_anchor=(0.5, -0.04),
            ncol=5,
            frameon=False,
            fontsize=8,
            handlelength=1.5,
            handleheight=1,
            columnspacing=1.0,
        )

    def _compute_local_snr_threshold(
        self,
        signal: np.ndarray,
        min_snr: float,
        window_size: int = 5
    ) -> np.ndarray:
        """
        Compute local SNR threshold signal for visualization.

        Computes a rolling local SNR estimate across the signal, then determines
        the adaptive threshold signal based on local noise characteristics.

        Parameters
        ----------
        signal : np.ndarray
            Raw signal array
        min_snr : float
            SNR threshold parameter (standard deviations above minimum SNR)
        window_size : int, optional
            Size of local window for noise estimation (default: 5)

        Returns
        -------
        np.ndarray
            Threshold signal - points below this are considered below SNR threshold

        Notes
        -----
        Algorithm:
        1. For each point, compute local noise from surrounding windows
        2. For each point, compute local baseline (median of surrounding)
        3. Compute local SNR for each point
        4. Compute adaptive threshold: min(SNRs) + min_snr * std(SNRs)
        5. Convert SNR threshold back to signal threshold using local characteristics

        This follows the peak detector logic but applies it across the entire signal
        for visualization purposes.
        """
        n = len(signal)
        background = np.percentile(signal, 10)

        # Compute local characteristics for each point
        local_snrs = np.zeros(n)

        for i in range(n):
            # Define noise windows (regions AROUND point, not including it)
            left_start = max(0, i - window_size)
            left_end = max(0, i)
            right_start = min(n, i + 1)
            right_end = min(n, i + 1 + window_size)

            # Collect surrounding points
            left_noise = signal[left_start:left_end] if left_end > left_start else np.array([])
            right_noise = signal[right_start:right_end] if right_end > right_start else np.array([])
            surrounding = np.concatenate([left_noise, right_noise])

            if len(surrounding) < 3:
                # Not enough data - use global background
                local_baseline = background
                local_noise = np.sqrt(background + 1.0)
            else:
                local_baseline = np.median(surrounding)
                local_noise = np.std(surrounding)
                if local_noise < 1e-6:
                    local_noise = np.sqrt(background + 1.0)

            # Compute local SNR for this point
            signal_amplitude = signal[i] - local_baseline
            local_snrs[i] = signal_amplitude / local_noise

        # Compute adaptive SNR threshold from all local SNRs
        snr_min = np.min(local_snrs)
        snr_std = np.std(local_snrs)
        snr_threshold = snr_min + (min_snr * snr_std)

        # Convert SNR threshold back to signal threshold for each point
        threshold_signal = np.zeros(n)

        for i in range(n):
            # Recompute local characteristics
            left_start = max(0, i - window_size)
            left_end = max(0, i)
            right_start = min(n, i + 1)
            right_end = min(n, i + 1 + window_size)

            left_noise = signal[left_start:left_end] if left_end > left_start else np.array([])
            right_noise = signal[right_start:right_end] if right_end > right_start else np.array([])
            surrounding = np.concatenate([left_noise, right_noise])

            if len(surrounding) < 3:
                local_baseline = background
                local_noise = np.sqrt(background + 1.0)
            else:
                local_baseline = np.median(surrounding)
                local_noise = np.std(surrounding)
                if local_noise < 1e-6:
                    local_noise = np.sqrt(background + 1.0)

            # Threshold signal at this point: baseline + (SNR_threshold * noise)
            threshold_signal[i] = local_baseline + (snr_threshold * local_noise)

        return threshold_signal
