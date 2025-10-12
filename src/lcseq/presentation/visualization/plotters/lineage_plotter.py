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
from ....domain.services import LineageFinderService
from ....config import (
    VisualizationConfig,
    PeakDetectionConfig,
    DEFAULT_TRUNCATION_MARGIN,
)


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
        # Secondary sort by canonical sequence to group chemically identical compounds
        level_attr = "monomer_level" if hierarchy_mode == HierarchyMode.MONOMER else "level"
        lineage_sorted = sorted(lineage, key=lambda c: (-getattr(c, level_attr), c.residue_sequence))

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

        # Calculate offsets with extra spacing between canonical sequence groups
        offsets = []
        prev_canonical = None
        cumulative_offset = 0.0

        for compound in lineage_sorted:
            canonical = compound.residue_sequence

            # Add extra space when transitioning to a new canonical sequence group
            if prev_canonical is not None and canonical != prev_canonical:
                cumulative_offset += VisualizationConfig.GROUP_SPACING_EXTRA

            offsets.append(cumulative_offset)
            cumulative_offset += VisualizationConfig.OFFSET_SPACING
            prev_canonical = canonical

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
        Assign colors to chemically identical compounds (same canonical sequence).

        All positional isomers (same canonical/residue sequence) get the same color.
        Colors are assigned using a rainbow colormap, organized by level groups so that
        compounds at the same level get similar colors from the same region of the
        color spectrum. Reference and null compounds are always black.

        Color mode controlled by VisualizationConfig.USE_COLORMAP:
        - True: Rainbow colormap by canonical sequence, organized by level
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
            Mapping from positional_sequence to hex color
        """
        # Check if colormap is enabled
        if not VisualizationConfig.USE_COLORMAP:
            # All black mode
            return {compound.positional_sequence: "black" for compound in compounds}

        # Determine level attribute based on hierarchy mode
        level_attr = "monomer_level" if hierarchy_mode == HierarchyMode.MONOMER else "level"

        # Group by level, then by canonical sequence within each level
        level_groups: Dict[int, Dict[str, List[Compound]]] = {}

        for compound in compounds:
            level = getattr(compound, level_attr)
            canonical = compound.residue_sequence

            if level not in level_groups:
                level_groups[level] = {}
            if canonical not in level_groups[level]:
                level_groups[level][canonical] = []
            level_groups[level][canonical].append(compound)

        # Sort levels (descending: maximal to minimal)
        sorted_levels = sorted(level_groups.keys(), reverse=True)

        # Generate colors using rainbow colormap
        # Divide spectrum into chunks, one per level
        n_levels = len(sorted_levels)
        cmap = plt.colormaps.get_cmap("rainbow")

        # Build color mapping
        color_map: Dict[str, str] = {}
        reference_canonical = reference.residue_sequence if reference else None

        for level_idx, level in enumerate(sorted_levels):
            canonical_seqs = sorted(level_groups[level].keys())
            n_canonical_at_level = len(canonical_seqs)

            # Calculate color range for this level
            # Each level gets an equal chunk of the rainbow spectrum
            level_start = level_idx / n_levels
            level_end = (level_idx + 1) / n_levels

            # Generate colors for this level's canonical sequences
            if n_canonical_at_level == 1:
                # Single canonical sequence - use midpoint of level's color range
                level_colors = [cmap((level_start + level_end) / 2)]
            else:
                # Multiple sequences - spread across level's color range
                level_colors = [cmap(level_start + (i / (n_canonical_at_level - 1)) * (level_end - level_start))
                               for i in range(n_canonical_at_level)]

            # Assign colors to canonical sequences at this level
            for canonical_idx, canonical in enumerate(canonical_seqs):
                is_reference_canonical = canonical == reference_canonical
                is_null = canonical == ""

                if is_reference_canonical or is_null:
                    color = "black"
                else:
                    # Convert RGBA to hex
                    rgba = level_colors[canonical_idx]
                    color = f"#{int(rgba[0]*255):02x}{int(rgba[1]*255):02x}{int(rgba[2]*255):02x}"

                # Assign same color to all positional variants
                for compound in level_groups[level][canonical]:
                    color_map[compound.positional_sequence] = color

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
        min_baseline_sds: float = 1.0,
        min_snr: float = 0.5,
        truncation_margin: float = 0.02
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

        # Get signals (using raw signal - no baseline correction per THEORY.md)
        time = compound.chromatogram.time_points
        raw_signal = compound.chromatogram.get_signal("raw")

        # Calculate global baseline threshold (for dotted/solid split)
        signal_min = np.min(raw_signal)
        signal_std = np.std(raw_signal)
        baseline_threshold = signal_min + (min_baseline_sds * signal_std)

        # NOTE: Local SNR threshold visualization is disabled because it doesn't match
        # the peak detector's threshold (detector uses SNR distribution of peak candidates,
        # not all points). Visualizing it would be misleading.
        # TODO: Fix by either computing SNR only for peak candidates or using a different approach
        local_snr_threshold_signal = np.zeros_like(raw_signal)  # Disabled for now

        # Normalize signal for display (presentation concern)
        signal_max = np.max(raw_signal)
        if signal_max > 0:
            normalized = raw_signal / signal_max
            normalized_baseline_threshold = baseline_threshold / signal_max
            normalized_snr_threshold = local_snr_threshold_signal / signal_max
        else:
            normalized = raw_signal
            normalized_baseline_threshold = 0.0
            normalized_snr_threshold = np.zeros_like(raw_signal)

        # Get color and styling
        color = color_map.get(compound.positional_sequence, "gray")
        is_reference = compound == reference
        linewidth = VisualizationConfig.LINEWIDTH_REFERENCE if is_reference else VisualizationConfig.LINEWIDTH_DEFAULT

        # Interpolate signal to double resolution for smooth boundary transitions
        # This creates intermediate points between each pair of consecutive time points
        # allowing segment boundaries to meet smoothly without visual gaps
        time_original = time  # Keep original for peak lookups
        time_interp = np.zeros(2 * len(time) - 1)
        time_interp[::2] = time  # Original points at even indices
        time_interp[1::2] = (time[:-1] + time[1:]) / 2  # Midpoints at odd indices

        # Interpolate normalized signal to match new time points
        normalized_interp = np.interp(time_interp, time, normalized)
        normalized_baseline_threshold_interp = np.interp(time_interp, time,
                                                         np.full_like(time, normalized_baseline_threshold))
        normalized_snr_threshold_interp = np.interp(time_interp, time, normalized_snr_threshold)

        # Use interpolated data for plotting traces
        time_minutes = time_interp / VisualizationConfig.SECONDS_PER_MINUTE
        normalized = normalized_interp
        normalized_baseline_threshold = normalized_baseline_threshold_interp[0]  # Scalar
        normalized_snr_threshold = normalized_snr_threshold_interp

        # Reverse z-order: maximal (index 0) on top, minimal (index n-1) on bottom
        trace_zorder = 5 + (n_compounds - index)
        # Peaks should be above their own trace but below the next trace
        peak_zorder = trace_zorder + 0.5

        # Determine truncation time for dashed line region
        max_truncation_time = None
        if hierarchy is not None:
            # All compounds have the minimal (NULL, L0) compound as a truncation product
            # Descendants' product peaks show where truncation products appear
            truncation_times = []

            # Get all descendants (truncation products) of this compound
            descendants = hierarchy.get_descendants(compound)
            if descendants:
                # Find the latest (maximum) truncation product retention time
                for desc in descendants:
                    desc_peaks = peaks_dict.get(desc, [])
                    # Get product peaks (PUTATIVE_PRODUCT for most, NULL for minimal compound)
                    product_peaks = [p for p in desc_peaks
                                   if p.peak_type in (PeakType.PUTATIVE_PRODUCT, PeakType.NULL)]
                    if product_peaks:
                        # Use the product peak position as the truncation time
                        truncation_times.extend([p.position for p in product_peaks])

                # Use the maximum truncation time (includes NULL at ~645s and all other descendants)
                if truncation_times:
                    max_truncation_time_base = max(truncation_times)
                    # Apply truncation margin to extend the dashed region
                    # This shows the buffer zone where product peaks cannot be selected
                    # Margin is absolute time in seconds
                    max_truncation_time = max_truncation_time_base + truncation_margin

        # Plot signal with appropriate linestyle based on thresholds and truncation region
        # Linestyle rules:
        #   - Below any threshold (global OR local SNR): dotted (":") or dot-dash ("-.")
        #   - Dotted (":") if NOT in truncation region
        #   - Dot-dash ("-.") if in truncation region AND below any threshold
        #   - Dashed ("--") if in truncation region AND above all thresholds
        #   - Solid ("-") if after truncation region AND above all thresholds

        # Create masks for different regions
        # Point is above ALL thresholds only if it exceeds both global baseline AND local SNR
        above_global_baseline = normalized >= normalized_baseline_threshold
        above_local_snr = normalized >= normalized_snr_threshold
        above_all_thresholds = above_global_baseline & above_local_snr

        if max_truncation_time is not None:
            max_trunc_minutes = max_truncation_time / VisualizationConfig.SECONDS_PER_MINUTE
            in_truncation_mask = time_minutes <= max_trunc_minutes
            # Find boundary index for segment extension capping
            boundary_idx = np.searchsorted(time_minutes, max_trunc_minutes, side='right') - 1
        else:
            in_truncation_mask = np.zeros(len(time_minutes), dtype=bool)
            boundary_idx = None

        # Plot segments efficiently with boundary-respecting extension for continuity
        def plot_segments(indices, linestyle, is_truncation_region):
            """Helper to plot continuous segments with boundary-respecting overlap and white outline.

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
                end = min(len(time_minutes), seg[-1] + 2)

                # Cap extension at truncation boundary to prevent overlap
                if boundary_idx is not None:
                    if is_truncation_region:
                        # Truncation segments: don't extend past boundary
                        end = min(end, boundary_idx + 1)
                    else:
                        # Post-truncation segments: don't extend before boundary
                        start = max(start, boundary_idx + 1)

                # Draw white outline first (thicker, behind)
                ax.plot(
                    time_minutes[start:end],
                    normalized[start:end] + offset,
                    color='white',
                    linewidth=linewidth + 1.5,  # Thicker for outline effect
                    linestyle=linestyle,
                    alpha=1.0,
                    zorder=trace_zorder - 0.1,  # Behind the colored line
                    solid_capstyle='round',
                    solid_joinstyle='round',
                )

                # Draw colored line on top
                ax.plot(
                    time_minutes[start:end],
                    normalized[start:end] + offset,
                    color=color,
                    linewidth=linewidth,
                    linestyle=linestyle,
                    alpha=VisualizationConfig.ALPHA_TRACE,
                    zorder=trace_zorder,
                )

        # Plot each region with appropriate linestyle
        if max_truncation_time is not None:
            # Case 1: In truncation region, below any threshold → dot-dash ("-.")
            trunc_below_idx = np.where(in_truncation_mask & ~above_all_thresholds)[0]
            plot_segments(trunc_below_idx, "-.", is_truncation_region=True)

            # Case 2: In truncation region, above all thresholds → dashed ("--")
            trunc_above_idx = np.where(in_truncation_mask & above_all_thresholds)[0]
            plot_segments(trunc_above_idx, "--", is_truncation_region=True)

            # Case 3: After truncation, below any threshold → dotted (":")
            post_trunc_below_idx = np.where(~in_truncation_mask & ~above_all_thresholds)[0]
            plot_segments(post_trunc_below_idx, ":", is_truncation_region=False)

            # Case 4: After truncation, above all thresholds → solid ("-")
            post_trunc_above_idx = np.where(~in_truncation_mask & above_all_thresholds)[0]
            plot_segments(post_trunc_above_idx, "-", is_truncation_region=False)
        else:
            # No truncation region: simpler logic
            # Below any threshold → dotted (":")
            below_idx = np.where(~above_all_thresholds)[0]
            plot_segments(below_idx, ":", is_truncation_region=False)

            # Above all thresholds → solid ("-")
            above_idx = np.where(above_all_thresholds)[0]
            plot_segments(above_idx, "-", is_truncation_region=False)

        # Plot peaks with different markers based on peak type
        peaks = peaks_dict.get(compound, [])
        for peak in peaks:
            # Find peak position in interpolated data
            # Peaks were detected on original time points, so find closest interpolated point
            idx = np.argmin(np.abs(time_interp - peak.position))
            peak_time = peak.position / VisualizationConfig.SECONDS_PER_MINUTE

            # Get marker shape and size based on peak type
            marker, marker_size = self._get_peak_marker(peak.peak_type)

            ax.plot(
                peak_time,
                normalized[idx] + offset,
                marker,
                color=color,
                markersize=marker_size,
                markeredgecolor="white",
                markeredgewidth=0.5,
                zorder=peak_zorder,
            )

        # Add validation indicator on left side
        # Green check if putative product found, red X if not
        has_product = (compound.selected_peak is not None and
                      compound.selected_peak.peak_type == PeakType.PUTATIVE_PRODUCT)

        indicator = "✓" if has_product else "✗"
        indicator_color = "green" if has_product else "red"

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

    def _get_peak_marker(self, peak_type: PeakType) -> tuple[str, float]:
        """
        Get matplotlib marker shape and size for peak type.

        Parameters
        ----------
        peak_type : PeakType
            Peak classification type

        Returns
        -------
        tuple[str, float]
            (marker_code, marker_size)

        Notes
        -----
        Marker mapping:
        - NULL: 's' (square) - DNA tag only (L₀)
        - TRUNCATION: 'D' (diamond) - incomplete synthesis
        - PUTATIVE_PRODUCT: 'o' (large circle) - expected product
        - UNKNOWN: '^' (triangle) - unclassified
        """
        marker_map = {
            PeakType.NULL: ('s', VisualizationConfig.MARKER_SIZE),           # Square
            PeakType.TRUNCATION: ('D', VisualizationConfig.MARKER_SIZE),     # Diamond
            PeakType.PUTATIVE_PRODUCT: ('o', VisualizationConfig.MARKER_SIZE * 1.5),  # Large circle
            PeakType.UNKNOWN: ('^', VisualizationConfig.MARKER_SIZE),        # Triangle
        }
        return marker_map.get(peak_type, ('o', VisualizationConfig.MARKER_SIZE))  # Default to circle

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

        # Perform MSA-style alignment
        aligned_seq = self._align_to_reference(
            compound, alignment_reference, hierarchy_mode
        )

        label = f"{aligned_seq}  (L{level})"
        return label

    def _align_to_reference(
        self,
        compound: Compound,
        reference: Compound,
        hierarchy_mode: HierarchyMode
    ) -> str:
        """
        Align compound sequence to reference using MSA-style gaps.

        Building-block mode: Position-based alignment (trivial)
        Monomer mode: Subsequence alignment to reference

        Parameters
        ----------
        compound : Compound
            Compound to align
        reference : Compound
            Reference compound (alignment template)
        hierarchy_mode : HierarchyMode
            Alignment mode

        Returns
        -------
        str
            Aligned sequence with "----" gaps for nulls/missing monomers
        """
        if hierarchy_mode == HierarchyMode.BUILDING_BLOCK:
            return self._align_building_blocks(compound, reference)
        else:  # MONOMER
            return self._align_monomers(compound, reference)

    def _align_building_blocks(self, compound: Compound, reference: Compound) -> str:
        """
        Align building blocks by synthesis position.

        Each position in compound aligns to same position in reference.
        AgxNull becomes gap (dashes matching reference residue length).
        No spaces in output - positions joined with "-".

        Parameters
        ----------
        compound : Compound
            Compound to align
        reference : Compound
            Reference compound (defines gap lengths)

        Returns
        -------
        str
            Position-aligned sequence with gaps (e.g., "---DNvl-DPhe")

        Examples
        --------
        >>> # Reference: Phe-DNvl-DPhe
        >>> # Compound:  AgxNull-DNvl-DPhe
        >>> aligned = _align_building_blocks(compound, reference)
        >>> aligned
        '---DNvl-DPhe'  # "---" matches "Phe" length (3 chars)
        """
        # Split sequences into building blocks
        compound_blocks = compound.positional_sequence.split("-")
        reference_blocks = reference.positional_sequence.split("-")

        # Align each position
        aligned = []
        for cpd_bb, ref_bb in zip(compound_blocks, reference_blocks):
            if "AgxNull" in cpd_bb or cpd_bb == "Null":
                # Replace null with dashes matching reference residue length
                aligned.append("-" * len(ref_bb))
            else:
                # Keep non-null residue as-is
                aligned.append(cpd_bb)

        return "-".join(aligned)

    def _align_monomers(self, compound: Compound, reference: Compound) -> str:
        """
        Align monomers using subsequence mapping to reference.

        Uses LineageFinderService to get position mappings via greedy subsequence
        matching. Inserts gaps (dashes) where monomers are missing.

        Parameters
        ----------
        compound : Compound
            Compound to align (subsequence of reference)
        reference : Compound
            Reference compound (full sequence template)

        Returns
        -------
        str
            Subsequence-aligned sequence with gaps (no spaces, joined by "-")

        Examples
        --------
        >>> # Reference: Leu-LA03-Pro-Leu-DPro (5 monomers, positions 0-4)
        >>> # Compound:  Leu-Pro-DPro (3 monomers)
        >>> # Mapping:   [0, 2, 4] (matches at ref positions 0, 2, 4)
        >>> aligned = _align_monomers(compound, reference)
        >>> aligned
        'Leu-----Pro---DPro'  # Gaps at pos 1 (LA03=5 chars) and 3 (Leu=3 chars)
        """
        # Use domain service to get subsequence alignment mapping
        lineage_finder = LineageFinderService()
        mapping = lineage_finder.get_monomer_alignment_mapping(compound, reference)

        # Get reference monomers for gap sizing
        reference_monomers = []
        for bb in reversed(reference.building_blocks):
            reference_monomers.extend(bb.decompose_to_monomers())

        # Get candidate monomers
        candidate_monomers = []
        for bb in reversed(compound.building_blocks):
            candidate_monomers.extend(bb.decompose_to_monomers())

        # Handle all-null case
        if not mapping:
            # All gaps (all positions in reference)
            return "-".join("-" * len(m) for m in reference_monomers)

        # Build aligned sequence using mapping
        aligned = []
        cand_idx = 0

        for ref_idx, ref_monomer in enumerate(reference_monomers):
            if cand_idx < len(mapping) and mapping[cand_idx] == ref_idx:
                # Candidate has monomer at this reference position
                aligned.append(candidate_monomers[cand_idx])
                cand_idx += 1
            else:
                # Gap: candidate missing monomer at this position
                # Use dashes matching reference monomer length
                aligned.append("-" * len(ref_monomer))

        return "-".join(aligned)

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
        reference_label = reference.residue_sequence if reference else lineage_sorted[-1].residue_sequence
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
        # Create custom legend handles for peak markers
        marker_handles = [
            Line2D([0], [0], marker='s', color='black', linestyle='None',
                   markersize=6, markeredgecolor='white', markeredgewidth=0.5,
                   label='NULL (L₀)'),
            Line2D([0], [0], marker='D', color='black', linestyle='None',
                   markersize=6, markeredgecolor='white', markeredgewidth=0.5,
                   label='Truncation'),
            Line2D([0], [0], marker='o', color='black', linestyle='None',
                   markersize=9, markeredgecolor='white', markeredgewidth=0.5,
                   label='Product'),
            Line2D([0], [0], marker='^', color='black', linestyle='None',
                   markersize=6, markeredgecolor='white', markeredgewidth=0.5,
                   label='Unknown'),
        ]

        # Create custom legend handles for line styles
        line_handles = [
            Line2D([0], [0], color='black', linestyle='-', linewidth=2,
                   label='Valid region'),
            Line2D([0], [0], color='black', linestyle='--', linewidth=2,
                   label='Before latest truncation'),
            Line2D([0], [0], color='black', linestyle='-.', linewidth=2,
                   label='Before latest truncation and below threshold'),
            Line2D([0], [0], color='black', linestyle=':', linewidth=2,
                   label='Below thresholds'),
        ]

        # Combine all handles and labels (markers first row, line styles second row)
        all_handles = marker_handles + line_handles
        all_labels = [h.get_label() for h in all_handles]

        # Add legend at bottom center with 2-row layout (no frame)
        ax.legend(
            all_handles,
            all_labels,
            loc='upper center',
            bbox_to_anchor=(0.5, -0.04),  # Centered horizontally, closer to x-axis
            ncol=4,  # 4 columns: wraps to 2 rows (markers above, line styles below)
            frameon=False,  # No box around legend
            fontsize=9,
            handlelength=2,
            handleheight=1,
            columnspacing=1.5,
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
