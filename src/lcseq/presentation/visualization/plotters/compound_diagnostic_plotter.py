"""
Single compound diagnostic plotter for detailed peak visualization.

Creates diagnostic plots for individual compounds showing:
- Raw signal (actual counts) vs time (minutes)
- All detected peaks with type-specific markers
- Peak boundaries as shaded regions
- % area labels and origin info for truncation peaks
- Truncation region shading and boundary line

Used for validating peak detection and classification across many compounds.
"""

from pathlib import Path
from typing import List, Optional, Tuple
import numpy as np

# Use Agg backend for thread-safe plotting (required for parallel generation on macOS)
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.lines import Line2D

# Optional: adjustText for automatic annotation repositioning
try:
    from adjustText import adjust_text
    HAS_ADJUST_TEXT = True
except ImportError:
    HAS_ADJUST_TEXT = False

from .base_plotter import BasePlotter
from ....domain.entities import Compound, Peak, RejectionReason
from ....domain.entities.peak import PeakType
from ....domain.models import CompoundHierarchy, HierarchyMode
from ....config import VisualizationConfig, PeakAppearanceConfig, DEFAULT_TRUNCATION_MARGIN, PeakDetectionConfig
from ....domain.services import SequenceAligner
from ....domain.services.baseline_estimator import BaselineEstimatorService
from ..utils.signal_styler import SignalStyler


class CompoundDiagnosticPlotter(BasePlotter):
    """
    Single-compound diagnostic plot with full peak annotations.

    Shows:
    - Raw signal (actual counts) vs time (minutes)
    - All detected peaks with type markers
    - Peak boundaries as shaded regions
    - % area labels near each peak
    - Origin info for truncation peaks
    - Truncation region shading + boundary line

    Example usage::

        plotter = CompoundDiagnosticPlotter()
        fig = plotter.plot(compound, hierarchy=hierarchy)
        plotter.save(fig, Path("output/compound_diagnostic.png"))

    For batch generation::

        generate_diagnostic_plots(compounds, hierarchy, Path("output/diagnostics/"))
    """

    def __init__(self, **kwargs):
        """Initialize plotter with default figure size for diagnostic plots."""
        super().__init__(
            figsize=(14, 6),
            dpi=100,  # Reduced from 150 for faster rendering
            **kwargs
        )
        # Reuse SequenceAligner across all plots (single source of truth)
        self._aligner = SequenceAligner()

    def plot(
        self,
        compound: Compound,
        hierarchy: Optional[CompoundHierarchy] = None,
        truncation_margin: float = DEFAULT_TRUNCATION_MARGIN,
        output_path: Optional[Path] = None,
        show_height: bool = True,
    ) -> Figure:
        """
        Generate diagnostic plot for a single compound.

        Parameters
        ----------
        compound : Compound
            Compound to visualize
        hierarchy : CompoundHierarchy, optional
            Hierarchy for computing truncation boundary from descendants.
            If None, truncation region is not shown.
        truncation_margin : float, optional
            Margin beyond max truncation position (in seconds).
            Default from config.
        output_path : Path, optional
            Path to save the plot. If None, figure is returned without saving.
        show_height : bool, optional
            Whether to show peak height values in annotations. Default True.

        Returns
        -------
        Figure
            Matplotlib figure object
        """
        # Build compound alias mapping for descendants
        alias_map = {}  # sequence -> alias letter
        descendants_ordered = []  # (compound, alias) tuples, minimal→maximal

        if hierarchy is not None:
            descendants = hierarchy.get_descendants(compound)
            if descendants:
                # Sort minimal at top (ascending level), then by sequence
                descendants_sorted = sorted(
                    descendants,
                    key=lambda c: (c.level, c.block_support_sequence)
                )
                descendants_ordered = descendants_sorted

                # Assign aliases A, B, C, ...
                for i, desc in enumerate(descendants_sorted):
                    alias = chr(ord('A') + i) if i < 26 else f"Z{i-25}"
                    alias_map[desc.block_support_sequence] = alias

        # Create figure - wider to accommodate compound key on right
        fig, ax = plt.subplots(figsize=(16, 6))

        # Get signal data - use baseline-corrected signal
        time_seconds = compound.chromatogram.time_points
        time_minutes = time_seconds / VisualizationConfig.SECONDS_PER_MINUTE
        signal = compound.chromatogram.get_signal("corrected")

        # Compute truncation boundary using cached descendants (used for line styling)
        truncation_boundary = None
        if hierarchy is not None and descendants_ordered:
            truncation_boundary = self._compute_truncation_boundary_from_descendants(
                descendants_ordered, truncation_margin
            )

        # Plot peak boundaries (shaded regions behind signal) - only for accepted peaks
        peaks = compound.detected_peaks or []
        for peak in peaks:
            # Skip boundary regions for rejected peaks - they aren't "true peaks"
            if peak.is_rejected:
                continue
            color = PeakAppearanceConfig.get_color(peak.peak_type)
            left_min = peak.left_base / VisualizationConfig.SECONDS_PER_MINUTE
            right_min = peak.right_base / VisualizationConfig.SECONDS_PER_MINUTE
            ax.axvspan(left_min, right_min, alpha=0.2, color=color, zorder=3)

        # Plot main signal with line style variations
        self._plot_signal_with_styles(
            ax, time_minutes, signal, truncation_boundary
        )

        # Calculate total area for percentage calculation
        total_area = self._calculate_total_area(peaks, signal, time_seconds)

        # Plot peak markers and collect annotation data
        annotation_data = []  # List of (text, (x, y), color)

        for peak in peaks:
            self._plot_peak_marker(ax, peak, signal, time_seconds, time_minutes)
            result = self._create_peak_annotation(
                ax, peak, signal, time_seconds,
                total_area, show_height, alias_map
            )
            if result[0] is not None:
                annotation_data.append(result)

        # Create annotations with calculated positions to avoid overlaps
        if annotation_data:
            self._place_annotations(ax, annotation_data, time_minutes, signal)

        # Format plot
        self._format_plot(ax, compound, truncation_boundary)

        # Add compound key on the right
        if descendants_ordered:
            self._add_compound_key(ax, compound, descendants_ordered, alias_map)

        # Add legend at bottom
        self._add_legend(ax, peaks, has_truncation_region=(truncation_boundary is not None))

        plt.tight_layout(rect=[0, 0.12, 1, 1])  # Leave room at bottom for legend

        # Save if path provided
        if output_path:
            self.save(fig, output_path)

        return fig

    def _plot_signal_with_styles(
        self,
        ax,
        time_minutes: np.ndarray,
        signal: np.ndarray,
        truncation_boundary: Optional[float],
        min_baseline_sds: float = PeakDetectionConfig.MIN_BASELINE_SDS,
    ) -> None:
        """
        Plot signal with different line styles based on region and threshold.

        Uses SignalStyler utility for consistent multi-style signal plotting.

        Line styles:
        - Solid (-): Valid region (after truncation, above threshold)
        - Dashed (--): In truncation region, above threshold
        - Dot-dash (-.): In truncation region, below threshold
        - Dotted (:): After truncation, below threshold

        Parameters
        ----------
        ax : Axes
            Matplotlib axes
        time_minutes : np.ndarray
            Time array in minutes
        signal : np.ndarray
            Signal array (raw counts)
        truncation_boundary : float or None
            Truncation boundary in seconds (None = no truncation region)
        min_baseline_sds : float
            Baseline threshold in standard deviations
        """
        # Calculate baseline threshold using same method as peak detector
        # Uses sigma-clipping to exclude peaks, then: baseline + min_baseline_sds * noise_std
        baseline_estimator = BaselineEstimatorService()
        background, noise_std = baseline_estimator.estimate_with_noise(signal)
        baseline_threshold = background + (min_baseline_sds * noise_std)

        # Convert truncation boundary from seconds to minutes (same units as time_minutes)
        if truncation_boundary is not None:
            trunc_minutes = truncation_boundary / VisualizationConfig.SECONDS_PER_MINUTE
        else:
            trunc_minutes = None

        # Use SignalStyler for analysis and plotting
        styler = SignalStyler()
        styled = styler.analyze(
            time=time_minutes,
            signal=signal,
            truncation_boundary=trunc_minutes,
            baseline_threshold=baseline_threshold,
        )

        # Plot with black color (like stacked plot)
        styler.plot_styled(
            ax=ax,
            styled=styled,
            color='black',
            linewidth=VisualizationConfig.LINEWIDTH_DEFAULT,
            offset=0.0,
            alpha=VisualizationConfig.ALPHA_TRACE,
            zorder=4,
            white_outline=False,
        )

    def _compute_truncation_boundary_from_descendants(
        self,
        descendants: List[Compound],
        margin: float
    ) -> Optional[float]:
        """
        Calculate truncation boundary from pre-fetched descendants.

        The boundary is max(descendant product positions) + margin.
        All peaks before this boundary could be truncation products.

        Parameters
        ----------
        descendants : List[Compound]
            Pre-fetched descendants (avoids redundant hierarchy traversal)
        margin : float
            Margin to add beyond max truncation position (seconds)

        Returns
        -------
        float or None
            Truncation boundary in seconds, or None if no descendants with peaks
        """
        if not descendants:
            return None

        truncation_times = []
        for desc in descendants:
            if desc.selected_peak is not None:
                truncation_times.append(desc.selected_peak.position)

        if not truncation_times:
            return None

        return max(truncation_times) + margin

    def _calculate_total_area(
        self,
        peaks: List[Peak],
        signal: np.ndarray,
        time_seconds: np.ndarray
    ) -> float:
        """Calculate total area under accepted peaks for percentage calculation."""
        total = 0.0
        for peak in peaks:
            # Only include accepted peaks in total area (matches purity calculation)
            if not peak.is_accepted:
                continue
            area = self._get_peak_area(peak, signal, time_seconds)
            total += area
        return total if total > 0 else 1.0  # Avoid division by zero

    def _get_peak_area(
        self,
        peak: Peak,
        signal: np.ndarray,
        time_seconds: np.ndarray
    ) -> float:
        """
        Get area for a peak, either from stored value or by integration.

        Uses trapezoidal integration within peak boundaries.
        """
        # Try to use stored area if available
        if hasattr(peak, 'area') and peak.area is not None and peak.area > 0:
            return peak.area

        # Otherwise integrate within boundaries
        mask = (time_seconds >= peak.left_base) & (time_seconds <= peak.right_base)
        if not np.any(mask):
            return 0.0

        return np.trapz(signal[mask], time_seconds[mask])

    def _plot_peak_marker(
        self,
        ax,
        peak: Peak,
        signal: np.ndarray,
        time_seconds: np.ndarray,
        time_minutes: np.ndarray
    ) -> None:
        """Plot peak marker at peak position."""
        # Find closest index to peak position
        idx = np.argmin(np.abs(time_seconds - peak.position))
        peak_time_min = peak.position / VisualizationConfig.SECONDS_PER_MINUTE
        height = signal[idx]

        marker, size = PeakAppearanceConfig.get_marker(peak.peak_type)
        color = PeakAppearanceConfig.get_color(peak.peak_type)

        # Use hollow markers for rejected peaks
        if peak.is_rejected:
            ax.plot(
                peak_time_min, height,
                marker,
                color='gray',
                markersize=size * 0.8,
                markerfacecolor='none',
                markeredgecolor='gray',
                markeredgewidth=1.5,
                zorder=9  # Below accepted peaks
            )
        else:
            ax.plot(
                peak_time_min, height,
                marker,
                color=color,
                markersize=size,
                markeredgecolor='white',
                markeredgewidth=1.0,
                zorder=10
            )

    def _create_peak_annotation(
        self,
        ax,
        peak: Peak,
        signal: np.ndarray,
        time_seconds: np.ndarray,
        total_area: float,
        show_height: bool,
        alias_map: Optional[dict] = None
    ) -> Tuple[Optional[str], Tuple[float, float], str]:
        """
        Create text annotation for a peak (without final positioning).

        Parameters
        ----------
        ax : Axes
            Matplotlib axes
        peak : Peak
            Peak to annotate
        signal : np.ndarray
            Signal array
        time_seconds : np.ndarray
            Time array in seconds
        total_area : float
            Total area for percentage calculation
        show_height : bool
            Whether to include height in annotation
        alias_map : dict, optional
            Mapping from compound sequence to alias letter (e.g., "Leu-Pro" -> "A")

        Returns
        -------
        Tuple[Optional[plt.Text], Optional[Tuple[float, float]]]
            (text_object, (x, y) position) or (None, None) if no annotation
        """
        # Find peak height
        idx = np.argmin(np.abs(time_seconds - peak.position))
        height = signal[idx]
        peak_time_min = peak.position / VisualizationConfig.SECONDS_PER_MINUTE

        # Calculate percentage
        peak_area = self._get_peak_area(peak, signal, time_seconds)
        pct = 100.0 * peak_area / total_area

        # Build compact annotation lines
        lines = []

        # Peak type label - unified notation
        # Format: (compound)suffix where:
        #   compound: * = self, A/B/C = descendant alias
        #   suffix: (none) = accepted product, ? = accepted unknown, - = rejected (significance), etc.
        type_label = self._build_peak_label(peak, alias_map)
        lines.append(type_label)

        # Retention time - compact (just value + m)
        rt_min = peak.position / VisualizationConfig.SECONDS_PER_MINUTE
        lines.append(f"{rt_min:.2f}m")

        # Percentage - only for accepted peaks (rejected peaks excluded from purity)
        if peak.is_accepted:
            lines.append(f"{pct:.1f}%")

        # Height (optional) - compact
        if show_height:
            lines.append(f"{height:.0f}")

        # Create annotation text
        annotation_text = '\n'.join(lines)

        # Use gray for rejected peaks, otherwise use peak type color
        color = 'gray' if peak.is_rejected else PeakAppearanceConfig.get_color(peak.peak_type)

        # Return annotation data (text will be created after positions are calculated)
        return annotation_text, (peak_time_min, height), color

    def _place_annotations(
        self,
        ax,
        annotation_data: List[Tuple[str, Tuple[float, float], str]],
        time_minutes: np.ndarray,
        signal: np.ndarray,
        headroom_fraction: float = 0.35
    ) -> None:
        """
        Place annotations above the local signal without overlapping.

        Each annotation is placed above the local signal in its region,
        with tier offsets to avoid overlapping nearby annotations.

        Parameters
        ----------
        ax : Axes
            Matplotlib axes
        annotation_data : List[Tuple[str, Tuple[float, float], str]]
            List of (text, (x, y), color) tuples
        time_minutes : np.ndarray
            Time array in minutes
        signal : np.ndarray
            Signal array
        headroom_fraction : float
            Fraction of signal range to add as headroom (default 0.35)
        """
        if not annotation_data:
            return

        y_min = np.min(signal)
        y_max = np.max(signal)
        y_range = y_max - y_min

        # Add headroom to y-axis for annotations
        headroom = y_range * headroom_fraction
        ax.set_ylim(y_min - y_range * 0.05, y_max + headroom)

        # Annotation box height estimate (as fraction of y_range)
        box_height = y_range * 0.12
        # Minimum vertical gap between annotation tiers
        tier_gap = y_range * 0.08
        # Base offset above local signal
        base_offset = y_range * 0.05

        # Window size for finding local signal max (in data points)
        x_range = time_minutes[-1] - time_minutes[0]
        window_fraction = 0.05  # 5% of x-range
        window_size = max(3, int(len(time_minutes) * window_fraction))

        # Sort annotations by x position
        sorted_data = sorted(enumerate(annotation_data), key=lambda x: x[1][1][0])

        # Track placed annotations: list of (x, y_bottom, y_top) for collision detection
        placed = []
        min_x_gap = x_range * 0.06  # Minimum x gap to consider annotations as non-overlapping

        annotation_placements = []

        for orig_idx, (text, (x, y), color) in sorted_data:
            # Find local signal max in window around this x position
            x_idx = np.argmin(np.abs(time_minutes - x))
            window_start = max(0, x_idx - window_size // 2)
            window_end = min(len(signal), x_idx + window_size // 2 + 1)
            local_max = np.max(signal[window_start:window_end])

            # Start placement above local signal max
            text_y = local_max + base_offset

            # Check for collisions with already-placed annotations
            # Only check annotations that are close in x
            for px, py_bottom, py_top in placed:
                if abs(x - px) < min_x_gap:
                    # These annotations might overlap horizontally
                    # Check if our y position would overlap
                    our_bottom = text_y
                    our_top = text_y + box_height

                    if not (our_top < py_bottom or our_bottom > py_top):
                        # Collision! Move above the conflicting annotation
                        text_y = py_top + tier_gap

            # Record this placement
            placed.append((x, text_y, text_y + box_height))
            annotation_placements.append((orig_idx, text, x, y, color, text_y))

        # Create annotations with proper arrows
        for orig_idx, text, x, y, color, text_y in annotation_placements:
            ax.annotate(
                text,
                xy=(x, y),  # Arrow points to peak position
                xytext=(x, text_y),  # Text placed above local signal
                ha='center',
                va='bottom',
                fontsize=8,
                bbox=dict(
                    boxstyle='round,pad=0.3',
                    facecolor='white',
                    edgecolor=color,
                    alpha=0.9
                ),
                arrowprops=dict(
                    arrowstyle='-',
                    color='gray',
                    alpha=0.6,
                    lw=0.8,
                    connectionstyle='arc3,rad=0'
                ),
                zorder=15
            )

    def _get_alias(self, peak: Peak, alias_map: Optional[dict]) -> Optional[str]:
        """Get alias letter for a truncation peak's matched compound."""
        if not alias_map:
            return None
        if hasattr(peak, 'matched_compound_sequence') and peak.matched_compound_sequence:
            return alias_map.get(peak.matched_compound_sequence)
        return None

    def _build_peak_label(self, peak: Peak, alias_map: Optional[dict]) -> str:
        """
        Build unified peak label using new notation system.

        Notation: (compound)suffix
        - compound: * = this compound, A/B/C = descendant alias
        - suffix:
            (none) = accepted, matched to product
            ? = accepted, matched to non-product (unknown)
            - = rejected (significance)
            ~ = rejected (prominence)
            _ = rejected (baseline)
            ! = rejected (SNR)
            x = rejected (not maximum)

        Examples:
        - (*) = this compound's product (accepted)
        - (*)? = this compound's unknown peak (accepted)
        - (*)- = would be this compound's peak, rejected by significance test
        - (A) = truncation from A's product (accepted)
        - (A)? = truncation from A's unknown (accepted)
        - (A)- = would match A, rejected by significance test
        """
        # Determine the compound identifier
        alias = self._get_alias(peak, alias_map)

        # Determine if this is a truncation (has matched compound) or self
        if peak.peak_type in (PeakType.TRUNCATION, PeakType.TRUNCATION_UNKNOWN):
            compound_id = alias if alias else "?"
        elif peak.peak_type == PeakType.NULL:
            compound_id = "0"  # L0 / null
        else:
            compound_id = "*"  # Self (product or unknown)

        # Determine the suffix based on rejection status and peak type
        if peak.is_rejected:
            # Rejected peak - suffix indicates reason
            suffix = self._get_rejection_suffix(peak.rejection_reason)
        elif peak.peak_type in (PeakType.TRUNCATION_UNKNOWN, PeakType.UNKNOWN):
            # Accepted but unknown/non-product match
            suffix = "?"
        else:
            # Accepted product match (TRUNCATION, PUTATIVE_PRODUCT, NULL)
            suffix = ""

        return f"({compound_id}){suffix}"

    def _get_rejection_suffix(self, reason: RejectionReason) -> str:
        """Get suffix character for rejection reason."""
        suffix_map = {
            RejectionReason.NONE: "",
            RejectionReason.SIGNIFICANCE: "-",
            RejectionReason.PROMINENCE: "~",
            RejectionReason.BASELINE: "_",
            RejectionReason.SNR: "!",
            RejectionReason.NOT_MAXIMUM: "/",
        }
        return suffix_map.get(reason, "?")

    def _get_type_label(self, peak_type: PeakType) -> str:
        """Get short label for peak type."""
        labels = {
            PeakType.NULL: "NL",
            PeakType.TRUNCATION: "TR",
            PeakType.TRUNCATION_UNKNOWN: "UN",
            PeakType.PUTATIVE_PRODUCT: "PR",
            PeakType.UNKNOWN: "UN",
        }
        return labels.get(peak_type, "?")

    def _format_plot(
        self,
        ax,
        compound: Compound,
        truncation_boundary: Optional[float]
    ) -> None:
        """Format plot axes and labels."""
        ax.set_xlabel("Time (minutes)", fontsize=12)
        ax.set_ylabel("Signal (counts)", fontsize=12)

        # No grid
        ax.grid(False)

        # Remove top and right spines
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

    def _add_compound_key(
        self,
        ax,
        reference: Compound,
        descendants: List[Compound],
        alias_map: dict
    ) -> None:
        """
        Add compound key on the right side of the plot.

        Shows descendants with MSA-style aligned sequences and alias letters.
        Minimal compounds at top, maximal (reference) at bottom.

        Parameters
        ----------
        ax : Axes
            Matplotlib axes
        reference : Compound
            Reference compound (shown at bottom)
        descendants : List[Compound]
            Descendants sorted minimal→maximal
        alias_map : dict
            Mapping from sequence to alias letter
        """
        # Get x position to the right of the plot
        xlim = ax.get_xlim()
        x_pos = xlim[1] + (xlim[1] - xlim[0]) * 0.02  # Just past right edge

        # Get y range for positioning
        ylim = ax.get_ylim()
        y_range = ylim[1] - ylim[0]

        # Total items: descendants + reference
        n_items = len(descendants) + 1
        y_spacing = min(y_range * 0.08, y_range / (n_items + 2))

        # Start from top for minimal compounds
        y_start = ylim[1] - y_spacing

        # Plot descendants (minimal at top)
        for i, desc in enumerate(descendants):
            alias = alias_map.get(desc.block_support_sequence, "?")
            y_pos = y_start - (i * y_spacing)

            # Get aligned sequence using SequenceAligner (single source of truth)
            aligned_seq = self._aligner.align_to_reference(desc, reference, HierarchyMode.BUILDING_BLOCK)

            # Format: "A: Leu---Pro"
            label = f"{alias}: {aligned_seq}  (L{desc.level})"

            ax.text(
                x_pos, y_pos, label,
                fontsize=8,
                fontfamily='monospace',
                verticalalignment='center',
                horizontalalignment='left',
                clip_on=False
            )

        # Plot reference at bottom (maximal)
        y_pos = y_start - (len(descendants) * y_spacing)
        ref_seq = reference.block_support_sequence
        label = f"*: {ref_seq}  (L{reference.level})"

        ax.text(
            x_pos, y_pos, label,
            fontsize=8,
            fontfamily='monospace',
            fontweight='bold',
            verticalalignment='center',
            horizontalalignment='left',
            clip_on=False
        )

    def _add_legend(self, ax, peaks: List[Peak], has_truncation_region: bool = False) -> None:
        """Add legend at bottom of plot, horizontal, outside the axes."""
        # Check if we have rejected peaks
        has_rejected = any(p.is_rejected for p in peaks)

        # Create legend handles using PeakAppearanceConfig
        marker_handles = PeakAppearanceConfig.create_legend_handles(include_rejected=has_rejected)
        line_handles = PeakAppearanceConfig.create_linestyle_handles(include_truncation=has_truncation_region)

        all_handles = marker_handles + line_handles

        if all_handles:
            # Legend at bottom center, below x-axis label
            ax.legend(
                handles=all_handles,
                loc='upper center',
                bbox_to_anchor=(0.5, -0.15),
                ncol=min(len(all_handles), 6),
                frameon=False,
                fontsize=8,
                handlelength=1.5,
                handleheight=1,
                columnspacing=1.0,
            )


def generate_diagnostic_plots(
    compounds: List[Compound],
    hierarchy: CompoundHierarchy,
    output_dir: Path,
    filename_attr: str = "block_support_sequence",
    truncation_margin: float = DEFAULT_TRUNCATION_MARGIN,
    show_height: bool = True,
) -> int:
    """
    Generate diagnostic PNGs for multiple compounds.

    Parameters
    ----------
    compounds : List[Compound]
        Compounds to generate plots for
    hierarchy : CompoundHierarchy
        Hierarchy for truncation boundary calculation
    output_dir : Path
        Directory to save PNG files
    filename_attr : str, optional
        Compound attribute to use for filename.
        Default "block_support_sequence", could also be "positional_block_sequence"
    truncation_margin : float, optional
        Margin beyond truncation positions (seconds)
    show_height : bool, optional
        Whether to show peak heights in annotations

    Returns
    -------
    int
        Number of plots generated
    """
    plotter = CompoundDiagnosticPlotter()
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    count = 0
    for compound in compounds:
        # Get filename from attribute
        name = getattr(compound, filename_attr, None)
        if not name:
            name = compound.positional_block_sequence or f"compound_{count}"

        # Sanitize filename (replace problematic characters)
        safe_name = name.replace("-", "_").replace(" ", "_").replace("/", "_")
        output_path = output_dir / f"{safe_name}.png"

        try:
            fig = plotter.plot(
                compound,
                hierarchy=hierarchy,
                truncation_margin=truncation_margin,
                output_path=output_path,
                show_height=show_height,
            )
            plt.close(fig)  # Free memory
            count += 1
        except Exception as e:
            print(f"Warning: Failed to generate plot for {name}: {e}")

    print(f"Generated {count} diagnostic plots in {output_dir}")
    return count
