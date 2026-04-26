"""
Lineage Peak Distribution Heatmap Plotter.

Visualizes how each compound's chromatogram signal distributes across
its ancestors (truncation products) as a compound × compound matrix.

Implementation based on plan: Lineage Peak Distribution Heatmap
"""

from dataclasses import dataclass
from typing import List, Optional, Tuple
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.colors import LinearSegmentedColormap

from ....domain.entities.compound import Compound
from ....domain.entities.peak import Peak, PeakType
from ....domain.models.compound_hierarchy import CompoundHierarchy


@dataclass
class LineagePeakMatrix:
    """
    Peak distribution matrix for a lineage.

    Represents how each compound's signal distributes across its ancestors.

    Attributes
    ----------
    compounds : List[Compound]
        Compounds in the lineage, ordered by level then sequence
    matrix : np.ndarray
        [n_compounds, n_compounds] raw peak areas
    matrix_pct : np.ndarray
        Same as matrix, but as percentages of row total
    labels : List[str]
        Compound labels for axes
    ancestry_mask : np.ndarray
        Boolean mask where True = valid cell (col is ancestor of row or diagonal)
    """
    compounds: List[Compound]
    matrix: np.ndarray
    matrix_pct: np.ndarray
    labels: List[str]
    ancestry_mask: np.ndarray


class LineageHeatmapPlotter:
    """
    Plots lineage peak distribution as a compound × compound heatmap.

    Shows how each compound's chromatogram signal distributes across
    its ancestors (truncation products), revealing synthesis quality
    and truncation patterns.

    The matrix has:
    - Rows: All compounds in lineage
    - Columns: Same compounds
    - Cell [row, col]: Area in row's chromatogram at col's RT (as % of row total)
    - Diagonal: Each compound's own product peak (purity)
    - Off-diagonal: Truncation peaks (respecting DAG ancestry)
    """

    def __init__(self, rt_tolerance: float = 0.02):
        """
        Initialize the plotter.

        Parameters
        ----------
        rt_tolerance : float
            Relative tolerance for matching peaks by retention time
        """
        self.rt_tolerance = rt_tolerance

    def build_matrix(
        self,
        compounds: List[Compound],
        hierarchy: CompoundHierarchy
    ) -> LineagePeakMatrix:
        """
        Build the compound × compound peak distribution matrix.

        Parameters
        ----------
        compounds : List[Compound]
            All compounds in the lineage
        hierarchy : CompoundHierarchy
            Hierarchy defining ancestry relationships

        Returns
        -------
        LineagePeakMatrix
            The peak distribution matrix
        """
        # Sort compounds by level (ascending), then by block_support_sequence (descending)
        # N at top, maximal at bottom (synthesis order)
        # Reverse within-level order to match stacked chromatogram when flipped
        sorted_compounds = sorted(
            compounds,
            key=lambda c: (-c.level, c.block_support_sequence)
        )
        sorted_compounds = list(reversed(sorted_compounds))

        n = len(sorted_compounds)
        matrix = np.zeros((n, n))
        ancestry_mask = np.zeros((n, n), dtype=bool)

        # Build labels
        labels = []
        for c in sorted_compounds:
            if c.level == 0 or c.is_null_compound:
                labels.append("N")
            else:
                labels.append(c.block_support_sequence)

        # Build matrix using peak classification (matched_compound_sequence)
        # This uses the same attribution as the diagnostic plotter
        for i, row_compound in enumerate(sorted_compounds):
            # Get descendants (truncations) of this compound
            descendants = hierarchy.get_descendants(row_compound)
            descendant_set = set(descendants)

            # Build lookup from sequence to column index
            seq_to_col = {c.block_support_sequence: j for j, c in enumerate(sorted_compounds)}

            for j, col_compound in enumerate(sorted_compounds):
                # Valid cell if col is a descendant (truncation) of row, or same compound
                is_self = (col_compound == row_compound)
                is_descendant = (col_compound in descendant_set)

                if is_self or is_descendant:
                    ancestry_mask[i, j] = True

            # Sum all accepted peaks by their matched_compound_sequence
            for peak in row_compound.detected_peaks:
                if not peak.is_accepted:
                    continue
                if peak.area is None or peak.area <= 0:
                    continue

                # Determine which column this peak belongs to
                if peak.peak_type == PeakType.PUTATIVE_PRODUCT:
                    # Product peak belongs to self (diagonal)
                    col_idx = i
                elif peak.peak_type == PeakType.NULL:
                    # NULL peak - find the null compound column
                    for j, c in enumerate(sorted_compounds):
                        if c.is_null_compound or c.level == 0:
                            col_idx = j
                            break
                    else:
                        continue
                elif peak.matched_compound_sequence:
                    # Truncation peak - use matched_compound_sequence
                    col_idx = seq_to_col.get(peak.matched_compound_sequence)
                    if col_idx is None:
                        continue
                else:
                    # Unknown peak without match - skip
                    continue

                # Add to matrix if this is a valid cell
                if ancestry_mask[i, col_idx]:
                    matrix[i, col_idx] += peak.area

        # Convert to percentages using total accepted peak area (same as diagnostic)
        # This ensures purity matches the diagnostic plotter's calculation
        row_totals = np.zeros(n)
        for i, compound in enumerate(sorted_compounds):
            # Sum all accepted peak areas (same calculation as diagnostic plotter)
            total = sum(
                p.area for p in compound.detected_peaks
                if p.is_accepted and p.area is not None and p.area > 0
            )
            row_totals[i] = total if total > 0 else 1.0

        with np.errstate(divide='ignore', invalid='ignore'):
            matrix_pct = (matrix / row_totals[:, np.newaxis]) * 100

        return LineagePeakMatrix(
            compounds=sorted_compounds,
            matrix=matrix,
            matrix_pct=matrix_pct,
            labels=labels,
            ancestry_mask=ancestry_mask
        )

    def _get_reference_rt(self, compound: Compound) -> Optional[float]:
        """
        Get the reference retention time for a compound.

        Uses selected_peak if available, otherwise looks for
        the product peak in detected_peaks.
        """
        if compound.selected_peak is not None:
            return compound.selected_peak.position

        # Fallback: find product peak
        for peak in compound.detected_peaks:
            if peak.peak_type == PeakType.PUTATIVE_PRODUCT:
                return peak.position
            if compound.is_null_compound and peak.peak_type == PeakType.NULL:
                return peak.position

        return None

    def _find_peak_at_position(
        self,
        peaks: List[Peak],
        target_rt: float,
        accepted_only: bool = False
    ) -> float:
        """
        Find peak area at a given retention time.

        Parameters
        ----------
        peaks : List[Peak]
            Peaks to search
        target_rt : float
            Target retention time
        accepted_only : bool
            If True, only consider accepted peaks (matches diagnostic plotter)

        Returns
        -------
        float
            Peak area if found, 0 otherwise
        """
        best_match = None
        best_distance = float('inf')

        for peak in peaks:
            # Skip rejected peaks if accepted_only is True
            if accepted_only and not peak.is_accepted:
                continue

            distance = abs(peak.position - target_rt)
            rel_distance = distance / max(target_rt, 1.0)

            if rel_distance < best_distance and rel_distance <= self.rt_tolerance:
                best_distance = rel_distance
                best_match = peak

        if best_match is None:
            return 0.0
        return best_match.area if best_match.area is not None else 0.0

    def plot(
        self,
        matrix: LineagePeakMatrix,
        title: Optional[str] = None,
        figsize: Tuple[int, int] = (12, 10),
        annotate: bool = True,
        annotation_threshold: float = 0.0
    ) -> Figure:
        """
        Plot the lineage peak distribution heatmap.

        Uses different colormaps for diagonal vs off-diagonal:
        - Diagonal (purity): RdYlGn (high=green=good)
        - Off-diagonal (truncation): RdYlGn_r (high=red=bad)

        Parameters
        ----------
        matrix : LineagePeakMatrix
            The matrix to plot (from build_matrix)
        title : str, optional
            Plot title
        figsize : Tuple[int, int]
            Figure size in inches
        annotate : bool
            Whether to annotate cells with percentages
        annotation_threshold : float
            Only annotate cells with value >= threshold (%)

        Returns
        -------
        Figure
            Matplotlib figure
        """
        # Use GridSpec for deterministic spacing
        # Columns: [heatmap, gap, purity_cbar, trunc_cbar]
        fig = plt.figure(figsize=figsize)
        gs = fig.add_gridspec(1, 4, width_ratios=[20, 2, 1, 1], wspace=0.02)
        ax = fig.add_subplot(gs[0])
        # gs[1] is empty spacer
        cax_purity = fig.add_subplot(gs[2])
        cax_trunc = fig.add_subplot(gs[3])
        n = len(matrix.compounds)

        # Create diagonal mask
        diagonal_mask = np.eye(n, dtype=bool)

        # Off-diagonal: valid cells that are not on diagonal
        off_diag_mask = matrix.ancestry_mask & ~diagonal_mask
        off_diag_data = np.ma.masked_where(~off_diag_mask, matrix.matrix_pct)

        # Diagonal: only diagonal cells that are valid
        diag_mask = matrix.ancestry_mask & diagonal_mask
        diag_data = np.ma.masked_where(~diag_mask, matrix.matrix_pct)

        # Plot off-diagonal first (reversed colormap: high truncation = red = bad)
        im_off = ax.imshow(
            off_diag_data,
            cmap='RdYlGn_r',
            vmin=0,
            vmax=100,
            aspect='auto'
        )

        # Plot diagonal on top (normal colormap: high purity = green = good)
        im_diag = ax.imshow(
            diag_data,
            cmap='RdYlGn',
            vmin=0,
            vmax=100,
            aspect='auto'
        )

        # Add annotations
        if annotate:
            for i in range(len(matrix.compounds)):
                for j in range(len(matrix.compounds)):
                    if matrix.ancestry_mask[i, j]:
                        val = matrix.matrix_pct[i, j]
                        if val >= annotation_threshold:
                            is_diagonal = (i == j)
                            # Text color based on background brightness
                            # RdYlGn: 0%=dark red, 50%=yellow, 100%=dark green
                            # RdYlGn_r: 0%=dark green, 50%=yellow, 100%=dark red
                            if is_diagonal:
                                # Dark backgrounds (red <30%, green >70%) need white text
                                text_color = 'black' if 30 <= val <= 70 else 'white'
                            else:
                                # Reversed colormap: dark backgrounds at opposite ends
                                text_color = 'black' if 30 <= val <= 70 else 'white'
                            ax.text(
                                j, i,
                                f'{val:.0f}%',
                                ha='center',
                                va='center',
                                fontsize=8,
                                color=text_color,
                                fontweight='bold' if is_diagonal else 'normal'
                            )

        # Highlight diagonal
        for i in range(len(matrix.compounds)):
            if matrix.ancestry_mask[i, i]:
                ax.add_patch(plt.Rectangle(
                    (i - 0.5, i - 0.5), 1, 1,
                    fill=False,
                    edgecolor='black',
                    linewidth=2
                ))

        # Labels
        ax.set_xticks(range(len(matrix.labels)))
        ax.set_xticklabels(matrix.labels, rotation=45, ha='right', fontsize=9)
        ax.set_yticks(range(len(matrix.labels)))
        ax.set_yticklabels(matrix.labels, fontsize=9)

        ax.set_xlabel('Ancestor Compound (truncation)', fontsize=11)
        ax.set_ylabel('Compound (chromatogram source)', fontsize=11)

        if title:
            ax.set_title(title, fontsize=13, fontweight='bold')
        else:
            ax.set_title('Lineage Peak Distribution', fontsize=13, fontweight='bold')

        # Colorbars - axes already created via GridSpec
        # Purity colorbar (left) - diagonal: high = green = good
        cbar_purity = plt.colorbar(im_diag, cax=cax_purity)
        cbar_purity.outline.set_visible(False)
        cax_purity.yaxis.set_ticks_position('left')
        cax_purity.yaxis.set_label_position('left')
        cax_purity.set_ylabel('Purity %', fontsize=9, rotation=90, va='center', ha='center')

        # Truncation colorbar (right) - off-diagonal: high = red = bad
        cbar_trunc = plt.colorbar(im_off, cax=cax_trunc)
        cbar_trunc.outline.set_visible(False)
        cax_trunc.yaxis.set_ticks_position('right')
        cax_trunc.yaxis.set_label_position('right')
        cax_trunc.set_ylabel('Truncation %', fontsize=9, rotation=270, va='center', ha='center')
        return fig

    def plot_from_compounds(
        self,
        compounds: List[Compound],
        hierarchy: CompoundHierarchy,
        **kwargs
    ) -> Figure:
        """
        Convenience method to build matrix and plot in one step.

        Parameters
        ----------
        compounds : List[Compound]
            All compounds in the lineage
        hierarchy : CompoundHierarchy
            Hierarchy defining ancestry relationships
        **kwargs
            Additional arguments passed to plot()

        Returns
        -------
        Figure
            Matplotlib figure
        """
        matrix = self.build_matrix(compounds, hierarchy)
        return self.plot(matrix, **kwargs)
