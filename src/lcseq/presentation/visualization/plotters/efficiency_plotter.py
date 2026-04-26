"""
Efficiency visualization plotter.

Visualizes per-compound and library-wide coupling efficiency.
"""

from pathlib import Path
from typing import Optional, Tuple, List
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.axes import Axes
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.patches as mpatches

from .base_plotter import BasePlotter
from ....domain.models.coupling_efficiency import CompoundEfficiency, CycleEfficiency


class EfficiencyPlotter(BasePlotter):
    """
    Plotter for coupling efficiency visualization.

    Provides visualizations for:
    - Per-compound cycle efficiency (bar chart)
    - Per-compound efficiency breakdown (stacked areas)

    Examples
    --------
    >>> plotter = EfficiencyPlotter()
    >>> fig = plotter.plot_compound_efficiency(compound_efficiency)
    >>> plotter.save(fig, Path("efficiency.png"))
    """

    def __init__(
        self,
        figsize: Tuple[float, float] = (10, 6),
        dpi: int = 300,
        style: str = "seaborn-v0_8-paper",
    ):
        super().__init__(figsize=figsize, dpi=dpi, style=style)

        # Color scheme for efficiency
        self.cmap = LinearSegmentedColormap.from_list(
            "efficiency",
            [(0.8, 0.2, 0.2), (1.0, 0.8, 0.2), (0.2, 0.7, 0.3)],  # red -> yellow -> green
            N=100
        )

    def plot(self, compound_efficiency: CompoundEfficiency, **kwargs) -> Figure:
        """
        Default plot method - delegates to plot_compound_efficiency.

        Parameters
        ----------
        compound_efficiency : CompoundEfficiency
            Efficiency data for a compound

        Returns
        -------
        Figure
            Matplotlib figure
        """
        return self.plot_compound_efficiency(compound_efficiency, **kwargs)

    def plot_compound_efficiency(
        self,
        compound_efficiency: CompoundEfficiency,
        title: Optional[str] = None,
        show_areas: bool = True,
        figsize: Optional[Tuple[float, float]] = None
    ) -> Figure:
        """
        Plot per-cycle efficiency for a single compound.

        Creates a horizontal bar chart showing efficiency at each cycle,
        color-coded from red (low) to green (high).

        Parameters
        ----------
        compound_efficiency : CompoundEfficiency
            Efficiency data for the compound
        title : str, optional
            Custom title (default: uses compound sequence)
        show_areas : bool, optional
            Show area values on bars
        figsize : Tuple[float, float], optional
            Figure size override

        Returns
        -------
        Figure
            Matplotlib figure with efficiency bar chart
        """
        ce = compound_efficiency
        cycles = ce.cycle_efficiencies

        if not cycles:
            # No cycles to plot
            fig, ax = self.create_figure(figsize or (8, 4))
            ax.text(0.5, 0.5, "No cycle data available",
                    ha='center', va='center', transform=ax.transAxes)
            return fig

        # Create figure
        n_cycles = len(cycles)
        fig_height = max(4, 1.5 + n_cycles * 0.8)
        fig, ax = self.create_figure(figsize or (10, fig_height))

        # Prepare data
        labels = [f"Cycle {c.cycle}: {c.transition}" for c in cycles]
        efficiencies = [c.efficiency if c.is_valid else 0 for c in cycles]
        colors = [self.cmap(e) for e in efficiencies]

        # Create horizontal bar chart
        y_pos = np.arange(n_cycles)
        bars = ax.barh(y_pos, efficiencies, color=colors, edgecolor='black', linewidth=0.5)

        # Add efficiency labels on bars
        for i, (bar, cycle) in enumerate(zip(bars, cycles)):
            width = bar.get_width()
            if cycle.is_valid:
                label = f"{width:.1%}"
                if show_areas:
                    label += f"\n({cycle.area_passed:.0f}/{cycle.area_reached:.0f})"
            else:
                label = "N/A"

            # Position label inside or outside bar
            if width > 0.3:
                ax.text(width - 0.02, bar.get_y() + bar.get_height()/2,
                       label, ha='right', va='center', fontsize=9,
                       color='white', fontweight='bold')
            else:
                ax.text(width + 0.02, bar.get_y() + bar.get_height()/2,
                       label, ha='left', va='center', fontsize=9,
                       color='black')

        # Styling
        ax.set_yticks(y_pos)
        ax.set_yticklabels(labels)
        ax.set_xlim(0, 1.05)
        ax.set_xlabel("Coupling Efficiency")

        # Add reference lines
        ax.axvline(x=0.9, color='green', linestyle='--', alpha=0.5, label='90%')
        ax.axvline(x=0.8, color='orange', linestyle='--', alpha=0.5, label='80%')
        ax.axvline(x=0.7, color='red', linestyle='--', alpha=0.5, label='70%')

        # Title
        plot_title = title or f"Coupling Efficiency: {ce.sequence}"
        overall = ce.overall_efficiency
        if not np.isnan(overall):
            plot_title += f" (Overall: {overall:.1%})"
        ax.set_title(plot_title, fontsize=12, fontweight='bold')

        # Add legend for reference lines
        ax.legend(loc='lower right', fontsize=8)

        # Invert y-axis so cycle 0 is at top
        ax.invert_yaxis()

        fig.tight_layout()
        return fig

    def plot_area_breakdown(
        self,
        compound_efficiency: CompoundEfficiency,
        title: Optional[str] = None,
        figsize: Optional[Tuple[float, float]] = None
    ) -> Figure:
        """
        Plot area breakdown by truncation level.

        Creates a pie or bar chart showing how total signal is distributed
        across product and various truncation levels.

        Parameters
        ----------
        compound_efficiency : CompoundEfficiency
            Efficiency data for the compound
        title : str, optional
            Custom title
        figsize : Tuple[float, float], optional
            Figure size override

        Returns
        -------
        Figure
            Matplotlib figure with area breakdown
        """
        ce = compound_efficiency
        areas = ce.areas_by_level

        if not areas:
            fig, ax = self.create_figure(figsize or (8, 6))
            ax.text(0.5, 0.5, "No area data available",
                    ha='center', va='center', transform=ax.transAxes)
            return fig

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize or (12, 5))

        # Sort levels
        sorted_levels = sorted(areas.keys())
        labels = []
        values = []
        colors = []

        for lvl in sorted_levels:
            if lvl == 0:
                labels.append("L0 (Null)")
                colors.append('#d62728')  # red
            elif lvl == ce.level:
                labels.append(f"L{lvl} (Product)")
                colors.append('#2ca02c')  # green
            else:
                labels.append(f"L{lvl} (Truncation)")
                colors.append('#ff7f0e')  # orange
            values.append(areas[lvl])

        # Pie chart
        wedges, texts, autotexts = ax1.pie(
            values, labels=labels, autopct='%1.1f%%',
            colors=colors, startangle=90
        )
        ax1.set_title("Area Distribution by Level", fontsize=11, fontweight='bold')

        # Bar chart
        x_pos = np.arange(len(labels))
        bars = ax2.bar(x_pos, values, color=colors, edgecolor='black', linewidth=0.5)

        # Add value labels on bars
        for bar, val in zip(bars, values):
            height = bar.get_height()
            ax2.text(bar.get_x() + bar.get_width()/2, height,
                    f'{val:.0f}', ha='center', va='bottom', fontsize=9)

        ax2.set_xticks(x_pos)
        ax2.set_xticklabels(labels, rotation=45, ha='right')
        ax2.set_ylabel("Peak Area")
        ax2.set_title("Peak Areas by Level", fontsize=11, fontweight='bold')

        # Overall title
        plot_title = title or f"Area Breakdown: {ce.sequence}"
        fig.suptitle(plot_title, fontsize=13, fontweight='bold', y=1.02)

        fig.tight_layout()
        return fig

    def plot_pairwise_heatmap(
        self,
        compound_efficiency: CompoundEfficiency,
        title: Optional[str] = None,
        figsize: Optional[Tuple[float, float]] = None,
        annot: bool = True
    ) -> Figure:
        """
        Plot pairwise coupling efficiency as a heatmap for a single compound.

        Creates a matrix where rows are previous blocks and columns are next blocks.
        Each cell shows the efficiency of that transition.

        Parameters
        ----------
        compound_efficiency : CompoundEfficiency
            Efficiency data for the compound
        title : str, optional
            Custom title
        figsize : Tuple[float, float], optional
            Figure size override
        annot : bool, optional
            Show efficiency values in cells

        Returns
        -------
        Figure
            Matplotlib figure with pairwise heatmap
        """
        ce = compound_efficiency
        cycles = ce.cycle_efficiencies

        if not cycles:
            fig, ax = self.create_figure(figsize or (8, 6))
            ax.text(0.5, 0.5, "No cycle data available",
                    ha='center', va='center', transform=ax.transAxes)
            return fig

        # Collect all unique blocks
        all_blocks = set()
        all_blocks.add("N")  # Always include null/start
        for cycle in cycles:
            all_blocks.add(cycle.prev_block)
            all_blocks.add(cycle.next_block)

        # Sort blocks: N first, then alphabetically
        block_list = ["N"] + sorted([b for b in all_blocks if b != "N"])
        n_blocks = len(block_list)
        block_to_idx = {b: i for i, b in enumerate(block_list)}

        # Build efficiency matrix
        matrix = np.full((n_blocks, n_blocks), np.nan)
        for cycle in cycles:
            i = block_to_idx[cycle.prev_block]
            j = block_to_idx[cycle.next_block]
            matrix[i, j] = cycle.efficiency

        # Create figure
        fig, ax = self.create_figure(figsize or (8, 7))

        # Plot heatmap
        im = ax.imshow(matrix, cmap=self.cmap, aspect='equal', vmin=0, vmax=1)

        # Axis labels
        ax.set_xticks(np.arange(n_blocks))
        ax.set_yticks(np.arange(n_blocks))
        ax.set_xticklabels(block_list)
        ax.set_yticklabels(block_list)

        # Rotate x labels
        plt.setp(ax.get_xticklabels(), rotation=45, ha='right', rotation_mode='anchor')

        # Add text annotations
        if annot:
            for i in range(n_blocks):
                for j in range(n_blocks):
                    val = matrix[i, j]
                    if not np.isnan(val):
                        text_color = 'white' if val < 0.5 else 'black'
                        ax.text(j, i, f"{val:.0%}", ha='center', va='center',
                               fontsize=12, fontweight='bold', color=text_color)

        # Add colorbar
        cbar = fig.colorbar(im, ax=ax, shrink=0.8)
        cbar.set_label("Coupling Efficiency", fontsize=11)

        # Labels
        ax.set_xlabel("Next Block (coupled)", fontsize=12)
        ax.set_ylabel("Previous Block", fontsize=12)

        # Title
        overall = ce.overall_efficiency
        plot_title = title or f"Pairwise Coupling Efficiency: {ce.sequence}"
        if not np.isnan(overall):
            plot_title += f"\n(Overall: {overall:.1%})"
        ax.set_title(plot_title, fontsize=13, fontweight='bold')

        fig.tight_layout()
        return fig

    def plot_aggregated_heatmap(
        self,
        compound_efficiencies: List[CompoundEfficiency],
        title: Optional[str] = None,
        figsize: Optional[Tuple[float, float]] = None,
        show_counts: bool = True,
        min_count: int = 1,
        annotate: bool = True,
        all_blocks: Optional[List[str]] = None,
        row_blocks: Optional[List[str]] = None,
        col_blocks: Optional[List[str]] = None
    ) -> Figure:
        """
        Plot aggregated pairwise coupling efficiency heatmap across multiple compounds.

        Aggregates all transitions from all compounds into a single matrix.

        Parameters
        ----------
        compound_efficiencies : List[CompoundEfficiency]
            List of efficiency results to aggregate
        title : str, optional
            Custom title
        figsize : Tuple[float, float], optional
            Figure size override
        show_counts : bool, optional
            Annotate cells with observation count
        min_count : int, optional
            Minimum observations to show a cell (default 1)
        annotate : bool, optional
            Whether to show text annotations in cells (default True)
        all_blocks : List[str], optional
            Complete list of all blocks to include in matrix (used for both axes).
            If None, only blocks with valid transitions are shown.
        row_blocks : List[str], optional
            Ordered list of blocks for rows (prev blocks). Overrides all_blocks for rows.
        col_blocks : List[str], optional
            Ordered list of blocks for columns (next blocks). Overrides all_blocks for cols.

        Returns
        -------
        Figure
            Matplotlib figure with aggregated heatmap
        """
        if not compound_efficiencies:
            fig, ax = self.create_figure(figsize or (10, 8))
            ax.text(0.5, 0.5, "No data available",
                    ha='center', va='center', transform=ax.transAxes)
            return fig

        # Collect all transitions with their efficiencies and weights
        from collections import defaultdict
        transitions = defaultdict(list)  # (prev, next) -> [(efficiency, weight), ...]

        blocks_from_data = set()
        blocks_from_data.add("N")

        for ce in compound_efficiencies:
            for cycle in ce.cycle_efficiencies:
                if cycle.is_valid:
                    key = (cycle.prev_block, cycle.next_block)
                    transitions[key].append((cycle.efficiency, cycle.area_reached))
                    blocks_from_data.add(cycle.prev_block)
                    blocks_from_data.add(cycle.next_block)

        # Determine row ordering (prev blocks)
        if row_blocks is not None:
            row_list = list(row_blocks)
            if "N" not in row_list:
                row_list = ["N"] + row_list
        elif all_blocks is not None:
            block_set = set(all_blocks)
            block_set.add("N")
            row_list = ["N"] + sorted([b for b in block_set if b != "N"])
        else:
            row_list = ["N"] + sorted([b for b in blocks_from_data if b != "N"])

        # Determine column ordering (next blocks)
        if col_blocks is not None:
            col_list = list(col_blocks)
            if "N" not in col_list:
                col_list = ["N"] + col_list
        elif all_blocks is not None:
            block_set = set(all_blocks)
            block_set.add("N")
            col_list = ["N"] + sorted([b for b in block_set if b != "N"])
        else:
            col_list = ["N"] + sorted([b for b in blocks_from_data if b != "N"])

        n_rows = len(row_list)
        n_cols = len(col_list)
        row_to_idx = {b: i for i, b in enumerate(row_list)}
        col_to_idx = {b: i for i, b in enumerate(col_list)}

        # Build aggregated matrices
        eff_matrix = np.full((n_rows, n_cols), np.nan)
        count_matrix = np.zeros((n_rows, n_cols), dtype=int)

        for (prev_block, next_block), measurements in transitions.items():
            if prev_block not in row_to_idx or next_block not in col_to_idx:
                continue
            i = row_to_idx[prev_block]
            j = col_to_idx[next_block]
            count_matrix[i, j] = len(measurements)

            if len(measurements) >= min_count:
                # Weighted mean
                efficiencies = [e for e, w in measurements]
                weights = [w for e, w in measurements]
                total_weight = sum(weights)
                if total_weight > 0:
                    weighted_mean = sum(e * w for e, w in measurements) / total_weight
                    eff_matrix[i, j] = weighted_mean

        # Create figure
        fig, ax = self.create_figure(figsize or (16, 14))

        # Plot heatmap
        im = ax.imshow(eff_matrix, cmap=self.cmap, aspect='equal', vmin=0, vmax=1)

        # Axis labels
        ax.set_xticks(np.arange(n_cols))
        ax.set_yticks(np.arange(n_rows))
        ax.set_xticklabels(col_list, rotation=90)
        ax.set_yticklabels(row_list)

        # Colorbar
        cbar = fig.colorbar(im, ax=ax)
        cbar.set_label("Coupling Efficiency")

        # Labels and title
        ax.set_xlabel("Next Block (coupled)")
        ax.set_ylabel("Previous Block")
        plot_title = title or f"Coupling Efficiency (n={len(compound_efficiencies)} compounds)"
        ax.set_title(plot_title)

        fig.tight_layout()
        return fig

    def plot_efficiency_summary(
        self,
        compound_efficiencies: List[CompoundEfficiency],
        title: Optional[str] = None,
        figsize: Optional[Tuple[float, float]] = None,
        sort_by: str = "overall"
    ) -> Figure:
        """
        Plot efficiency summary for multiple compounds.

        Creates a heatmap-style visualization showing per-cycle efficiency
        for each compound.

        Parameters
        ----------
        compound_efficiencies : List[CompoundEfficiency]
            List of compound efficiency results
        title : str, optional
            Custom title
        figsize : Tuple[float, float], optional
            Figure size override
        sort_by : str, optional
            Sort compounds by: "overall", "sequence", "level", or "none"

        Returns
        -------
        Figure
            Matplotlib figure with efficiency summary heatmap
        """
        if not compound_efficiencies:
            fig, ax = self.create_figure(figsize or (10, 6))
            ax.text(0.5, 0.5, "No data available",
                    ha='center', va='center', transform=ax.transAxes)
            return fig

        # Sort compounds
        if sort_by == "overall":
            compound_efficiencies = sorted(
                compound_efficiencies,
                key=lambda x: x.overall_efficiency if not np.isnan(x.overall_efficiency) else 0,
                reverse=True
            )
        elif sort_by == "sequence":
            compound_efficiencies = sorted(compound_efficiencies, key=lambda x: x.sequence)
        elif sort_by == "level":
            compound_efficiencies = sorted(compound_efficiencies, key=lambda x: x.level, reverse=True)

        # Find max cycles across all compounds
        max_cycles = max(len(ce.cycle_efficiencies) for ce in compound_efficiencies)

        # Build matrix
        n_compounds = len(compound_efficiencies)
        matrix = np.full((n_compounds, max_cycles + 1), np.nan)  # +1 for overall column

        compound_labels = []
        for i, ce in enumerate(compound_efficiencies):
            compound_labels.append(f"{ce.sequence[:20]}..." if len(ce.sequence) > 20 else ce.sequence)
            for cycle_eff in ce.cycle_efficiencies:
                if cycle_eff.is_valid:
                    matrix[i, cycle_eff.cycle] = cycle_eff.efficiency
            # Add overall efficiency in last column
            if not np.isnan(ce.overall_efficiency):
                matrix[i, max_cycles] = ce.overall_efficiency

        # Create figure
        fig_height = max(6, 2 + n_compounds * 0.4)
        fig, ax = self.create_figure(figsize or (12, fig_height))

        # Plot heatmap
        im = ax.imshow(matrix, cmap=self.cmap, aspect='auto', vmin=0, vmax=1)

        # Axis labels
        cycle_labels = [f"Cycle {i}" for i in range(max_cycles)] + ["Overall"]
        ax.set_xticks(np.arange(max_cycles + 1))
        ax.set_xticklabels(cycle_labels, rotation=45, ha='right')
        ax.set_yticks(np.arange(n_compounds))
        ax.set_yticklabels(compound_labels)

        # Add colorbar
        cbar = fig.colorbar(im, ax=ax, shrink=0.8)
        cbar.set_label("Efficiency", fontsize=10)

        # Add text annotations
        for i in range(n_compounds):
            for j in range(max_cycles + 1):
                val = matrix[i, j]
                if not np.isnan(val):
                    text_color = 'white' if val < 0.5 else 'black'
                    ax.text(j, i, f"{val:.0%}", ha='center', va='center',
                           fontsize=8, color=text_color)

        # Title
        plot_title = title or "Coupling Efficiency Summary"
        ax.set_title(plot_title, fontsize=13, fontweight='bold')

        fig.tight_layout()
        return fig
