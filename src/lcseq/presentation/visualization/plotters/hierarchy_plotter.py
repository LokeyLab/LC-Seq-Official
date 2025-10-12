"""
Hierarchy plotter for visualizing compound hierarchy DAG structures.

Implementation based on THEORY.md Section 3.3, 4.2.
"""

from pathlib import Path
from typing import Optional, Dict
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
import networkx as nx

from .base_plotter import BasePlotter
from ....domain.models.compound_hierarchy import CompoundHierarchy


class HierarchyPlotter(BasePlotter):
    """
    Visualize CompoundHierarchy as directed acyclic graph (DAG).

    Shows:
    - Nodes (compounds) colored by validation status
    - Edges (truncation relationships)
    - Level-based layout (building-block or monomer mode)
    """

    # Color scheme for validation status
    STATUS_COLORS = {
        "valid": "#90EE90",  # Light green
        "invalid": "#FFB6C1",  # Light red
        "uncertain": "#FFD700",  # Gold
        "unvalidated": "#D3D3D3",  # Light gray
    }

    def plot(
        self,
        hierarchy: CompoundHierarchy,
        layout: str = "hierarchical",
        node_size: int = 500,
        title: Optional[str] = None,
    ) -> Figure:
        """
        Plot compound hierarchy as DAG.

        Parameters
        ----------
        hierarchy : CompoundHierarchy
            Hierarchy to visualize
        layout : str, optional
            Layout algorithm: "hierarchical", "spring", or "circular"
        node_size : int, optional
            Node size
        title : str, optional
            Plot title

        Returns
        -------
        Figure
            Matplotlib figure
        """
        fig, ax = self.create_figure(figsize=(12, 8))

        # Get NetworkX graph
        G = hierarchy.graph

        # Compute layout
        if layout == "hierarchical":
            pos = self._hierarchical_layout(G, hierarchy)
        elif layout == "spring":
            pos = nx.spring_layout(G, k=0.5, iterations=50)
        elif layout == "circular":
            pos = nx.circular_layout(G)
        else:
            pos = nx.kamada_kawai_layout(G)

        # Get node colors based on validation status (if available)
        node_colors = self._get_node_colors(G)

        # Draw graph
        nx.draw_networkx_nodes(
            G, pos, ax=ax, node_color=node_colors, node_size=node_size, alpha=0.9, edgecolors="black", linewidths=1.5
        )

        nx.draw_networkx_edges(G, pos, ax=ax, edge_color="gray", arrows=True, arrowsize=15, width=1.5, alpha=0.6)

        # Add labels
        labels = {node: node[:10] for node in G.nodes()}  # Truncate long IDs
        nx.draw_networkx_labels(G, pos, labels, ax=ax, font_size=8)

        # Style
        plot_title = title or f"Hierarchy DAG ({len(G.nodes)} compounds, {len(G.edges)} edges)"
        ax.set_title(plot_title, fontsize=14, fontweight="bold")
        ax.axis("off")

        # Add legend
        self._add_status_legend(ax)

        return fig

    def _hierarchical_layout(self, G: nx.DiGraph, hierarchy: CompoundHierarchy) -> Dict:
        """Compute hierarchical layout based on compound levels."""
        pos = {}

        # Get compounds at each level
        levels = hierarchy.get_level_distribution()
        max_level = max(levels.keys()) if levels else 0

        y_spacing = 1.0 / (max_level + 1) if max_level > 0 else 1.0

        for level, compounds in levels.items():
            n_compounds = len(compounds)
            x_spacing = 1.0 / (n_compounds + 1) if n_compounds > 0 else 0.5

            for i, compound_id in enumerate(sorted(compounds)):
                x = (i + 1) * x_spacing
                y = 1.0 - level * y_spacing
                pos[compound_id] = (x, y)

        return pos

    def _get_node_colors(self, G: nx.DiGraph) -> list:
        """Get node colors based on validation status."""
        colors = []

        for node in G.nodes():
            # Try to get validation status from node attributes
            node_data = G.nodes[node]
            status = node_data.get("validation_status", "unvalidated")
            colors.append(self.STATUS_COLORS.get(status, self.STATUS_COLORS["unvalidated"]))

        return colors

    def _add_status_legend(self, ax) -> None:
        """Add legend for validation status colors."""
        from matplotlib.patches import Patch

        legend_elements = [
            Patch(facecolor=color, edgecolor="black", label=status.capitalize())
            for status, color in self.STATUS_COLORS.items()
        ]

        ax.legend(handles=legend_elements, loc="upper right", framealpha=0.9, fontsize=10)
