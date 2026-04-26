"""
Base plotter class with template method pattern.

Provides common styling and functionality for all plotters.
"""

from abc import ABC, abstractmethod
from pathlib import Path
from typing import Optional, Tuple

# Use Agg backend for thread-safe plotting (required for parallel generation on macOS)
# Must be set before importing pyplot
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.axes import Axes


class BasePlotter(ABC):
    """
    Abstract base class for all plotters.

    Provides template method pattern for consistent plotting workflow:
    1. Create figure
    2. Plot data (subclass implements)
    3. Style plot
    4. Save or show

    Attributes
    ----------
    figsize : Tuple[float, float]
        Default figure size (width, height) in inches
    dpi : int
        Resolution for saved figures
    style : str
        Matplotlib style to use
    """

    def __init__(
        self,
        figsize: Tuple[float, float] = (10, 6),
        dpi: int = 300,
        style: str = "seaborn-v0_8-paper",
    ):
        """
        Initialize plotter with styling parameters.

        Parameters
        ----------
        figsize : Tuple[float, float], optional
            Figure size (width, height) in inches
        dpi : int, optional
            Resolution for saved figures
        style : str, optional
            Matplotlib style name
        """
        self.figsize = figsize
        self.dpi = dpi
        self.style = style

        # Apply matplotlib style
        plt.style.use(self.style)

    def create_figure(self, figsize: Optional[Tuple[float, float]] = None) -> Tuple[Figure, Axes]:
        """
        Create a matplotlib figure and axes.

        Parameters
        ----------
        figsize : Tuple[float, float], optional
            Override default figure size

        Returns
        -------
        Tuple[Figure, Axes]
            Created figure and axes
        """
        size = figsize or self.figsize
        fig, ax = plt.subplots(figsize=size)
        return fig, ax

    @abstractmethod
    def plot(self, *args, **kwargs) -> Figure:
        """
        Plot data. Must be implemented by subclasses.

        Returns
        -------
        Figure
            Matplotlib figure object
        """
        pass

    def save(self, fig: Figure, output_path: Path, tight: bool = True) -> None:
        """
        Save figure to file.

        Parameters
        ----------
        fig : Figure
            Matplotlib figure to save
        output_path : Path
            Output file path
        tight : bool, optional
            Use tight layout
        """
        output_path.parent.mkdir(parents=True, exist_ok=True)

        if tight:
            fig.tight_layout()

        fig.savefig(output_path, dpi=self.dpi, bbox_inches="tight")
        plt.close(fig)

    def show(self, fig: Figure) -> None:
        """
        Display figure interactively.

        Parameters
        ----------
        fig : Figure
            Matplotlib figure to show
        """
        fig.tight_layout()
        plt.show()

    def apply_common_styling(
        self,
        ax: Axes,
        title: Optional[str] = None,
        xlabel: Optional[str] = None,
        ylabel: Optional[str] = None,
        grid: bool = True,
    ) -> None:
        """
        Apply common styling to axes.

        Parameters
        ----------
        ax : Axes
            Matplotlib axes to style
        title : str, optional
            Plot title
        xlabel : str, optional
            X-axis label
        ylabel : str, optional
            Y-axis label
        grid : bool, optional
            Show grid
        """
        if title:
            ax.set_title(title, fontsize=14, fontweight="bold")

        if xlabel:
            ax.set_xlabel(xlabel, fontsize=12)

        if ylabel:
            ax.set_ylabel(ylabel, fontsize=12)

        if grid:
            ax.grid(True, alpha=0.3)

        # Improve tick label readability
        ax.tick_params(labelsize=10)
