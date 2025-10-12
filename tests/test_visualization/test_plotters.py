"""Basic tests for visualization plotters."""

import pytest
import numpy as np

# Skip if visualization dependencies not installed
pytest.importorskip("matplotlib")
pytest.importorskip("networkx")

from lcseq.domain.entities.chromatogram import Chromatogram
from lcseq.domain.entities.compound import Compound
from lcseq.domain.entities.building_block import BuildingBlock
from lcseq.domain.entities.peak import Peak, PeakType
from lcseq.presentation.visualization.plotters import (
    ChromatogramPlotter,
    BasePlotter,
)


class TestChromatogramPlotter:
    """Test chromatogram plotting."""

    def test_create_plotter(self):
        """Test plotter instantiation."""
        plotter = ChromatogramPlotter()
        assert plotter is not None
        assert isinstance(plotter, BasePlotter)

    def test_plot_simple_chromatogram(self):
        """Test plotting a simple chromatogram."""
        # Create test chromatogram
        time = np.array([0.0, 1.0, 2.0, 3.0, 4.0], dtype=np.float64)
        counts = np.array([10.0, 50.0, 100.0, 50.0, 10.0], dtype=np.float64)
        chrom = Chromatogram(time_points=time, counts=counts)

        # Plot
        plotter = ChromatogramPlotter()
        fig = plotter.plot(chrom, title="Test Chromatogram")

        assert fig is not None
        assert len(fig.axes) == 1

    def test_plot_with_peaks(self):
        """Test plotting chromatogram with peaks."""
        # Create test data
        time = np.array([0.0, 1.0, 2.0, 3.0, 4.0], dtype=np.float64)
        counts = np.array([10.0, 50.0, 100.0, 50.0, 10.0], dtype=np.float64)
        chrom = Chromatogram(time_points=time, counts=counts)

        # Create test peak
        peak = Peak(
            position=2.0,
            left_base=1.0,
            right_base=3.0,
            height=100.0,
            area=200.0,
            peak_type=PeakType.PUTATIVE_PRODUCT,
        )

        # Plot
        plotter = ChromatogramPlotter()
        fig = plotter.plot(chrom, peaks=[peak], title="Test with Peaks")

        assert fig is not None
        assert len(fig.axes) == 1


