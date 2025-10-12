"""
Tests for PurityCalculator service.

Tests the single-source-of-truth purity calculation implementation used
by all validation services.
"""

import pytest
import numpy as np

from lcseq.domain.entities.building_block import BuildingBlock
from lcseq.domain.entities.chromatogram import Chromatogram
from lcseq.domain.entities.compound import Compound
from lcseq.domain.entities.peak import Peak, PeakType
from lcseq.domain.services.purity_calculator import PurityCalculator


class TestPurityCalculator:
    """Test PurityCalculator service."""

    def test_calculate_mixed_peaks(self):
        """Test purity calculation with mixed peak types."""
        time_points = np.linspace(0, 100, 200)
        counts = np.ones(200)
        chrom = Chromatogram(time_points=time_points, counts=counts)

        compound = Compound([BuildingBlock.from_code(0, "Pro")], chrom)

        # Add peaks: product=100, truncation=50, null=30 → purity = 100/180 ≈ 0.556
        compound.detected_peaks = [
            Peak(position=10.0, left_base=5, right_base=15, height=20, area=30, peak_type=PeakType.NULL),
            Peak(position=20.0, left_base=15, right_base=25, height=30, area=50, peak_type=PeakType.TRUNCATION),
            Peak(position=40.0, left_base=35, right_base=45, height=50, area=100, peak_type=PeakType.PUTATIVE_PRODUCT),
        ]

        purity = PurityCalculator.calculate(compound)

        assert purity == pytest.approx(100 / 180, rel=0.01)

    def test_calculate_no_peaks(self):
        """Test purity calculation with no peaks."""
        time_points = np.linspace(0, 100, 200)
        counts = np.ones(200)
        chrom = Chromatogram(time_points=time_points, counts=counts)

        compound = Compound([BuildingBlock.from_code(0, "Pro")], chrom)
        compound.detected_peaks = []

        purity = PurityCalculator.calculate(compound)

        assert purity == 0.0

    def test_calculate_perfect_purity(self):
        """Test purity calculation with only product peak (purity = 1.0)."""
        time_points = np.linspace(0, 100, 200)
        counts = np.ones(200)
        chrom = Chromatogram(time_points=time_points, counts=counts)

        compound = Compound([BuildingBlock.from_code(0, "Pro")], chrom)

        # Only product peak → purity = 1.0
        compound.detected_peaks = [
            Peak(position=40.0, left_base=35, right_base=45, height=50, area=100, peak_type=PeakType.PUTATIVE_PRODUCT),
        ]

        purity = PurityCalculator.calculate(compound)

        assert purity == 1.0

    def test_calculate_zero_purity(self):
        """Test purity calculation with no product peaks (purity = 0.0)."""
        time_points = np.linspace(0, 100, 200)
        counts = np.ones(200)
        chrom = Chromatogram(time_points=time_points, counts=counts)

        compound = Compound([BuildingBlock.from_code(0, "Pro")], chrom)

        # Only impurity peaks → purity = 0.0
        compound.detected_peaks = [
            Peak(position=10.0, left_base=5, right_base=15, height=20, area=30, peak_type=PeakType.NULL),
            Peak(position=20.0, left_base=15, right_base=25, height=30, area=50, peak_type=PeakType.TRUNCATION),
        ]

        purity = PurityCalculator.calculate(compound)

        assert purity == 0.0

    def test_calculate_multiple_product_peaks(self):
        """Test purity calculation with multiple product peaks."""
        time_points = np.linspace(0, 100, 200)
        counts = np.ones(200)
        chrom = Chromatogram(time_points=time_points, counts=counts)

        compound = Compound([BuildingBlock.from_code(0, "Pro")], chrom)

        # Two product peaks + one truncation
        compound.detected_peaks = [
            Peak(position=20.0, left_base=15, right_base=25, height=30, area=50, peak_type=PeakType.TRUNCATION),
            Peak(position=40.0, left_base=35, right_base=45, height=50, area=100, peak_type=PeakType.PUTATIVE_PRODUCT),
            Peak(position=60.0, left_base=55, right_base=65, height=40, area=80, peak_type=PeakType.PUTATIVE_PRODUCT),
        ]

        purity = PurityCalculator.calculate(compound)

        # Product = 100 + 80 = 180, Total = 50 + 100 + 80 = 230, Purity = 180/230
        assert purity == pytest.approx(180 / 230, rel=0.01)

    def test_calculate_clipping(self):
        """Test that purity is clipped to [0, 1] range."""
        time_points = np.linspace(0, 100, 200)
        counts = np.ones(200)
        chrom = Chromatogram(time_points=time_points, counts=counts)

        compound = Compound([BuildingBlock.from_code(0, "Pro")], chrom)

        # Only product peak → purity = 1.0
        compound.detected_peaks = [
            Peak(position=40.0, left_base=35, right_base=45, height=50, area=100, peak_type=PeakType.PUTATIVE_PRODUCT),
        ]

        purity = PurityCalculator.calculate(compound)

        # Should be clipped to exactly 1.0
        assert purity >= 0.0
        assert purity <= 1.0
        assert purity == 1.0

    def test_calculate_zero_area(self):
        """Test purity calculation with zero total area."""
        time_points = np.linspace(0, 100, 200)
        counts = np.ones(200)
        chrom = Chromatogram(time_points=time_points, counts=counts)

        compound = Compound([BuildingBlock.from_code(0, "Pro")], chrom)

        # Peaks with zero area
        compound.detected_peaks = [
            Peak(position=10.0, left_base=5, right_base=15, height=20, area=0, peak_type=PeakType.NULL),
            Peak(position=40.0, left_base=35, right_base=45, height=50, area=0, peak_type=PeakType.PUTATIVE_PRODUCT),
        ]

        purity = PurityCalculator.calculate(compound)

        # Should return 0.0 when total area is 0
        assert purity == 0.0
