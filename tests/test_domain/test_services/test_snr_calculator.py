"""
Tests for SNRCalculator service.

Tests the single-source-of-truth SNR calculation implementation used
by all validation services.
"""

import pytest
import numpy as np

from lcseq.domain.entities.building_block import BuildingBlock
from lcseq.domain.entities.chromatogram import Chromatogram
from lcseq.domain.entities.compound import Compound
from lcseq.domain.entities.peak import Peak, PeakType
from lcseq.domain.services.snr_calculator import SNRCalculator


class TestSNRCalculator:
    """Test SNRCalculator service."""

    def test_calculate_normal_snr(self):
        """Test SNR calculation with normal values."""
        time_points = np.linspace(0, 100, 200)
        counts = np.ones(200)
        chrom = Chromatogram(time_points=time_points, counts=counts)

        compound = Compound([BuildingBlock.from_code(0, "Pro")], chrom)

        # Peak with height 100, background 10 → SNR = 10
        peak = Peak(position=40.0, left_base=35, right_base=45, height=100, area=200, peak_type=PeakType.PUTATIVE_PRODUCT)
        compound.selected_peak = peak

        snr = SNRCalculator.calculate(compound, background=10.0)

        assert snr == pytest.approx(10.0)

    def test_calculate_high_snr(self):
        """Test SNR calculation with high confidence signal."""
        time_points = np.linspace(0, 100, 200)
        counts = np.ones(200)
        chrom = Chromatogram(time_points=time_points, counts=counts)

        compound = Compound([BuildingBlock.from_code(0, "Pro")], chrom)

        # High SNR: peak=500, background=10 → SNR = 50
        peak = Peak(position=40.0, left_base=35, right_base=45, height=500, area=1000, peak_type=PeakType.PUTATIVE_PRODUCT)
        compound.selected_peak = peak

        snr = SNRCalculator.calculate(compound, background=10.0)

        assert snr == pytest.approx(50.0)
        assert snr > 10.0  # High confidence per THEORY.md Section 6.5

    def test_calculate_low_snr(self):
        """Test SNR calculation with low confidence signal (near noise floor)."""
        time_points = np.linspace(0, 100, 200)
        counts = np.ones(200)
        chrom = Chromatogram(time_points=time_points, counts=counts)

        compound = Compound([BuildingBlock.from_code(0, "Pro")], chrom)

        # Low SNR: peak=25, background=10 → SNR = 2.5
        peak = Peak(position=40.0, left_base=35, right_base=45, height=25, area=50, peak_type=PeakType.PUTATIVE_PRODUCT)
        compound.selected_peak = peak

        snr = SNRCalculator.calculate(compound, background=10.0)

        assert snr == pytest.approx(2.5)
        assert snr < 3.0  # Near noise floor per THEORY.md Section 6.5

    def test_calculate_no_selected_peak(self):
        """Test SNR calculation with no selected peak."""
        time_points = np.linspace(0, 100, 200)
        counts = np.ones(200)
        chrom = Chromatogram(time_points=time_points, counts=counts)

        compound = Compound([BuildingBlock.from_code(0, "Pro")], chrom)
        compound.selected_peak = None

        snr = SNRCalculator.calculate(compound, background=10.0)

        assert snr == 0.0

    def test_calculate_zero_background(self):
        """Test SNR calculation with zero background (edge case)."""
        time_points = np.linspace(0, 100, 200)
        counts = np.ones(200)
        chrom = Chromatogram(time_points=time_points, counts=counts)

        compound = Compound([BuildingBlock.from_code(0, "Pro")], chrom)

        peak = Peak(position=40.0, left_base=35, right_base=45, height=100, area=200, peak_type=PeakType.PUTATIVE_PRODUCT)
        compound.selected_peak = peak

        snr = SNRCalculator.calculate(compound, background=0.0)

        # Should return 0.0 for invalid background
        assert snr == 0.0

    def test_calculate_negative_background(self):
        """Test SNR calculation with negative background (edge case)."""
        time_points = np.linspace(0, 100, 200)
        counts = np.ones(200)
        chrom = Chromatogram(time_points=time_points, counts=counts)

        compound = Compound([BuildingBlock.from_code(0, "Pro")], chrom)

        peak = Peak(position=40.0, left_base=35, right_base=45, height=100, area=200, peak_type=PeakType.PUTATIVE_PRODUCT)
        compound.selected_peak = peak

        snr = SNRCalculator.calculate(compound, background=-5.0)

        # Should return 0.0 for invalid background
        assert snr == 0.0

    def test_calculate_moderate_snr(self):
        """Test SNR calculation with moderate confidence signal."""
        time_points = np.linspace(0, 100, 200)
        counts = np.ones(200)
        chrom = Chromatogram(time_points=time_points, counts=counts)

        compound = Compound([BuildingBlock.from_code(0, "Pro")], chrom)

        # Moderate SNR: peak=60, background=10 → SNR = 6.0
        peak = Peak(position=40.0, left_base=35, right_base=45, height=60, area=120, peak_type=PeakType.PUTATIVE_PRODUCT)
        compound.selected_peak = peak

        snr = SNRCalculator.calculate(compound, background=10.0)

        assert snr == pytest.approx(6.0)
        assert 3.0 <= snr <= 10.0  # Moderate confidence per THEORY.md Section 6.5

    def test_calculate_float_precision(self):
        """Test SNR calculation with non-integer values."""
        time_points = np.linspace(0, 100, 200)
        counts = np.ones(200)
        chrom = Chromatogram(time_points=time_points, counts=counts)

        compound = Compound([BuildingBlock.from_code(0, "Pro")], chrom)

        # Float values: peak=123.45, background=12.34 → SNR ≈ 10.004
        peak = Peak(position=40.0, left_base=35, right_base=45, height=123.45, area=246.9, peak_type=PeakType.PUTATIVE_PRODUCT)
        compound.selected_peak = peak

        snr = SNRCalculator.calculate(compound, background=12.34)

        assert snr == pytest.approx(10.004, rel=0.001)
        assert isinstance(snr, float)
