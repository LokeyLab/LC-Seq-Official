"""
Tests for PeakIntegrator service.

Tests peak area integration and boundary detection.
"""

import pytest
import numpy as np

from lcseq.domain.entities.chromatogram import Chromatogram
from lcseq.domain.services.peak_integrator import PeakIntegrator


class TestPeakIntegrator:
    """Test PeakIntegrator service."""

    def test_integrate_peak_auto_boundaries(self):
        """Test peak integration with automatic boundary detection."""
        integrator = PeakIntegrator()

        # Create chromatogram with single peak
        time_points = np.linspace(0, 100, 200)
        signal = 20 * np.exp(-((time_points - 50) ** 2) / 25)

        chromatogram = Chromatogram(time_points=time_points, counts=signal)
        chromatogram.add_signal_variant("corrected", signal)

        # Integrate at peak position
        left_base, right_base, area = integrator.integrate_peak(
            chromatogram, peak_position=50.0
        )

        # Boundaries should bracket peak
        assert left_base < 50.0
        assert right_base > 50.0

        # Area should be positive
        assert area > 0

    def test_integrate_peak_fixed_boundaries(self):
        """Test integration with provided boundaries."""
        integrator = PeakIntegrator()

        time_points = np.linspace(0, 100, 200)
        signal = 15 * np.exp(-((time_points - 50) ** 2) / 20)

        chromatogram = Chromatogram(time_points=time_points, counts=signal)
        chromatogram.add_signal_variant("corrected", signal)

        # Integrate with fixed boundaries
        left_base, right_base, area = integrator.integrate_peak(
            chromatogram, peak_position=50.0,
            left_boundary=40.0, right_boundary=60.0
        )

        # Should return provided boundaries
        assert left_base == pytest.approx(40.0, abs=1.0)
        assert right_base == pytest.approx(60.0, abs=1.0)

        # Area should be positive
        assert area > 0

    def test_integrate_multiple_positions(self):
        """Test integration at multiple peak positions."""
        integrator = PeakIntegrator()

        time_points = np.linspace(0, 150, 300)
        signal = (
            10 * np.exp(-((time_points - 40) ** 2) / 15) +
            15 * np.exp(-((time_points - 100) ** 2) / 20)
        )

        chromatogram = Chromatogram(time_points=time_points, counts=signal)
        chromatogram.add_signal_variant("corrected", signal)

        # Integrate first peak
        left1, right1, area1 = integrator.integrate_peak(chromatogram, peak_position=40.0)

        # Integrate second peak
        left2, right2, area2 = integrator.integrate_peak(chromatogram, peak_position=100.0)

        # Peaks should not overlap
        assert right1 < left2

        # Both areas should be positive
        assert area1 > 0
        assert area2 > 0

        # Second peak is larger
        assert area2 > area1

    def test_missing_signal_variant_raises_error(self):
        """Test that missing signal variant raises error."""
        integrator = PeakIntegrator()

        time_points = np.linspace(0, 50, 100)
        signal = np.ones(100)

        chromatogram = Chromatogram(time_points=time_points, counts=signal)

        with pytest.raises(ValueError, match="not found"):
            integrator.integrate_peak(chromatogram, peak_position=25.0, signal_variant="nonexistent")


class TestPeakIntegratorBoundaryDetection:
    """Test boundary detection methods."""

    def test_find_valley_boundaries(self):
        """Test valley boundary detection."""
        integrator = PeakIntegrator()

        # Signal with clear valleys
        signal = np.array([5, 4, 3, 2, 10, 15, 20, 15, 10, 2, 3, 4, 5])
        peak_idx = 6

        left_idx, right_idx = integrator._find_valley_boundaries(signal, peak_idx)

        # Should find valleys
        assert left_idx < peak_idx
        assert right_idx > peak_idx

        # Valleys should be around indices 3 and 9
        assert left_idx <= 4
        assert right_idx >= 8

    def test_threshold_boundary_detection(self):
        """Test 5% threshold boundary detection."""
        integrator = PeakIntegrator()

        # Peak without clear valleys
        signal = np.array([0, 2, 5, 10, 20, 10, 5, 2, 0])
        peak_idx = 4

        left_idx, right_idx = integrator._find_valley_boundaries(
            signal, peak_idx, threshold_fraction=0.05
        )

        # Should find boundaries based on threshold (5% of 20 = 1.0)
        assert left_idx < peak_idx
        assert right_idx > peak_idx

    def test_integrate_area(self):
        """Test area integration."""
        integrator = PeakIntegrator()

        signal = np.array([1.0, 2.0, 5.0, 10.0, 5.0, 2.0, 1.0])

        # Integrate from index 1 to 5
        area = integrator._integrate_area(signal, left_idx=1, right_idx=5)

        # Should sum values: 2 + 5 + 10 + 5 + 2
        assert area == pytest.approx(24.0)

    def test_integrate_empty_region(self):
        """Test integration of empty region."""
        integrator = PeakIntegrator()

        signal = np.array([1, 2, 3, 4, 5])

        # Empty region (left >= right)
        area = integrator._integrate_area(signal, left_idx=3, right_idx=2)
        assert area == 0.0

    def test_find_nearest_index(self):
        """Test finding nearest time point index."""
        integrator = PeakIntegrator()

        time_points = np.array([0.0, 10.0, 20.0, 30.0, 40.0])

        # Exact match
        idx = integrator._find_nearest_index(time_points, 20.0)
        assert idx == 2

        # Closest to 20.0
        idx = integrator._find_nearest_index(time_points, 18.0)
        assert idx == 2

        # Closest to 0.0
        idx = integrator._find_nearest_index(time_points, 2.0)
        assert idx == 0


class TestPeakIntegratorEdgeCases:
    """Test edge cases."""

    def test_integrate_at_signal_edge(self):
        """Test integration near signal boundaries."""
        integrator = PeakIntegrator()

        time_points = np.linspace(0, 50, 100)
        signal = 10 * np.exp(-((time_points - 5) ** 2) / 5)  # Peak near start

        chromatogram = Chromatogram(time_points=time_points, counts=signal)
        chromatogram.add_signal_variant("corrected", signal)

        # Should handle edge case
        left_base, right_base, area = integrator.integrate_peak(chromatogram, peak_position=5.0)

        assert left_base >= 0
        assert area > 0

    def test_integrate_flat_signal(self):
        """Test integration of flat signal."""
        integrator = PeakIntegrator()

        time_points = np.linspace(0, 50, 100)
        signal = np.ones(100) * 5.0

        chromatogram = Chromatogram(time_points=time_points, counts=signal)
        chromatogram.add_signal_variant("corrected", signal)

        left_base, right_base, area = integrator.integrate_peak(chromatogram, peak_position=25.0)

        # Should return some boundaries and area
        assert left_base < 25.0
        assert right_base > 25.0
        assert area >= 0

    def test_negative_area_handling(self):
        """Test that negative areas are clipped to zero."""
        integrator = PeakIntegrator()

        # Negative signal (shouldn't happen with proper baseline correction)
        signal = np.array([-5, -3, -1, 0, -1, -3, -5])

        area = integrator._integrate_area(signal, left_idx=1, right_idx=5)

        # Should be non-negative
        assert area >= 0
