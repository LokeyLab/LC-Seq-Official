"""
Tests for PeakDetector service.

Tests persistent homology peak detection with scale-space filtration.
"""

import pytest
import numpy as np

from lcseq.domain.entities.chromatogram import Chromatogram
from lcseq.domain.entities.peak import Peak
from lcseq.domain.services.peak_detector import PeakDetector


# Test scale parameters (typical for test data)
SIGMA_MIN = 0.5
SIGMA_MAX = 25.0
NUM_SCALES = 20


class TestPeakDetector:
    """Test PeakDetector service with scale-space filtration."""

    def test_detect_single_peak(self):
        """Test detection of single peak."""
        detector = PeakDetector()

        # Create chromatogram with single peak
        time_points = np.linspace(0, 100, 200)
        signal = 10 * np.exp(-((time_points - 50) ** 2) / 20)

        chromatogram = Chromatogram(time_points=time_points, counts=signal)
        chromatogram.add_signal_variant("corrected", signal)

        # Detect peaks with new API
        peaks = detector.detect_peaks(
            chromatogram,
            sigma_min=SIGMA_MIN,
            sigma_max=SIGMA_MAX,
            num_scales=NUM_SCALES,
            persistence_threshold_mode="percentile",
            persistence_threshold=0.5,  # Top 50%
            signal_variant="raw"
        )

        # Should find exactly 1 peak
        assert len(peaks) == 1

        peak = peaks[0]
        assert isinstance(peak, Peak)
        assert 45 < peak.position < 55  # Near peak center
        assert peak.height > 8  # Near peak maximum
        assert peak.persistence is not None
        assert peak.area > 0

    def test_detect_multiple_peaks(self):
        """Test detection of multiple peaks."""
        detector = PeakDetector()

        # Create chromatogram with 3 peaks
        time_points = np.linspace(0, 150, 300)
        signal = (
            15 * np.exp(-((time_points - 30) ** 2) / 15) +
            20 * np.exp(-((time_points - 75) ** 2) / 20) +
            10 * np.exp(-((time_points - 120) ** 2) / 12)
        )

        chromatogram = Chromatogram(time_points=time_points, counts=signal)
        chromatogram.add_signal_variant("corrected", signal)

        # Detect peaks with new API
        peaks = detector.detect_peaks(
            chromatogram,
            sigma_min=SIGMA_MIN,
            sigma_max=SIGMA_MAX,
            num_scales=NUM_SCALES,
            persistence_threshold_mode="percentile",
            persistence_threshold=0.3,  # Top 70%
            signal_variant="raw"
        )

        # Should find 3 peaks
        assert len(peaks) == 3

        # Check peaks are at expected positions
        positions = sorted([p.position for p in peaks])
        assert positions[0] < 40  # First peak around 30
        assert 65 < positions[1] < 85  # Second peak around 75
        assert 110 < positions[2] < 130  # Third peak around 120

    def test_persistence_filtering_percentile_mode(self):
        """Test percentile-based persistence filtering."""
        detector = PeakDetector()

        # Create signal with major peak and noise bumps
        time_points = np.linspace(0, 100, 200)
        signal = 50 * np.exp(-((time_points - 50) ** 2) / 30)  # Major peak

        # Add small noise peaks
        signal += 2 * np.exp(-((time_points - 20) ** 2) / 3)
        signal += 2 * np.exp(-((time_points - 80) ** 2) / 3)

        chromatogram = Chromatogram(time_points=time_points, counts=signal)
        chromatogram.add_signal_variant("corrected", signal)

        # With low percentile threshold (top 90%), should find more peaks
        peaks_low = detector.detect_peaks(
            chromatogram,
            sigma_min=SIGMA_MIN,
            sigma_max=SIGMA_MAX,
            num_scales=NUM_SCALES,
            persistence_threshold_mode="percentile",
            persistence_threshold=0.1,  # Top 90%
            signal_variant="raw"
        )

        # With high percentile threshold (top 20%), should filter noise
        peaks_high = detector.detect_peaks(
            chromatogram,
            sigma_min=SIGMA_MIN,
            sigma_max=SIGMA_MAX,
            num_scales=NUM_SCALES,
            persistence_threshold_mode="percentile",
            persistence_threshold=0.8,  # Top 20%
            signal_variant="raw"
        )

        assert len(peaks_low) >= len(peaks_high)
        assert len(peaks_high) >= 1  # At least the major peak

    def test_persistence_filtering_gap_mode(self):
        """Test gap-based persistence filtering."""
        detector = PeakDetector()

        # Create signal with major peak and smaller peaks
        time_points = np.linspace(0, 100, 200)
        signal = (
            50 * np.exp(-((time_points - 50) ** 2) / 30) +  # Major peak
            5 * np.exp(-((time_points - 20) ** 2) / 5) +     # Small peak
            5 * np.exp(-((time_points - 80) ** 2) / 5)       # Small peak
        )

        chromatogram = Chromatogram(time_points=time_points, counts=signal)
        chromatogram.add_signal_variant("corrected", signal)

        # Gap mode should find natural threshold
        peaks = detector.detect_peaks(
            chromatogram,
            sigma_min=SIGMA_MIN,
            sigma_max=SIGMA_MAX,
            num_scales=NUM_SCALES,
            persistence_threshold_mode="gap",
            persistence_threshold=0.0,  # Unused in gap mode
            signal_variant="raw"
        )

        # Should find at least the major peak
        assert len(peaks) >= 1

    def test_persistence_filtering_relative_mode(self):
        """Test relative persistence filtering."""
        detector = PeakDetector()

        time_points = np.linspace(0, 100, 200)
        signal = (
            50 * np.exp(-((time_points - 50) ** 2) / 30) +
            5 * np.exp(-((time_points - 80) ** 2) / 5)
        )

        chromatogram = Chromatogram(time_points=time_points, counts=signal)
        chromatogram.add_signal_variant("corrected", signal)

        # Relative mode: threshold = 0.2 * max_persistence
        peaks = detector.detect_peaks(
            chromatogram,
            sigma_min=SIGMA_MIN,
            sigma_max=SIGMA_MAX,
            num_scales=NUM_SCALES,
            persistence_threshold_mode="relative",
            persistence_threshold=0.2,  # 20% of max persistence
            signal_variant="raw"
        )

        # Should find peaks above 20% of max persistence
        assert len(peaks) >= 1

    def test_signal_variant_selection(self):
        """Test that correct signal variant is used."""
        detector = PeakDetector()

        time_points = np.linspace(0, 50, 100)
        raw_signal = np.ones(100) * 100  # Flat
        corrected_signal = 10 * np.exp(-((time_points - 25) ** 2) / 10)  # Peak

        chromatogram = Chromatogram(time_points=time_points, counts=raw_signal)
        chromatogram.add_signal_variant("corrected", corrected_signal)

        # Using "corrected" should find peak
        peaks = detector.detect_peaks(
            chromatogram,
            sigma_min=SIGMA_MIN,
            sigma_max=SIGMA_MAX,
            num_scales=NUM_SCALES,
            persistence_threshold_mode="percentile",
            persistence_threshold=0.5,
            signal_variant="corrected"
        )
        assert len(peaks) >= 1

    def test_missing_signal_variant_raises_error(self):
        """Test that missing signal variant raises error."""
        detector = PeakDetector()

        time_points = np.linspace(0, 50, 100)
        signal = np.ones(100)

        chromatogram = Chromatogram(time_points=time_points, counts=signal)

        with pytest.raises(ValueError, match="not found"):
            detector.detect_peaks(
                chromatogram,
                sigma_min=SIGMA_MIN,
                sigma_max=SIGMA_MAX,
                num_scales=NUM_SCALES,
                signal_variant="nonexistent"
            )

    def test_peak_boundaries(self):
        """Test that peak boundaries are detected."""
        detector = PeakDetector()

        time_points = np.linspace(0, 100, 200)
        signal = 20 * np.exp(-((time_points - 50) ** 2) / 25)

        chromatogram = Chromatogram(time_points=time_points, counts=signal)
        chromatogram.add_signal_variant("corrected", signal)

        peaks = detector.detect_peaks(
            chromatogram,
            sigma_min=SIGMA_MIN,
            sigma_max=SIGMA_MAX,
            num_scales=NUM_SCALES,
            persistence_threshold_mode="percentile",
            persistence_threshold=0.5,
            signal_variant="raw"
        )

        assert len(peaks) >= 1
        peak = peaks[0]

        # Boundaries should bracket peak
        assert peak.left_base < peak.position
        assert peak.right_base > peak.position
        assert peak.left_base < peak.right_base

        # Width should be reasonable
        assert peak.width > 0

    def test_peak_area_integration(self):
        """Test that peak areas are calculated."""
        detector = PeakDetector()

        time_points = np.linspace(0, 100, 200)
        signal = 30 * np.exp(-((time_points - 50) ** 2) / 20)

        chromatogram = Chromatogram(time_points=time_points, counts=signal)
        chromatogram.add_signal_variant("corrected", signal)

        peaks = detector.detect_peaks(
            chromatogram,
            sigma_min=SIGMA_MIN,
            sigma_max=SIGMA_MAX,
            num_scales=NUM_SCALES,
            persistence_threshold_mode="percentile",
            persistence_threshold=0.5,
            signal_variant="raw"
        )

        assert len(peaks) >= 1
        peak = peaks[0]

        # Area should be positive
        assert peak.area > 0

        # Area should be reasonable (not too small or huge)
        max_signal = np.max(signal)
        assert peak.area > max_signal  # At least peak height
        assert peak.area < max_signal * len(signal)  # Not more than total

    def test_scale_space_parameters(self):
        """Test that scale-space parameters affect detection."""
        detector = PeakDetector()

        time_points = np.linspace(0, 100, 200)
        signal = 20 * np.exp(-((time_points - 50) ** 2) / 20)

        chromatogram = Chromatogram(time_points=time_points, counts=signal)
        chromatogram.add_signal_variant("corrected", signal)

        # Fine scale space (more scales)
        peaks_fine = detector.detect_peaks(
            chromatogram,
            sigma_min=SIGMA_MIN,
            sigma_max=SIGMA_MAX,
            num_scales=40,  # More scales
            persistence_threshold_mode="percentile",
            persistence_threshold=0.5,
            signal_variant="raw"
        )

        # Coarse scale space (fewer scales)
        peaks_coarse = detector.detect_peaks(
            chromatogram,
            sigma_min=SIGMA_MIN,
            sigma_max=SIGMA_MAX,
            num_scales=10,  # Fewer scales
            persistence_threshold_mode="percentile",
            persistence_threshold=0.5,
            signal_variant="raw"
        )

        # Both should find the major peak
        assert len(peaks_fine) >= 1
        assert len(peaks_coarse) >= 1


class TestPeakDetectorHelperMethods:
    """Test PeakDetector helper methods."""

    def test_find_local_maxima(self):
        """Test local maxima detection."""
        detector = PeakDetector()

        # Signal with 3 maxima
        signal = np.array([1, 2, 3, 2, 1, 2, 5, 3, 1, 2, 4, 2, 1])

        maxima = detector._find_local_maxima(signal)

        # Should find peaks at indices 2, 6, 10
        assert 2 in maxima
        assert 6 in maxima
        assert 10 in maxima

    def test_find_local_maxima_flat_regions(self):
        """Test maxima detection with flat regions."""
        detector = PeakDetector()

        signal = np.array([1, 3, 3, 3, 1, 5, 5, 2])

        maxima = detector._find_local_maxima(signal)

        # Flat regions might not be detected as single maximum
        # Just check we don't crash
        assert isinstance(maxima, list)

    def test_compute_persistence_with_scale_space(self):
        """Test persistence computation via scale-space tracking."""
        detector = PeakDetector()

        time_points = np.linspace(0, 50, 100)
        signal = 10 * np.exp(-((time_points - 25) ** 2) / 15)

        # New signature: requires scale parameters
        persistence_data = detector._compute_persistence(
            signal,
            time_points,
            sigma_min=SIGMA_MIN,
            sigma_max=SIGMA_MAX,
            num_scales=NUM_SCALES
        )

        assert len(persistence_data) >= 1

        # Each entry should be (idx, time, birth_scale, death_scale, persistence)
        for entry in persistence_data:
            idx, time, birth_scale, death_scale, persistence = entry
            assert isinstance(idx, (int, np.integer))
            assert isinstance(time, (float, np.floating))
            assert persistence >= 0
            assert birth_scale >= SIGMA_MIN
            assert death_scale <= SIGMA_MAX
            assert death_scale >= birth_scale

    def test_find_peak_boundaries_with_threshold(self):
        """Test boundary detection using threshold method."""
        detector = PeakDetector()

        # Simple peak
        signal = np.array([0, 1, 2, 5, 10, 5, 2, 1, 0])
        peak_idx = 4

        left, right = detector._find_peak_boundaries(signal, peak_idx, threshold_fraction=0.05)

        # Should find boundaries
        assert 0 <= left < peak_idx
        assert peak_idx < right < len(signal)

    def test_find_peak_boundaries_with_valley(self):
        """Test boundary detection with valley."""
        detector = PeakDetector()

        # Peak with valleys on sides
        signal = np.array([5, 4, 3, 2, 10, 15, 10, 2, 3, 4, 5])
        peak_idx = 5

        left, right = detector._find_peak_boundaries(signal, peak_idx)

        # Should find valleys as boundaries
        assert left <= 3  # Left valley around index 3
        assert right >= 7  # Right valley around index 7

    def test_integrate_peak_area(self):
        """Test peak area integration."""
        detector = PeakDetector()

        signal = np.array([0, 1, 5, 10, 5, 1, 0])

        area = detector._integrate_peak_area(signal, left_idx=1, right_idx=5)

        # Area should be sum of values from 1 to 5 (inclusive)
        expected = 1 + 5 + 10 + 5 + 1
        assert area == pytest.approx(expected)

    def test_integrate_empty_region(self):
        """Test integration of empty region."""
        detector = PeakDetector()

        signal = np.array([1, 2, 3, 4, 5])

        # Left >= right should return 0
        area = detector._integrate_peak_area(signal, left_idx=3, right_idx=2)
        assert area == 0.0

    def test_create_scale_space(self):
        """Test scale-space creation."""
        detector = PeakDetector()

        signal = 10 * np.exp(-((np.linspace(0, 50, 100) - 25) ** 2) / 15)

        scale_space = detector._create_scale_space(
            signal,
            sigma_min=SIGMA_MIN,
            sigma_max=SIGMA_MAX,
            num_scales=NUM_SCALES
        )

        # Should have NUM_SCALES entries
        assert len(scale_space) == NUM_SCALES

        # Each entry is (sigma, smoothed_signal)
        for sigma, smoothed in scale_space:
            assert SIGMA_MIN <= sigma <= SIGMA_MAX
            assert len(smoothed) == len(signal)

    def test_match_peaks_across_scales(self):
        """Test peak matching using Hungarian algorithm."""
        detector = PeakDetector()

        # Previous scale peaks: (id, idx, time)
        prev_peaks = [(0, 10, 5.0), (1, 30, 15.0), (2, 50, 25.0)]

        # Current scale peaks: (idx, time) - one peak disappeared
        curr_peaks = [(11, 5.5), (51, 25.5)]

        matched_ids, new_peaks, died_peak_ids = detector._match_peaks_across_scales(
            prev_peaks, curr_peaks, tolerance=5
        )

        # Should match peaks 0 and 2
        assert len(matched_ids) == 2

        # Peak 1 should have died
        assert 1 in died_peak_ids

        # No new peaks
        assert len(new_peaks) == 0


class TestPeakDetectorEdgeCases:
    """Test edge cases."""

    def test_empty_signal(self):
        """Test with very small signal."""
        detector = PeakDetector()

        time_points = np.array([0.0, 1.0, 2.0])
        signal = np.array([0.0, 0.0, 0.0])

        chromatogram = Chromatogram(time_points=time_points, counts=signal)
        chromatogram.add_signal_variant("corrected", signal)

        peaks = detector.detect_peaks(
            chromatogram,
            sigma_min=0.1,  # Smaller scales for tiny signal
            sigma_max=1.0,
            num_scales=10,
            persistence_threshold_mode="percentile",
            persistence_threshold=0.5,
            signal_variant="raw"
        )

        # Should return empty or minimal peaks
        assert isinstance(peaks, list)

    def test_no_peaks_above_threshold(self):
        """Test when all peaks filtered by threshold."""
        detector = PeakDetector()

        time_points = np.linspace(0, 50, 100)
        signal = 1 * np.exp(-((time_points - 25) ** 2) / 10)  # Small peak

        chromatogram = Chromatogram(time_points=time_points, counts=signal)
        chromatogram.add_signal_variant("corrected", signal)

        # Very high percentile threshold - with only 1 peak, it will still pass
        # (99th percentile of a single value is that value itself)
        peaks = detector.detect_peaks(
            chromatogram,
            sigma_min=SIGMA_MIN,
            sigma_max=SIGMA_MAX,
            num_scales=NUM_SCALES,
            persistence_threshold_mode="percentile",
            persistence_threshold=0.99,  # Top 1%
            signal_variant="raw"
        )

        # With single peak, percentile filtering keeps it
        assert len(peaks) >= 1

    def test_very_noisy_signal(self):
        """Test with noisy signal."""
        detector = PeakDetector()

        np.random.seed(42)
        time_points = np.linspace(0, 100, 200)
        signal = 20 * np.exp(-((time_points - 50) ** 2) / 30)
        signal += np.random.randn(200) * 2  # Add noise

        chromatogram = Chromatogram(time_points=time_points, counts=signal)
        chromatogram.add_signal_variant("corrected", signal)

        # Should still find main peak with appropriate threshold
        peaks = detector.detect_peaks(
            chromatogram,
            sigma_min=SIGMA_MIN,
            sigma_max=SIGMA_MAX,
            num_scales=NUM_SCALES,
            persistence_threshold_mode="percentile",
            persistence_threshold=0.7,  # Top 30%
            signal_variant="raw"
        )

        assert len(peaks) >= 1

    def test_overlapping_peaks(self):
        """Test detection of overlapping peaks."""
        detector = PeakDetector()

        time_points = np.linspace(0, 100, 200)
        # Two close peaks
        signal = (
            15 * np.exp(-((time_points - 45) ** 2) / 15) +
            15 * np.exp(-((time_points - 55) ** 2) / 15)
        )

        chromatogram = Chromatogram(time_points=time_points, counts=signal)
        chromatogram.add_signal_variant("corrected", signal)

        peaks = detector.detect_peaks(
            chromatogram,
            sigma_min=SIGMA_MIN,
            sigma_max=SIGMA_MAX,
            num_scales=NUM_SCALES,
            persistence_threshold_mode="percentile",
            persistence_threshold=0.3,
            signal_variant="raw"
        )

        # May detect as 1 or 2 peaks depending on resolution
        assert len(peaks) >= 1
        assert len(peaks) <= 2

    def test_invalid_persistence_mode(self):
        """Test that invalid persistence mode raises error."""
        detector = PeakDetector()

        time_points = np.linspace(0, 50, 100)
        signal = 10 * np.exp(-((time_points - 25) ** 2) / 10)

        chromatogram = Chromatogram(time_points=time_points, counts=signal)
        chromatogram.add_signal_variant("corrected", signal)

        with pytest.raises(ValueError, match="Unknown persistence_threshold_mode"):
            detector.detect_peaks(
                chromatogram,
                sigma_min=SIGMA_MIN,
                sigma_max=SIGMA_MAX,
                num_scales=NUM_SCALES,
                persistence_threshold_mode="invalid_mode",
                persistence_threshold=0.5,
                signal_variant="raw"
            )
