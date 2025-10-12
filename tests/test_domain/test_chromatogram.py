"""
Tests for Chromatogram entity.

Tests implementation against THEORY.md Section 2.1, 5.0.1, 2.3.1, 2.3.2 specifications.
"""

import pytest
import numpy as np
from lcseq.domain.entities.chromatogram import Chromatogram


class TestChromatogramCreation:
    """Test Chromatogram instantiation and validation."""

    def test_create_simple_chromatogram(self):
        """Basic chromatogram with time and counts."""
        time_points = np.array([0.0, 1.0, 2.0], dtype=np.float64)
        counts = np.array([10.0, 20.0, 15.0], dtype=np.float64)

        chrom = Chromatogram(time_points=time_points, counts=counts)

        assert len(chrom.time_points) == 3
        assert len(chrom.counts) == 3
        np.testing.assert_array_equal(chrom.time_points, time_points)
        np.testing.assert_array_equal(chrom.counts, counts)

    def test_create_with_signal_variants(self):
        """Chromatogram with additional signal variants."""
        time_points = np.array([0.0, 1.0, 2.0], dtype=np.float64)
        counts = np.array([10.0, 20.0, 15.0], dtype=np.float64)
        corrected = np.array([5.0, 15.0, 10.0], dtype=np.float64)

        chrom = Chromatogram(
            time_points=time_points,
            counts=counts,
            signal_variants={"corrected": corrected}
        )

        assert "corrected" in chrom.signal_variants
        np.testing.assert_array_equal(chrom.signal_variants["corrected"], corrected)

    def test_create_from_list_converts_to_array(self):
        """Lists are automatically converted to numpy arrays."""
        chrom = Chromatogram(
            time_points=[0.0, 1.0, 2.0],
            counts=[10.0, 20.0, 15.0]
        )

        assert isinstance(chrom.time_points, np.ndarray)
        assert isinstance(chrom.counts, np.ndarray)
        assert chrom.time_points.dtype == np.float64
        assert chrom.counts.dtype == np.float64

    def test_create_with_absolute_time_values(self):
        """Time points can start at arbitrary values (THEORY.md Section 2.3.1)."""
        # Time points starting at 600 seconds (10 minutes)
        time_points = np.array([600.0, 660.0, 720.0], dtype=np.float64)
        counts = np.array([10.0, 20.0, 15.0], dtype=np.float64)

        chrom = Chromatogram(time_points=time_points, counts=counts)

        assert chrom.time_points[0] == 600.0
        assert chrom.time_range == (600.0, 720.0)


class TestChromatogramValidation:
    """Test Chromatogram validation rules."""

    def test_empty_time_points_raises_error(self):
        """Chromatogram must have at least one time point."""
        with pytest.raises(ValueError, match="at least one time point"):
            Chromatogram(
                time_points=np.array([], dtype=np.float64),
                counts=np.array([], dtype=np.float64)
            )

    def test_length_mismatch_raises_error(self):
        """Time points and counts must have same length."""
        with pytest.raises(ValueError, match="must have same length"):
            Chromatogram(
                time_points=np.array([0.0, 1.0, 2.0], dtype=np.float64),
                counts=np.array([10.0, 20.0], dtype=np.float64)  # Wrong length
            )

    def test_non_increasing_time_raises_error(self):
        """Time points must be strictly increasing (THEORY.md Section 2.3.1)."""
        with pytest.raises(ValueError, match="strictly increasing"):
            Chromatogram(
                time_points=np.array([0.0, 2.0, 1.0], dtype=np.float64),  # Not increasing
                counts=np.array([10.0, 20.0, 15.0], dtype=np.float64)
            )

    def test_duplicate_time_points_raises_error(self):
        """Duplicate time points violate strictly increasing requirement."""
        with pytest.raises(ValueError, match="strictly increasing"):
            Chromatogram(
                time_points=np.array([0.0, 1.0, 1.0], dtype=np.float64),  # Duplicate
                counts=np.array([10.0, 20.0, 15.0], dtype=np.float64)
            )

    def test_signal_variant_length_mismatch_raises_error(self):
        """Signal variants must match time_points length."""
        with pytest.raises(ValueError, match="expected 3"):
            Chromatogram(
                time_points=np.array([0.0, 1.0, 2.0], dtype=np.float64),
                counts=np.array([10.0, 20.0, 15.0], dtype=np.float64),
                signal_variants={"corrected": np.array([5.0, 15.0], dtype=np.float64)}  # Wrong length
            )


class TestChromatogramProperties:
    """Test Chromatogram computed properties."""

    def test_duration_single_point(self):
        """Duration is 0 for single time point."""
        chrom = Chromatogram(
            time_points=np.array([5.0], dtype=np.float64),
            counts=np.array([10.0], dtype=np.float64)
        )

        assert chrom.duration == 0.0

    def test_duration_multiple_points(self):
        """Duration is last - first time point."""
        chrom = Chromatogram(
            time_points=np.array([0.0, 30.0, 60.0], dtype=np.float64),
            counts=np.array([10.0, 20.0, 15.0], dtype=np.float64)
        )

        assert chrom.duration == 60.0

    def test_duration_absolute_time(self):
        """Duration works with absolute time values."""
        chrom = Chromatogram(
            time_points=np.array([600.0, 900.0, 1200.0], dtype=np.float64),
            counts=np.array([10.0, 20.0, 15.0], dtype=np.float64)
        )

        assert chrom.duration == 600.0

    def test_time_range(self):
        """Time range returns (start, end) tuple."""
        chrom = Chromatogram(
            time_points=np.array([10.0, 20.0, 30.0], dtype=np.float64),
            counts=np.array([10.0, 20.0, 15.0], dtype=np.float64)
        )

        assert chrom.time_range == (10.0, 30.0)

    def test_time_range_variable_boundaries(self):
        """Different chromatograms can have different time ranges (THEORY.md Section 2.3.2)."""
        chrom1 = Chromatogram(
            time_points=np.array([0.0, 60.0], dtype=np.float64),
            counts=np.array([10.0, 20.0], dtype=np.float64)
        )

        chrom2 = Chromatogram(
            time_points=np.array([600.0, 1200.0], dtype=np.float64),
            counts=np.array([15.0, 25.0], dtype=np.float64)
        )

        assert chrom1.time_range != chrom2.time_range
        assert chrom1.time_range[0] < chrom2.time_range[0]


class TestSignalVariantManagement:
    """Test signal variant management methods."""

    def test_get_raw_signal(self):
        """get_signal('raw') returns original counts."""
        chrom = Chromatogram(
            time_points=np.array([0.0, 1.0, 2.0], dtype=np.float64),
            counts=np.array([10.0, 20.0, 15.0], dtype=np.float64)
        )

        signal = chrom.get_signal("raw")
        np.testing.assert_array_equal(signal, chrom.counts)

    def test_get_variant_signal(self):
        """get_signal returns requested variant."""
        corrected = np.array([5.0, 15.0, 10.0], dtype=np.float64)
        chrom = Chromatogram(
            time_points=np.array([0.0, 1.0, 2.0], dtype=np.float64),
            counts=np.array([10.0, 20.0, 15.0], dtype=np.float64),
            signal_variants={"corrected": corrected}
        )

        signal = chrom.get_signal("corrected")
        np.testing.assert_array_equal(signal, corrected)

    def test_get_nonexistent_variant_raises_error(self):
        """Requesting non-existent variant raises KeyError."""
        chrom = Chromatogram(
            time_points=np.array([0.0, 1.0, 2.0], dtype=np.float64),
            counts=np.array([10.0, 20.0, 15.0], dtype=np.float64)
        )

        with pytest.raises(KeyError, match="not found"):
            chrom.get_signal("nonexistent")

    def test_add_signal_variant(self):
        """Add new signal variant to chromatogram."""
        chrom = Chromatogram(
            time_points=np.array([0.0, 1.0, 2.0], dtype=np.float64),
            counts=np.array([10.0, 20.0, 15.0], dtype=np.float64)
        )

        derivative = np.array([10.0, -5.0, -10.0], dtype=np.float64)
        chrom.add_signal_variant("derivative", derivative)

        assert "derivative" in chrom.signal_variants
        np.testing.assert_array_equal(chrom.signal_variants["derivative"], derivative)

    def test_add_multiple_signal_variants(self):
        """Add multiple signal variants (THEORY.md Section 5.0.1)."""
        chrom = Chromatogram(
            time_points=np.array([0.0, 1.0, 2.0], dtype=np.float64),
            counts=np.array([10.0, 20.0, 15.0], dtype=np.float64)
        )

        chrom.add_signal_variant("corrected", np.array([5.0, 15.0, 10.0], dtype=np.float64))
        chrom.add_signal_variant("derivative", np.array([10.0, -5.0, -10.0], dtype=np.float64))
        chrom.add_signal_variant("derivative_2", np.array([-15.0, -5.0, 5.0], dtype=np.float64))

        assert len(chrom.signal_variants) == 3
        assert "corrected" in chrom.signal_variants
        assert "derivative" in chrom.signal_variants
        assert "derivative_2" in chrom.signal_variants

    def test_add_variant_wrong_length_raises_error(self):
        """Signal variant must match time_points length."""
        chrom = Chromatogram(
            time_points=np.array([0.0, 1.0, 2.0], dtype=np.float64),
            counts=np.array([10.0, 20.0, 15.0], dtype=np.float64)
        )

        with pytest.raises(ValueError, match="must match"):
            chrom.add_signal_variant("bad", np.array([1.0, 2.0], dtype=np.float64))

    def test_has_signal_variant_raw(self):
        """'raw' variant always exists."""
        chrom = Chromatogram(
            time_points=np.array([0.0, 1.0], dtype=np.float64),
            counts=np.array([10.0, 20.0], dtype=np.float64)
        )

        assert chrom.has_signal_variant("raw") is True

    def test_has_signal_variant_exists(self):
        """Check if variant exists."""
        chrom = Chromatogram(
            time_points=np.array([0.0, 1.0], dtype=np.float64),
            counts=np.array([10.0, 20.0], dtype=np.float64),
            signal_variants={"corrected": np.array([5.0, 15.0], dtype=np.float64)}
        )

        assert chrom.has_signal_variant("corrected") is True

    def test_has_signal_variant_not_exists(self):
        """Check returns False for non-existent variant."""
        chrom = Chromatogram(
            time_points=np.array([0.0, 1.0], dtype=np.float64),
            counts=np.array([10.0, 20.0], dtype=np.float64)
        )

        assert chrom.has_signal_variant("corrected") is False


class TestChromatogramSpecialMethods:
    """Test special methods (__len__, __repr__)."""

    def test_len_returns_number_of_points(self):
        """len() returns number of time points."""
        chrom = Chromatogram(
            time_points=np.array([0.0, 1.0, 2.0, 3.0], dtype=np.float64),
            counts=np.array([10.0, 20.0, 15.0, 12.0], dtype=np.float64)
        )

        assert len(chrom) == 4

    def test_len_single_point(self):
        """len() works for single time point."""
        chrom = Chromatogram(
            time_points=np.array([0.0], dtype=np.float64),
            counts=np.array([10.0], dtype=np.float64)
        )

        assert len(chrom) == 1

    def test_repr_shows_key_info(self):
        """repr() shows n_points, time_range, and variants."""
        chrom = Chromatogram(
            time_points=np.array([0.0, 60.0], dtype=np.float64),
            counts=np.array([10.0, 20.0], dtype=np.float64),
            signal_variants={"corrected": np.array([5.0, 15.0], dtype=np.float64)}
        )

        repr_str = repr(chrom)
        assert "n_points=2" in repr_str
        assert "time_range=(0.0, 60.0)" in repr_str
        assert "variants=" in repr_str
        assert "raw" in repr_str
        assert "corrected" in repr_str


class TestVariableSignalBoundaries:
    """Test variable signal boundaries (THEORY.md Section 2.3.2)."""

    def test_signals_can_have_different_start_times(self):
        """Different chromatograms can start at different times."""
        chrom1 = Chromatogram(
            time_points=np.array([0.0, 30.0, 60.0], dtype=np.float64),
            counts=np.array([10.0, 20.0, 15.0], dtype=np.float64)
        )

        chrom2 = Chromatogram(
            time_points=np.array([300.0, 330.0, 360.0], dtype=np.float64),
            counts=np.array([10.0, 20.0, 15.0], dtype=np.float64)
        )

        assert chrom1.time_range[0] == 0.0
        assert chrom2.time_range[0] == 300.0
        assert chrom1.time_range[0] != chrom2.time_range[0]

    def test_signals_can_have_different_end_times(self):
        """Different chromatograms can end at different times."""
        chrom1 = Chromatogram(
            time_points=np.array([0.0, 30.0], dtype=np.float64),
            counts=np.array([10.0, 20.0], dtype=np.float64)
        )

        chrom2 = Chromatogram(
            time_points=np.array([0.0, 30.0, 60.0, 90.0], dtype=np.float64),
            counts=np.array([10.0, 20.0, 15.0, 12.0], dtype=np.float64)
        )

        assert chrom1.time_range[1] == 30.0
        assert chrom2.time_range[1] == 90.0
        assert chrom1.duration < chrom2.duration


class TestEdgeCases:
    """Test edge cases and boundary conditions."""

    def test_single_time_point_chromatogram(self):
        """Single time point is valid."""
        chrom = Chromatogram(
            time_points=np.array([42.0], dtype=np.float64),
            counts=np.array([100.0], dtype=np.float64)
        )

        assert len(chrom) == 1
        assert chrom.duration == 0.0
        assert chrom.time_range == (42.0, 42.0)

    def test_very_large_chromatogram(self):
        """Chromatogram with many points."""
        n_points = 10000
        time_points = np.linspace(0, 1000, n_points, dtype=np.float64)
        counts = np.random.rand(n_points).astype(np.float64) * 100

        chrom = Chromatogram(time_points=time_points, counts=counts)

        assert len(chrom) == n_points
        assert chrom.duration == 1000.0

    def test_zero_counts_valid(self):
        """Zero counts are valid."""
        chrom = Chromatogram(
            time_points=np.array([0.0, 1.0, 2.0], dtype=np.float64),
            counts=np.array([0.0, 0.0, 0.0], dtype=np.float64)
        )

        assert len(chrom) == 3
        np.testing.assert_array_equal(chrom.counts, np.zeros(3))

    def test_very_small_time_increments(self):
        """Very small time increments are valid."""
        chrom = Chromatogram(
            time_points=np.array([0.0, 0.001, 0.002, 0.003], dtype=np.float64),
            counts=np.array([10.0, 20.0, 15.0, 12.0], dtype=np.float64)
        )

        assert chrom.duration == 0.003
        assert len(chrom) == 4


class TestSignalVariantTypes:
    """Test different signal variant types (THEORY.md Section 5.0.1)."""

    def test_baseline_corrected_variant(self):
        """Baseline-corrected signal variant."""
        chrom = Chromatogram(
            time_points=np.array([0.0, 1.0, 2.0], dtype=np.float64),
            counts=np.array([15.0, 25.0, 20.0], dtype=np.float64)  # With baseline
        )

        # Subtract baseline of 5
        corrected = chrom.counts - 5.0
        chrom.add_signal_variant("corrected", corrected)

        assert chrom.has_signal_variant("corrected")
        np.testing.assert_array_equal(chrom.get_signal("corrected"), [10.0, 20.0, 15.0])

    def test_derivative_variant(self):
        """First derivative signal variant."""
        chrom = Chromatogram(
            time_points=np.array([0.0, 1.0, 2.0, 3.0], dtype=np.float64),
            counts=np.array([10.0, 15.0, 25.0, 20.0], dtype=np.float64)
        )

        # Simple forward difference derivative
        derivative = np.diff(chrom.counts, prepend=chrom.counts[0])
        chrom.add_signal_variant("derivative", derivative)

        assert chrom.has_signal_variant("derivative")
        assert len(chrom.get_signal("derivative")) == len(chrom)

    def test_second_derivative_variant(self):
        """Second derivative signal variant."""
        chrom = Chromatogram(
            time_points=np.array([0.0, 1.0, 2.0, 3.0, 4.0], dtype=np.float64),
            counts=np.array([10.0, 15.0, 25.0, 20.0, 18.0], dtype=np.float64)
        )

        # Compute second derivative (simplified)
        first_deriv = np.diff(chrom.counts, prepend=chrom.counts[0])
        second_deriv = np.diff(first_deriv, prepend=first_deriv[0])

        chrom.add_signal_variant("derivative_2", second_deriv)

        assert chrom.has_signal_variant("derivative_2")
        assert len(chrom.get_signal("derivative_2")) == len(chrom)

    def test_custom_processing_variant(self):
        """Custom processing can add arbitrary variants."""
        chrom = Chromatogram(
            time_points=np.array([0.0, 1.0, 2.0], dtype=np.float64),
            counts=np.array([10.0, 20.0, 15.0], dtype=np.float64)
        )

        # Custom smoothed signal
        smoothed = np.array([12.0, 17.5, 15.0], dtype=np.float64)
        chrom.add_signal_variant("smoothed", smoothed)

        # Custom normalized signal
        normalized = chrom.counts / np.max(chrom.counts)
        chrom.add_signal_variant("normalized", normalized)

        assert chrom.has_signal_variant("smoothed")
        assert chrom.has_signal_variant("normalized")
        assert len(chrom.signal_variants) == 2
