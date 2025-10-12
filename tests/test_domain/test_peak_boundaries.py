"""
Tests for PeakBoundaries value object.

Tests implementation against THEORY.md Section 5.2.4 specifications.
"""

import pytest
from lcseq.domain.value_objects.peak_boundaries import PeakBoundaries


class TestPeakBoundariesCreation:
    """Test PeakBoundaries instantiation and factory methods."""

    def test_create_with_valleys(self):
        """Create boundaries with valley information."""
        bounds = PeakBoundaries(
            left_base=44.0,
            right_base=46.0,
            left_valley=44.2,
            right_valley=45.8,
        )

        assert bounds.left_base == 44.0
        assert bounds.right_base == 46.0
        assert bounds.left_valley == 44.2
        assert bounds.right_valley == 45.8

    def test_create_without_valleys(self):
        """Create boundaries without valleys (threshold-based)."""
        bounds = PeakBoundaries(
            left_base=44.0,
            right_base=46.0,
            left_valley=None,
            right_valley=None,
        )

        assert bounds.left_base == 44.0
        assert bounds.right_base == 46.0
        assert bounds.left_valley is None
        assert bounds.right_valley is None

    def test_from_valleys_factory(self):
        """Factory method from valley positions."""
        bounds = PeakBoundaries.from_valleys(44.2, 45.8)

        assert bounds.left_base == 44.2
        assert bounds.right_base == 45.8
        assert bounds.left_valley == 44.2
        assert bounds.right_valley == 45.8

    def test_from_threshold_factory(self):
        """Factory method for threshold-based boundaries (THEORY.md Section 5.2.4)."""
        bounds = PeakBoundaries.from_threshold(44.0, 46.0)

        assert bounds.left_base == 44.0
        assert bounds.right_base == 46.0
        assert bounds.left_valley is None
        assert bounds.right_valley is None

    def test_partial_valleys_left_only(self):
        """Boundaries with only left valley."""
        bounds = PeakBoundaries(
            left_base=44.0,
            right_base=46.0,
            left_valley=44.2,
            right_valley=None,
        )

        assert bounds.left_valley == 44.2
        assert bounds.right_valley is None

    def test_partial_valleys_right_only(self):
        """Boundaries with only right valley."""
        bounds = PeakBoundaries(
            left_base=44.0,
            right_base=46.0,
            left_valley=None,
            right_valley=45.8,
        )

        assert bounds.left_valley is None
        assert bounds.right_valley == 45.8

    def test_immutable(self):
        """PeakBoundaries is immutable (frozen dataclass)."""
        bounds = PeakBoundaries.from_threshold(44.0, 46.0)

        with pytest.raises(AttributeError):
            bounds.left_base = 45.0


class TestPeakBoundariesValidation:
    """Test PeakBoundaries validation rules."""

    def test_left_must_be_less_than_right(self):
        """left_base must be < right_base."""
        with pytest.raises(ValueError, match="left_base.*must be < right_base"):
            PeakBoundaries(left_base=46.0, right_base=44.0)

    def test_equal_boundaries_invalid(self):
        """Equal left and right boundaries are invalid."""
        with pytest.raises(ValueError, match="left_base.*must be < right_base"):
            PeakBoundaries(left_base=45.0, right_base=45.0)

    def test_left_valley_must_be_within_boundaries(self):
        """left_valley must be within [left_base, right_base]."""
        # Left valley before left_base
        with pytest.raises(ValueError, match="left_valley.*must be within"):
            PeakBoundaries(
                left_base=44.0,
                right_base=46.0,
                left_valley=43.0,  # Before left_base
            )

        # Left valley after right_base
        with pytest.raises(ValueError, match="left_valley.*must be within"):
            PeakBoundaries(
                left_base=44.0,
                right_base=46.0,
                left_valley=47.0,  # After right_base
            )

    def test_right_valley_must_be_within_boundaries(self):
        """right_valley must be within [left_base, right_base]."""
        # Right valley before left_base
        with pytest.raises(ValueError, match="right_valley.*must be within"):
            PeakBoundaries(
                left_base=44.0,
                right_base=46.0,
                right_valley=43.0,
            )

        # Right valley after right_base
        with pytest.raises(ValueError, match="right_valley.*must be within"):
            PeakBoundaries(
                left_base=44.0,
                right_base=46.0,
                right_valley=47.0,
            )

    def test_left_valley_must_be_less_than_right_valley(self):
        """If both valleys present, left_valley < right_valley."""
        with pytest.raises(ValueError, match="left_valley.*must be <.*right_valley"):
            PeakBoundaries(
                left_base=44.0,
                right_base=46.0,
                left_valley=45.5,
                right_valley=44.5,
            )

    def test_equal_valleys_invalid(self):
        """Equal valley positions are invalid."""
        with pytest.raises(ValueError, match="left_valley.*must be <.*right_valley"):
            PeakBoundaries(
                left_base=44.0,
                right_base=46.0,
                left_valley=45.0,
                right_valley=45.0,
            )


class TestBoundaryProperties:
    """Test boundary property calculations."""

    def test_width(self):
        """Width = right_base - left_base."""
        bounds = PeakBoundaries(left_base=44.0, right_base=46.0)

        assert bounds.width() == 2.0

    def test_width_small(self):
        """Width calculation for narrow peak."""
        bounds = PeakBoundaries(left_base=45.0, right_base=45.5)

        assert bounds.width() == 0.5

    def test_width_large(self):
        """Width calculation for broad peak."""
        bounds = PeakBoundaries(left_base=40.0, right_base=50.0)

        assert bounds.width() == 10.0

    def test_valley_width_with_valleys(self):
        """Valley width when both valleys present."""
        bounds = PeakBoundaries.from_valleys(44.2, 45.8)

        assert abs(bounds.valley_width() - 1.6) < 1e-9  # Floating point comparison

    def test_valley_width_without_valleys(self):
        """Valley width is None when no valleys."""
        bounds = PeakBoundaries.from_threshold(44.0, 46.0)

        assert bounds.valley_width() is None

    def test_valley_width_partial_valleys(self):
        """Valley width is None when only one valley."""
        bounds = PeakBoundaries(
            left_base=44.0,
            right_base=46.0,
            left_valley=44.2,
        )

        assert bounds.valley_width() is None


class TestContainment:
    """Test time point containment checks."""

    def test_contains_within_boundaries(self):
        """Point within boundaries returns True."""
        bounds = PeakBoundaries(left_base=44.0, right_base=46.0)

        assert bounds.contains(45.0) is True

    def test_contains_at_left_boundary(self):
        """Point at left boundary is contained (inclusive)."""
        bounds = PeakBoundaries(left_base=44.0, right_base=46.0)

        assert bounds.contains(44.0) is True

    def test_contains_at_right_boundary(self):
        """Point at right boundary is contained (inclusive)."""
        bounds = PeakBoundaries(left_base=44.0, right_base=46.0)

        assert bounds.contains(46.0) is True

    def test_contains_before_boundaries(self):
        """Point before boundaries returns False."""
        bounds = PeakBoundaries(left_base=44.0, right_base=46.0)

        assert bounds.contains(43.0) is False

    def test_contains_after_boundaries(self):
        """Point after boundaries returns False."""
        bounds = PeakBoundaries(left_base=44.0, right_base=46.0)

        assert bounds.contains(47.0) is False

    def test_contains_just_outside_left(self):
        """Point just outside left boundary."""
        bounds = PeakBoundaries(left_base=44.0, right_base=46.0)

        assert bounds.contains(43.999) is False

    def test_contains_just_outside_right(self):
        """Point just outside right boundary."""
        bounds = PeakBoundaries(left_base=44.0, right_base=46.0)

        assert bounds.contains(46.001) is False


class TestValleyDetection:
    """Test valley detection status methods."""

    def test_has_valleys_both_present(self):
        """has_valleys() returns True when both valleys present."""
        bounds = PeakBoundaries.from_valleys(44.0, 46.0)

        assert bounds.has_valleys() is True

    def test_has_valleys_none_present(self):
        """has_valleys() returns False when no valleys."""
        bounds = PeakBoundaries.from_threshold(44.0, 46.0)

        assert bounds.has_valleys() is False

    def test_has_valleys_only_left(self):
        """has_valleys() returns False when only left valley."""
        bounds = PeakBoundaries(
            left_base=44.0,
            right_base=46.0,
            left_valley=44.2,
        )

        assert bounds.has_valleys() is False

    def test_has_valleys_only_right(self):
        """has_valleys() returns False when only right valley."""
        bounds = PeakBoundaries(
            left_base=44.0,
            right_base=46.0,
            right_valley=45.8,
        )

        assert bounds.has_valleys() is False

    def test_has_partial_valleys_left_only(self):
        """has_partial_valleys() returns True when only left valley."""
        bounds = PeakBoundaries(
            left_base=44.0,
            right_base=46.0,
            left_valley=44.2,
        )

        assert bounds.has_partial_valleys() is True

    def test_has_partial_valleys_right_only(self):
        """has_partial_valleys() returns True when only right valley."""
        bounds = PeakBoundaries(
            left_base=44.0,
            right_base=46.0,
            right_valley=45.8,
        )

        assert bounds.has_partial_valleys() is True

    def test_has_partial_valleys_both_present(self):
        """has_partial_valleys() returns False when both present."""
        bounds = PeakBoundaries.from_valleys(44.0, 46.0)

        assert bounds.has_partial_valleys() is False

    def test_has_partial_valleys_none_present(self):
        """has_partial_valleys() returns False when none present."""
        bounds = PeakBoundaries.from_threshold(44.0, 46.0)

        assert bounds.has_partial_valleys() is False


class TestEquality:
    """Test equality and hashing."""

    def test_equality_same_values(self):
        """Boundaries with same values are equal."""
        bounds1 = PeakBoundaries(left_base=44.0, right_base=46.0)
        bounds2 = PeakBoundaries(left_base=44.0, right_base=46.0)

        assert bounds1 == bounds2

    def test_equality_with_valleys(self):
        """Boundaries with same valleys are equal."""
        bounds1 = PeakBoundaries.from_valleys(44.2, 45.8)
        bounds2 = PeakBoundaries.from_valleys(44.2, 45.8)

        assert bounds1 == bounds2

    def test_inequality_different_left(self):
        """Different left_base → not equal."""
        bounds1 = PeakBoundaries(left_base=44.0, right_base=46.0)
        bounds2 = PeakBoundaries(left_base=44.5, right_base=46.0)

        assert bounds1 != bounds2

    def test_inequality_different_right(self):
        """Different right_base → not equal."""
        bounds1 = PeakBoundaries(left_base=44.0, right_base=46.0)
        bounds2 = PeakBoundaries(left_base=44.0, right_base=46.5)

        assert bounds1 != bounds2

    def test_inequality_different_valleys(self):
        """Different valleys → not equal."""
        bounds1 = PeakBoundaries.from_valleys(44.2, 45.8)
        bounds2 = PeakBoundaries.from_threshold(44.2, 45.8)

        assert bounds1 != bounds2

    def test_hashable(self):
        """PeakBoundaries is hashable (can use in sets/dicts)."""
        bounds1 = PeakBoundaries.from_threshold(44.0, 46.0)
        bounds2 = PeakBoundaries.from_valleys(44.2, 45.8)
        bounds3 = PeakBoundaries.from_threshold(44.0, 46.0)  # Same as bounds1

        boundaries = {bounds1, bounds2, bounds3}

        # bounds1 and bounds3 are same, so set has 2 elements
        assert len(boundaries) == 2


class TestStringRepresentation:
    """Test string conversion methods."""

    def test_str_without_valleys(self):
        """str() shows boundaries without valleys."""
        bounds = PeakBoundaries.from_threshold(44.0, 46.0)

        assert str(bounds) == "[44.00, 46.00]"

    def test_str_with_valleys(self):
        """str() shows boundaries with valley information."""
        bounds = PeakBoundaries.from_valleys(44.2, 45.8)

        result = str(bounds)
        assert "44.20, 45.80" in result
        assert "valleys:" in result

    def test_repr_shows_all_fields(self):
        """repr() shows all fields."""
        bounds = PeakBoundaries(
            left_base=44.0,
            right_base=46.0,
            left_valley=44.2,
            right_valley=45.8,
        )

        repr_str = repr(bounds)
        assert "left_base=44.0" in repr_str
        assert "right_base=46.0" in repr_str
        assert "left_valley=44.2" in repr_str
        assert "right_valley=45.8" in repr_str


class TestTheoryAlignment:
    """Test alignment with THEORY.md Section 5.2.4 specifications."""

    def test_valley_detection_scenario(self):
        """
        Test valley detection scenario from THEORY.md.

        "If reach valley (local minimum where f''(t) > 0):
         t_start = t_valley"
        """
        # Valleys detected on both sides
        bounds = PeakBoundaries.from_valleys(44.2, 45.8)

        assert bounds.left_base == 44.2  # Valley becomes base
        assert bounds.right_base == 45.8  # Valley becomes base
        assert bounds.has_valleys() is True

    def test_threshold_detection_scenario(self):
        """
        Test threshold detection scenario from THEORY.md.

        "If signal(t_i) < threshold_fraction × height_peak:
         t_start = t_i
         threshold_fraction = 0.05 (5% of peak height)"
        """
        # Threshold-based boundaries (no valleys detected)
        bounds = PeakBoundaries.from_threshold(44.0, 46.0)

        assert bounds.has_valleys() is False
        assert bounds.left_base == 44.0
        assert bounds.right_base == 46.0

    def test_boundary_at_signal_edge(self):
        """
        Test boundary at signal edge scenario from THEORY.md.

        "If reach signal start: t_start = t_min
         If reach signal end: t_end = t_max"
        """
        # Peak at signal boundary (e.g., only right valley detected)
        bounds = PeakBoundaries(
            left_base=0.0,  # Signal start
            right_base=46.0,
            left_valley=None,  # No left valley (at boundary)
            right_valley=45.8,
        )

        assert bounds.has_partial_valleys() is True
        assert bounds.left_valley is None


class TestEdgeCases:
    """Test edge cases and boundary conditions."""

    def test_very_narrow_peak(self):
        """Very narrow peak boundaries."""
        bounds = PeakBoundaries(left_base=45.0, right_base=45.01)

        assert abs(bounds.width() - 0.01) < 1e-9  # Floating point comparison

    def test_very_wide_peak(self):
        """Very wide peak boundaries."""
        bounds = PeakBoundaries(left_base=0.0, right_base=100.0)

        assert bounds.width() == 100.0

    def test_valleys_at_boundaries(self):
        """Valleys can be exactly at boundaries."""
        bounds = PeakBoundaries(
            left_base=44.0,
            right_base=46.0,
            left_valley=44.0,  # Exactly at left_base
            right_valley=46.0,  # Exactly at right_base
        )

        assert bounds.has_valleys() is True

    def test_fractional_times(self):
        """Fractional time values work correctly."""
        bounds = PeakBoundaries(
            left_base=44.123456,
            right_base=45.987654,
            left_valley=44.234567,
            right_valley=45.876543,
        )

        assert bounds.left_base == 44.123456
        assert bounds.width() == pytest.approx(45.987654 - 44.123456)

    def test_negative_times_allowed(self):
        """Negative times are allowed (unusual but valid)."""
        bounds = PeakBoundaries(left_base=-1.0, right_base=1.0)

        assert bounds.left_base == -1.0
        assert bounds.width() == 2.0

    def test_large_time_values(self):
        """Large time values work correctly."""
        bounds = PeakBoundaries(left_base=3000.0, right_base=3600.0)

        assert bounds.width() == 600.0
