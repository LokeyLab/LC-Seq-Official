"""
Tests for RetentionTime value object.

Tests implementation against THEORY.md Section 2.3.1 specifications.
"""

import pytest
from lcseq.domain.value_objects.retention_time import RetentionTime, TimeUnit


class TestRetentionTimeCreation:
    """Test RetentionTime instantiation and factory methods."""

    def test_create_seconds(self):
        """Create retention time in seconds."""
        rt = RetentionTime(value=45.2, unit=TimeUnit.SECONDS)

        assert rt.value == 45.2
        assert rt.unit == TimeUnit.SECONDS

    def test_create_minutes(self):
        """Create retention time in minutes."""
        rt = RetentionTime(value=0.75, unit=TimeUnit.MINUTES)

        assert rt.value == 0.75
        assert rt.unit == TimeUnit.MINUTES

    def test_from_seconds_factory(self):
        """Factory method from seconds."""
        rt = RetentionTime.from_seconds(45.2)

        assert rt.value == 45.2
        assert rt.unit == TimeUnit.SECONDS

    def test_from_minutes_factory(self):
        """Factory method from minutes."""
        rt = RetentionTime.from_minutes(0.75)

        assert rt.value == 0.75
        assert rt.unit == TimeUnit.MINUTES

    def test_zero_retention_time(self):
        """Zero retention time is valid (injection time)."""
        rt = RetentionTime.from_seconds(0.0)

        assert rt.value == 0.0

    def test_immutable(self):
        """RetentionTime is immutable (frozen dataclass)."""
        rt = RetentionTime.from_seconds(45.0)

        with pytest.raises(AttributeError):
            rt.value = 50.0


class TestRetentionTimeValidation:
    """Test RetentionTime validation rules."""

    def test_negative_time_raises_error(self):
        """Retention time must be non-negative (THEORY.md Section 2.3.1)."""
        with pytest.raises(ValueError, match="non-negative"):
            RetentionTime.from_seconds(-10.0)

    def test_negative_minutes_raises_error(self):
        """Negative minutes also raises error."""
        with pytest.raises(ValueError, match="non-negative"):
            RetentionTime.from_minutes(-0.5)

    def test_unreasonably_large_time_raises_error(self):
        """Unreasonably large times raise error (>10 hours)."""
        # 11 hours in seconds
        with pytest.raises(ValueError, match="reasonable limit"):
            RetentionTime.from_seconds(39600.0)

    def test_unreasonably_large_minutes_raises_error(self):
        """Unreasonably large times in minutes also raise error."""
        # 11 hours in minutes
        with pytest.raises(ValueError, match="reasonable limit"):
            RetentionTime.from_minutes(660.0)

    def test_maximum_valid_time(self):
        """10 hours is maximum valid time."""
        # 10 hours = 600 minutes = 36000 seconds
        rt_s = RetentionTime.from_seconds(36000.0)
        rt_m = RetentionTime.from_minutes(600.0)

        assert rt_s.value == 36000.0
        assert rt_m.value == 600.0


class TestUnitConversion:
    """Test unit conversion methods."""

    def test_in_seconds_from_seconds(self):
        """in_seconds() when already in seconds."""
        rt = RetentionTime.from_seconds(45.2)

        assert rt.in_seconds() == 45.2

    def test_in_seconds_from_minutes(self):
        """in_seconds() converts from minutes."""
        rt = RetentionTime.from_minutes(1.5)

        assert rt.in_seconds() == 90.0

    def test_in_minutes_from_minutes(self):
        """in_minutes() when already in minutes."""
        rt = RetentionTime.from_minutes(1.5)

        assert rt.in_minutes() == 1.5

    def test_in_minutes_from_seconds(self):
        """in_minutes() converts from seconds."""
        rt = RetentionTime.from_seconds(90.0)

        assert rt.in_minutes() == 1.5

    def test_to_unit_seconds_to_minutes(self):
        """Convert from seconds to minutes."""
        rt = RetentionTime.from_seconds(90.0)
        rt_min = rt.to_unit(TimeUnit.MINUTES)

        assert rt_min.value == 1.5
        assert rt_min.unit == TimeUnit.MINUTES

    def test_to_unit_minutes_to_seconds(self):
        """Convert from minutes to seconds."""
        rt = RetentionTime.from_minutes(1.5)
        rt_sec = rt.to_unit(TimeUnit.SECONDS)

        assert rt_sec.value == 90.0
        assert rt_sec.unit == TimeUnit.SECONDS

    def test_to_unit_same_unit_returns_self(self):
        """Converting to same unit returns same instance."""
        rt = RetentionTime.from_seconds(45.0)
        rt_same = rt.to_unit(TimeUnit.SECONDS)

        assert rt_same is rt

    def test_conversion_precision(self):
        """Unit conversion maintains precision."""
        rt = RetentionTime.from_seconds(45.5)

        # Convert to minutes and back
        rt_min = rt.to_unit(TimeUnit.MINUTES)
        rt_sec = rt_min.to_unit(TimeUnit.SECONDS)

        assert abs(rt_sec.in_seconds() - 45.5) < 1e-9


class TestComparisons:
    """Test comparison operators (THEORY.md Section 2.3.1)."""

    def test_equality_same_unit(self):
        """Equal times in same unit."""
        rt1 = RetentionTime.from_seconds(45.0)
        rt2 = RetentionTime.from_seconds(45.0)

        assert rt1 == rt2

    def test_equality_with_non_retention_time(self):
        """Comparing with non-RetentionTime returns NotImplemented."""
        rt = RetentionTime.from_seconds(45.0)

        assert rt.__eq__(45.0) == NotImplemented
        assert rt.__eq__("45.0s") == NotImplemented

    def test_equality_different_units(self):
        """Equal times in different units (automatic conversion)."""
        rt1 = RetentionTime.from_seconds(45.0)
        rt2 = RetentionTime.from_minutes(0.75)

        assert rt1 == rt2

    def test_inequality(self):
        """Different times are not equal."""
        rt1 = RetentionTime.from_seconds(45.0)
        rt2 = RetentionTime.from_seconds(50.0)

        assert rt1 != rt2

    def test_less_than_same_unit(self):
        """Less than comparison in same unit."""
        rt1 = RetentionTime.from_seconds(45.0)
        rt2 = RetentionTime.from_seconds(50.0)

        assert rt1 < rt2
        assert not rt2 < rt1

    def test_less_than_different_units(self):
        """Less than comparison across units."""
        rt1 = RetentionTime.from_seconds(45.0)  # 45 seconds
        rt2 = RetentionTime.from_minutes(1.0)   # 60 seconds

        assert rt1 < rt2

    def test_less_than_or_equal(self):
        """Less than or equal comparison."""
        rt1 = RetentionTime.from_seconds(45.0)
        rt2 = RetentionTime.from_seconds(45.0)
        rt3 = RetentionTime.from_seconds(50.0)

        assert rt1 <= rt2  # Equal
        assert rt1 <= rt3  # Less than

    def test_greater_than(self):
        """Greater than comparison."""
        rt1 = RetentionTime.from_seconds(50.0)
        rt2 = RetentionTime.from_seconds(45.0)

        assert rt1 > rt2
        assert not rt2 > rt1

    def test_greater_than_or_equal(self):
        """Greater than or equal comparison."""
        rt1 = RetentionTime.from_seconds(50.0)
        rt2 = RetentionTime.from_seconds(50.0)
        rt3 = RetentionTime.from_seconds(45.0)

        assert rt1 >= rt2  # Equal
        assert rt1 >= rt3  # Greater than

    def test_sorting(self):
        """Can sort retention times."""
        times = [
            RetentionTime.from_seconds(50.0),
            RetentionTime.from_minutes(0.5),  # 30 seconds
            RetentionTime.from_seconds(45.0),
        ]

        sorted_times = sorted(times)

        assert sorted_times[0].in_seconds() == 30.0
        assert sorted_times[1].in_seconds() == 45.0
        assert sorted_times[2].in_seconds() == 50.0


class TestMatching:
    """Test retention time matching with tolerance (THEORY.md Section 2.3.1)."""

    def test_matches_within_tolerance(self):
        """Times match if difference <= tolerance."""
        rt1 = RetentionTime.from_seconds(45.0)
        rt2 = RetentionTime.from_seconds(45.2)

        # |45.2 - 45.0| = 0.2s <= 0.5s tolerance
        assert rt1.matches(rt2, tolerance=0.5)

    def test_matches_exact(self):
        """Exact match always within tolerance."""
        rt1 = RetentionTime.from_seconds(45.0)
        rt2 = RetentionTime.from_seconds(45.0)

        assert rt1.matches(rt2, tolerance=0.1)

    def test_does_not_match_outside_tolerance(self):
        """Times don't match if difference > tolerance."""
        rt1 = RetentionTime.from_seconds(45.0)
        rt2 = RetentionTime.from_seconds(45.2)

        # |45.2 - 45.0| = 0.2s > 0.1s tolerance
        assert not rt1.matches(rt2, tolerance=0.1)

    def test_matches_different_units(self):
        """Matching works across units."""
        rt1 = RetentionTime.from_seconds(45.0)
        rt2 = RetentionTime.from_minutes(0.75)  # 45 seconds

        assert rt1.matches(rt2, tolerance=0.1)

    def test_matches_symmetric(self):
        """Matching is symmetric."""
        rt1 = RetentionTime.from_seconds(45.0)
        rt2 = RetentionTime.from_seconds(45.2)

        assert rt1.matches(rt2, tolerance=0.5) == rt2.matches(rt1, tolerance=0.5)

    def test_matches_example_from_theory(self):
        """
        Test example from THEORY.md Section 2.3.1.

        "Does peak at 45.2s match expected 45.0s? → |45.2 - 45.0| = 0.2s"
        """
        expected = RetentionTime.from_seconds(45.0)
        observed = RetentionTime.from_seconds(45.2)

        diff = abs(observed.in_seconds() - expected.in_seconds())
        assert abs(diff - 0.2) < 1e-9  # Floating point comparison

        # Should match with 0.5s tolerance
        assert observed.matches(expected, tolerance=0.5)


class TestHashing:
    """Test hashing for use in sets/dicts."""

    def test_hashable(self):
        """RetentionTime is hashable."""
        rt1 = RetentionTime.from_seconds(45.0)
        rt2 = RetentionTime.from_minutes(0.75)  # Same time
        rt3 = RetentionTime.from_seconds(50.0)

        times = {rt1, rt2, rt3}

        # rt1 and rt2 are same time, so set has 2 elements
        assert len(times) == 2

    def test_hash_same_for_equal_times(self):
        """Equal times have same hash."""
        rt1 = RetentionTime.from_seconds(45.0)
        rt2 = RetentionTime.from_minutes(0.75)

        assert hash(rt1) == hash(rt2)

    def test_can_use_as_dict_key(self):
        """Can use as dictionary key."""
        rt1 = RetentionTime.from_seconds(45.0)
        rt2 = RetentionTime.from_seconds(50.0)

        peak_data = {
            rt1: "peak_1",
            rt2: "peak_2",
        }

        assert peak_data[rt1] == "peak_1"


class TestStringRepresentation:
    """Test string conversion methods."""

    def test_str_seconds(self):
        """str() shows value with unit (seconds)."""
        rt = RetentionTime.from_seconds(45.2)

        assert str(rt) == "45.2s"

    def test_str_minutes(self):
        """str() shows value with unit (minutes)."""
        rt = RetentionTime.from_minutes(1.5)

        assert str(rt) == "1.5min"

    def test_repr_seconds(self):
        """repr() shows detailed representation."""
        rt = RetentionTime.from_seconds(45.2)

        repr_str = repr(rt)
        assert "45.2" in repr_str
        assert "SECONDS" in repr_str

    def test_repr_minutes(self):
        """repr() shows unit name for minutes."""
        rt = RetentionTime.from_minutes(1.5)

        repr_str = repr(rt)
        assert "1.5" in repr_str
        assert "MINUTES" in repr_str


class TestPhysicalMeaning:
    """Test that retention times are physically meaningful (THEORY.md Section 2.3.1)."""

    def test_absolute_time_not_index(self):
        """
        Retention time is absolute time, not array index.

        THEORY.md: "Retention times are physically meaningful (not array indices)"
        """
        rt = RetentionTime.from_seconds(45.2)

        # This is a time value, not an index
        assert rt.in_seconds() == 45.2
        assert isinstance(rt.in_seconds(), float)

    def test_scalar_arithmetic(self):
        """
        Peak position comparisons use direct scalar arithmetic.

        THEORY.md: "Peak position comparisons: direct scalar arithmetic"
        """
        rt1 = RetentionTime.from_seconds(45.2)
        rt2 = RetentionTime.from_seconds(45.0)

        # Direct arithmetic on scalars
        diff = abs(rt1.in_seconds() - rt2.in_seconds())
        assert abs(diff - 0.2) < 1e-9  # Floating point comparison


class TestEdgeCases:
    """Test edge cases and boundary conditions."""

    def test_very_small_time(self):
        """Very small but positive time is valid."""
        rt = RetentionTime.from_seconds(0.001)

        assert rt.value == 0.001

    def test_fractional_seconds(self):
        """Fractional seconds work correctly."""
        rt = RetentionTime.from_seconds(45.123456)

        assert rt.in_seconds() == 45.123456

    def test_fractional_minutes(self):
        """Fractional minutes work correctly."""
        rt = RetentionTime.from_minutes(1.123456)

        assert rt.in_minutes() == 1.123456

    def test_conversion_edge_case(self):
        """Edge case: very small minutes converts correctly."""
        rt = RetentionTime.from_minutes(0.001)  # 0.06 seconds

        assert abs(rt.in_seconds() - 0.06) < 1e-9

    def test_large_valid_time(self):
        """Large but valid time (near 10 hour limit)."""
        rt = RetentionTime.from_seconds(35999.0)  # Just under 10 hours

        assert rt.in_seconds() == 35999.0

    def test_typical_lc_range(self):
        """Typical LC retention times (0-60 minutes)."""
        rt_start = RetentionTime.from_minutes(0.0)
        rt_end = RetentionTime.from_minutes(60.0)

        assert rt_start.in_seconds() == 0.0
        assert rt_end.in_seconds() == 3600.0
