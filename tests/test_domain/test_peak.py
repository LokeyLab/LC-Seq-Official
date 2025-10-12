"""
Tests for Peak entity.

Tests implementation against THEORY.md Section 2.1, 5.3, 6.10, 2.3.1 specifications.
"""

import pytest
from lcseq.domain.entities.peak import Peak, PeakType, ValidationStatus


class TestPeakTypeEnum:
    """Test PeakType enum values and behavior."""

    def test_peak_type_values(self):
        """PeakType enum has correct values."""
        assert PeakType.NULL.value == "NULL"
        assert PeakType.TRUNCATION.value == "TRUNCATION"
        assert PeakType.PUTATIVE_PRODUCT.value == "PUTATIVE_PRODUCT"
        assert PeakType.UNKNOWN.value == "UNKNOWN"

    def test_peak_type_membership(self):
        """All expected peak types are present."""
        peak_types = list(PeakType)
        assert len(peak_types) == 4
        assert PeakType.NULL in peak_types
        assert PeakType.TRUNCATION in peak_types
        assert PeakType.PUTATIVE_PRODUCT in peak_types
        assert PeakType.UNKNOWN in peak_types


class TestValidationStatusEnum:
    """Test ValidationStatus enum values and behavior."""

    def test_validation_status_values(self):
        """ValidationStatus enum has correct values."""
        assert ValidationStatus.VALIDATED.value == "VALIDATED"
        assert ValidationStatus.LIKELY_SUCCESS.value == "LIKELY_SUCCESS"
        assert ValidationStatus.UNCERTAIN.value == "UNCERTAIN"
        assert ValidationStatus.LIKELY_FAILURE.value == "LIKELY_FAILURE"
        assert ValidationStatus.FAILED.value == "FAILED"
        assert ValidationStatus.NOT_VALIDATED.value == "NOT_VALIDATED"

    def test_validation_status_membership(self):
        """All expected validation statuses are present."""
        statuses = list(ValidationStatus)
        assert len(statuses) == 6


class TestPeakCreation:
    """Test Peak instantiation and validation."""

    def test_create_basic_peak(self):
        """Basic peak with required fields."""
        peak = Peak(
            position=50.0,
            left_base=45.0,
            right_base=55.0,
            height=100.0,
            area=500.0
        )

        assert peak.position == 50.0
        assert peak.left_base == 45.0
        assert peak.right_base == 55.0
        assert peak.height == 100.0
        assert peak.area == 500.0
        assert peak.peak_type == PeakType.UNKNOWN
        assert peak.validation_status == ValidationStatus.NOT_VALIDATED

    def test_create_peak_with_classification(self):
        """Peak with peak_type classification."""
        peak = Peak(
            position=50.0,
            left_base=45.0,
            right_base=55.0,
            height=100.0,
            area=500.0,
            peak_type=PeakType.PUTATIVE_PRODUCT
        )

        assert peak.peak_type == PeakType.PUTATIVE_PRODUCT

    def test_create_peak_with_validation(self):
        """Peak with validation status."""
        peak = Peak(
            position=50.0,
            left_base=45.0,
            right_base=55.0,
            height=100.0,
            area=500.0,
            validation_status=ValidationStatus.VALIDATED
        )

        assert peak.validation_status == ValidationStatus.VALIDATED

    def test_create_peak_with_valleys(self):
        """Peak with valley positions."""
        peak = Peak(
            position=50.0,
            left_base=45.0,
            right_base=55.0,
            height=100.0,
            area=500.0,
            left_valley=46.0,
            right_valley=54.0
        )

        assert peak.left_valley == 46.0
        assert peak.right_valley == 54.0

    def test_create_peak_with_persistence(self):
        """Peak with topological persistence value."""
        peak = Peak(
            position=50.0,
            left_base=45.0,
            right_base=55.0,
            height=100.0,
            area=500.0,
            persistence=80.0
        )

        assert peak.persistence == 80.0


class TestPeakValidation:
    """Test Peak validation rules."""

    def test_negative_height_raises_error(self):
        """Peak height must be non-negative."""
        with pytest.raises(ValueError, match="non-negative"):
            Peak(
                position=50.0,
                left_base=45.0,
                right_base=55.0,
                height=-10.0,
                area=500.0
            )

    def test_negative_area_raises_error(self):
        """Peak area must be non-negative."""
        with pytest.raises(ValueError, match="non-negative"):
            Peak(
                position=50.0,
                left_base=45.0,
                right_base=55.0,
                height=100.0,
                area=-100.0
            )

    def test_left_base_greater_than_right_raises_error(self):
        """Left base must be less than right base."""
        with pytest.raises(ValueError, match="Left base.*must be.*right base"):
            Peak(
                position=50.0,
                left_base=60.0,  # Greater than right_base
                right_base=55.0,
                height=100.0,
                area=500.0
            )

    def test_left_base_equal_to_right_raises_error(self):
        """Left base cannot equal right base."""
        with pytest.raises(ValueError, match="Left base.*must be.*right base"):
            Peak(
                position=50.0,
                left_base=55.0,
                right_base=55.0,  # Equal
                height=100.0,
                area=500.0
            )

    def test_position_less_than_left_base_raises_error(self):
        """Peak position must be within boundaries."""
        with pytest.raises(ValueError, match="must be within boundaries"):
            Peak(
                position=40.0,  # Less than left_base
                left_base=45.0,
                right_base=55.0,
                height=100.0,
                area=500.0
            )

    def test_position_greater_than_right_base_raises_error(self):
        """Peak position must be within boundaries."""
        with pytest.raises(ValueError, match="must be within boundaries"):
            Peak(
                position=60.0,  # Greater than right_base
                left_base=45.0,
                right_base=55.0,
                height=100.0,
                area=500.0
            )

    def test_position_at_left_boundary_valid(self):
        """Position can be at left boundary."""
        peak = Peak(
            position=45.0,  # At left_base
            left_base=45.0,
            right_base=55.0,
            height=100.0,
            area=500.0
        )

        assert peak.position == 45.0

    def test_position_at_right_boundary_valid(self):
        """Position can be at right boundary."""
        peak = Peak(
            position=55.0,  # At right_base
            left_base=45.0,
            right_base=55.0,
            height=100.0,
            area=500.0
        )

        assert peak.position == 55.0

    def test_negative_persistence_raises_error(self):
        """Persistence must be non-negative if provided."""
        with pytest.raises(ValueError, match="non-negative"):
            Peak(
                position=50.0,
                left_base=45.0,
                right_base=55.0,
                height=100.0,
                area=500.0,
                persistence=-10.0
            )


class TestPeakProperties:
    """Test Peak computed properties."""

    def test_width_calculation(self):
        """Width is right_base - left_base."""
        peak = Peak(
            position=50.0,
            left_base=45.0,
            right_base=55.0,
            height=100.0,
            area=500.0
        )

        assert peak.width == 10.0

    def test_width_narrow_peak(self):
        """Width calculation for narrow peak."""
        peak = Peak(
            position=50.0,
            left_base=49.5,
            right_base=50.5,
            height=100.0,
            area=50.0
        )

        assert peak.width == 1.0

    def test_width_wide_peak(self):
        """Width calculation for wide peak."""
        peak = Peak(
            position=50.0,
            left_base=30.0,
            right_base=70.0,
            height=100.0,
            area=2000.0
        )

        assert peak.width == 40.0

    def test_is_product_peak_true(self):
        """is_product_peak returns True for PUTATIVE_PRODUCT."""
        peak = Peak(
            position=50.0,
            left_base=45.0,
            right_base=55.0,
            height=100.0,
            area=500.0,
            peak_type=PeakType.PUTATIVE_PRODUCT
        )

        assert peak.is_product_peak is True

    def test_is_product_peak_false(self):
        """is_product_peak returns False for non-PUTATIVE_PRODUCT types."""
        test_cases = [
            PeakType.NULL,
            PeakType.TRUNCATION,
            PeakType.UNKNOWN
        ]

        for peak_type in test_cases:
            peak = Peak(
                position=50.0,
                left_base=45.0,
                right_base=55.0,
                height=100.0,
                area=500.0,
                peak_type=peak_type
            )
            assert peak.is_product_peak is False, f"Failed for {peak_type}"

    def test_is_validated_true(self):
        """is_validated returns True for VALIDATED status."""
        peak = Peak(
            position=50.0,
            left_base=45.0,
            right_base=55.0,
            height=100.0,
            area=500.0,
            validation_status=ValidationStatus.VALIDATED
        )

        assert peak.is_validated is True

    def test_is_validated_false(self):
        """is_validated returns False for non-VALIDATED statuses."""
        test_cases = [
            ValidationStatus.NOT_VALIDATED,
            ValidationStatus.LIKELY_SUCCESS,
            ValidationStatus.UNCERTAIN,
            ValidationStatus.LIKELY_FAILURE,
            ValidationStatus.FAILED
        ]

        for status in test_cases:
            peak = Peak(
                position=50.0,
                left_base=45.0,
                right_base=55.0,
                height=100.0,
                area=500.0,
                validation_status=status
            )
            assert peak.is_validated is False, f"Failed for {status}"


class TestPeakClassification:
    """Test peak classification scenarios (THEORY.md Section 5.3)."""

    def test_null_peak_classification(self):
        """NULL peak at L₀ retention time."""
        peak = Peak(
            position=30.0,  # L₀ time
            left_base=28.0,
            right_base=32.0,
            height=50.0,
            area=200.0,
            peak_type=PeakType.NULL
        )

        assert peak.peak_type == PeakType.NULL
        assert peak.is_product_peak is False

    def test_truncation_peak_classification(self):
        """TRUNCATION peak at ancestor retention time."""
        peak = Peak(
            position=45.0,
            left_base=42.0,
            right_base=48.0,
            height=80.0,
            area=400.0,
            peak_type=PeakType.TRUNCATION
        )

        assert peak.peak_type == PeakType.TRUNCATION
        assert peak.is_product_peak is False

    def test_putative_product_classification(self):
        """PUTATIVE_PRODUCT positionally consistent with expected elution."""
        peak = Peak(
            position=60.0,
            left_base=55.0,
            right_base=65.0,
            height=120.0,
            area=800.0,
            peak_type=PeakType.PUTATIVE_PRODUCT
        )

        assert peak.peak_type == PeakType.PUTATIVE_PRODUCT
        assert peak.is_product_peak is True

    def test_unknown_peak_classification(self):
        """UNKNOWN peak that cannot be classified."""
        peak = Peak(
            position=90.0,
            left_base=85.0,
            right_base=95.0,
            height=40.0,
            area=200.0,
            peak_type=PeakType.UNKNOWN
        )

        assert peak.peak_type == PeakType.UNKNOWN
        assert peak.is_product_peak is False


class TestValidationStatuses:
    """Test validation status scenarios (THEORY.md Section 6.10)."""

    def test_validated_status(self):
        """VALIDATED: Synthesis succeeded with high confidence."""
        peak = Peak(
            position=60.0,
            left_base=55.0,
            right_base=65.0,
            height=120.0,
            area=800.0,
            peak_type=PeakType.PUTATIVE_PRODUCT,
            validation_status=ValidationStatus.VALIDATED
        )

        assert peak.validation_status == ValidationStatus.VALIDATED
        assert peak.is_validated is True

    def test_likely_success_status(self):
        """LIKELY_SUCCESS: Probably succeeded with moderate confidence."""
        peak = Peak(
            position=60.0,
            left_base=55.0,
            right_base=65.0,
            height=100.0,
            area=600.0,
            peak_type=PeakType.PUTATIVE_PRODUCT,
            validation_status=ValidationStatus.LIKELY_SUCCESS
        )

        assert peak.validation_status == ValidationStatus.LIKELY_SUCCESS
        assert peak.is_validated is False

    def test_uncertain_status(self):
        """UNCERTAIN: Ambiguous result."""
        peak = Peak(
            position=60.0,
            left_base=55.0,
            right_base=65.0,
            height=80.0,
            area=400.0,
            validation_status=ValidationStatus.UNCERTAIN
        )

        assert peak.validation_status == ValidationStatus.UNCERTAIN
        assert peak.is_validated is False

    def test_likely_failure_status(self):
        """LIKELY_FAILURE: Probably failed with moderate confidence."""
        peak = Peak(
            position=30.0,
            left_base=28.0,
            right_base=32.0,
            height=30.0,
            area=150.0,
            validation_status=ValidationStatus.LIKELY_FAILURE
        )

        assert peak.validation_status == ValidationStatus.LIKELY_FAILURE

    def test_failed_status(self):
        """FAILED: Synthesis failed with high confidence."""
        peak = Peak(
            position=30.0,
            left_base=28.0,
            right_base=32.0,
            height=20.0,
            area=100.0,
            peak_type=PeakType.NULL,
            validation_status=ValidationStatus.FAILED
        )

        assert peak.validation_status == ValidationStatus.FAILED

    def test_not_validated_default(self):
        """NOT_VALIDATED: Default status."""
        peak = Peak(
            position=50.0,
            left_base=45.0,
            right_base=55.0,
            height=100.0,
            area=500.0
        )

        assert peak.validation_status == ValidationStatus.NOT_VALIDATED


class TestClassificationVsValidation:
    """Test separation of classification and validation (THEORY.md Section 6.13)."""

    def test_putative_product_not_validated(self):
        """PUTATIVE_PRODUCT can be NOT_VALIDATED (positional hypothesis only)."""
        peak = Peak(
            position=60.0,
            left_base=55.0,
            right_base=65.0,
            height=120.0,
            area=800.0,
            peak_type=PeakType.PUTATIVE_PRODUCT,
            validation_status=ValidationStatus.NOT_VALIDATED
        )

        assert peak.is_product_peak is True
        assert peak.is_validated is False

    def test_putative_product_failed_validation(self):
        """PUTATIVE_PRODUCT can be FAILED (positional match but synthesis failed)."""
        peak = Peak(
            position=60.0,
            left_base=55.0,
            right_base=65.0,
            height=40.0,
            area=200.0,
            peak_type=PeakType.PUTATIVE_PRODUCT,
            validation_status=ValidationStatus.FAILED
        )

        assert peak.is_product_peak is True
        assert peak.is_validated is False

    def test_unknown_can_be_validated(self):
        """Even UNKNOWN peaks can have validation status."""
        peak = Peak(
            position=90.0,
            left_base=85.0,
            right_base=95.0,
            height=150.0,
            area=1000.0,
            peak_type=PeakType.UNKNOWN,
            validation_status=ValidationStatus.LIKELY_SUCCESS
        )

        assert peak.peak_type == PeakType.UNKNOWN
        assert peak.validation_status == ValidationStatus.LIKELY_SUCCESS


class TestAbsoluteTimeRepresentation:
    """Test absolute time representation (THEORY.md Section 2.3.1)."""

    def test_peak_position_absolute_time(self):
        """Position is absolute time (not array index)."""
        # Peak at 10 minutes (600 seconds)
        peak = Peak(
            position=600.0,
            left_base=580.0,
            right_base=620.0,
            height=100.0,
            area=800.0
        )

        assert peak.position == 600.0
        assert peak.width == 40.0

    def test_peak_at_late_retention_time(self):
        """Peak can be at any absolute time value."""
        peak = Peak(
            position=1800.0,  # 30 minutes
            left_base=1750.0,
            right_base=1850.0,
            height=80.0,
            area=600.0
        )

        assert peak.position == 1800.0


class TestPeakStringRepresentation:
    """Test string representation methods."""

    def test_repr_shows_key_info(self):
        """repr() shows position, height, area, type, and validation."""
        peak = Peak(
            position=50.25,
            left_base=45.0,
            right_base=55.0,
            height=123.4,
            area=678.9,
            peak_type=PeakType.PUTATIVE_PRODUCT,
            validation_status=ValidationStatus.VALIDATED
        )

        repr_str = repr(peak)
        assert "position=50.25" in repr_str
        assert "height=123.4" in repr_str
        assert "area=678.9" in repr_str
        assert "PUTATIVE_PRODUCT" in repr_str
        assert "VALIDATED" in repr_str


class TestEdgeCases:
    """Test edge cases and boundary conditions."""

    def test_zero_height_peak(self):
        """Zero height is valid (though unlikely in practice)."""
        peak = Peak(
            position=50.0,
            left_base=45.0,
            right_base=55.0,
            height=0.0,
            area=0.0
        )

        assert peak.height == 0.0
        assert peak.area == 0.0

    def test_zero_area_peak(self):
        """Zero area is valid."""
        peak = Peak(
            position=50.0,
            left_base=45.0,
            right_base=55.0,
            height=100.0,
            area=0.0
        )

        assert peak.area == 0.0

    def test_very_narrow_peak(self):
        """Very narrow peak is valid."""
        peak = Peak(
            position=50.0,
            left_base=49.99,
            right_base=50.01,
            height=100.0,
            area=1.0
        )

        assert peak.width == pytest.approx(0.02, abs=1e-10)

    def test_very_wide_peak(self):
        """Very wide peak is valid."""
        peak = Peak(
            position=500.0,
            left_base=100.0,
            right_base=900.0,
            height=50.0,
            area=20000.0
        )

        assert peak.width == 800.0

    def test_zero_persistence_valid(self):
        """Zero persistence is valid."""
        peak = Peak(
            position=50.0,
            left_base=45.0,
            right_base=55.0,
            height=100.0,
            area=500.0,
            persistence=0.0
        )

        assert peak.persistence == 0.0

    def test_none_persistence_valid(self):
        """None persistence is valid (not yet computed)."""
        peak = Peak(
            position=50.0,
            left_base=45.0,
            right_base=55.0,
            height=100.0,
            area=500.0,
            persistence=None
        )

        assert peak.persistence is None


class TestValleyPositions:
    """Test valley position handling."""

    def test_valleys_within_boundaries(self):
        """Valleys typically within peak boundaries."""
        peak = Peak(
            position=50.0,
            left_base=45.0,
            right_base=55.0,
            height=100.0,
            area=500.0,
            left_valley=46.0,
            right_valley=54.0
        )

        assert peak.left_valley > peak.left_base
        assert peak.right_valley < peak.right_base

    def test_valleys_at_boundaries(self):
        """Valleys can be at boundaries."""
        peak = Peak(
            position=50.0,
            left_base=45.0,
            right_base=55.0,
            height=100.0,
            area=500.0,
            left_valley=45.0,
            right_valley=55.0
        )

        assert peak.left_valley == peak.left_base
        assert peak.right_valley == peak.right_base

    def test_none_valleys_valid(self):
        """Valleys are optional."""
        peak = Peak(
            position=50.0,
            left_base=45.0,
            right_base=55.0,
            height=100.0,
            area=500.0,
            left_valley=None,
            right_valley=None
        )

        assert peak.left_valley is None
        assert peak.right_valley is None
