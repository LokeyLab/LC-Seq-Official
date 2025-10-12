"""
Comprehensive tests for PeakClassification model.

Tests peak classification result structure.
"""

import pytest
import numpy as np
from lcseq.domain.models.peak_classification import PeakClassification
from lcseq.domain.entities.compound import Compound
from lcseq.domain.entities.peak import Peak, PeakType
from lcseq.domain.entities.building_block import BuildingBlock
from lcseq.domain.entities.chromatogram import Chromatogram


@pytest.fixture
def sample_chromatogram():
    """Create sample chromatogram for testing."""
    return Chromatogram(
        time_points=np.array([1.0, 2.0, 3.0, 4.0, 5.0]),
        counts=np.array([100.0, 200.0, 300.0, 200.0, 100.0]),
    )


@pytest.fixture
def sample_compound(sample_chromatogram):
    """Create sample compound for testing."""
    bb0 = BuildingBlock.from_code(0, "Pro")
    bb1 = BuildingBlock.from_code(1, "Leu")
    return Compound([bb0, bb1], sample_chromatogram)


@pytest.fixture
def sample_peak():
    """Create sample peak for testing."""
    return Peak(
        position=2.5,
        left_base=2.0,
        right_base=3.0,
        height=300.0,
        area=450.0,
        peak_type=PeakType.PUTATIVE_PRODUCT,
    )


class TestPeakClassificationCreation:
    """Test peak classification creation."""

    def test_create_classification_without_confidence(self, sample_compound, sample_peak):
        """Test creating classification without confidence score."""
        classification = PeakClassification(
            compound=sample_compound,
            peak=sample_peak,
            classification=PeakType.PUTATIVE_PRODUCT,
        )

        assert classification.compound == sample_compound
        assert classification.peak == sample_peak
        assert classification.classification == PeakType.PUTATIVE_PRODUCT
        assert classification.confidence is None

    def test_create_classification_with_confidence(self, sample_compound, sample_peak):
        """Test creating classification with confidence score."""
        classification = PeakClassification(
            compound=sample_compound,
            peak=sample_peak,
            classification=PeakType.PUTATIVE_PRODUCT,
            confidence=0.95,
        )

        assert classification.confidence == 0.95

    def test_create_truncation_classification(self, sample_compound, sample_peak):
        """Test creating truncation classification."""
        peak = Peak(
            position=1.8,
            left_base=1.5,
            right_base=2.1,
            height=150.0,
            area=200.0,
            peak_type=PeakType.TRUNCATION,
        )

        classification = PeakClassification(
            compound=sample_compound,
            peak=peak,
            classification=PeakType.TRUNCATION,
            confidence=0.88,
        )

        assert classification.classification == PeakType.TRUNCATION
        assert classification.is_truncation
        assert not classification.is_product

    def test_create_null_classification(self, sample_compound, sample_peak):
        """Test creating NULL classification."""
        peak = Peak(
            position=1.2,
            left_base=1.0,
            right_base=1.4,
            height=80.0,
            area=100.0,
            peak_type=PeakType.NULL,
        )

        classification = PeakClassification(
            compound=sample_compound,
            peak=peak,
            classification=PeakType.NULL,
        )

        assert classification.classification == PeakType.NULL
        assert classification.is_null
        assert not classification.is_product

    def test_create_unknown_classification(self, sample_compound, sample_peak):
        """Test creating UNKNOWN classification."""
        peak = Peak(
            position=4.5,
            left_base=4.0,
            right_base=5.0,
            height=120.0,
            area=150.0,
            peak_type=PeakType.UNKNOWN,
        )

        classification = PeakClassification(
            compound=sample_compound,
            peak=peak,
            classification=PeakType.UNKNOWN,
        )

        assert classification.classification == PeakType.UNKNOWN
        assert classification.is_unknown
        assert not classification.is_product


class TestConfidenceValidation:
    """Test confidence validation."""

    def test_confidence_in_valid_range(self, sample_compound, sample_peak):
        """Test confidence values in valid range [0, 1]."""
        # Test boundary values
        c0 = PeakClassification(
            compound=sample_compound,
            peak=sample_peak,
            classification=PeakType.PUTATIVE_PRODUCT,
            confidence=0.0,
        )
        assert c0.confidence == 0.0

        c1 = PeakClassification(
            compound=sample_compound,
            peak=sample_peak,
            classification=PeakType.PUTATIVE_PRODUCT,
            confidence=1.0,
        )
        assert c1.confidence == 1.0

        # Test mid-range
        c_mid = PeakClassification(
            compound=sample_compound,
            peak=sample_peak,
            classification=PeakType.PUTATIVE_PRODUCT,
            confidence=0.75,
        )
        assert c_mid.confidence == 0.75

    def test_confidence_below_zero_raises_error(self, sample_compound, sample_peak):
        """Test confidence below 0 raises error."""
        with pytest.raises(ValueError, match="Confidence must be in"):
            PeakClassification(
                compound=sample_compound,
                peak=sample_peak,
                classification=PeakType.PUTATIVE_PRODUCT,
                confidence=-0.1,
            )

    def test_confidence_above_one_raises_error(self, sample_compound, sample_peak):
        """Test confidence above 1 raises error."""
        with pytest.raises(ValueError, match="Confidence must be in"):
            PeakClassification(
                compound=sample_compound,
                peak=sample_peak,
                classification=PeakType.PUTATIVE_PRODUCT,
                confidence=1.1,
            )


class TestClassificationTypeChecks:
    """Test classification type boolean properties."""

    def test_is_product_true(self, sample_compound, sample_peak):
        """Test is_product returns True for PUTATIVE_PRODUCT."""
        classification = PeakClassification(
            compound=sample_compound,
            peak=sample_peak,
            classification=PeakType.PUTATIVE_PRODUCT,
        )

        assert classification.is_product
        assert not classification.is_truncation
        assert not classification.is_null
        assert not classification.is_unknown

    def test_is_truncation_true(self, sample_compound, sample_peak):
        """Test is_truncation returns True for TRUNCATION."""
        classification = PeakClassification(
            compound=sample_compound,
            peak=sample_peak,
            classification=PeakType.TRUNCATION,
        )

        assert classification.is_truncation
        assert not classification.is_product
        assert not classification.is_null
        assert not classification.is_unknown

    def test_is_null_true(self, sample_compound, sample_peak):
        """Test is_null returns True for NULL."""
        classification = PeakClassification(
            compound=sample_compound,
            peak=sample_peak,
            classification=PeakType.NULL,
        )

        assert classification.is_null
        assert not classification.is_product
        assert not classification.is_truncation
        assert not classification.is_unknown

    def test_is_unknown_true(self, sample_compound, sample_peak):
        """Test is_unknown returns True for UNKNOWN."""
        classification = PeakClassification(
            compound=sample_compound,
            peak=sample_peak,
            classification=PeakType.UNKNOWN,
        )

        assert classification.is_unknown
        assert not classification.is_product
        assert not classification.is_truncation
        assert not classification.is_null


class TestImmutability:
    """Test that PeakClassification is immutable (frozen)."""

    def test_cannot_modify_compound(self, sample_compound, sample_peak):
        """Test cannot modify compound after creation."""
        classification = PeakClassification(
            compound=sample_compound,
            peak=sample_peak,
            classification=PeakType.PUTATIVE_PRODUCT,
        )

        with pytest.raises(Exception):  # FrozenInstanceError
            classification.compound = None

    def test_cannot_modify_peak(self, sample_compound, sample_peak):
        """Test cannot modify peak after creation."""
        classification = PeakClassification(
            compound=sample_compound,
            peak=sample_peak,
            classification=PeakType.PUTATIVE_PRODUCT,
        )

        with pytest.raises(Exception):  # FrozenInstanceError
            classification.peak = None

    def test_cannot_modify_classification(self, sample_compound, sample_peak):
        """Test cannot modify classification after creation."""
        classification = PeakClassification(
            compound=sample_compound,
            peak=sample_peak,
            classification=PeakType.PUTATIVE_PRODUCT,
        )

        with pytest.raises(Exception):  # FrozenInstanceError
            classification.classification = PeakType.TRUNCATION

    def test_cannot_modify_confidence(self, sample_compound, sample_peak):
        """Test cannot modify confidence after creation."""
        classification = PeakClassification(
            compound=sample_compound,
            peak=sample_peak,
            classification=PeakType.PUTATIVE_PRODUCT,
            confidence=0.95,
        )

        with pytest.raises(Exception):  # FrozenInstanceError
            classification.confidence = 0.5


class TestStringRepresentations:
    """Test string representations."""

    def test_repr_without_confidence(self, sample_compound, sample_peak):
        """Test repr without confidence."""
        classification = PeakClassification(
            compound=sample_compound,
            peak=sample_peak,
            classification=PeakType.PUTATIVE_PRODUCT,
        )

        repr_str = repr(classification)
        assert "PeakClassification" in repr_str
        assert "Leu-Pro" in repr_str
        assert "PUTATIVE_PRODUCT" in repr_str
        assert "2.50" in repr_str  # peak position

    def test_repr_with_confidence(self, sample_compound, sample_peak):
        """Test repr with confidence."""
        classification = PeakClassification(
            compound=sample_compound,
            peak=sample_peak,
            classification=PeakType.PUTATIVE_PRODUCT,
            confidence=0.95,
        )

        repr_str = repr(classification)
        assert "confidence=0.950" in repr_str

    def test_str_without_confidence(self, sample_compound, sample_peak):
        """Test str representation without confidence."""
        classification = PeakClassification(
            compound=sample_compound,
            peak=sample_peak,
            classification=PeakType.PUTATIVE_PRODUCT,
        )

        str_repr = str(classification)
        assert "Leu-Pro" in str_repr
        assert "PUTATIVE_PRODUCT" in str_repr
        assert "2.50" in str_repr

    def test_str_with_confidence(self, sample_compound, sample_peak):
        """Test str representation with confidence."""
        classification = PeakClassification(
            compound=sample_compound,
            peak=sample_peak,
            classification=PeakType.TRUNCATION,
            confidence=0.88,
        )

        str_repr = str(classification)
        assert "conf=0.88" in str_repr


class TestDifferentCompoundsAndPeaks:
    """Test classification with various compounds and peaks."""

    def test_classification_preserves_compound_info(self, sample_chromatogram):
        """Test that compound information is preserved."""
        bb0 = BuildingBlock.from_code(0, "Ala")
        bb1 = BuildingBlock.from_code(1, "Val")
        bb2 = BuildingBlock.from_code(2, "Pro")
        compound = Compound([bb0, bb1, bb2], sample_chromatogram)

        peak = Peak(
            position=3.0,
            left_base=2.5,
            right_base=3.5,
            height=400.0,
            area=600.0,
        )

        classification = PeakClassification(
            compound=compound,
            peak=peak,
            classification=PeakType.PUTATIVE_PRODUCT,
        )

        assert classification.compound.positional_sequence == "Pro-Val-Ala"
        assert classification.compound.level == 3

    def test_classification_preserves_peak_info(self, sample_compound):
        """Test that peak information is preserved."""
        peak = Peak(
            position=2.75,
            left_base=2.25,
            right_base=3.25,
            height=250.0,
            area=375.0,
            persistence=0.15,
        )

        classification = PeakClassification(
            compound=sample_compound,
            peak=peak,
            classification=PeakType.PUTATIVE_PRODUCT,
        )

        assert classification.peak.position == 2.75
        assert classification.peak.height == 250.0
        assert classification.peak.area == 375.0
        assert classification.peak.persistence == 0.15


class TestAllPeakTypes:
    """Test all peak type classifications."""

    def test_all_peak_types_representable(self, sample_compound, sample_peak):
        """Test that all PeakType values can be classified."""
        for peak_type in PeakType:
            classification = PeakClassification(
                compound=sample_compound,
                peak=sample_peak,
                classification=peak_type,
            )

            assert classification.classification == peak_type
