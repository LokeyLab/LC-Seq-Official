"""
Comprehensive tests for AnalysisResult model.

Tests complete analysis output structure.
"""

import pytest
import numpy as np
from lcseq.domain.models.analysis_result import AnalysisResult
from lcseq.domain.models.compound_hierarchy import CompoundHierarchy, HierarchyMode
from lcseq.domain.models.equivalence_class import EquivalenceClass
from lcseq.domain.models.peak_classification import PeakClassification
from lcseq.domain.entities.compound import Compound
from lcseq.domain.entities.peak import Peak, PeakType, ValidationStatus
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
def sample_hierarchy():
    """Create sample hierarchy for testing."""
    return CompoundHierarchy(mode=HierarchyMode.BUILDING_BLOCK)


@pytest.fixture
def sample_compound(sample_chromatogram):
    """Create sample compound for testing."""
    bb0 = BuildingBlock.from_code(0, "Pro")
    bb1 = BuildingBlock.from_code(1, "Leu")
    return Compound([bb0, bb1], sample_chromatogram)


class TestAnalysisResultCreation:
    """Test analysis result creation."""

    def test_create_empty_result(self, sample_hierarchy):
        """Test creating empty analysis result."""
        result = AnalysisResult(hierarchy=sample_hierarchy)

        assert result.hierarchy == sample_hierarchy
        assert len(result.equivalence_classes) == 0
        assert len(result.peak_classifications) == 0
        assert len(result.validation_results) == 0
        assert len(result.metadata) == 0

    def test_create_with_equivalence_classes(self, sample_hierarchy):
        """Test creating with equivalence classes."""
        eq_class = EquivalenceClass(residue_sequence="Val-Pro")

        result = AnalysisResult(
            hierarchy=sample_hierarchy, equivalence_classes=[eq_class]
        )

        assert len(result.equivalence_classes) == 1
        assert eq_class in result.equivalence_classes

    def test_create_with_metadata(self, sample_hierarchy):
        """Test creating with metadata."""
        metadata = {"timestamp": "2025-01-01", "version": "1.0"}

        result = AnalysisResult(hierarchy=sample_hierarchy, metadata=metadata)

        assert result.metadata["timestamp"] == "2025-01-01"
        assert result.metadata["version"] == "1.0"


class TestAddPeakClassification:
    """Test adding peak classifications."""

    def test_add_single_classification(self, sample_hierarchy, sample_compound, sample_chromatogram):
        """Test adding single peak classification."""
        result = AnalysisResult(hierarchy=sample_hierarchy)

        peak = Peak(
            position=2.5,
            left_base=2.0,
            right_base=3.0,
            height=300.0,
            area=450.0,
        )
        classification = PeakClassification(
            compound=sample_compound,
            peak=peak,
            classification=PeakType.PUTATIVE_PRODUCT,
        )

        result.add_peak_classification(sample_compound, classification)

        assert sample_compound in result.peak_classifications
        assert len(result.peak_classifications[sample_compound]) == 1
        assert classification in result.peak_classifications[sample_compound]

    def test_add_multiple_classifications_same_compound(
        self, sample_hierarchy, sample_compound, sample_chromatogram
    ):
        """Test adding multiple classifications for same compound."""
        result = AnalysisResult(hierarchy=sample_hierarchy)

        # Product peak
        peak1 = Peak(
            position=2.5,
            left_base=2.0,
            right_base=3.0,
            height=300.0,
            area=450.0,
        )
        classification1 = PeakClassification(
            compound=sample_compound,
            peak=peak1,
            classification=PeakType.PUTATIVE_PRODUCT,
        )

        # Truncation peak
        peak2 = Peak(
            position=1.8,
            left_base=1.5,
            right_base=2.1,
            height=150.0,
            area=200.0,
        )
        classification2 = PeakClassification(
            compound=sample_compound,
            peak=peak2,
            classification=PeakType.TRUNCATION,
        )

        result.add_peak_classification(sample_compound, classification1)
        result.add_peak_classification(sample_compound, classification2)

        assert len(result.peak_classifications[sample_compound]) == 2


class TestSetValidationStatus:
    """Test setting validation status."""

    def test_set_validation_status(self, sample_hierarchy, sample_compound):
        """Test setting validation status for compound."""
        result = AnalysisResult(hierarchy=sample_hierarchy)

        result.set_validation_status(sample_compound, ValidationStatus.VALIDATED)

        assert sample_compound in result.validation_results
        assert result.validation_results[sample_compound] == ValidationStatus.VALIDATED

    def test_set_validation_status_overwrites_previous(
        self, sample_hierarchy, sample_compound
    ):
        """Test that setting validation status overwrites previous."""
        result = AnalysisResult(hierarchy=sample_hierarchy)

        result.set_validation_status(sample_compound, ValidationStatus.UNCERTAIN)
        result.set_validation_status(sample_compound, ValidationStatus.VALIDATED)

        assert result.validation_results[sample_compound] == ValidationStatus.VALIDATED


class TestGetProductPeak:
    """Test getting product peak classification."""

    def test_get_product_peak_single_product(
        self, sample_hierarchy, sample_compound, sample_chromatogram
    ):
        """Test getting product peak when one exists."""
        result = AnalysisResult(hierarchy=sample_hierarchy)

        peak = Peak(
            position=2.5,
            left_base=2.0,
            right_base=3.0,
            height=300.0,
            area=450.0,
        )
        classification = PeakClassification(
            compound=sample_compound,
            peak=peak,
            classification=PeakType.PUTATIVE_PRODUCT,
        )

        result.add_peak_classification(sample_compound, classification)

        product = result.get_product_peak(sample_compound)
        assert product == classification

    def test_get_product_peak_with_truncations(
        self, sample_hierarchy, sample_compound, sample_chromatogram
    ):
        """Test getting product peak when truncations also exist."""
        result = AnalysisResult(hierarchy=sample_hierarchy)

        # Add truncation
        peak1 = Peak(
            position=1.8,
            left_base=1.5,
            right_base=2.1,
            height=150.0,
            area=200.0,
        )
        truncation = PeakClassification(
            compound=sample_compound,
            peak=peak1,
            classification=PeakType.TRUNCATION,
        )
        result.add_peak_classification(sample_compound, truncation)

        # Add product
        peak2 = Peak(
            position=2.5,
            left_base=2.0,
            right_base=3.0,
            height=300.0,
            area=450.0,
        )
        product_class = PeakClassification(
            compound=sample_compound,
            peak=peak2,
            classification=PeakType.PUTATIVE_PRODUCT,
        )
        result.add_peak_classification(sample_compound, product_class)

        product = result.get_product_peak(sample_compound)
        assert product == product_class

    def test_get_product_peak_no_product(self, sample_hierarchy, sample_compound):
        """Test getting product peak when none exists."""
        result = AnalysisResult(hierarchy=sample_hierarchy)

        product = result.get_product_peak(sample_compound)
        assert product is None

    def test_get_product_peak_compound_not_analyzed(self, sample_hierarchy, sample_compound):
        """Test getting product peak for compound not analyzed."""
        result = AnalysisResult(hierarchy=sample_hierarchy)

        product = result.get_product_peak(sample_compound)
        assert product is None


class TestGetTruncationPeaks:
    """Test getting truncation peak classifications."""

    def test_get_truncation_peaks_single_truncation(
        self, sample_hierarchy, sample_compound, sample_chromatogram
    ):
        """Test getting truncation peaks when one exists."""
        result = AnalysisResult(hierarchy=sample_hierarchy)

        peak = Peak(
            position=1.8,
            left_base=1.5,
            right_base=2.1,
            height=150.0,
            area=200.0,
        )
        classification = PeakClassification(
            compound=sample_compound,
            peak=peak,
            classification=PeakType.TRUNCATION,
        )

        result.add_peak_classification(sample_compound, classification)

        truncations = result.get_truncation_peaks(sample_compound)
        assert len(truncations) == 1
        assert classification in truncations

    def test_get_truncation_peaks_multiple_truncations(
        self, sample_hierarchy, sample_compound, sample_chromatogram
    ):
        """Test getting multiple truncation peaks."""
        result = AnalysisResult(hierarchy=sample_hierarchy)

        # Add two truncations
        peak1 = Peak(
            position=1.5,
            left_base=1.2,
            right_base=1.8,
            height=100.0,
            area=150.0,
        )
        trunc1 = PeakClassification(
            compound=sample_compound,
            peak=peak1,
            classification=PeakType.TRUNCATION,
        )

        peak2 = Peak(
            position=2.0,
            left_base=1.8,
            right_base=2.2,
            height=120.0,
            area=180.0,
        )
        trunc2 = PeakClassification(
            compound=sample_compound,
            peak=peak2,
            classification=PeakType.TRUNCATION,
        )

        result.add_peak_classification(sample_compound, trunc1)
        result.add_peak_classification(sample_compound, trunc2)

        truncations = result.get_truncation_peaks(sample_compound)
        assert len(truncations) == 2
        assert trunc1 in truncations
        assert trunc2 in truncations

    def test_get_truncation_peaks_no_truncations(
        self, sample_hierarchy, sample_compound, sample_chromatogram
    ):
        """Test getting truncation peaks when none exist."""
        result = AnalysisResult(hierarchy=sample_hierarchy)

        # Add only product peak
        peak = Peak(
            position=2.5,
            left_base=2.0,
            right_base=3.0,
            height=300.0,
            area=450.0,
        )
        product = PeakClassification(
            compound=sample_compound,
            peak=peak,
            classification=PeakType.PUTATIVE_PRODUCT,
        )
        result.add_peak_classification(sample_compound, product)

        truncations = result.get_truncation_peaks(sample_compound)
        assert len(truncations) == 0


class TestGetValidatedCompounds:
    """Test getting validated compounds."""

    def test_get_validated_compounds(
        self, sample_hierarchy, sample_compound, sample_chromatogram
    ):
        """Test getting validated compounds."""
        result = AnalysisResult(hierarchy=sample_hierarchy)

        bb0 = BuildingBlock.from_code(0, "Val")
        bb1 = BuildingBlock.from_code(1, "Pro")
        compound2 = Compound([bb0, bb1], sample_chromatogram)

        result.set_validation_status(sample_compound, ValidationStatus.VALIDATED)
        result.set_validation_status(compound2, ValidationStatus.FAILED)

        validated = result.get_validated_compounds()
        assert len(validated) == 1
        assert sample_compound in validated
        assert compound2 not in validated

    def test_get_validated_compounds_empty(self, sample_hierarchy, sample_compound):
        """Test getting validated compounds when none exist."""
        result = AnalysisResult(hierarchy=sample_hierarchy)

        result.set_validation_status(sample_compound, ValidationStatus.FAILED)

        validated = result.get_validated_compounds()
        assert len(validated) == 0


class TestGetFailedCompounds:
    """Test getting failed compounds."""

    def test_get_failed_compounds(
        self, sample_hierarchy, sample_compound, sample_chromatogram
    ):
        """Test getting failed compounds."""
        result = AnalysisResult(hierarchy=sample_hierarchy)

        bb0 = BuildingBlock.from_code(0, "Val")
        bb1 = BuildingBlock.from_code(1, "Pro")
        compound2 = Compound([bb0, bb1], sample_chromatogram)

        result.set_validation_status(sample_compound, ValidationStatus.VALIDATED)
        result.set_validation_status(compound2, ValidationStatus.FAILED)

        failed = result.get_failed_compounds()
        assert len(failed) == 1
        assert compound2 in failed
        assert sample_compound not in failed

    def test_get_failed_compounds_empty(self, sample_hierarchy, sample_compound):
        """Test getting failed compounds when none exist."""
        result = AnalysisResult(hierarchy=sample_hierarchy)

        result.set_validation_status(sample_compound, ValidationStatus.VALIDATED)

        failed = result.get_failed_compounds()
        assert len(failed) == 0


class TestGetValidationSummary:
    """Test getting validation summary."""

    def test_get_validation_summary(
        self, sample_hierarchy, sample_compound, sample_chromatogram
    ):
        """Test getting validation summary."""
        result = AnalysisResult(hierarchy=sample_hierarchy)

        bb0 = BuildingBlock.from_code(0, "Val")
        bb1 = BuildingBlock.from_code(1, "Pro")
        bb2 = BuildingBlock.from_code(2, "Ala")

        compound2 = Compound([bb0, bb1], sample_chromatogram)
        compound3 = Compound([bb0, bb2], sample_chromatogram)
        compound4 = Compound([bb1, bb2], sample_chromatogram)

        result.set_validation_status(sample_compound, ValidationStatus.VALIDATED)
        result.set_validation_status(compound2, ValidationStatus.VALIDATED)
        result.set_validation_status(compound3, ValidationStatus.FAILED)
        result.set_validation_status(compound4, ValidationStatus.UNCERTAIN)

        summary = result.get_validation_summary()

        assert summary[ValidationStatus.VALIDATED] == 2
        assert summary[ValidationStatus.FAILED] == 1
        assert summary[ValidationStatus.UNCERTAIN] == 1
        assert ValidationStatus.LIKELY_SUCCESS not in summary

    def test_get_validation_summary_empty(self, sample_hierarchy):
        """Test getting validation summary when empty."""
        result = AnalysisResult(hierarchy=sample_hierarchy)

        summary = result.get_validation_summary()
        assert len(summary) == 0


class TestGetEquivalenceClass:
    """Test getting equivalence class by residue sequence."""

    def test_get_equivalence_class_exists(self, sample_hierarchy):
        """Test getting equivalence class that exists."""
        eq_class = EquivalenceClass(residue_sequence="Val-Pro")
        result = AnalysisResult(
            hierarchy=sample_hierarchy, equivalence_classes=[eq_class]
        )

        found = result.get_equivalence_class("Val-Pro")
        assert found == eq_class

    def test_get_equivalence_class_not_exists(self, sample_hierarchy):
        """Test getting equivalence class that doesn't exist."""
        eq_class = EquivalenceClass(residue_sequence="Val-Pro")
        result = AnalysisResult(
            hierarchy=sample_hierarchy, equivalence_classes=[eq_class]
        )

        found = result.get_equivalence_class("Leu-Ala")
        assert found is None

    def test_get_equivalence_class_empty_list(self, sample_hierarchy):
        """Test getting equivalence class from empty list."""
        result = AnalysisResult(hierarchy=sample_hierarchy)

        found = result.get_equivalence_class("Val-Pro")
        assert found is None


class TestTotalCompounds:
    """Test total compounds count."""

    def test_total_compounds(self, sample_hierarchy, sample_compound, sample_chromatogram):
        """Test getting total compounds from hierarchy."""
        sample_hierarchy.add_compound(sample_compound)

        bb0 = BuildingBlock.from_code(0, "Val")
        bb1 = BuildingBlock.from_code(1, "Pro")
        compound2 = Compound([bb0, bb1], sample_chromatogram)
        sample_hierarchy.add_compound(compound2)

        result = AnalysisResult(hierarchy=sample_hierarchy)

        assert result.total_compounds() == 2

    def test_total_compounds_empty(self, sample_hierarchy):
        """Test total compounds when hierarchy is empty."""
        result = AnalysisResult(hierarchy=sample_hierarchy)

        assert result.total_compounds() == 0


class TestTotalPeaksDetected:
    """Test total peaks detected count."""

    def test_total_peaks_detected(
        self, sample_hierarchy, sample_compound, sample_chromatogram
    ):
        """Test getting total peaks detected."""
        result = AnalysisResult(hierarchy=sample_hierarchy)

        # Add 2 peaks for compound 1
        peak1 = Peak(
            position=2.5,
            left_base=2.0,
            right_base=3.0,
            height=300.0,
            area=450.0,
        )
        c1 = PeakClassification(
            compound=sample_compound,
            peak=peak1,
            classification=PeakType.PUTATIVE_PRODUCT,
        )

        peak2 = Peak(
            position=1.8,
            left_base=1.5,
            right_base=2.1,
            height=150.0,
            area=200.0,
        )
        c2 = PeakClassification(
            compound=sample_compound,
            peak=peak2,
            classification=PeakType.TRUNCATION,
        )

        result.add_peak_classification(sample_compound, c1)
        result.add_peak_classification(sample_compound, c2)

        # Add 1 peak for compound 2
        bb0 = BuildingBlock.from_code(0, "Val")
        bb1 = BuildingBlock.from_code(1, "Pro")
        compound2 = Compound([bb0, bb1], sample_chromatogram)

        peak3 = Peak(
            position=3.0,
            left_base=2.5,
            right_base=3.5,
            height=400.0,
            area=600.0,
        )
        c3 = PeakClassification(
            compound=compound2,
            peak=peak3,
            classification=PeakType.PUTATIVE_PRODUCT,
        )
        result.add_peak_classification(compound2, c3)

        assert result.total_peaks_detected() == 3

    def test_total_peaks_detected_empty(self, sample_hierarchy):
        """Test total peaks detected when none exist."""
        result = AnalysisResult(hierarchy=sample_hierarchy)

        assert result.total_peaks_detected() == 0


class TestStringRepresentation:
    """Test string representation."""

    def test_repr(self, sample_hierarchy, sample_compound, sample_chromatogram):
        """Test repr representation."""
        sample_hierarchy.add_compound(sample_compound)
        result = AnalysisResult(hierarchy=sample_hierarchy)

        result.set_validation_status(sample_compound, ValidationStatus.VALIDATED)

        peak = Peak(
            position=2.5,
            left_base=2.0,
            right_base=3.0,
            height=300.0,
            area=450.0,
        )
        classification = PeakClassification(
            compound=sample_compound,
            peak=peak,
            classification=PeakType.PUTATIVE_PRODUCT,
        )
        result.add_peak_classification(sample_compound, classification)

        repr_str = repr(result)
        assert "AnalysisResult" in repr_str
        assert "compounds=1" in repr_str
        assert "peaks=1" in repr_str
        assert "validated=1" in repr_str
