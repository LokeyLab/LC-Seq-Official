"""
Tests for PeakClassifier service.

Tests peak type classification via DAG constraints.
"""

import pytest
import numpy as np

from lcseq.domain.entities.building_block import BuildingBlock
from lcseq.domain.entities.chromatogram import Chromatogram
from lcseq.domain.entities.compound import Compound
from lcseq.domain.entities.peak import Peak, PeakType
from lcseq.domain.models.compound_hierarchy import CompoundHierarchy, HierarchyMode
from lcseq.domain.services.peak_classifier import PeakClassifier


class TestPeakClassifier:
    """Test PeakClassifier service."""

    @pytest.fixture
    def simple_hierarchy(self):
        """Create simple hierarchy for testing."""
        hierarchy = CompoundHierarchy(mode=HierarchyMode.BUILDING_BLOCK)

        # Create chromatogram
        time_points = np.linspace(0, 100, 200)
        counts = np.ones(200)
        chrom = Chromatogram(time_points=time_points, counts=counts)

        # L0 (null compound)
        l0 = Compound(
            [BuildingBlock.from_code(0, "Null"),
             BuildingBlock.from_code(1, "Null")],
            chrom
        )
        l0_peak = Peak(position=10.0, left_base=5, right_base=15, height=50, area=100)
        l0.detected_peaks = [l0_peak]
        l0.selected_peak = l0_peak

        # Level 1 compound
        l1 = Compound(
            [BuildingBlock.from_code(0, "Pro"),
             BuildingBlock.from_code(1, "Null")],
            chrom
        )
        l1_peak = Peak(position=30.0, left_base=25, right_base=35, height=60, area=120)
        l1.detected_peaks = [l1_peak]
        l1.selected_peak = l1_peak

        # Level 2 compound
        l2 = Compound(
            [BuildingBlock.from_code(0, "Pro"),
             BuildingBlock.from_code(1, "Leu")],
            chrom
        )

        # Add to hierarchy
        hierarchy.add_compound(l0)
        hierarchy.add_compound(l1)
        hierarchy.add_compound(l2)
        # Edges go from ancestor (longer) to descendant (truncation/shorter)
        hierarchy.add_edge(l1, l0)  # L1 -> L0 (L1 is ancestor, L0 is descendant)
        hierarchy.add_edge(l2, l1)  # L2 -> L1 (L2 is ancestor, L1 is descendant)
        hierarchy.add_edge(l2, l0)  # L2 -> L0 (L2 is ancestor, L0 is descendant)

        return hierarchy

    def test_classify_null_peak(self, simple_hierarchy):
        """Test classification of NULL peak."""
        classifier = PeakClassifier()

        # Get L2 compound
        l2 = [c for c in simple_hierarchy.compounds if c.level == 2][0]

        # Create peak at L0 position (10.0)
        peak = Peak(position=10.0, left_base=8, right_base=12, height=20, area=40)

        # Classify
        classified = classifier.classify_peak(
            l2, peak, simple_hierarchy, l0_retention_time=10.0, tolerance=0.01
        )

        assert classified.peak_type == PeakType.NULL

    def test_classify_truncation_peak(self, simple_hierarchy):
        """Test classification of TRUNCATION peak."""
        classifier = PeakClassifier()

        # Get L2 compound
        l2 = [c for c in simple_hierarchy.compounds if c.level == 2][0]

        # Create peak at L1 product position (30.0)
        peak = Peak(position=30.0, left_base=28, right_base=32, height=25, area=50)

        # Classify
        classified = classifier.classify_peak(
            l2, peak, simple_hierarchy, l0_retention_time=10.0, tolerance=0.05
        )

        assert classified.peak_type == PeakType.TRUNCATION

    def test_classify_putative_product_peak(self, simple_hierarchy):
        """Test classification of PUTATIVE_PRODUCT peak."""
        classifier = PeakClassifier()

        # Get L2 compound
        l2 = [c for c in simple_hierarchy.compounds if c.level == 2][0]

        # Create peak after truncation positions (>30.0)
        peak = Peak(position=50.0, left_base=45, right_base=55, height=40, area=80)

        # Classify
        classified = classifier.classify_peak(
            l2, peak, simple_hierarchy, l0_retention_time=10.0, tolerance=0.05
        )

        assert classified.peak_type == PeakType.PUTATIVE_PRODUCT

    def test_classify_unknown_peak(self, simple_hierarchy):
        """Test classification of UNKNOWN peak."""
        classifier = PeakClassifier()

        # Get L1 compound
        l1 = [c for c in simple_hierarchy.compounds if c.level == 1][0]

        # Create peak BEFORE truncation positions - should be UNKNOWN
        # L0 is at 10.0, so a peak at 5.0 is before all truncations
        peak = Peak(position=5.0, left_base=3, right_base=7, height=15, area=30)

        # Classify
        classified = classifier.classify_peak(
            l1, peak, simple_hierarchy, l0_retention_time=10.0, tolerance=0.01
        )

        # Should be UNKNOWN (doesn't match any expected position and before truncations)
        assert classified.peak_type == PeakType.UNKNOWN

    def test_classify_all_peaks(self, simple_hierarchy):
        """Test classification of all peaks for a compound."""
        classifier = PeakClassifier()

        # Get L2 compound
        l2 = [c for c in simple_hierarchy.compounds if c.level == 2][0]

        # Add multiple peaks
        l2.detected_peaks = [
            Peak(position=10.0, left_base=8, right_base=12, height=20, area=40),   # NULL
            Peak(position=30.0, left_base=28, right_base=32, height=25, area=50),  # TRUNCATION
            Peak(position=50.0, left_base=45, right_base=55, height=40, area=80),  # PUTATIVE_PRODUCT
            Peak(position=70.0, left_base=65, right_base=75, height=15, area=30),  # UNKNOWN
        ]

        # Classify all
        classifier.classify_all_peaks(l2, simple_hierarchy, l0_retention_time=10.0, tolerance=0.05)

        # Check classifications
        assert l2.detected_peaks[0].peak_type == PeakType.NULL
        assert l2.detected_peaks[1].peak_type == PeakType.TRUNCATION
        assert l2.detected_peaks[2].peak_type == PeakType.PUTATIVE_PRODUCT
        assert l2.detected_peaks[3].peak_type == PeakType.UNKNOWN

        # Selected peak should be the PUTATIVE_PRODUCT
        assert l2.selected_peak == l2.detected_peaks[2]


class TestPeakClassifierTolerances:
    """Test tolerance handling."""

    def test_strict_tolerance(self):
        """Test that strict tolerance requires close match."""
        classifier = PeakClassifier()

        peak = Peak(position=10.5, left_base=9, right_base=12, height=20, area=40)

        # Should match with loose tolerance
        assert classifier._is_null_peak(peak, l0_retention_time=10.0, tolerance=0.1)

        # Should not match with strict tolerance
        assert not classifier._is_null_peak(peak, l0_retention_time=10.0, tolerance=0.01)

    def test_relative_tolerance(self):
        """Test that tolerance is relative to position."""
        classifier = PeakClassifier()

        # At position 100, tolerance of 0.1 allows ±10
        peak1 = Peak(position=105, left_base=100, right_base=110, height=20, area=40)
        assert classifier._is_null_peak(peak1, l0_retention_time=100.0, tolerance=0.1)

        # At position 10, tolerance of 0.1 allows ±1
        peak2 = Peak(position=11.5, left_base=10, right_base=13, height=20, area=40)
        assert not classifier._is_null_peak(peak2, l0_retention_time=10.0, tolerance=0.1)


class TestPeakClassifierHelperMethods:
    """Test helper methods."""

    def test_is_null_peak(self):
        """Test null peak detection."""
        classifier = PeakClassifier()

        peak_match = Peak(position=10.0, left_base=9, right_base=11, height=20, area=40)
        peak_close = Peak(position=10.05, left_base=9, right_base=11, height=20, area=40)
        peak_far = Peak(position=15.0, left_base=14, right_base=16, height=20, area=40)

        assert classifier._is_null_peak(peak_match, 10.0, tolerance=0.01)
        assert classifier._is_null_peak(peak_close, 10.0, tolerance=0.01)
        assert not classifier._is_null_peak(peak_far, 10.0, tolerance=0.01)

    def test_is_truncation_peak(self):
        """Test truncation peak detection."""
        classifier = PeakClassifier()

        ancestor_positions = [20.0, 30.0]

        peak_match = Peak(position=30.0, left_base=28, right_base=32, height=20, area=40)
        peak_far = Peak(position=50.0, left_base=48, right_base=52, height=20, area=40)

        assert classifier._is_truncation_peak(
            peak_match, ancestor_positions, l0_retention_time=10.0, tolerance=0.05
        )
        assert not classifier._is_truncation_peak(
            peak_far, ancestor_positions, l0_retention_time=10.0, tolerance=0.05
        )

    def test_find_l0_retention_time(self):
        """Test finding L0 retention time from hierarchy."""
        classifier = PeakClassifier()

        hierarchy = CompoundHierarchy(mode=HierarchyMode.BUILDING_BLOCK)
        time_points = np.linspace(0, 100, 200)
        counts = np.ones(200)
        chrom = Chromatogram(time_points=time_points, counts=counts)

        # Create L0 with selected peak
        l0 = Compound(
            [BuildingBlock.from_code(0, "Null")],
            chrom
        )
        l0_peak = Peak(position=12.5, left_base=10, right_base=15, height=50, area=100)
        l0.selected_peak = l0_peak
        hierarchy.add_compound(l0)

        rt = classifier._find_l0_retention_time(hierarchy)
        assert rt == 12.5

    def test_get_product_positions(self):
        """Test extracting product positions from compounds."""
        classifier = PeakClassifier()

        time_points = np.linspace(0, 100, 200)
        counts = np.ones(200)
        chrom = Chromatogram(time_points=time_points, counts=counts)

        compounds = [
            Compound([BuildingBlock.from_code(0, "Pro")], chrom),
            Compound([BuildingBlock.from_code(0, "Leu")], chrom),
        ]

        compounds[0].selected_peak = Peak(position=20.0, left_base=18, right_base=22, height=30, area=60)
        compounds[1].selected_peak = Peak(position=35.0, left_base=33, right_base=37, height=40, area=80)

        positions = classifier._get_product_positions(compounds)

        assert positions == [20.0, 35.0]


class TestPeakClassifierEdgeCases:
    """Test edge cases."""

    def test_classify_with_no_ancestors(self):
        """Test classification when compound has no ancestors."""
        classifier = PeakClassifier()

        hierarchy = CompoundHierarchy(mode=HierarchyMode.BUILDING_BLOCK)
        time_points = np.linspace(0, 100, 200)
        counts = np.ones(200)
        chrom = Chromatogram(time_points=time_points, counts=counts)

        # Single compound (no ancestors)
        compound = Compound([BuildingBlock.from_code(0, "Pro")], chrom)
        hierarchy.add_compound(compound)

        peak = Peak(position=40.0, left_base=35, right_base=45, height=30, area=60)

        # Should classify as PUTATIVE_PRODUCT (first peak after no truncations)
        classified = classifier.classify_peak(
            compound, peak, hierarchy, l0_retention_time=None, tolerance=0.05
        )

        assert classified.peak_type == PeakType.PUTATIVE_PRODUCT

    def test_classify_with_no_peaks(self):
        """Test classifying compound with no detected peaks."""
        classifier = PeakClassifier()

        hierarchy = CompoundHierarchy(mode=HierarchyMode.BUILDING_BLOCK)
        time_points = np.linspace(0, 100, 200)
        counts = np.ones(200)
        chrom = Chromatogram(time_points=time_points, counts=counts)

        compound = Compound([BuildingBlock.from_code(0, "Pro")], chrom)
        hierarchy.add_compound(compound)

        # No peaks
        compound.detected_peaks = []

        # Should handle gracefully
        classifier.classify_all_peaks(compound, hierarchy)

        assert len(compound.detected_peaks) == 0
