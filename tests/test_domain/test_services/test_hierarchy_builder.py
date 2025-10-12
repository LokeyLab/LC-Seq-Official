"""
Simplified tests for HierarchyBuilder service focused on key functionality.
"""

import pytest
import numpy as np
from lcseq.domain.entities.building_block import BuildingBlock
from lcseq.domain.entities.compound import Compound
from lcseq.domain.entities.chromatogram import Chromatogram
from lcseq.domain.models.compound_hierarchy import HierarchyMode
from lcseq.domain.services.hierarchy_builder import HierarchyBuilder


@pytest.fixture
def chromatogram():
    """Simple chromatogram for testing."""
    return Chromatogram(
        time_points=np.array([1.0, 2.0, 3.0]),
        counts=np.array([100.0, 200.0, 150.0])
    )


@pytest.fixture
def builder():
    """HierarchyBuilder instance."""
    return HierarchyBuilder()


def make_compound(codes, chromatogram):
    """Helper to make compound with codes in order.

    codes: list of (cycle, code) tuples
    """
    # Sort by cycle to ensure proper order
    sorted_codes = sorted(codes, key=lambda x: x[0])
    blocks = [BuildingBlock.from_code(cycle, code) for cycle, code in sorted_codes]
    return Compound(blocks, chromatogram)


class TestBasicFunctionality:
    """Test basic hierarchy building."""

    def test_empty_hierarchy(self, builder):
        """Test building empty hierarchy."""
        hierarchy = builder.build([], HierarchyMode.BUILDING_BLOCK)
        assert len(hierarchy.compounds) == 0
        assert hierarchy.edge_count() == 0

    def test_single_compound(self, builder, chromatogram):
        """Test hierarchy with one compound."""
        compound = make_compound([(0, "Leu")], chromatogram)
        hierarchy = builder.build([compound], HierarchyMode.BUILDING_BLOCK)
        assert len(hierarchy.compounds) == 1
        assert hierarchy.edge_count() == 0

    def test_simple_truncation(self, builder, chromatogram):
        """Test two compounds with truncation relationship."""
        # Leu-Pro (2 blocks)
        maximal = make_compound([(0, "Leu"), (1, "Pro")], chromatogram)
        # Leu-Null (1 block)
        truncation = make_compound([(0, "Leu"), (1, "Null")], chromatogram)

        hierarchy = builder.build([maximal, truncation], HierarchyMode.BUILDING_BLOCK)

        assert len(hierarchy.compounds) == 2
        assert hierarchy.edge_count() == 1
        assert truncation in hierarchy.get_descendants(maximal)

    def test_linear_chain(self, builder, chromatogram):
        """Test linear truncation chain."""
        # Leu-Pro-Val (3 blocks)
        c1 = make_compound([(0, "Leu"), (1, "Pro"), (2, "Val")], chromatogram)
        # Leu-Pro-Null (2 blocks)
        c2 = make_compound([(0, "Leu"), (1, "Pro"), (2, "Null")], chromatogram)
        # Leu-Null-Null (1 block)
        c3 = make_compound([(0, "Leu"), (1, "Null"), (2, "Null")], chromatogram)
        # All null (0 blocks)
        c4 = make_compound([(0, "Null"), (1, "Null"), (2, "Null")], chromatogram)

        hierarchy = builder.build([c1, c2, c3, c4], HierarchyMode.BUILDING_BLOCK)

        assert len(hierarchy.compounds) == 4
        # 3 edges in linear chain
        assert hierarchy.edge_count() == 3

        # Verify descendants
        all_desc = hierarchy.get_descendants(c1)
        assert c2 in all_desc
        assert c3 in all_desc
        assert c4 in all_desc

    def test_no_truncation_relationship(self, builder, chromatogram):
        """Test compounds with no truncation relationship."""
        # Leu-Pro
        c1 = make_compound([(0, "Leu"), (1, "Pro")], chromatogram)
        # Val-Ala (different blocks, no truncation)
        c2 = make_compound([(0, "Val"), (1, "Ala")], chromatogram)

        hierarchy = builder.build([c1, c2], HierarchyMode.BUILDING_BLOCK)

        assert len(hierarchy.compounds) == 2
        # No edges - no truncation
        assert hierarchy.edge_count() == 0


class TestEdgeDetection:
    """Test truncation edge detection logic."""

    def test_is_building_block_truncation_true(self, builder, chromatogram):
        """Test detection of valid truncation."""
        ancestor = make_compound([(0, "Leu"), (1, "Pro"), (2, "Val")], chromatogram)
        descendant = make_compound([(0, "Leu"), (1, "Null"), (2, "Val")], chromatogram)

        is_trunc = builder._is_building_block_truncation(ancestor, descendant)
        assert is_trunc is True

    def test_is_building_block_truncation_false_same(self, builder, chromatogram):
        """Test same compound is not truncation."""
        compound = make_compound([(0, "Leu"), (1, "Pro")], chromatogram)

        is_trunc = builder._is_building_block_truncation(compound, compound)
        assert is_trunc is False

    def test_is_building_block_truncation_false_two_removed(self, builder, chromatogram):
        """Test two removals is not direct truncation."""
        ancestor = make_compound([(0, "Leu"), (1, "Pro"), (2, "Val")], chromatogram)
        descendant = make_compound([(0, "Null"), (1, "Null"), (2, "Val")], chromatogram)

        is_trunc = builder._is_building_block_truncation(ancestor, descendant)
        assert is_trunc is False


class TestMonomerMode:
    """Test monomer-level hierarchy building."""

    def test_monomer_truncation_detection(self, builder, chromatogram):
        """Test monomer-level truncation detection."""
        # Leu-Pro (2 monomers)
        ancestor = make_compound([(0, "Leu"), (1, "Pro")], chromatogram)
        # Leu-Null (1 monomer)
        descendant = make_compound([(0, "Leu"), (1, "Null")], chromatogram)

        is_trunc = builder._is_monomer_truncation(ancestor, descendant)
        assert is_trunc is True
