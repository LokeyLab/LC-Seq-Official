"""Tests for LevelAnalyzer service."""

import pytest
import numpy as np
from lcseq.domain.entities.building_block import BuildingBlock
from lcseq.domain.entities.compound import Compound
from lcseq.domain.entities.chromatogram import Chromatogram
from lcseq.domain.models.compound_hierarchy import CompoundHierarchy, HierarchyMode
from lcseq.domain.services.level_analyzer import LevelAnalyzer


@pytest.fixture
def chromatogram():
    return Chromatogram(
        time_points=np.array([1.0, 2.0, 3.0]),
        counts=np.array([100.0, 200.0, 150.0])
    )


@pytest.fixture
def analyzer():
    return LevelAnalyzer()


def make_compound(codes, chromatogram):
    sorted_codes = sorted(codes, key=lambda x: x[0])
    blocks = [BuildingBlock.from_code(cycle, code) for cycle, code in sorted_codes]
    return Compound(blocks, chromatogram)


class TestLevelAnalysis:
    """Test level analysis functionality."""

    def test_get_level_distribution(self, analyzer, chromatogram):
        """Test getting level distribution."""
        hierarchy = CompoundHierarchy(mode=HierarchyMode.BUILDING_BLOCK)

        c1 = make_compound([(0, "Leu"), (1, "Pro")], chromatogram)
        c2 = make_compound([(0, "Leu"), (1, "Null")], chromatogram)
        c3 = make_compound([(0, "Null"), (1, "Null")], chromatogram)

        hierarchy.add_compound(c1)
        hierarchy.add_compound(c2)
        hierarchy.add_compound(c3)

        distribution = analyzer.get_level_distribution(hierarchy)

        assert distribution[2] == 1  # c1 has level 2
        assert distribution[1] == 1  # c2 has level 1
        assert distribution[0] == 1  # c3 has level 0

    def test_get_max_level(self, analyzer, chromatogram):
        """Test getting max level."""
        hierarchy = CompoundHierarchy(mode=HierarchyMode.BUILDING_BLOCK)

        c1 = make_compound([(0, "Leu"), (1, "Pro")], chromatogram)
        c2 = make_compound([(0, "Leu"), (1, "Null")], chromatogram)

        hierarchy.add_compound(c1)
        hierarchy.add_compound(c2)

        max_level = analyzer.get_max_level(hierarchy)

        assert max_level == 2

    def test_get_compounds_at_level(self, analyzer, chromatogram):
        """Test getting compounds at specific level."""
        hierarchy = CompoundHierarchy(mode=HierarchyMode.BUILDING_BLOCK)

        c1 = make_compound([(0, "Leu"), (1, "Pro")], chromatogram)
        c2 = make_compound([(0, "Val"), (1, "Ala")], chromatogram)
        c3 = make_compound([(0, "Leu"), (1, "Null")], chromatogram)

        hierarchy.add_compound(c1)
        hierarchy.add_compound(c2)
        hierarchy.add_compound(c3)

        level_2 = analyzer.get_compounds_at_level(hierarchy, 2)

        assert len(level_2) == 2
        assert c1 in level_2
        assert c2 in level_2
