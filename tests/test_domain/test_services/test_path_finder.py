"""Tests for PathFinder service."""

import pytest
import numpy as np
from lcseq.domain.entities.building_block import BuildingBlock
from lcseq.domain.entities.compound import Compound
from lcseq.domain.entities.chromatogram import Chromatogram
from lcseq.domain.models.compound_hierarchy import CompoundHierarchy, HierarchyMode
from lcseq.domain.services.path_finder import PathFinder


@pytest.fixture
def chromatogram():
    return Chromatogram(
        time_points=np.array([1.0, 2.0, 3.0]),
        counts=np.array([100.0, 200.0, 150.0])
    )


@pytest.fixture
def finder():
    return PathFinder()


def make_compound(codes, chromatogram):
    sorted_codes = sorted(codes, key=lambda x: x[0])
    blocks = [BuildingBlock.from_code(cycle, code) for cycle, code in sorted_codes]
    return Compound(blocks, chromatogram)


class TestPathFinding:
    """Test path finding functionality."""

    def test_find_path_exists(self, finder, chromatogram):
        """Test finding path when it exists."""
        hierarchy = CompoundHierarchy(mode=HierarchyMode.BUILDING_BLOCK)

        c1 = make_compound([(0, "Leu"), (1, "Pro")], chromatogram)
        c2 = make_compound([(0, "Leu"), (1, "Null")], chromatogram)

        hierarchy.add_compound(c1)
        hierarchy.add_compound(c2)
        hierarchy.add_edge(c1, c2)

        path = finder.find_path(hierarchy, c1, c2)

        assert path is not None
        assert path == [c1, c2]

    def test_find_path_not_exists(self, finder, chromatogram):
        """Test when no path exists."""
        hierarchy = CompoundHierarchy(mode=HierarchyMode.BUILDING_BLOCK)

        c1 = make_compound([(0, "Leu"), (1, "Pro")], chromatogram)
        c2 = make_compound([(0, "Val"), (1, "Ala")], chromatogram)

        hierarchy.add_compound(c1)
        hierarchy.add_compound(c2)

        path = finder.find_path(hierarchy, c1, c2)

        assert path is None

    def test_shortest_path(self, finder, chromatogram):
        """Test shortest path finding."""
        hierarchy = CompoundHierarchy(mode=HierarchyMode.BUILDING_BLOCK)

        c1 = make_compound([(0, "Leu"), (1, "Pro"), (2, "Val")], chromatogram)
        c2 = make_compound([(0, "Leu"), (1, "Pro"), (2, "Null")], chromatogram)
        c3 = make_compound([(0, "Leu"), (1, "Null"), (2, "Null")], chromatogram)

        hierarchy.add_compound(c1)
        hierarchy.add_compound(c2)
        hierarchy.add_compound(c3)
        hierarchy.add_edge(c1, c2)
        hierarchy.add_edge(c2, c3)

        path = finder.shortest_path(hierarchy, c1, c3)

        assert len(path) == 3
        assert path == [c1, c2, c3]
