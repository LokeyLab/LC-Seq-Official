"""Tests for ValidationChecker service."""

import pytest
import numpy as np
from lcseq.domain.entities.building_block import BuildingBlock
from lcseq.domain.entities.compound import Compound
from lcseq.domain.entities.chromatogram import Chromatogram
from lcseq.domain.models.compound_hierarchy import CompoundHierarchy, HierarchyMode
from lcseq.domain.services.validation_checker import ValidationChecker


@pytest.fixture
def chromatogram():
    return Chromatogram(
        time_points=np.array([1.0, 2.0, 3.0]),
        counts=np.array([100.0, 200.0, 150.0])
    )


@pytest.fixture
def checker():
    return ValidationChecker()


def make_compound(codes, chromatogram):
    sorted_codes = sorted(codes, key=lambda x: x[0])
    blocks = [BuildingBlock.from_code(cycle, code) for cycle, code in sorted_codes]
    return Compound(blocks, chromatogram)


class TestValidation:
    """Test validation functionality."""

    def test_is_valid_dag_true(self, checker, chromatogram):
        """Test valid DAG."""
        hierarchy = CompoundHierarchy(mode=HierarchyMode.BUILDING_BLOCK)

        c1 = make_compound([(0, "Leu"), (1, "Pro")], chromatogram)
        c2 = make_compound([(0, "Leu"), (1, "Null")], chromatogram)

        hierarchy.add_compound(c1)
        hierarchy.add_compound(c2)
        hierarchy.add_edge(c1, c2)

        assert checker.is_valid_dag(hierarchy) is True

    def test_has_cycles_false(self, checker, chromatogram):
        """Test no cycles in valid hierarchy."""
        hierarchy = CompoundHierarchy(mode=HierarchyMode.BUILDING_BLOCK)

        c1 = make_compound([(0, "Leu"), (1, "Pro")], chromatogram)
        c2 = make_compound([(0, "Leu"), (1, "Null")], chromatogram)

        hierarchy.add_compound(c1)
        hierarchy.add_compound(c2)
        hierarchy.add_edge(c1, c2)

        assert checker.has_cycles(hierarchy) is False

    def test_validate_level_ordering(self, checker, chromatogram):
        """Test level ordering validation."""
        hierarchy = CompoundHierarchy(mode=HierarchyMode.BUILDING_BLOCK)

        c1 = make_compound([(0, "Leu"), (1, "Pro")], chromatogram)
        c2 = make_compound([(0, "Leu"), (1, "Null")], chromatogram)

        hierarchy.add_compound(c1)
        hierarchy.add_compound(c2)
        hierarchy.add_edge(c1, c2)

        assert checker.validate_level_ordering(hierarchy) is True
