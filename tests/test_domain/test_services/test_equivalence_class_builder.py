"""
Tests for EquivalenceClassBuilder service.
"""

import pytest
import numpy as np
from lcseq.domain.entities.building_block import BuildingBlock
from lcseq.domain.entities.compound import Compound
from lcseq.domain.entities.chromatogram import Chromatogram
from lcseq.domain.services.equivalence_class_builder import EquivalenceClassBuilder


@pytest.fixture
def chromatogram():
    """Simple chromatogram for testing."""
    return Chromatogram(
        time_points=np.array([1.0, 2.0, 3.0]),
        counts=np.array([100.0, 200.0, 150.0])
    )


@pytest.fixture
def builder():
    """EquivalenceClassBuilder instance."""
    return EquivalenceClassBuilder()


def make_compound(codes, chromatogram):
    """Helper to make compound with codes in order."""
    sorted_codes = sorted(codes, key=lambda x: x[0])
    blocks = [BuildingBlock.from_code(cycle, code) for cycle, code in sorted_codes]
    return Compound(blocks, chromatogram)


class TestBasicFunctionality:
    """Test basic equivalence class building."""

    def test_empty_list(self, builder):
        """Test with empty compound list."""
        classes = builder.build([])
        assert len(classes) == 0

    def test_single_compound(self, builder, chromatogram):
        """Test with single compound."""
        compound = make_compound([(0, "Leu")], chromatogram)
        classes = builder.build([compound])

        assert len(classes) == 1
        assert len(classes[0].compounds) == 1
        assert classes[0].residue_sequence == "Leu"

    def test_same_residue_sequence(self, builder, chromatogram):
        """Test compounds with same residue sequence grouped together."""
        # Leu-Pro (different null positions)
        c1 = make_compound([(0, "Leu"), (1, "Pro"), (2, "Null")], chromatogram)
        c2 = make_compound([(0, "Leu"), (1, "Null"), (2, "Pro")], chromatogram)
        c3 = make_compound([(0, "Null"), (1, "Leu"), (2, "Pro")], chromatogram)

        classes = builder.build([c1, c2, c3])

        # All have same residue sequence "Leu-Pro"
        assert len(classes) == 1
        assert len(classes[0].compounds) == 3
        assert classes[0].residue_sequence == "Pro-Leu"

    def test_different_residue_sequences(self, builder, chromatogram):
        """Test compounds with different residue sequences."""
        c1 = make_compound([(0, "Leu"), (1, "Pro")], chromatogram)
        c2 = make_compound([(0, "Val"), (1, "Ala")], chromatogram)

        classes = builder.build([c1, c2])

        assert len(classes) == 2
        # Sorted by residue sequence
        assert classes[0].residue_sequence in ["Pro-Leu", "Ala-Val"]
        assert classes[1].residue_sequence in ["Pro-Leu", "Ala-Val"]

    def test_all_null_compound(self, builder, chromatogram):
        """Test all-null compound."""
        c1 = make_compound([(0, "Null"), (1, "Null")], chromatogram)

        classes = builder.build([c1])

        assert len(classes) == 1
        assert classes[0].residue_sequence == ""


class TestGrouping:
    """Test grouping logic."""

    def test_get_residue_sequence(self, builder, chromatogram):
        """Test residue sequence extraction."""
        compound = make_compound([(0, "Leu"), (1, "Null"), (2, "Pro")], chromatogram)

        residue_seq = builder._get_residue_sequence(compound)

        assert residue_seq == "Pro-Leu"

    def test_multiple_groups(self, builder, chromatogram):
        """Test creating multiple equivalence classes."""
        # 3 different residue sequences
        c1 = make_compound([(0, "Leu")], chromatogram)
        c2 = make_compound([(0, "Pro")], chromatogram)
        c3 = make_compound([(0, "Val")], chromatogram)

        classes = builder.build([c1, c2, c3])

        assert len(classes) == 3
        # Each class has one compound
        assert all(len(c.compounds) == 1 for c in classes)
