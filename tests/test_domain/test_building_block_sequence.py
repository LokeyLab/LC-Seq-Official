"""
Tests for BuildingBlockSequence value object.

Tests implementation against THEORY.md Section 2.2, 1.5.1 specifications.
"""

import pytest
from lcseq.domain.entities.building_block import BuildingBlock
from lcseq.domain.value_objects.building_block_sequence import BuildingBlockSequence


class TestBuildingBlockSequenceCreation:
    """Test BuildingBlockSequence instantiation and factory methods."""

    def test_create_from_blocks_dict(self):
        """Direct creation with blocks dictionary."""
        bb0 = BuildingBlock.from_code(0, "Pro")
        bb1 = BuildingBlock.from_code(1, "Leu")

        seq = BuildingBlockSequence(blocks={0: bb0, 1: bb1})

        assert len(seq) == 2
        assert seq.get_block_at(0) == bb0
        assert seq.get_block_at(1) == bb1

    def test_from_blocks_list(self):
        """Factory method from ordered block list."""
        bb0 = BuildingBlock.from_code(0, "Pro")
        bb1 = BuildingBlock.from_code(1, "Leu")
        bb2 = BuildingBlock.from_code(2, "Ala")

        seq = BuildingBlockSequence.from_blocks([bb0, bb1, bb2])

        assert len(seq) == 3
        assert seq.get_block_at(0).code == "Pro"
        assert seq.get_block_at(1).code == "Leu"
        assert seq.get_block_at(2).code == "Ala"

    def test_from_codes_simple(self):
        """Factory method from code strings."""
        seq = BuildingBlockSequence.from_codes(["Pro", "Leu", "Ala"])

        assert len(seq) == 3
        assert seq.get_block_at(0).code == "Pro"
        assert seq.get_block_at(1).code == "Leu"
        assert seq.get_block_at(2).code == "Ala"

    def test_from_codes_with_nulls(self):
        """Factory method with null blocks."""
        seq = BuildingBlockSequence.from_codes(["Pro", "Null", "Leu"])

        assert len(seq) == 3
        assert seq.get_block_at(0).code == "Pro"
        assert seq.get_block_at(1).is_null is True
        assert seq.get_block_at(2).code == "Leu"

    def test_from_codes_auto_detects_nulls(self):
        """Null detection is automatic from codes."""
        seq = BuildingBlockSequence.from_codes(["Pro", "null", "Leu"])

        assert seq.get_block_at(1).is_null is True

    def test_single_block_sequence(self):
        """Single building block sequence."""
        seq = BuildingBlockSequence.from_codes(["Pro"])

        assert len(seq) == 1
        assert seq.get_block_at(0).code == "Pro"

    def test_immutable(self):
        """BuildingBlockSequence is immutable (frozen dataclass)."""
        seq = BuildingBlockSequence.from_codes(["Pro", "Leu"])

        with pytest.raises(AttributeError):
            seq.blocks = {}


class TestBuildingBlockSequenceValidation:
    """Test BuildingBlockSequence validation rules."""

    def test_empty_sequence_raises_error(self):
        """Sequence cannot be empty."""
        with pytest.raises(ValueError, match="cannot be empty"):
            BuildingBlockSequence(blocks={})

    def test_negative_position_raises_error(self):
        """Positions must be non-negative (THEORY.md Section 1.5.1)."""
        # BuildingBlock validates cycle, so create with valid cycle first
        bb0 = BuildingBlock.from_code(0, "Leu")
        bb1 = BuildingBlock.from_code(1, "Ala")

        # Try to create sequence with negative position key
        with pytest.raises(ValueError, match="non-negative"):
            BuildingBlockSequence(blocks={-1: bb0, 0: bb1})

    def test_non_contiguous_positions_raises_error(self):
        """Positions must be contiguous 0, 1, 2, ..., n-1."""
        bb0 = BuildingBlock.from_code(0, "Pro")
        bb2 = BuildingBlock.from_code(2, "Leu")  # Missing position 1

        with pytest.raises(ValueError, match="contiguous"):
            BuildingBlockSequence(blocks={0: bb0, 2: bb2})

    def test_mismatched_cycle_raises_error(self):
        """Block cycle must match its position."""
        bb0 = BuildingBlock.from_code(0, "Pro")
        bb1_wrong = BuildingBlock.from_code(5, "Leu")  # Says cycle 5

        with pytest.raises(ValueError, match="mismatched cycle"):
            BuildingBlockSequence(blocks={0: bb0, 1: bb1_wrong})  # But at position 1

    def test_positions_must_start_at_zero(self):
        """Positions must start at 0 (C-terminus)."""
        bb1 = BuildingBlock.from_code(1, "Leu")
        bb2 = BuildingBlock.from_code(2, "Ala")

        with pytest.raises(ValueError, match="contiguous"):
            BuildingBlockSequence(blocks={1: bb1, 2: bb2})  # Missing position 0


class TestPositionAccess:
    """Test accessing blocks by position."""

    def test_get_block_at_valid_position(self):
        """Get block at valid position."""
        seq = BuildingBlockSequence.from_codes(["Pro", "Leu", "Ala"])

        block = seq.get_block_at(1)
        assert block.code == "Leu"

    def test_get_block_at_c_terminus(self):
        """Position 0 is C-terminus (THEORY.md Section 1.5.1)."""
        seq = BuildingBlockSequence.from_codes(["Pro", "Leu", "Ala"])

        c_terminus = seq.get_block_at(0)
        assert c_terminus.code == "Pro"

    def test_get_block_at_n_terminus(self):
        """Highest position is N-terminus."""
        seq = BuildingBlockSequence.from_codes(["Pro", "Leu", "Ala"])

        n_terminus = seq.get_block_at(2)
        assert n_terminus.code == "Ala"

    def test_get_block_at_invalid_position_raises_error(self):
        """Invalid position raises KeyError."""
        seq = BuildingBlockSequence.from_codes(["Pro", "Leu"])

        with pytest.raises(KeyError):
            seq.get_block_at(5)

    def test_get_non_null_blocks(self):
        """Get only non-null blocks in position order."""
        seq = BuildingBlockSequence.from_codes(["Pro", "Null", "Leu", "Ala", "Null"])

        non_null = seq.get_non_null_blocks()

        assert len(non_null) == 3
        assert non_null[0].code == "Pro"  # Position 0
        assert non_null[1].code == "Leu"  # Position 2
        assert non_null[2].code == "Ala"  # Position 3

    def test_get_non_null_blocks_all_null(self):
        """All null sequence returns empty list."""
        seq = BuildingBlockSequence.from_codes(["Null", "Null"])

        non_null = seq.get_non_null_blocks()
        assert non_null == []

    def test_get_non_null_blocks_no_nulls(self):
        """No null sequence returns all blocks."""
        seq = BuildingBlockSequence.from_codes(["Pro", "Leu", "Ala"])

        non_null = seq.get_non_null_blocks()
        assert len(non_null) == 3


class TestSequenceConversions:
    """Test string conversion methods."""

    def test_to_positional_string_no_nulls(self):
        """Positional string with no nulls (N→C display)."""
        seq = BuildingBlockSequence.from_codes(["Pro", "Leu", "Ala"])

        # N→C: Ala-Leu-Pro (pos 2 - pos 1 - pos 0)
        assert seq.to_positional_string() == "Ala-Leu-Pro"

    def test_to_positional_string_with_nulls(self):
        """Positional string includes null positions."""
        seq = BuildingBlockSequence.from_codes(["Pro", "Null", "Leu"])

        # N→C: Leu-Null-Pro
        assert seq.to_positional_string() == "Leu-Null-Pro"

    def test_to_residue_string_removes_nulls(self):
        """Residue string excludes null blocks (THEORY.md Section 2.2)."""
        seq = BuildingBlockSequence.from_codes(["Pro", "Null", "Leu"])

        # N→C non-null only: Leu-Pro
        assert seq.to_residue_string() == "Leu-Pro"

    def test_to_residue_string_no_nulls(self):
        """Residue string same as positional when no nulls."""
        seq = BuildingBlockSequence.from_codes(["Pro", "Leu", "Ala"])

        assert seq.to_residue_string() == "Ala-Leu-Pro"

    def test_to_residue_string_all_nulls(self):
        """All null sequence produces empty residue string."""
        seq = BuildingBlockSequence.from_codes(["Null", "Null"])

        assert seq.to_residue_string() == ""

    def test_str_returns_positional_string(self):
        """str() returns positional sequence."""
        seq = BuildingBlockSequence.from_codes(["Pro", "Null", "Leu"])

        assert str(seq) == "Leu-Null-Pro"

    def test_repr_shows_position_mapping(self):
        """repr() shows position: code mapping."""
        seq = BuildingBlockSequence.from_codes(["Pro", "Leu"])

        repr_str = repr(seq)
        assert "0: Pro" in repr_str
        assert "1: Leu" in repr_str


class TestNSequenceConvention:
    """Test N→C sequence convention (THEORY.md Section 1.5.1)."""

    def test_position_0_is_c_terminus(self):
        """Position 0 = C-terminus (synthesized first)."""
        seq = BuildingBlockSequence.from_codes(["Pro", "Leu", "Ala"])

        # Pro at position 0 = C-terminus
        assert seq.get_block_at(0).code == "Pro"

    def test_highest_position_is_n_terminus(self):
        """Highest position = N-terminus (synthesized last)."""
        seq = BuildingBlockSequence.from_codes(["Pro", "Leu", "Ala"])

        # Ala at position 2 = N-terminus
        assert seq.get_block_at(2).code == "Ala"

    def test_string_displays_n_to_c(self):
        """String representation is N→C (left to right) for readability."""
        seq = BuildingBlockSequence.from_codes(["Pro", "Leu", "Ala"])

        # Display: Ala-Leu-Pro (N→C, left to right)
        # Storage: {0: Pro, 1: Leu, 2: Ala} (C→N indexing)
        assert str(seq) == "Ala-Leu-Pro"


class TestPositionalVariants:
    """Test that different positions create different sequences (THEORY.md Section 1.2)."""

    def test_same_residues_different_positions_are_different(self):
        """Same chemical identity, different positions → different sequences."""
        # Val at position 0, nulls elsewhere
        seq1 = BuildingBlockSequence.from_codes(["Val", "Null", "Null"])

        # Val at position 1, nulls elsewhere
        seq2 = BuildingBlockSequence.from_codes(["Null", "Val", "Null"])

        # Different positional sequences
        assert seq1.to_positional_string() != seq2.to_positional_string()
        assert seq1 != seq2

    def test_same_residue_sequence_different_positions(self):
        """Multiple positional sequences can have same residue sequence."""
        seq1 = BuildingBlockSequence.from_codes(["Val", "Null", "Null"])
        seq2 = BuildingBlockSequence.from_codes(["Null", "Val", "Null"])
        seq3 = BuildingBlockSequence.from_codes(["Null", "Null", "Val"])

        # All have same residue sequence (just "Val")
        assert seq1.to_residue_string() == "Val"
        assert seq2.to_residue_string() == "Val"
        assert seq3.to_residue_string() == "Val"

        # But different positional sequences
        assert seq1 != seq2 != seq3


class TestEquality:
    """Test equality and hashing."""

    def test_equality_same_values(self):
        """Sequences with same blocks are equal."""
        seq1 = BuildingBlockSequence.from_codes(["Pro", "Leu"])
        seq2 = BuildingBlockSequence.from_codes(["Pro", "Leu"])

        assert seq1 == seq2

    def test_inequality_different_codes(self):
        """Different codes → not equal."""
        seq1 = BuildingBlockSequence.from_codes(["Pro", "Leu"])
        seq2 = BuildingBlockSequence.from_codes(["Pro", "Ala"])

        assert seq1 != seq2

    def test_inequality_different_null_positions(self):
        """Different null positions → not equal."""
        seq1 = BuildingBlockSequence.from_codes(["Pro", "Null"])
        seq2 = BuildingBlockSequence.from_codes(["Null", "Pro"])

        assert seq1 != seq2

    def test_not_hashable_due_to_dict(self):
        """BuildingBlockSequence is not hashable (contains dict)."""
        seq1 = BuildingBlockSequence.from_codes(["Pro", "Leu"])

        # Cannot use in set because dict is unhashable
        with pytest.raises(TypeError, match="unhashable"):
            {seq1}


class TestCompositeBlocks:
    """Test sequences with composite building blocks."""

    def test_composite_block_in_sequence(self):
        """Composite blocks work in sequences."""
        seq = BuildingBlockSequence.from_codes(["Pro", "Leu-Ala-Val", "Gly"])

        assert len(seq) == 3
        assert seq.get_block_at(1).code == "Leu-Ala-Val"

    def test_positional_string_preserves_composite(self):
        """Positional string preserves composite block notation."""
        seq = BuildingBlockSequence.from_codes(["Pro", "Leu-Ala-Val", "Gly"])

        # N→C: Gly - (Leu-Ala-Val) - Pro
        assert seq.to_positional_string() == "Gly-Leu-Ala-Val-Pro"

    def test_residue_string_preserves_composite(self):
        """Residue string preserves composite blocks."""
        seq = BuildingBlockSequence.from_codes(["Pro", "Null", "Leu-Ala-Val"])

        # Non-null only: (Leu-Ala-Val) - Pro
        assert seq.to_residue_string() == "Leu-Ala-Val-Pro"


class TestEdgeCases:
    """Test edge cases and boundary conditions."""

    def test_long_sequence(self):
        """Handle long sequences."""
        codes = [f"BB{i}" for i in range(20)]
        seq = BuildingBlockSequence.from_codes(codes)

        assert len(seq) == 20
        assert seq.get_block_at(0).code == "BB0"
        assert seq.get_block_at(19).code == "BB19"

    def test_all_null_sequence(self):
        """All null blocks is valid."""
        seq = BuildingBlockSequence.from_codes(["Null", "Null", "Null"])

        assert len(seq) == 3
        assert seq.get_non_null_blocks() == []
        assert seq.to_residue_string() == ""

    def test_modified_residues_in_sequence(self):
        """Modified residues work correctly."""
        seq = BuildingBlockSequence.from_codes(["DPro", "NMeLeu", "Ala"])

        assert seq.get_block_at(0).code == "DPro"
        assert seq.get_block_at(1).code == "NMeLeu"
        assert seq.to_positional_string() == "Ala-NMeLeu-DPro"
