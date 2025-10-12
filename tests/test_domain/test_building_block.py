"""
Tests for BuildingBlock entity.

Tests implementation against THEORY.md Section 2.1, 1.5.3 specifications.
"""

import pytest
from lcseq.domain.entities.building_block import BuildingBlock


class TestBuildingBlockCreation:
    """Test BuildingBlock instantiation and validation."""

    def test_create_monomer_block(self):
        """Single monomer building block."""
        bb = BuildingBlock(cycle=0, code="Leu", is_null=False)
        
        assert bb.cycle == 0
        assert bb.code == "Leu"
        assert bb.is_null is False

    def test_create_composite_block(self):
        """Composite building block (trimeric)."""
        bb = BuildingBlock(cycle=1, code="Leu-Ala-Val", is_null=False)
        
        assert bb.cycle == 1
        assert bb.code == "Leu-Ala-Val"
        assert bb.is_null is False

    def test_create_null_block(self):
        """Null building block."""
        bb = BuildingBlock(cycle=2, code="Null", is_null=True)
        
        assert bb.cycle == 2
        assert bb.code == "Null"
        assert bb.is_null is True

    def test_from_code_monomer(self):
        """Factory method with monomer code."""
        bb = BuildingBlock.from_code(0, "Leu")
        
        assert bb.cycle == 0
        assert bb.code == "Leu"
        assert bb.is_null is False

    def test_from_code_null_standard(self):
        """Factory method with standard Null code."""
        bb = BuildingBlock.from_code(1, "Null")
        
        assert bb.cycle == 1
        assert bb.code == "Null"
        assert bb.is_null is True

    def test_from_code_null_case_insensitive(self):
        """Null detection is case-insensitive (THEORY.md Section 2.1)."""
        test_cases = ["null", "NULL", "Null", "null_variant", "NullBlock"]
        
        for code in test_cases:
            bb = BuildingBlock.from_code(0, code)
            assert bb.is_null is True, f"Failed for code: {code}"

    def test_immutable(self):
        """BuildingBlock is immutable (frozen dataclass)."""
        bb = BuildingBlock.from_code(0, "Leu")
        
        with pytest.raises(AttributeError):
            bb.cycle = 1
        
        with pytest.raises(AttributeError):
            bb.code = "Ala"


class TestBuildingBlockValidation:
    """Test BuildingBlock validation rules."""

    def test_negative_cycle_raises_error(self):
        """Cycle must be non-negative."""
        with pytest.raises(ValueError, match="Cycle must be non-negative"):
            BuildingBlock(cycle=-1, code="Leu", is_null=False)

    def test_empty_code_raises_error(self):
        """Code cannot be empty."""
        with pytest.raises(ValueError, match="code cannot be empty"):
            BuildingBlock(cycle=0, code="", is_null=False)


class TestMonomerDecomposition:
    """Test decompose_to_monomers() - THEORY.md Section 1.5.3."""

    def test_decompose_single_monomer(self):
        """Single monomer returns single-element list."""
        bb = BuildingBlock.from_code(0, "Leu")
        
        monomers = bb.decompose_to_monomers()
        assert monomers == ["Leu"]

    def test_decompose_composite_dimeric(self):
        """Dimeric composite block splits into two monomers."""
        bb = BuildingBlock.from_code(1, "Leu-Ala")
        
        monomers = bb.decompose_to_monomers()
        assert monomers == ["Leu", "Ala"]

    def test_decompose_composite_trimeric(self):
        """Trimeric composite block splits into three monomers."""
        bb = BuildingBlock.from_code(1, "Leu-Ala-Val")
        
        monomers = bb.decompose_to_monomers()
        assert monomers == ["Leu", "Ala", "Val"]

    def test_decompose_null_block(self):
        """Null block returns empty list (THEORY.md Section 1.5.3)."""
        bb = BuildingBlock.from_code(2, "Null")
        
        monomers = bb.decompose_to_monomers()
        assert monomers == []

    def test_decompose_handles_whitespace(self):
        """Decomposition strips whitespace from monomers."""
        bb = BuildingBlock.from_code(1, "Leu - Ala - Val")
        
        monomers = bb.decompose_to_monomers()
        assert monomers == ["Leu", "Ala", "Val"]

    def test_decompose_modified_residues(self):
        """Works with modified residues (D-amino acids, methylated, etc.)."""
        bb = BuildingBlock.from_code(0, "DLeuMe-DPro-Leu")
        
        monomers = bb.decompose_to_monomers()
        assert monomers == ["DLeuMe", "DPro", "Leu"]


class TestBuildingBlockEquality:
    """Test equality and hashing (frozen dataclass behavior)."""

    def test_equality_same_values(self):
        """BuildingBlocks with same values are equal."""
        bb1 = BuildingBlock(cycle=0, code="Leu", is_null=False)
        bb2 = BuildingBlock(cycle=0, code="Leu", is_null=False)
        
        assert bb1 == bb2
        assert hash(bb1) == hash(bb2)

    def test_inequality_different_cycle(self):
        """Different cycle → not equal."""
        bb1 = BuildingBlock(cycle=0, code="Leu", is_null=False)
        bb2 = BuildingBlock(cycle=1, code="Leu", is_null=False)
        
        assert bb1 != bb2

    def test_inequality_different_code(self):
        """Different code → not equal."""
        bb1 = BuildingBlock(cycle=0, code="Leu", is_null=False)
        bb2 = BuildingBlock(cycle=0, code="Ala", is_null=False)
        
        assert bb1 != bb2

    def test_hashable_can_use_in_set(self):
        """BuildingBlock is hashable (can use in sets/dicts)."""
        bb1 = BuildingBlock.from_code(0, "Leu")
        bb2 = BuildingBlock.from_code(1, "Ala")
        bb3 = BuildingBlock.from_code(0, "Leu")  # Same as bb1
        
        building_blocks = {bb1, bb2, bb3}
        
        # bb1 and bb3 are same, so set has 2 elements
        assert len(building_blocks) == 2


class TestBuildingBlockStringRepresentation:
    """Test string conversion methods."""

    def test_str_returns_code(self):
        """str() returns the building block code."""
        bb = BuildingBlock.from_code(0, "Leu")
        
        assert str(bb) == "Leu"

    def test_str_composite_block(self):
        """str() preserves composite block notation."""
        bb = BuildingBlock.from_code(1, "Leu-Ala-Val")
        
        assert str(bb) == "Leu-Ala-Val"

    def test_repr_shows_all_fields(self):
        """repr() shows cycle, code, and is_null."""
        bb = BuildingBlock.from_code(0, "Leu")
        
        repr_str = repr(bb)
        assert "cycle=0" in repr_str
        assert "code='Leu'" in repr_str
        assert "is_null=False" in repr_str


class TestNSequenceConvention:
    """Test N→C sequence convention (THEORY.md Section 1.5.1)."""

    def test_position_0_is_c_terminus(self):
        """
        Position 0 = C-terminus (rightmost, synthesized first).
        
        This is documented behavior, not enforced by BuildingBlock,
        but important for library design.
        """
        # Position 0 should be C-terminus building block
        c_terminus_block = BuildingBlock.from_code(cycle=0, code="Pro")
        
        # Position 8 would be N-terminus (synthesized last)
        n_terminus_block = BuildingBlock.from_code(cycle=8, code="Leu")
        
        # Just verify they can be created with correct cycle assignments
        assert c_terminus_block.cycle == 0
        assert n_terminus_block.cycle == 8


class TestPropertyBasedDecomposition:
    """Property-based tests for decomposition invariants."""

    def test_decomposition_preserves_order(self):
        """
        Decomposed monomers maintain N→C order.
        
        For "Leu-Ala-Val", first monomer is N-terminus (Leu),
        last monomer is C-terminus (Val).
        """
        bb = BuildingBlock.from_code(1, "Leu-Ala-Val")
        monomers = bb.decompose_to_monomers()
        
        # Order preserved: N→C
        assert monomers[0] == "Leu"  # N-terminus
        assert monomers[-1] == "Val"  # C-terminus

    def test_null_block_always_empty_decomposition(self):
        """Any block with 'null' in code decomposes to empty list."""
        null_variants = ["Null", "null", "NULL_BLOCK", "null_v1"]
        
        for code in null_variants:
            bb = BuildingBlock.from_code(0, code)
            assert bb.decompose_to_monomers() == []
