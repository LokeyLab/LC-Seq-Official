"""
Comprehensive tests for EquivalenceClass model.

Tests grouping of positional variants by chemical identity.
"""

import pytest
import numpy as np
from lcseq.domain.models.equivalence_class import EquivalenceClass
from lcseq.domain.entities.compound import Compound
from lcseq.domain.entities.building_block import BuildingBlock
from lcseq.domain.entities.chromatogram import Chromatogram


@pytest.fixture
def sample_chromatogram():
    """Create sample chromatogram for testing."""
    return Chromatogram(
        time_points=np.array([1.0, 2.0, 3.0, 4.0, 5.0]),
        counts=np.array([100.0, 200.0, 300.0, 200.0, 100.0]),
    )


class TestEquivalenceClassCreation:
    """Test equivalence class creation and initialization."""

    def test_create_empty_equivalence_class(self):
        """Test creating empty equivalence class."""
        eq_class = EquivalenceClass(residue_sequence="Val-Pro")

        assert eq_class.residue_sequence == "Val-Pro"
        assert len(eq_class.compounds) == 0
        assert eq_class.is_empty()

    def test_create_with_compounds(self, sample_chromatogram):
        """Test creating equivalence class with initial compounds."""
        bb0 = BuildingBlock.from_code(0, "Val")
        bb1 = BuildingBlock.from_code(1, "Null")
        bb2 = BuildingBlock.from_code(2, "Null")
        compound = Compound([bb0, bb1, bb2], sample_chromatogram)

        eq_class = EquivalenceClass(residue_sequence="Val", compounds={compound})

        assert eq_class.residue_sequence == "Val"
        assert len(eq_class.compounds) == 1
        assert compound in eq_class.compounds

    def test_create_with_mismatched_residue_sequence_raises_error(self, sample_chromatogram):
        """Test creating with mismatched residue sequence raises error."""
        bb0 = BuildingBlock.from_code(0, "Val")
        bb1 = BuildingBlock.from_code(1, "Pro")
        compound = Compound([bb0, bb1], sample_chromatogram)

        # Compound has residue_sequence "Pro-Val" but class is "Leu-Pro"
        with pytest.raises(ValueError, match="does not match class residue sequence"):
            EquivalenceClass(residue_sequence="Leu-Pro", compounds={compound})


class TestAddCompound:
    """Test adding compounds to equivalence class."""

    def test_add_matching_compound(self, sample_chromatogram):
        """Test adding compound with matching residue sequence."""
        eq_class = EquivalenceClass(residue_sequence="Val")

        bb0 = BuildingBlock.from_code(0, "Val")
        bb1 = BuildingBlock.from_code(1, "Null")
        compound = Compound([bb0, bb1], sample_chromatogram)

        eq_class.add_compound(compound)

        assert compound in eq_class.compounds
        assert len(eq_class.compounds) == 1

    def test_add_multiple_positional_variants(self, sample_chromatogram):
        """Test adding multiple positional variants of same peptide."""
        eq_class = EquivalenceClass(residue_sequence="Val")

        # Three positional variants: [Val, Null, Null], [Null, Val, Null], [Null, Null, Val]
        bb_val_0 = BuildingBlock.from_code(0, "Val")
        bb_val_1 = BuildingBlock.from_code(1, "Val")
        bb_val_2 = BuildingBlock.from_code(2, "Val")
        bb_null_0 = BuildingBlock.from_code(0, "Null")
        bb_null_1 = BuildingBlock.from_code(1, "Null")
        bb_null_2 = BuildingBlock.from_code(2, "Null")

        variant1 = Compound([bb_val_0, bb_null_1, bb_null_2], sample_chromatogram)
        variant2 = Compound([bb_null_0, bb_val_1, bb_null_2], sample_chromatogram)
        variant3 = Compound([bb_null_0, bb_null_1, bb_val_2], sample_chromatogram)

        eq_class.add_compound(variant1)
        eq_class.add_compound(variant2)
        eq_class.add_compound(variant3)

        assert len(eq_class.compounds) == 3
        assert variant1 in eq_class.compounds
        assert variant2 in eq_class.compounds
        assert variant3 in eq_class.compounds

    def test_add_compound_with_mismatched_residue_sequence_raises_error(self, sample_chromatogram):
        """Test adding compound with wrong residue sequence raises error."""
        eq_class = EquivalenceClass(residue_sequence="Val")

        bb0 = BuildingBlock.from_code(0, "Leu")  # Different residue
        bb1 = BuildingBlock.from_code(1, "Null")
        compound = Compound([bb0, bb1], sample_chromatogram)

        with pytest.raises(ValueError, match="Cannot add compound with residue sequence"):
            eq_class.add_compound(compound)

    def test_add_compound_idempotent(self, sample_chromatogram):
        """Test that adding same compound multiple times is safe."""
        eq_class = EquivalenceClass(residue_sequence="Val")

        bb0 = BuildingBlock.from_code(0, "Val")
        bb1 = BuildingBlock.from_code(1, "Null")
        compound = Compound([bb0, bb1], sample_chromatogram)

        eq_class.add_compound(compound)
        eq_class.add_compound(compound)
        eq_class.add_compound(compound)

        assert len(eq_class.compounds) == 1


class TestGetRepresentative:
    """Test getting representative compound."""

    def test_get_representative_from_empty_class(self):
        """Test getting representative from empty class returns None."""
        eq_class = EquivalenceClass(residue_sequence="Val")

        representative = eq_class.get_representative()
        assert representative is None

    def test_get_representative_single_compound(self, sample_chromatogram):
        """Test getting representative with single compound."""
        eq_class = EquivalenceClass(residue_sequence="Val")

        bb0 = BuildingBlock.from_code(0, "Val")
        bb1 = BuildingBlock.from_code(1, "Null")
        compound = Compound([bb0, bb1], sample_chromatogram)
        eq_class.add_compound(compound)

        representative = eq_class.get_representative()
        assert representative == compound

    def test_get_representative_multiple_compounds_consistent(self, sample_chromatogram):
        """Test that representative is consistent across calls."""
        eq_class = EquivalenceClass(residue_sequence="Val")

        # Add variants in random order
        bb_val_0 = BuildingBlock.from_code(0, "Val")
        bb_val_1 = BuildingBlock.from_code(1, "Val")
        bb_val_2 = BuildingBlock.from_code(2, "Val")
        bb_null_0 = BuildingBlock.from_code(0, "Null")
        bb_null_1 = BuildingBlock.from_code(1, "Null")
        bb_null_2 = BuildingBlock.from_code(2, "Null")

        variant1 = Compound([bb_val_0, bb_null_1, bb_null_2], sample_chromatogram)
        variant2 = Compound([bb_null_0, bb_val_1, bb_null_2], sample_chromatogram)
        variant3 = Compound([bb_null_0, bb_null_1, bb_val_2], sample_chromatogram)

        eq_class.add_compound(variant2)
        eq_class.add_compound(variant1)
        eq_class.add_compound(variant3)

        # Should be consistent (lexicographically first)
        rep1 = eq_class.get_representative()
        rep2 = eq_class.get_representative()
        rep3 = eq_class.get_representative()

        assert rep1 == rep2 == rep3
        assert rep1.residue_sequence == "Val"

    def test_get_representative_lexicographic_ordering(self, sample_chromatogram):
        """Test that representative uses lexicographic ordering."""
        eq_class = EquivalenceClass(residue_sequence="Leu-Pro")

        # Create two variants with different positional sequences
        # variant1: [cycle0=Leu, cycle1=Pro, cycle2=Null] -> N→C: "Null-Pro-Leu" -> residue: "Pro-Leu"
        # variant2: [cycle0=Null, cycle1=Leu, cycle2=Pro] -> N→C: "Pro-Leu-Null" -> residue: "Pro-Leu"
        # We need to create "Leu-Pro" residue sequence instead

        bb_leu_0 = BuildingBlock.from_code(0, "Leu")
        bb_null_0 = BuildingBlock.from_code(0, "Null")
        bb_pro_1 = BuildingBlock.from_code(1, "Pro")
        bb_null_1 = BuildingBlock.from_code(1, "Null")
        bb_null_2 = BuildingBlock.from_code(2, "Null")

        # Positional: "Null-Pro-Leu" -> residue: "Pro-Leu" (WRONG, we need "Leu-Pro")
        # Let's fix to create proper "Leu-Pro" residue:
        # [cycle0=Pro, cycle1=Leu, cycle2=Null] -> N→C: "Null-Leu-Pro" -> residue: "Leu-Pro"
        # [cycle0=Null, cycle1=Pro, cycle2=Leu] -> N→C: "Leu-Pro-Null" -> residue: "Leu-Pro"

        bb_pro_0 = BuildingBlock.from_code(0, "Pro")
        bb_leu_1 = BuildingBlock.from_code(1, "Leu")
        bb_leu_2 = BuildingBlock.from_code(2, "Leu")
        bb_pro_2 = BuildingBlock.from_code(2, "Pro")

        # Positional: "Null-Leu-Pro" (N→C order)
        variant1 = Compound([bb_pro_0, bb_leu_1, bb_null_2], sample_chromatogram)
        # Positional: "Leu-Pro-Null" (N→C order)
        variant2 = Compound([bb_null_0, bb_pro_1, bb_leu_2], sample_chromatogram)

        eq_class.add_compound(variant1)
        eq_class.add_compound(variant2)

        rep = eq_class.get_representative()
        # "Leu-Pro-Null" comes before "Null-Leu-Pro" lexicographically
        assert rep.positional_sequence == "Leu-Pro-Null"


class TestSizeAndEmpty:
    """Test size and empty checks."""

    def test_size_empty_class(self):
        """Test size of empty class."""
        eq_class = EquivalenceClass(residue_sequence="Val")

        assert eq_class.size() == 0
        assert len(eq_class) == 0
        assert eq_class.is_empty()

    def test_size_with_compounds(self, sample_chromatogram):
        """Test size with multiple compounds."""
        eq_class = EquivalenceClass(residue_sequence="Val")

        bb_val_0 = BuildingBlock.from_code(0, "Val")
        bb_val_1 = BuildingBlock.from_code(1, "Val")
        bb_null_0 = BuildingBlock.from_code(0, "Null")
        bb_null_1 = BuildingBlock.from_code(1, "Null")

        variant1 = Compound([bb_val_0, bb_null_1], sample_chromatogram)
        variant2 = Compound([bb_null_0, bb_val_1], sample_chromatogram)

        eq_class.add_compound(variant1)
        eq_class.add_compound(variant2)

        assert eq_class.size() == 2
        assert len(eq_class) == 2
        assert not eq_class.is_empty()


class TestGetPositionalSequences:
    """Test getting positional sequences."""

    def test_get_positional_sequences_empty_class(self):
        """Test getting positional sequences from empty class."""
        eq_class = EquivalenceClass(residue_sequence="Val")

        sequences = eq_class.get_positional_sequences()
        assert len(sequences) == 0

    def test_get_positional_sequences_single_compound(self, sample_chromatogram):
        """Test getting positional sequences with single compound."""
        eq_class = EquivalenceClass(residue_sequence="Val")

        bb0 = BuildingBlock.from_code(0, "Val")
        bb1 = BuildingBlock.from_code(1, "Null")
        compound = Compound([bb0, bb1], sample_chromatogram)
        eq_class.add_compound(compound)

        sequences = eq_class.get_positional_sequences()
        assert len(sequences) == 1
        assert "Null-Val" in sequences

    def test_get_positional_sequences_multiple_variants(self, sample_chromatogram):
        """Test getting positional sequences with multiple variants."""
        eq_class = EquivalenceClass(residue_sequence="Val")

        bb_val_0 = BuildingBlock.from_code(0, "Val")
        bb_val_1 = BuildingBlock.from_code(1, "Val")
        bb_val_2 = BuildingBlock.from_code(2, "Val")
        bb_null_0 = BuildingBlock.from_code(0, "Null")
        bb_null_1 = BuildingBlock.from_code(1, "Null")
        bb_null_2 = BuildingBlock.from_code(2, "Null")

        variant1 = Compound([bb_val_0, bb_null_1, bb_null_2], sample_chromatogram)
        variant2 = Compound([bb_null_0, bb_val_1, bb_null_2], sample_chromatogram)
        variant3 = Compound([bb_null_0, bb_null_1, bb_val_2], sample_chromatogram)

        eq_class.add_compound(variant1)
        eq_class.add_compound(variant2)
        eq_class.add_compound(variant3)

        sequences = eq_class.get_positional_sequences()
        assert len(sequences) == 3
        assert "Null-Null-Val" in sequences
        assert "Null-Val-Null" in sequences
        assert "Val-Null-Null" in sequences


class TestStringRepresentations:
    """Test string representations."""

    def test_repr(self, sample_chromatogram):
        """Test repr representation."""
        eq_class = EquivalenceClass(residue_sequence="Leu-Pro")

        # [cycle0=Pro, cycle1=Leu] -> N→C: "Leu-Pro" -> residue: "Leu-Pro"
        bb_pro = BuildingBlock.from_code(0, "Pro")
        bb_leu = BuildingBlock.from_code(1, "Leu")
        compound = Compound([bb_pro, bb_leu], sample_chromatogram)
        eq_class.add_compound(compound)

        repr_str = repr(eq_class)
        assert "EquivalenceClass" in repr_str
        assert "Leu-Pro" in repr_str
        assert "variants=1" in repr_str

    def test_str(self, sample_chromatogram):
        """Test str representation."""
        eq_class = EquivalenceClass(residue_sequence="Leu-Pro")

        # [cycle0=Pro, cycle1=Leu, cycle2=Null] -> N→C: "Null-Leu-Pro" -> residue: "Leu-Pro"
        bb_pro_0 = BuildingBlock.from_code(0, "Pro")
        bb_leu_1 = BuildingBlock.from_code(1, "Leu")
        bb_null_2 = BuildingBlock.from_code(2, "Null")

        variant1 = Compound([bb_pro_0, bb_leu_1, bb_null_2], sample_chromatogram)

        # [cycle0=Null, cycle1=Pro, cycle2=Leu] -> N→C: "Leu-Pro-Null" -> residue: "Leu-Pro"
        bb_null_0 = BuildingBlock.from_code(0, "Null")
        bb_pro_1 = BuildingBlock.from_code(1, "Pro")
        bb_leu_2 = BuildingBlock.from_code(2, "Leu")
        variant2 = Compound([bb_null_0, bb_pro_1, bb_leu_2], sample_chromatogram)

        eq_class.add_compound(variant1)
        eq_class.add_compound(variant2)

        str_repr = str(eq_class)
        assert "[Leu-Pro]" in str_repr
        assert "2 variants" in str_repr


class TestEquivalenceRelationProperties:
    """Test that equivalence class satisfies equivalence relation properties."""

    def test_reflexive_property(self, sample_chromatogram):
        """Test reflexive: compound related to itself."""
        bb0 = BuildingBlock.from_code(0, "Val")
        bb1 = BuildingBlock.from_code(1, "Null")
        compound = Compound([bb0, bb1], sample_chromatogram)

        eq_class = EquivalenceClass(residue_sequence="Val")
        eq_class.add_compound(compound)

        # Compound should be in its own equivalence class
        assert compound in eq_class.compounds

    def test_symmetric_property(self, sample_chromatogram):
        """Test symmetric: if A relates to B, then B relates to A."""
        bb_val_0 = BuildingBlock.from_code(0, "Val")
        bb_val_1 = BuildingBlock.from_code(1, "Val")
        bb_null_0 = BuildingBlock.from_code(0, "Null")
        bb_null_1 = BuildingBlock.from_code(1, "Null")

        variant1 = Compound([bb_val_0, bb_null_1], sample_chromatogram)
        variant2 = Compound([bb_null_0, bb_val_1], sample_chromatogram)

        eq_class = EquivalenceClass(residue_sequence="Val")
        eq_class.add_compound(variant1)
        eq_class.add_compound(variant2)

        # Both should be in same class (symmetric)
        assert variant1 in eq_class.compounds
        assert variant2 in eq_class.compounds

    def test_transitive_property(self, sample_chromatogram):
        """Test transitive: if A~B and B~C, then A~C."""
        bb_val_0 = BuildingBlock.from_code(0, "Val")
        bb_val_1 = BuildingBlock.from_code(1, "Val")
        bb_val_2 = BuildingBlock.from_code(2, "Val")
        bb_null_0 = BuildingBlock.from_code(0, "Null")
        bb_null_1 = BuildingBlock.from_code(1, "Null")
        bb_null_2 = BuildingBlock.from_code(2, "Null")

        variantA = Compound([bb_val_0, bb_null_1, bb_null_2], sample_chromatogram)
        variantB = Compound([bb_null_0, bb_val_1, bb_null_2], sample_chromatogram)
        variantC = Compound([bb_null_0, bb_null_1, bb_val_2], sample_chromatogram)

        eq_class = EquivalenceClass(residue_sequence="Val")
        eq_class.add_compound(variantA)
        eq_class.add_compound(variantB)
        eq_class.add_compound(variantC)

        # All three should be in same class (transitive)
        assert variantA in eq_class.compounds
        assert variantB in eq_class.compounds
        assert variantC in eq_class.compounds
