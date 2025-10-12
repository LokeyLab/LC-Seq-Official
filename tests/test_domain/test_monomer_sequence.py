"""
Tests for MonomerSequence value object.

Tests implementation against THEORY.md Section 2.2, 1.5.3, 1.2 specifications.
"""

import pytest
from lcseq.domain.value_objects.monomer_sequence import MonomerSequence


class TestMonomerSequenceCreation:
    """Test MonomerSequence instantiation and factory methods."""

    def test_create_from_tuple(self):
        """Direct creation with monomer tuple."""
        seq = MonomerSequence(monomers=("Leu", "Ala", "Pro"))

        assert len(seq) == 3
        assert seq.monomers == ("Leu", "Ala", "Pro")

    def test_from_list(self):
        """Factory method from list."""
        seq = MonomerSequence.from_list(["Leu", "Ala", "Pro"])

        assert len(seq) == 3
        assert seq.monomers == ("Leu", "Ala", "Pro")

    def test_from_string_simple(self):
        """Factory method from hyphen-separated string."""
        seq = MonomerSequence.from_string("Leu-Ala-Pro")

        assert len(seq) == 3
        assert seq.to_string() == "Leu-Ala-Pro"

    def test_from_string_single_monomer(self):
        """Single monomer string."""
        seq = MonomerSequence.from_string("Leu")

        assert len(seq) == 1
        assert seq.to_string() == "Leu"

    def test_from_string_empty_produces_empty_sequence(self):
        """Empty string produces empty sequence."""
        seq = MonomerSequence.from_string("")

        assert seq.is_empty()
        assert len(seq) == 0

    def test_from_string_whitespace_produces_empty_sequence(self):
        """Whitespace-only string produces empty sequence."""
        seq = MonomerSequence.from_string("   ")

        assert seq.is_empty()

    def test_from_string_strips_whitespace(self):
        """Whitespace around monomers is stripped."""
        seq = MonomerSequence.from_string(" Leu - Ala - Pro ")

        assert seq.to_string() == "Leu-Ala-Pro"

    def test_empty_sequence(self):
        """Empty sequence (e.g., from all null blocks)."""
        seq = MonomerSequence.from_list([])

        assert len(seq) == 0
        assert seq.is_empty()

    def test_immutable(self):
        """MonomerSequence is immutable (frozen dataclass)."""
        seq = MonomerSequence.from_list(["Leu", "Ala"])

        with pytest.raises(AttributeError):
            seq.monomers = ("Pro",)


class TestMonomerSequenceValidation:
    """Test MonomerSequence validation rules."""

    def test_empty_monomer_code_raises_error(self):
        """Monomer code cannot be empty."""
        with pytest.raises(ValueError, match="cannot be empty"):
            MonomerSequence(monomers=("Leu", "", "Pro"))

    def test_whitespace_only_monomer_raises_error(self):
        """Monomer code cannot be whitespace only."""
        with pytest.raises(ValueError, match="cannot be empty or whitespace"):
            MonomerSequence(monomers=("Leu", "   ", "Pro"))


class TestSequenceAccess:
    """Test accessing monomers in sequence."""

    def test_get_n_terminus(self):
        """N-terminus is first monomer (THEORY.md N→C convention)."""
        seq = MonomerSequence.from_list(["Leu", "Ala", "Pro"])

        assert seq.get_n_terminus() == "Leu"

    def test_get_c_terminus(self):
        """C-terminus is last monomer."""
        seq = MonomerSequence.from_list(["Leu", "Ala", "Pro"])

        assert seq.get_c_terminus() == "Pro"

    def test_get_n_terminus_empty_raises_error(self):
        """Empty sequence has no N-terminus."""
        seq = MonomerSequence.from_list([])

        with pytest.raises(IndexError):
            seq.get_n_terminus()

    def test_get_c_terminus_empty_raises_error(self):
        """Empty sequence has no C-terminus."""
        seq = MonomerSequence.from_list([])

        with pytest.raises(IndexError):
            seq.get_c_terminus()

    def test_get_monomer_at_valid_position(self):
        """Get monomer at valid position (0-indexed from N-terminus)."""
        seq = MonomerSequence.from_list(["Leu", "Ala", "Pro"])

        assert seq.get_monomer_at(0) == "Leu"  # N-terminus
        assert seq.get_monomer_at(1) == "Ala"
        assert seq.get_monomer_at(2) == "Pro"  # C-terminus

    def test_get_monomer_at_negative_index(self):
        """Negative indexing works (Python style)."""
        seq = MonomerSequence.from_list(["Leu", "Ala", "Pro"])

        assert seq.get_monomer_at(-1) == "Pro"  # C-terminus
        assert seq.get_monomer_at(-2) == "Ala"
        assert seq.get_monomer_at(-3) == "Leu"  # N-terminus

    def test_get_monomer_at_invalid_position_raises_error(self):
        """Invalid position raises IndexError."""
        seq = MonomerSequence.from_list(["Leu", "Ala"])

        with pytest.raises(IndexError):
            seq.get_monomer_at(5)

    def test_iteration(self):
        """Can iterate over monomers (N→C order)."""
        seq = MonomerSequence.from_list(["Leu", "Ala", "Pro"])

        monomers = list(seq)
        assert monomers == ["Leu", "Ala", "Pro"]

    def test_indexing(self):
        """Can index directly into sequence."""
        seq = MonomerSequence.from_list(["Leu", "Ala", "Pro"])

        assert seq[0] == "Leu"
        assert seq[1] == "Ala"
        assert seq[2] == "Pro"

    def test_slicing(self):
        """Can slice sequence."""
        seq = MonomerSequence.from_list(["Leu", "Ala", "Pro", "Gly"])

        sub = seq[1:3]
        assert isinstance(sub, MonomerSequence)
        assert sub.to_string() == "Ala-Pro"


class TestSequenceOperations:
    """Test sequence operations and transformations."""

    def test_to_string(self):
        """Convert to hyphen-separated string."""
        seq = MonomerSequence.from_list(["Leu", "Ala", "Pro"])

        assert seq.to_string() == "Leu-Ala-Pro"

    def test_to_string_empty(self):
        """Empty sequence produces empty string."""
        seq = MonomerSequence.from_list([])

        assert seq.to_string() == ""

    def test_is_empty_true(self):
        """Empty sequence returns True."""
        seq = MonomerSequence.from_list([])

        assert seq.is_empty() is True

    def test_is_empty_false(self):
        """Non-empty sequence returns False."""
        seq = MonomerSequence.from_list(["Leu"])

        assert seq.is_empty() is False

    def test_reverse(self):
        """Reverse creates new sequence with reversed order."""
        seq = MonomerSequence.from_list(["Leu", "Ala", "Pro"])

        rev = seq.reverse()

        assert rev.to_string() == "Pro-Ala-Leu"
        # Original unchanged
        assert seq.to_string() == "Leu-Ala-Pro"

    def test_reverse_empty(self):
        """Reversing empty sequence produces empty sequence."""
        seq = MonomerSequence.from_list([])

        rev = seq.reverse()
        assert rev.is_empty()

    def test_reverse_single(self):
        """Reversing single monomer produces same sequence."""
        seq = MonomerSequence.from_list(["Leu"])

        rev = seq.reverse()
        assert rev.to_string() == "Leu"


class TestEquality:
    """Test equality and hashing (THEORY.md Section 1.2: chemical identity)."""

    def test_equality_same_monomers(self):
        """Sequences with same monomers are equal."""
        seq1 = MonomerSequence.from_list(["Leu", "Ala", "Pro"])
        seq2 = MonomerSequence.from_list(["Leu", "Ala", "Pro"])

        assert seq1 == seq2
        assert hash(seq1) == hash(seq2)

    def test_inequality_different_monomers(self):
        """Different monomers → not equal."""
        seq1 = MonomerSequence.from_list(["Leu", "Ala", "Pro"])
        seq2 = MonomerSequence.from_list(["Leu", "Ala", "Gly"])

        assert seq1 != seq2

    def test_inequality_different_order(self):
        """Different order → not equal (N→C matters)."""
        seq1 = MonomerSequence.from_list(["Leu", "Ala", "Pro"])
        seq2 = MonomerSequence.from_list(["Pro", "Ala", "Leu"])

        assert seq1 != seq2

    def test_inequality_different_length(self):
        """Different length → not equal."""
        seq1 = MonomerSequence.from_list(["Leu", "Ala", "Pro"])
        seq2 = MonomerSequence.from_list(["Leu", "Ala"])

        assert seq1 != seq2

    def test_empty_sequences_equal(self):
        """Empty sequences are equal."""
        seq1 = MonomerSequence.from_list([])
        seq2 = MonomerSequence.from_list([])

        assert seq1 == seq2

    def test_hashable(self):
        """MonomerSequence is hashable (can use in sets/dicts)."""
        seq1 = MonomerSequence.from_list(["Leu", "Ala"])
        seq2 = MonomerSequence.from_list(["Leu", "Pro"])
        seq3 = MonomerSequence.from_list(["Leu", "Ala"])  # Same as seq1

        sequences = {seq1, seq2, seq3}

        # seq1 and seq3 are same, so set has 2 elements
        assert len(sequences) == 2


class TestChemicalIdentity:
    """Test chemical identity concept (THEORY.md Section 1.2)."""

    def test_same_chemical_identity(self):
        """
        Same monomer sequence = same chemical identity.

        Different positional sequences can produce same monomer sequence
        (same molecule).
        """
        # These would come from different positional sequences:
        # [Val, Null, Null] → "Val"
        # [Null, Val, Null] → "Val"
        # [Null, Null, Val] → "Val"

        seq1 = MonomerSequence.from_list(["Val"])
        seq2 = MonomerSequence.from_list(["Val"])

        # Same chemical identity
        assert seq1 == seq2

    def test_different_chemical_identity(self):
        """Different monomer sequences = different molecules."""
        seq1 = MonomerSequence.from_list(["Leu", "Ala"])
        seq2 = MonomerSequence.from_list(["Ala", "Leu"])

        # Different chemical identities (different molecules)
        assert seq1 != seq2


class TestStringRepresentation:
    """Test string conversion methods."""

    def test_str_returns_sequence_string(self):
        """str() returns hyphen-separated sequence."""
        seq = MonomerSequence.from_list(["Leu", "Ala", "Pro"])

        assert str(seq) == "Leu-Ala-Pro"

    def test_repr_shows_sequence(self):
        """repr() shows sequence in quotes."""
        seq = MonomerSequence.from_list(["Leu", "Ala", "Pro"])

        repr_str = repr(seq)
        assert "Leu-Ala-Pro" in repr_str
        assert "MonomerSequence" in repr_str

    def test_str_empty_sequence(self):
        """str() of empty sequence is empty string."""
        seq = MonomerSequence.from_list([])

        assert str(seq) == ""


class TestModifiedResidues:
    """Test sequences with modified residues."""

    def test_d_amino_acids(self):
        """D-amino acids work correctly."""
        seq = MonomerSequence.from_list(["DLeu", "DAla", "Pro"])

        assert seq.to_string() == "DLeu-DAla-Pro"
        assert seq.get_n_terminus() == "DLeu"

    def test_methylated_residues(self):
        """Methylated residues work correctly."""
        seq = MonomerSequence.from_list(["NMeLeu", "Ala", "Pro"])

        assert seq.to_string() == "NMeLeu-Ala-Pro"

    def test_mixed_modifications(self):
        """Mixed modifications work correctly."""
        seq = MonomerSequence.from_list(["DLeuMe", "DPro", "Leu"])

        assert len(seq) == 3
        assert seq.get_monomer_at(0) == "DLeuMe"


class TestLongSequences:
    """Test handling of long sequences."""

    def test_long_sequence(self):
        """Handle long sequences."""
        monomers = [f"M{i}" for i in range(100)]
        seq = MonomerSequence.from_list(monomers)

        assert len(seq) == 100
        assert seq.get_n_terminus() == "M0"
        assert seq.get_c_terminus() == "M99"

    def test_long_sequence_iteration(self):
        """Can iterate over long sequences."""
        monomers = [f"M{i}" for i in range(50)]
        seq = MonomerSequence.from_list(monomers)

        count = 0
        for monomer in seq:
            count += 1

        assert count == 50


class TestEdgeCases:
    """Test edge cases and boundary conditions."""

    def test_single_monomer_sequence(self):
        """Single monomer sequence."""
        seq = MonomerSequence.from_list(["Leu"])

        assert len(seq) == 1
        assert seq.get_n_terminus() == "Leu"
        assert seq.get_c_terminus() == "Leu"
        assert seq.to_string() == "Leu"

    def test_unusual_monomer_codes(self):
        """Unusual but valid monomer codes."""
        seq = MonomerSequence.from_list(["Leu123", "Ala-modified", "X"])

        assert len(seq) == 3
        assert seq.to_string() == "Leu123-Ala-modified-X"

    def test_case_sensitive(self):
        """Monomer codes are case-sensitive."""
        seq1 = MonomerSequence.from_list(["Leu", "ALA"])
        seq2 = MonomerSequence.from_list(["Leu", "ala"])

        assert seq1 != seq2

    def test_tuple_immutability(self):
        """Monomers stored as immutable tuple."""
        seq = MonomerSequence.from_list(["Leu", "Ala"])

        assert isinstance(seq.monomers, tuple)
