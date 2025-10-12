"""Tests for truncation generation and hierarchy."""

import numpy as np
import pytest

from lcseq.domain.entities.chromatogram import Chromatogram
from lcseq.domain.entities.compound import Compound
from lcseq.domain.entities.building_block import BuildingBlock


class TestBuildingBlockNullDetection:
    """Test automatic null marker detection."""

    def test_agxnull_marker_detected(self):
        """Test that 'AgxNull' is automatically detected as null."""
        bb = BuildingBlock.from_code(cycle=0, code="AgxNull")
        assert bb.is_null is True

    def test_agxnull_case_insensitive(self):
        """Test that null detection is case-insensitive."""
        bb1 = BuildingBlock.from_code(cycle=0, code="agxnull")
        bb2 = BuildingBlock.from_code(cycle=0, code="AGXNULL")
        bb3 = BuildingBlock.from_code(cycle=0, code="AgXnUlL")
        assert bb1.is_null is True
        assert bb2.is_null is True
        assert bb3.is_null is True

    def test_null_marker_detected(self):
        """Test that 'NULL' is also detected as null."""
        bb = BuildingBlock.from_code(cycle=0, code="NULL")
        assert bb.is_null is True

    def test_explicit_is_null_overrides(self):
        """Test that explicit is_null parameter overrides auto-detection."""
        bb = BuildingBlock(cycle=0, code="AgxNull", is_null=False)
        assert bb.is_null is False

    def test_null_factory_method(self):
        """Test creating null building blocks with factory method."""
        bb = BuildingBlock.null(2)
        assert bb.cycle == 2
        assert bb.code == "AgxNull"
        assert bb.is_null is True


class TestResidueSequenceFormatting:
    """Test residue sequence string formatting."""

    def test_full_length_residue_sequence(self):
        """Test formatting full-length compound."""
        time = np.array([1.0, 2.0, 3.0], dtype=np.float64)
        counts = np.array([10, 20, 30], dtype=np.float64)
        chrom = Chromatogram(time_points=time, counts=counts)

        blocks = [
            BuildingBlock.from_code(cycle=0, code="Val"),
            BuildingBlock.from_code(cycle=1, code="Nvl"),
            BuildingBlock.from_code(cycle=2, code="Leu"),
        ]

        compound = Compound(
            building_blocks=blocks, chromatogram=chrom, compound_id="SEQ001"
        )

        assert compound.residue_sequence == "Leu-Nvl-Val"  # N->C order

    def test_truncated_residue_sequence(self):
        """Test formatting truncated compound uses CANONICAL form (nulls filtered)."""
        time = np.array([1.0, 2.0, 3.0], dtype=np.float64)
        counts = np.array([10, 20, 30], dtype=np.float64)
        chrom = Chromatogram(time_points=time, counts=counts)

        blocks = [
            BuildingBlock.from_code(cycle=0, code="Val"),
            BuildingBlock.null(cycle=1),
            BuildingBlock.from_code(cycle=2, code="Leu"),
        ]

        compound = Compound(
            building_blocks=blocks, chromatogram=chrom, compound_id="SEQ002"
        )

        # Canonical representation filters out nulls (N->C order)
        assert compound.residue_sequence == "Leu-Val"

    def test_format_building_blocks_static(self):
        """Test residue sequence formatting."""
        blocks = [
            BuildingBlock.from_code(cycle=0, code="Ala"),
            BuildingBlock.from_code(cycle=1, code="Gly"),
        ]

        time = np.array([1.0, 2.0], dtype=np.float64)
        counts = np.array([10, 20], dtype=np.float64)
        chrom = Chromatogram(time_points=time, counts=counts)

        compound = Compound(building_blocks=blocks, chromatogram=chrom, compound_id="TEST")
        assert compound.residue_sequence == "Gly-Ala"  # N->C order


class TestTruncationRelationships:
    """Test truncation relationship detection."""

    @pytest.fixture
    def full_length_compound(self):
        """Fixture for full-length compound."""
        time = np.array([1.0, 2.0, 3.0], dtype=np.float64)
        counts = np.array([10, 20, 30], dtype=np.float64)
        chrom = Chromatogram(time_points=time, counts=counts)

        blocks = [
            BuildingBlock.from_code(cycle=0, code="Val"),
            BuildingBlock.from_code(cycle=1, code="Nvl"),
            BuildingBlock.from_code(cycle=2, code="Leu"),
        ]

        return Compound(
            building_blocks=blocks, chromatogram=chrom, compound_id="FULL"
        )

    @pytest.mark.skip(reason="Truncation relationship methods not implemented in current design")
    def test_truncation_relationship_detected(self, full_length_compound):
        """Test that truncation is correctly identified."""
        pass

    @pytest.mark.skip(reason="Truncation relationship methods not implemented in current design")
    def test_full_length_not_truncation_of_itself(self, full_length_compound):
        """Test that full-length is not a truncation of itself."""
        pass

    @pytest.mark.skip(reason="Truncation relationship methods not implemented in current design")
    def test_different_sequence_not_truncation(self, full_length_compound):
        """Test that different sequences are not truncations."""
        pass

    @pytest.mark.skip(reason="Truncation relationship methods not implemented in current design")
    def test_mismatched_non_null_positions(self, full_length_compound):
        """Test that mismatched positions are not truncations."""
        pass


class TestTruncationGeneration:
    """Test generating all truncation variants."""

    @pytest.mark.skip(reason="Truncation generation not implemented in current Compound design")
    def test_generate_all_truncations_3_cycles(self):
        """Test generating ALL POSITIONAL truncations from 3-cycle compound."""
        pass

    @pytest.mark.skip(reason="Truncation generation not implemented in current Compound design")
    def test_generate_truncations_include_self(self):
        """Test including self in truncations (all positional forms)."""
        pass

    @pytest.mark.skip(reason="Truncation generation not implemented in current Compound design")
    def test_cannot_generate_from_truncation(self):
        """Test that truncations cannot generate further truncations."""
        pass

    @pytest.mark.skip(reason="Truncation generation not implemented in current Compound design")
    def test_truncation_count_formula(self):
        """Test that truncation count follows 2^n - 2 formula."""
        pass

    @pytest.mark.skip(reason="Truncation generation not implemented in current Compound design")
    def test_truncation_with_all_null(self):
        """Test including the all-null compound (L0) in truncations."""
        pass


class TestHierarchyIntegration:
    """Test hierarchy with truncations."""

    @pytest.mark.skip(reason="Hierarchy methods not implemented in current Compound design")
    def test_build_hierarchy_from_truncations(self):
        """Test building ancestor-descendant relationships."""
        pass

    @pytest.mark.skip(reason="Truncation relationship methods not implemented in current Compound design")
    def test_verify_truncation_relationship_after_creation(self):
        """Test that created truncations validate correctly."""
        pass
