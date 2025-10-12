"""
Tests for Compound entity.

Tests implementation against THEORY.md Section 2.1, 2.2, 3.3, 1.5.1, 1.5.3 specifications.
"""

import pytest
import numpy as np
from lcseq.domain.entities.compound import Compound
from lcseq.domain.entities.building_block import BuildingBlock
from lcseq.domain.entities.chromatogram import Chromatogram
from lcseq.domain.entities.peak import Peak, PeakType, ValidationStatus


class TestCompoundCreation:
    """Test Compound instantiation and validation."""

    def test_create_simple_compound(self):
        """Basic compound with building blocks and chromatogram."""
        blocks = [
            BuildingBlock.from_code(0, "Pro"),
            BuildingBlock.from_code(1, "Leu")
        ]
        chrom = Chromatogram(
            time_points=np.array([0.0, 1.0, 2.0], dtype=np.float64),
            counts=np.array([10.0, 20.0, 15.0], dtype=np.float64)
        )

        compound = Compound(building_blocks=blocks, chromatogram=chrom)

        assert len(compound.building_blocks) == 2
        assert compound.chromatogram == chrom
        assert len(compound.detected_peaks) == 0
        assert compound.selected_peak is None

    def test_create_compound_with_peaks(self):
        """Compound with detected peaks."""
        blocks = [BuildingBlock.from_code(0, "Pro")]
        chrom = Chromatogram(
            time_points=np.array([0.0, 1.0], dtype=np.float64),
            counts=np.array([10.0, 20.0], dtype=np.float64)
        )
        peak = Peak(
            position=0.5,
            left_base=0.2,
            right_base=0.8,
            height=20.0,
            area=10.0
        )

        compound = Compound(
            building_blocks=blocks,
            chromatogram=chrom,
            detected_peaks=[peak]
        )

        assert len(compound.detected_peaks) == 1
        assert compound.detected_peaks[0] == peak

    def test_create_compound_with_selected_peak(self):
        """Compound with selected peak."""
        blocks = [BuildingBlock.from_code(0, "Pro")]
        chrom = Chromatogram(
            time_points=np.array([0.0, 1.0], dtype=np.float64),
            counts=np.array([10.0, 20.0], dtype=np.float64)
        )
        peak = Peak(
            position=0.5,
            left_base=0.2,
            right_base=0.8,
            height=20.0,
            area=10.0
        )

        compound = Compound(
            building_blocks=blocks,
            chromatogram=chrom,
            detected_peaks=[peak],
            selected_peak=peak
        )

        assert compound.selected_peak == peak

    def test_create_compound_with_id(self):
        """Compound with compound_id."""
        blocks = [BuildingBlock.from_code(0, "Pro")]
        chrom = Chromatogram(
            time_points=np.array([0.0, 1.0], dtype=np.float64),
            counts=np.array([10.0, 20.0], dtype=np.float64)
        )

        compound = Compound(
            building_blocks=blocks,
            chromatogram=chrom,
            compound_id="COMP_001"
        )

        assert compound.compound_id == "COMP_001"


class TestCompoundValidation:
    """Test Compound validation rules."""

    def test_empty_building_blocks_raises_error(self):
        """Compound must have at least one building block."""
        chrom = Chromatogram(
            time_points=np.array([0.0, 1.0], dtype=np.float64),
            counts=np.array([10.0, 20.0], dtype=np.float64)
        )

        with pytest.raises(ValueError, match="at least one building block"):
            Compound(building_blocks=[], chromatogram=chrom)

    def test_building_blocks_not_ordered_raises_error(self):
        """Building blocks must be ordered by cycle."""
        blocks = [
            BuildingBlock.from_code(1, "Leu"),  # Cycle 1
            BuildingBlock.from_code(0, "Pro")   # Cycle 0 - out of order!
        ]
        chrom = Chromatogram(
            time_points=np.array([0.0, 1.0], dtype=np.float64),
            counts=np.array([10.0, 20.0], dtype=np.float64)
        )

        with pytest.raises(ValueError, match="ordered by cycle"):
            Compound(building_blocks=blocks, chromatogram=chrom)


class TestPositionalSequence:
    """Test positional_sequence property (THEORY.md Section 2.2)."""

    def test_positional_sequence_simple(self):
        """Positional sequence for simple compound."""
        blocks = [
            BuildingBlock.from_code(0, "Pro"),
            BuildingBlock.from_code(1, "Leu")
        ]
        chrom = Chromatogram(
            time_points=np.array([0.0, 1.0], dtype=np.float64),
            counts=np.array([10.0, 20.0], dtype=np.float64)
        )

        compound = Compound(building_blocks=blocks, chromatogram=chrom)

        # N→C order: position 0 = C-terminus (rightmost)
        assert compound.positional_sequence == "Leu-Pro"

    def test_positional_sequence_with_null(self):
        """Positional sequence includes null positions."""
        blocks = [
            BuildingBlock.from_code(0, "Pro"),
            BuildingBlock.from_code(1, "Null"),
            BuildingBlock.from_code(2, "Leu")
        ]
        chrom = Chromatogram(
            time_points=np.array([0.0, 1.0], dtype=np.float64),
            counts=np.array([10.0, 20.0], dtype=np.float64)
        )

        compound = Compound(building_blocks=blocks, chromatogram=chrom)

        assert compound.positional_sequence == "Leu-Null-Pro"

    def test_positional_sequence_all_null(self):
        """All-null compound has positional sequence."""
        blocks = [
            BuildingBlock.from_code(0, "Null"),
            BuildingBlock.from_code(1, "Null")
        ]
        chrom = Chromatogram(
            time_points=np.array([0.0, 1.0], dtype=np.float64),
            counts=np.array([10.0, 20.0], dtype=np.float64)
        )

        compound = Compound(building_blocks=blocks, chromatogram=chrom)

        assert compound.positional_sequence == "Null-Null"

    def test_positional_sequence_n_to_c_order(self):
        """Positional sequence follows N→C convention (THEORY.md Section 1.5.1)."""
        blocks = [
            BuildingBlock.from_code(0, "Val"),  # C-terminus (synthesized first)
            BuildingBlock.from_code(1, "Ala"),
            BuildingBlock.from_code(2, "Leu")   # N-terminus (synthesized last)
        ]
        chrom = Chromatogram(
            time_points=np.array([0.0, 1.0], dtype=np.float64),
            counts=np.array([10.0, 20.0], dtype=np.float64)
        )

        compound = Compound(building_blocks=blocks, chromatogram=chrom)

        # N→C: Leu (N-terminus) - Ala - Val (C-terminus)
        assert compound.positional_sequence == "Leu-Ala-Val"


class TestResidueSequence:
    """Test residue_sequence property (THEORY.md Section 2.2)."""

    def test_residue_sequence_no_nulls(self):
        """Residue sequence without nulls."""
        blocks = [
            BuildingBlock.from_code(0, "Pro"),
            BuildingBlock.from_code(1, "Leu")
        ]
        chrom = Chromatogram(
            time_points=np.array([0.0, 1.0], dtype=np.float64),
            counts=np.array([10.0, 20.0], dtype=np.float64)
        )

        compound = Compound(building_blocks=blocks, chromatogram=chrom)

        assert compound.residue_sequence == "Leu-Pro"

    def test_residue_sequence_with_null(self):
        """Residue sequence excludes null blocks."""
        blocks = [
            BuildingBlock.from_code(0, "Pro"),
            BuildingBlock.from_code(1, "Null"),
            BuildingBlock.from_code(2, "Leu")
        ]
        chrom = Chromatogram(
            time_points=np.array([0.0, 1.0], dtype=np.float64),
            counts=np.array([10.0, 20.0], dtype=np.float64)
        )

        compound = Compound(building_blocks=blocks, chromatogram=chrom)

        # Null excluded
        assert compound.residue_sequence == "Leu-Pro"

    def test_residue_sequence_all_null(self):
        """All-null compound has empty residue sequence."""
        blocks = [
            BuildingBlock.from_code(0, "Null"),
            BuildingBlock.from_code(1, "Null")
        ]
        chrom = Chromatogram(
            time_points=np.array([0.0, 1.0], dtype=np.float64),
            counts=np.array([10.0, 20.0], dtype=np.float64)
        )

        compound = Compound(building_blocks=blocks, chromatogram=chrom)

        assert compound.residue_sequence == ""

    def test_positional_variants_same_residue_sequence(self):
        """Different positional sequences can have same residue sequence (THEORY.md Section 4.2.1)."""
        # Variant 1: [Leu, Null, Pro]
        blocks1 = [
            BuildingBlock.from_code(0, "Pro"),
            BuildingBlock.from_code(1, "Null"),
            BuildingBlock.from_code(2, "Leu")
        ]

        # Variant 2: [Null, Leu, Pro]
        blocks2 = [
            BuildingBlock.from_code(0, "Pro"),
            BuildingBlock.from_code(1, "Leu"),
            BuildingBlock.from_code(2, "Null")
        ]

        chrom = Chromatogram(
            time_points=np.array([0.0, 1.0], dtype=np.float64),
            counts=np.array([10.0, 20.0], dtype=np.float64)
        )

        compound1 = Compound(building_blocks=blocks1, chromatogram=chrom)
        compound2 = Compound(building_blocks=blocks2, chromatogram=chrom)

        # Different positional sequences
        assert compound1.positional_sequence != compound2.positional_sequence
        # Same residue sequence
        assert compound1.residue_sequence == compound2.residue_sequence == "Leu-Pro"


class TestMonomerSequence:
    """Test monomer_sequence property (THEORY.md Section 2.2, 1.5.3)."""

    def test_monomer_sequence_simple(self):
        """Monomer sequence for simple monomers."""
        blocks = [
            BuildingBlock.from_code(0, "Pro"),
            BuildingBlock.from_code(1, "Leu")
        ]
        chrom = Chromatogram(
            time_points=np.array([0.0, 1.0], dtype=np.float64),
            counts=np.array([10.0, 20.0], dtype=np.float64)
        )

        compound = Compound(building_blocks=blocks, chromatogram=chrom)

        assert compound.monomer_sequence == "Leu-Pro"

    def test_monomer_sequence_with_composite(self):
        """Monomer sequence decomposes composite blocks."""
        blocks = [
            BuildingBlock.from_code(0, "Pro"),
            BuildingBlock.from_code(1, "Leu-Ala-Val"),  # Composite (trimeric)
            BuildingBlock.from_code(2, "Leu")
        ]
        chrom = Chromatogram(
            time_points=np.array([0.0, 1.0], dtype=np.float64),
            counts=np.array([10.0, 20.0], dtype=np.float64)
        )

        compound = Compound(building_blocks=blocks, chromatogram=chrom)

        # Composite block decomposed
        assert compound.monomer_sequence == "Leu-Leu-Ala-Val-Pro"

    def test_monomer_sequence_with_null(self):
        """Null blocks contribute nothing to monomer sequence."""
        blocks = [
            BuildingBlock.from_code(0, "Pro"),
            BuildingBlock.from_code(1, "Null"),
            BuildingBlock.from_code(2, "Leu")
        ]
        chrom = Chromatogram(
            time_points=np.array([0.0, 1.0], dtype=np.float64),
            counts=np.array([10.0, 20.0], dtype=np.float64)
        )

        compound = Compound(building_blocks=blocks, chromatogram=chrom)

        assert compound.monomer_sequence == "Leu-Pro"

    def test_monomer_sequence_all_null(self):
        """All-null compound has empty monomer sequence."""
        blocks = [
            BuildingBlock.from_code(0, "Null"),
            BuildingBlock.from_code(1, "Null")
        ]
        chrom = Chromatogram(
            time_points=np.array([0.0, 1.0], dtype=np.float64),
            counts=np.array([10.0, 20.0], dtype=np.float64)
        )

        compound = Compound(building_blocks=blocks, chromatogram=chrom)

        assert compound.monomer_sequence == ""

    def test_monomer_sequence_multiple_composites(self):
        """Multiple composite blocks all decompose."""
        blocks = [
            BuildingBlock.from_code(0, "Pro-Gly"),        # Dimeric
            BuildingBlock.from_code(1, "Leu-Ala-Val"),    # Trimeric
        ]
        chrom = Chromatogram(
            time_points=np.array([0.0, 1.0], dtype=np.float64),
            counts=np.array([10.0, 20.0], dtype=np.float64)
        )

        compound = Compound(building_blocks=blocks, chromatogram=chrom)

        assert compound.monomer_sequence == "Leu-Ala-Val-Pro-Gly"


class TestLevelCalculations:
    """Test level and monomer_level properties (THEORY.md Section 3.3)."""

    def test_level_all_non_null(self):
        """Level counts non-null building blocks."""
        blocks = [
            BuildingBlock.from_code(0, "Pro"),
            BuildingBlock.from_code(1, "Leu"),
            BuildingBlock.from_code(2, "Ala")
        ]
        chrom = Chromatogram(
            time_points=np.array([0.0, 1.0], dtype=np.float64),
            counts=np.array([10.0, 20.0], dtype=np.float64)
        )

        compound = Compound(building_blocks=blocks, chromatogram=chrom)

        assert compound.level == 3

    def test_level_with_nulls(self):
        """Level excludes null blocks."""
        blocks = [
            BuildingBlock.from_code(0, "Pro"),
            BuildingBlock.from_code(1, "Null"),
            BuildingBlock.from_code(2, "Leu"),
            BuildingBlock.from_code(3, "Null")
        ]
        chrom = Chromatogram(
            time_points=np.array([0.0, 1.0], dtype=np.float64),
            counts=np.array([10.0, 20.0], dtype=np.float64)
        )

        compound = Compound(building_blocks=blocks, chromatogram=chrom)

        assert compound.level == 2

    def test_level_all_null(self):
        """L₀ (all-null) has level 0."""
        blocks = [
            BuildingBlock.from_code(0, "Null"),
            BuildingBlock.from_code(1, "Null")
        ]
        chrom = Chromatogram(
            time_points=np.array([0.0, 1.0], dtype=np.float64),
            counts=np.array([10.0, 20.0], dtype=np.float64)
        )

        compound = Compound(building_blocks=blocks, chromatogram=chrom)

        assert compound.level == 0

    def test_monomer_level_simple(self):
        """Monomer level counts individual monomers."""
        blocks = [
            BuildingBlock.from_code(0, "Pro"),
            BuildingBlock.from_code(1, "Leu"),
            BuildingBlock.from_code(2, "Ala")
        ]
        chrom = Chromatogram(
            time_points=np.array([0.0, 1.0], dtype=np.float64),
            counts=np.array([10.0, 20.0], dtype=np.float64)
        )

        compound = Compound(building_blocks=blocks, chromatogram=chrom)

        assert compound.monomer_level == 3

    def test_monomer_level_with_composite(self):
        """Monomer level counts monomers from composite blocks."""
        blocks = [
            BuildingBlock.from_code(0, "Pro"),              # 1 monomer
            BuildingBlock.from_code(1, "Leu-Ala-Val"),     # 3 monomers
            BuildingBlock.from_code(2, "Leu")              # 1 monomer
        ]
        chrom = Chromatogram(
            time_points=np.array([0.0, 1.0], dtype=np.float64),
            counts=np.array([10.0, 20.0], dtype=np.float64)
        )

        compound = Compound(building_blocks=blocks, chromatogram=chrom)

        assert compound.level == 3           # 3 building blocks
        assert compound.monomer_level == 5   # 5 total monomers

    def test_monomer_level_with_null(self):
        """Null blocks contribute 0 to monomer level."""
        blocks = [
            BuildingBlock.from_code(0, "Pro"),
            BuildingBlock.from_code(1, "Null"),
            BuildingBlock.from_code(2, "Leu-Ala")
        ]
        chrom = Chromatogram(
            time_points=np.array([0.0, 1.0], dtype=np.float64),
            counts=np.array([10.0, 20.0], dtype=np.float64)
        )

        compound = Compound(building_blocks=blocks, chromatogram=chrom)

        assert compound.level == 2           # 2 non-null building blocks
        assert compound.monomer_level == 3   # Pro(1) + Leu-Ala(2) = 3 monomers


class TestNullCompound:
    """Test null compound handling (L₀)."""

    def test_is_null_compound_all_null(self):
        """is_null_compound returns True for all-null."""
        blocks = [
            BuildingBlock.from_code(0, "Null"),
            BuildingBlock.from_code(1, "Null")
        ]
        chrom = Chromatogram(
            time_points=np.array([0.0, 1.0], dtype=np.float64),
            counts=np.array([10.0, 20.0], dtype=np.float64)
        )

        compound = Compound(building_blocks=blocks, chromatogram=chrom)

        assert compound.is_null_compound is True

    def test_is_null_compound_with_non_null(self):
        """is_null_compound returns False if any non-null."""
        blocks = [
            BuildingBlock.from_code(0, "Pro"),
            BuildingBlock.from_code(1, "Null")
        ]
        chrom = Chromatogram(
            time_points=np.array([0.0, 1.0], dtype=np.float64),
            counts=np.array([10.0, 20.0], dtype=np.float64)
        )

        compound = Compound(building_blocks=blocks, chromatogram=chrom)

        assert compound.is_null_compound is False

    def test_null_compound_properties(self):
        """Null compound has expected properties."""
        blocks = [
            BuildingBlock.from_code(0, "Null"),
            BuildingBlock.from_code(1, "Null"),
            BuildingBlock.from_code(2, "Null")
        ]
        chrom = Chromatogram(
            time_points=np.array([0.0, 1.0], dtype=np.float64),
            counts=np.array([10.0, 20.0], dtype=np.float64)
        )

        compound = Compound(building_blocks=blocks, chromatogram=chrom)

        assert compound.is_null_compound is True
        assert compound.level == 0
        assert compound.monomer_level == 0
        assert compound.residue_sequence == ""
        assert compound.monomer_sequence == ""
        assert compound.positional_sequence == "Null-Null-Null"


class TestStringRepresentation:
    """Test string representation methods."""

    def test_str_returns_positional_sequence(self):
        """str() returns positional sequence."""
        blocks = [
            BuildingBlock.from_code(0, "Pro"),
            BuildingBlock.from_code(1, "Leu")
        ]
        chrom = Chromatogram(
            time_points=np.array([0.0, 1.0], dtype=np.float64),
            counts=np.array([10.0, 20.0], dtype=np.float64)
        )

        compound = Compound(building_blocks=blocks, chromatogram=chrom)

        assert str(compound) == "Leu-Pro"

    def test_repr_shows_key_info(self):
        """repr() shows sequence, level, and n_peaks."""
        blocks = [
            BuildingBlock.from_code(0, "Pro"),
            BuildingBlock.from_code(1, "Leu")
        ]
        chrom = Chromatogram(
            time_points=np.array([0.0, 1.0], dtype=np.float64),
            counts=np.array([10.0, 20.0], dtype=np.float64)
        )
        peak = Peak(
            position=0.5,
            left_base=0.2,
            right_base=0.8,
            height=20.0,
            area=10.0
        )

        compound = Compound(
            building_blocks=blocks,
            chromatogram=chrom,
            detected_peaks=[peak]
        )

        repr_str = repr(compound)
        assert "Leu-Pro" in repr_str
        assert "level=2" in repr_str
        assert "n_peaks=1" in repr_str


class TestIntegrationWithBuildingBlock:
    """Test integration with BuildingBlock entity."""

    def test_composite_building_blocks(self):
        """Compound handles composite building blocks correctly."""
        blocks = [
            BuildingBlock.from_code(0, "Pro"),
            BuildingBlock.from_code(1, "Leu-Ala-Val"),
            BuildingBlock.from_code(2, "DLeuMe-DPro")
        ]
        chrom = Chromatogram(
            time_points=np.array([0.0, 1.0], dtype=np.float64),
            counts=np.array([10.0, 20.0], dtype=np.float64)
        )

        compound = Compound(building_blocks=blocks, chromatogram=chrom)

        # Positional includes composite as-is
        assert compound.positional_sequence == "DLeuMe-DPro-Leu-Ala-Val-Pro"
        # Residue also includes composite as-is
        assert compound.residue_sequence == "DLeuMe-DPro-Leu-Ala-Val-Pro"
        # Monomer decomposes all
        assert compound.monomer_sequence == "DLeuMe-DPro-Leu-Ala-Val-Pro"
        # Level counts building blocks
        assert compound.level == 3
        # Monomer level counts individual monomers
        assert compound.monomer_level == 6

    def test_null_detection_from_building_block(self):
        """Compound respects BuildingBlock null detection."""
        blocks = [
            BuildingBlock.from_code(0, "Pro"),
            BuildingBlock.from_code(1, "NULL"),         # Standard null
            BuildingBlock.from_code(2, "null_variant")  # Variant null
        ]
        chrom = Chromatogram(
            time_points=np.array([0.0, 1.0], dtype=np.float64),
            counts=np.array([10.0, 20.0], dtype=np.float64)
        )

        compound = Compound(building_blocks=blocks, chromatogram=chrom)

        assert compound.level == 1  # Only Pro is non-null
        assert compound.residue_sequence == "Pro"


class TestIntegrationWithChromatogram:
    """Test integration with Chromatogram entity."""

    def test_compound_with_chromatogram_variants(self):
        """Compound works with chromatogram signal variants."""
        blocks = [BuildingBlock.from_code(0, "Pro")]
        chrom = Chromatogram(
            time_points=np.array([0.0, 1.0, 2.0], dtype=np.float64),
            counts=np.array([10.0, 20.0, 15.0], dtype=np.float64)
        )
        chrom.add_signal_variant("corrected", np.array([5.0, 15.0, 10.0], dtype=np.float64))

        compound = Compound(building_blocks=blocks, chromatogram=chrom)

        assert compound.chromatogram.has_signal_variant("corrected")
        assert len(compound.chromatogram) == 3


class TestIntegrationWithPeak:
    """Test integration with Peak entity."""

    def test_compound_with_multiple_peaks(self):
        """Compound can have multiple detected peaks."""
        blocks = [BuildingBlock.from_code(0, "Pro")]
        chrom = Chromatogram(
            time_points=np.array([0.0, 1.0, 2.0], dtype=np.float64),
            counts=np.array([10.0, 20.0, 15.0], dtype=np.float64)
        )
        peak1 = Peak(
            position=0.5,
            left_base=0.2,
            right_base=0.8,
            height=20.0,
            area=10.0,
            peak_type=PeakType.NULL
        )
        peak2 = Peak(
            position=1.5,
            left_base=1.2,
            right_base=1.8,
            height=15.0,
            area=8.0,
            peak_type=PeakType.PUTATIVE_PRODUCT
        )

        compound = Compound(
            building_blocks=blocks,
            chromatogram=chrom,
            detected_peaks=[peak1, peak2]
        )

        assert len(compound.detected_peaks) == 2
        assert peak1 in compound.detected_peaks
        assert peak2 in compound.detected_peaks

    def test_selected_peak_with_validation(self):
        """Selected peak can have validation status."""
        blocks = [BuildingBlock.from_code(0, "Pro")]
        chrom = Chromatogram(
            time_points=np.array([0.0, 1.0], dtype=np.float64),
            counts=np.array([10.0, 20.0], dtype=np.float64)
        )
        peak = Peak(
            position=0.5,
            left_base=0.2,
            right_base=0.8,
            height=20.0,
            area=10.0,
            peak_type=PeakType.PUTATIVE_PRODUCT,
            validation_status=ValidationStatus.VALIDATED
        )

        compound = Compound(
            building_blocks=blocks,
            chromatogram=chrom,
            detected_peaks=[peak],
            selected_peak=peak
        )

        assert compound.selected_peak.is_product_peak is True
        assert compound.selected_peak.is_validated is True


class TestEdgeCases:
    """Test edge cases and boundary conditions."""

    def test_single_building_block(self):
        """Compound with single building block."""
        blocks = [BuildingBlock.from_code(0, "Pro")]
        chrom = Chromatogram(
            time_points=np.array([0.0, 1.0], dtype=np.float64),
            counts=np.array([10.0, 20.0], dtype=np.float64)
        )

        compound = Compound(building_blocks=blocks, chromatogram=chrom)

        assert compound.positional_sequence == "Pro"
        assert compound.level == 1

    def test_many_building_blocks(self):
        """Compound with many building blocks."""
        blocks = [BuildingBlock.from_code(i, f"BB{i}") for i in range(10)]
        chrom = Chromatogram(
            time_points=np.array([0.0, 1.0], dtype=np.float64),
            counts=np.array([10.0, 20.0], dtype=np.float64)
        )

        compound = Compound(building_blocks=blocks, chromatogram=chrom)

        assert compound.level == 10
        assert len(compound.building_blocks) == 10

    def test_compound_with_no_peaks(self):
        """Compound with empty detected_peaks."""
        blocks = [BuildingBlock.from_code(0, "Pro")]
        chrom = Chromatogram(
            time_points=np.array([0.0, 1.0], dtype=np.float64),
            counts=np.array([10.0, 20.0], dtype=np.float64)
        )

        compound = Compound(building_blocks=blocks, chromatogram=chrom)

        assert len(compound.detected_peaks) == 0
        assert compound.selected_peak is None


class TestSequenceConventions:
    """Test N→C sequence conventions (THEORY.md Section 1.5.1)."""

    def test_position_0_is_c_terminus(self):
        """Position 0 = C-terminus (rightmost, synthesized first)."""
        blocks = [
            BuildingBlock.from_code(0, "Val"),  # C-terminus
            BuildingBlock.from_code(1, "Ala"),
            BuildingBlock.from_code(2, "Leu")   # N-terminus
        ]
        chrom = Chromatogram(
            time_points=np.array([0.0, 1.0], dtype=np.float64),
            counts=np.array([10.0, 20.0], dtype=np.float64)
        )

        compound = Compound(building_blocks=blocks, chromatogram=chrom)

        # N→C order: Leu-Ala-Val
        # N-terminus (Leu) is leftmost
        # C-terminus (Val) is rightmost
        assert compound.positional_sequence == "Leu-Ala-Val"
        assert compound.positional_sequence.split("-")[0] == "Leu"   # N-terminus
        assert compound.positional_sequence.split("-")[-1] == "Val"  # C-terminus

    def test_synthesis_order_vs_display_order(self):
        """Synthesis order (0→N) is opposite of display order (N→C)."""
        blocks = [
            BuildingBlock.from_code(0, "A"),  # Synthesized first (C-terminus)
            BuildingBlock.from_code(1, "B"),
            BuildingBlock.from_code(2, "C"),  # Synthesized last (N-terminus)
        ]
        chrom = Chromatogram(
            time_points=np.array([0.0, 1.0], dtype=np.float64),
            counts=np.array([10.0, 20.0], dtype=np.float64)
        )

        compound = Compound(building_blocks=blocks, chromatogram=chrom)

        # Display in N→C order (reverse of synthesis)
        assert compound.positional_sequence == "C-B-A"


class TestCompositeBlocks:
    """Test composite building block handling (THEORY.md Section 1.5.3)."""

    def test_dimeric_composite_block(self):
        """Dimeric composite block."""
        blocks = [
            BuildingBlock.from_code(0, "Pro"),
            BuildingBlock.from_code(1, "Leu-Ala")  # Dimeric
        ]
        chrom = Chromatogram(
            time_points=np.array([0.0, 1.0], dtype=np.float64),
            counts=np.array([10.0, 20.0], dtype=np.float64)
        )

        compound = Compound(building_blocks=blocks, chromatogram=chrom)

        assert compound.positional_sequence == "Leu-Ala-Pro"
        assert compound.monomer_sequence == "Leu-Ala-Pro"
        assert compound.level == 2           # 2 building blocks
        assert compound.monomer_level == 3   # 3 monomers

    def test_trimeric_composite_block(self):
        """Trimeric composite block."""
        blocks = [
            BuildingBlock.from_code(0, "Pro"),
            BuildingBlock.from_code(1, "Leu-Ala-Val")  # Trimeric
        ]
        chrom = Chromatogram(
            time_points=np.array([0.0, 1.0], dtype=np.float64),
            counts=np.array([10.0, 20.0], dtype=np.float64)
        )

        compound = Compound(building_blocks=blocks, chromatogram=chrom)

        assert compound.positional_sequence == "Leu-Ala-Val-Pro"
        assert compound.monomer_sequence == "Leu-Ala-Val-Pro"
        assert compound.level == 2           # 2 building blocks
        assert compound.monomer_level == 4   # 4 monomers
