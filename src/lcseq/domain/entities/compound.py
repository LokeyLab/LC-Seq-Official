"""
Compound entity - represents a library member with sequence and chromatogram.

Implementation based on THEORY.md Section 2.1, 2.2, 3.3.
"""

from dataclasses import dataclass, field
from typing import List, Optional
from .building_block import BuildingBlock
from .chromatogram import Chromatogram
from .peak import Peak


@dataclass
class Compound:
    """
    Represents a DNA-encoded library member.

    A compound is defined by its building block sequence and associated
    chromatogram. It may have detected peaks, hierarchical relationships
    (ancestors/descendants), and multiple sequence representations.

    Attributes
    ----------
    building_blocks : List[BuildingBlock]
        Building blocks in N→C order (position 0 = C-terminus)
    chromatogram : Chromatogram
        Elution profile for this compound
    detected_peaks : List[Peak]
        All detected peaks in chromatogram
    selected_peak : Optional[Peak]
        The peak selected as putative product (if any)
    compound_id : Optional[str]
        Unique identifier for this compound

    Notes
    -----
    - Building blocks stored in N→C order (THEORY.md Section 1.5.1)
    - Position 0 = C-terminus (rightmost, synthesized first)
    - Position N = N-terminus (leftmost, synthesized last)
    - Sequence representations computed from building blocks
    - Hashable based on building block sequence for use in sets/dicts
    - Equality based on building block sequence (not chromatogram)

    References
    ----------
    THEORY.md Section 2.1: Core Entities
    THEORY.md Section 2.2: Sequence Representations
    THEORY.md Section 3.3: Hierarchy Properties
    THEORY.md Section 1.5.1: Peptide Sequence Convention
    """

    building_blocks: List[BuildingBlock]
    chromatogram: Chromatogram
    detected_peaks: List[Peak] = field(default_factory=list)
    selected_peak: Optional[Peak] = None
    compound_id: Optional[str] = None

    # cLPE validation reference data (loaded from external source)
    # LogK is computed from OBSERVED RT, not loaded from reference
    alogp: Optional[float] = None  # Calculated lipophilicity from structure
    scaffold_group: Optional[str] = None  # Stereochemistry grouping for cLPE regression

    def __post_init__(self) -> None:
        """Validate compound properties."""
        if not self.building_blocks:
            raise ValueError("Compound must have at least one building block")
        
        # Validate building blocks are in order by cycle
        cycles = [bb.cycle for bb in self.building_blocks]
        if cycles != sorted(cycles):
            raise ValueError(
                f"Building blocks must be ordered by cycle. Got cycles: {cycles}"
            )

    @property
    def positional_block_sequence(self) -> str:
        """
        Positional block sequence including null positions.

        Returns building block codes joined by hyphens, preserving all
        positions including nulls. Shows synthesis path at block granularity.

        Returns
        -------
        str
            Full positional block sequence (e.g., "Leu-Null-Pro", "Leu-Ala-Val-Pro")

        Examples
        --------
        >>> blocks = [
        ...     BuildingBlock.from_code(0, "Pro"),
        ...     BuildingBlock.from_code(1, "Null"),
        ...     BuildingBlock.from_code(2, "Leu")
        ... ]
        >>> compound = Compound(blocks, chromatogram)
        >>> compound.positional_block_sequence
        'Leu-Null-Pro'

        Notes
        -----
        Sequence is in N→C order (position 0 = C-terminus on right).
        This is the unique representation per compound (synthesis path).

        References
        ----------
        THEORY.md Section 2.2: Positional Block Sequence
        """
        # Blocks are stored in order by cycle, so reverse for N→C display
        return "-".join(bb.code for bb in reversed(self.building_blocks))

    @property
    def block_support_sequence(self) -> str:
        """
        Block support sequence (non-null building blocks only).

        Restriction to support (non-null positions). Multiple positional
        block sequences can have the same block support sequence (positional
        variants). This is the equivalence class identifier.

        Returns
        -------
        str
            Non-null building blocks joined by hyphens

        Examples
        --------
        >>> blocks = [
        ...     BuildingBlock.from_code(0, "Pro"),
        ...     BuildingBlock.from_code(1, "Null"),
        ...     BuildingBlock.from_code(2, "Leu")
        ... ]
        >>> compound = Compound(blocks, chromatogram)
        >>> compound.block_support_sequence
        'Leu-Pro'

        Notes
        -----
        Used for grouping positional variants into equivalence classes.
        Same block support sequence = same chemical molecule at block granularity.

        References
        ----------
        THEORY.md Section 2.2: Block Support Sequence
        THEORY.md Section 4.2.1: EquivalenceClass Definition
        """
        non_null_blocks = [bb for bb in reversed(self.building_blocks) if not bb.is_null]
        if not non_null_blocks:
            return ""  # All-null compound
        return "-".join(bb.code for bb in non_null_blocks)

    @property
    def monomer_support_sequence(self) -> str:
        """
        Monomer support sequence (fully decomposed, non-null).

        Composite building blocks expanded to individual monomers, with nulls
        removed (support only). This represents the actual chemical peptide
        sequence at monomer granularity.

        Returns
        -------
        str
            All monomers joined by hyphens

        Examples
        --------
        >>> blocks = [
        ...     BuildingBlock.from_code(0, "Pro"),
        ...     BuildingBlock.from_code(1, "Leu-Ala-Val"),  # Composite
        ...     BuildingBlock.from_code(2, "Leu")
        ... ]
        >>> compound = Compound(blocks, chromatogram)
        >>> compound.monomer_support_sequence
        'Leu-Leu-Ala-Val-Pro'

        Notes
        -----
        - Composite blocks decompose: "Leu-Ala-Val" → ["Leu", "Ala", "Val"]
        - Null blocks contribute nothing (support only)
        - Result is chemical identity at monomer granularity

        References
        ----------
        THEORY.md Section 2.2: Monomer Support Sequence
        THEORY.md Section 1.5.3: Monomer-Level Decomposition
        """
        all_monomers = []
        # Process in reverse order for N→C
        for bb in reversed(self.building_blocks):
            monomers = bb.decompose_to_monomers()
            all_monomers.extend(monomers)
        
        if not all_monomers:
            return ""  # All-null compound
        return "-".join(all_monomers)

    @property
    def level(self) -> int:
        """
        Truncation level (number of non-null building blocks).

        Returns
        -------
        int
            Count of non-null building blocks

        Examples
        --------
        >>> blocks = [
        ...     BuildingBlock.from_code(0, "Pro"),
        ...     BuildingBlock.from_code(1, "Null"),
        ...     BuildingBlock.from_code(2, "Leu")
        ... ]
        >>> compound = Compound(blocks, chromatogram)
        >>> compound.level
        2

        Notes
        -----
        - Level 0 = all null (L₀, complete truncation)
        - Level N = maximal in dataset (N non-null building blocks)

        References
        ----------
        THEORY.md Section 3.3: Truncation Level
        """
        return sum(1 for bb in self.building_blocks if not bb.is_null)

    @property
    def monomer_level(self) -> int:
        """
        Monomer-level (total number of monomers after decomposition).

        Returns
        -------
        int
            Total count of monomers from all building blocks

        Examples
        --------
        >>> blocks = [
        ...     BuildingBlock.from_code(0, "Pro"),          # 1 monomer
        ...     BuildingBlock.from_code(1, "Leu-Ala-Val"), # 3 monomers
        ...     BuildingBlock.from_code(2, "Leu")          # 1 monomer
        ... ]
        >>> compound = Compound(blocks, chromatogram)
        >>> compound.monomer_level
        5

        Notes
        -----
        Used in monomer-based hierarchy mode.
        Composite blocks contribute multiple monomers.

        References
        ----------
        THEORY.md Section 3.3: Monomer Level
        """
        return sum(len(bb.decompose_to_monomers()) for bb in self.building_blocks)

    @property
    def is_null_compound(self) -> bool:
        """
        Check if this is the L₀ (all-null) compound.

        Returns
        -------
        bool
            True if all building blocks are null
        """
        return all(bb.is_null for bb in self.building_blocks)

    def __hash__(self) -> int:
        """
        Hash based on building block sequence.

        Enables use of Compound as dict key or in sets.
        Two compounds with same building block sequence hash to same value.
        """
        return hash(tuple(self.building_blocks))

    def __eq__(self, other: object) -> bool:
        """
        Equality based on building block sequence.

        Two compounds are equal if they have the same building blocks
        (same sequence), regardless of chromatogram or peaks.
        """
        if not isinstance(other, Compound):
            return NotImplemented
        return self.building_blocks == other.building_blocks

    def __str__(self) -> str:
        """String representation using positional block sequence."""
        return self.positional_block_sequence

    def __repr__(self) -> str:
        """Detailed representation for debugging."""
        return (
            f"Compound(sequence='{self.positional_block_sequence}', "
            f"level={self.level}, "
            f"n_peaks={len(self.detected_peaks)})"
        )
