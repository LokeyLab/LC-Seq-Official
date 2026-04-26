"""Sequence alignment domain service for MSA-style compound alignment.

This module provides the single source of truth for MSA-style sequence alignment,
eliminating duplicated alignment logic across visualization components.

Alignment Strategy (THEORY.md Section 2.3.4):
    - Greedy left-to-right subsequence matching
    - Gaps inserted where monomers/blocks are missing
    - Building-block mode: Position-agnostic block support alignment
    - Monomer mode: Residue sequence alignment using lineage mapping

Domain Responsibility:
    - MSA-style alignment algorithm (single source of truth)
    - Gap insertion for missing sequence elements
    - Delegation to LineageFinderService for monomer mappings

Not Responsible For:
    - Visualization or label formatting (that's presentation layer)
    - Sorting compounds or color assignment
    - Computing truncation levels

References:
    THEORY.md Section 2.3.4: Sequence alignment for visualization
    THEORY.md Section 1.5.6: Subsequence matching
"""

from typing import Optional
from ..entities import Compound
from ..models import HierarchyMode
from .lineage_finder import LineageFinderService


class SequenceAligner:
    """
    Single source of truth for MSA-style sequence alignment.

    Aligns compound sequences to a reference using greedy left-to-right
    subsequence matching, inserting gaps where monomers/blocks are missing.

    This service eliminates duplicated alignment logic in visualization components.

    Examples
    --------
    >>> aligner = SequenceAligner()
    >>> aligned = aligner.align_to_reference(compound, reference, HierarchyMode.BUILDING_BLOCK)
    >>> aligned
    'Leu----------Leu-DLeuMe-DPro'  # Gap for missing LA03-Pro block

    >>> aligned = aligner.align_to_reference(compound, reference, HierarchyMode.MONOMER)
    >>> aligned
    'Leu-----Pro---DPro'  # Gaps at positions 1 (LA03) and 3 (Leu)
    """

    def __init__(self, lineage_finder: Optional[LineageFinderService] = None):
        """
        Initialize sequence aligner.

        Parameters
        ----------
        lineage_finder : LineageFinderService, optional
            Service for monomer alignment mappings. If None, creates new instance.
            Reusing a single instance improves performance for batch operations.
        """
        self._lineage_finder = lineage_finder or LineageFinderService()

    def align_to_reference(
        self,
        compound: Compound,
        reference: Compound,
        mode: HierarchyMode = HierarchyMode.BUILDING_BLOCK,
    ) -> str:
        """
        Align compound sequence to reference using MSA-style gaps.

        Building-block mode: Position-agnostic block support alignment
        Monomer mode: Residue sequence alignment using lineage mapping

        Parameters
        ----------
        compound : Compound
            Compound to align
        reference : Compound
            Reference compound (alignment template)
        mode : HierarchyMode, optional
            Alignment mode (default: BUILDING_BLOCK)

        Returns
        -------
        str
            Aligned sequence with "----" gaps for nulls/missing monomers

        Examples
        --------
        >>> # Building-block mode
        >>> aligner.align_to_reference(compound, reference, HierarchyMode.BUILDING_BLOCK)
        'Leu----------Leu-DLeuMe-DPro'  # Gap matches missing block length

        >>> # Monomer mode
        >>> aligner.align_to_reference(compound, reference, HierarchyMode.MONOMER)
        'Leu-----Pro---DPro'  # Gaps match individual monomer lengths

        References
        ----------
        THEORY.md Section 2.3.4: Sequence alignment for visualization
        """
        if mode == HierarchyMode.BUILDING_BLOCK:
            return self._align_building_blocks(compound, reference)
        else:  # MONOMER
            return self._align_monomers(compound, reference)

    def _align_building_blocks(self, compound: Compound, reference: Compound) -> str:
        """
        Align building blocks by block support sequence (subsequence matching).

        In building block mode, alignment is based on non-null blocks only.
        Uses greedy left-to-right subsequence matching to find which reference
        blocks are present in the compound. Missing blocks become gaps.

        Parameters
        ----------
        compound : Compound
            Compound to align
        reference : Compound
            Reference compound (defines gap lengths)

        Returns
        -------
        str
            Block-support-aligned sequence with gaps for missing blocks

        Examples
        --------
        >>> # Reference support: Leu-LA03-Pro-Leu-DLeuMe-DPro (3 blocks)
        >>> # Compound support:  Leu----------Leu-DLeuMe-DPro (2 blocks, missing LA03-Pro)
        >>> aligned = _align_building_blocks(compound, reference)
        >>> aligned
        'Leu----------Leu-DLeuMe-DPro'  # Gap matches "LA03-Pro" length (10 chars including hyphen)

        Notes
        -----
        This aligns by block support sequence (non-null blocks), NOT by synthesis
        position. Positions are ignored in building block mode.

        References
        ----------
        THEORY.md Section 2.3.4: Sequence alignment for visualization
        THEORY.md Section 1.5.6: Greedy subsequence matching
        """
        # Get block support sequences (non-null blocks only)
        compound_support_blocks = [bb.code for bb in reversed(compound.building_blocks) if not bb.is_null]
        reference_support_blocks = [bb.code for bb in reversed(reference.building_blocks) if not bb.is_null]

        # Handle all-null case
        if not compound_support_blocks:
            # All gaps - total length should match reference support
            total_length = sum(len(bb) for bb in reference_support_blocks) + len(reference_support_blocks) - 1
            return "-" * total_length

        # Use greedy left-to-right subsequence matching (same logic as lineage finder)
        # Find which reference blocks are present in compound
        ref_idx = 0
        compound_idx = 0
        aligned = []

        while ref_idx < len(reference_support_blocks):
            ref_block = reference_support_blocks[ref_idx]

            # Check if current compound block matches this reference block
            if compound_idx < len(compound_support_blocks) and compound_support_blocks[compound_idx] == ref_block:
                # Match found - add the block
                aligned.append(ref_block)
                compound_idx += 1
            else:
                # No match - add gap with length matching reference block
                aligned.append("-" * len(ref_block))

            ref_idx += 1

        return "-".join(aligned)

    def _align_monomers(self, compound: Compound, reference: Compound) -> str:
        """
        Align monomers using subsequence mapping to reference.

        Uses LineageFinderService to get position mappings via greedy subsequence
        matching. Inserts gaps (dashes) where monomers are missing.

        Parameters
        ----------
        compound : Compound
            Compound to align (subsequence of reference)
        reference : Compound
            Reference compound (full sequence template)

        Returns
        -------
        str
            Subsequence-aligned sequence with gaps (no spaces, joined by "-")

        Examples
        --------
        >>> # Reference: Leu-LA03-Pro-Leu-DPro (5 monomers, positions 0-4)
        >>> # Compound:  Leu-Pro-DPro (3 monomers)
        >>> # Mapping:   [0, 2, 4] (matches at ref positions 0, 2, 4)
        >>> aligned = _align_monomers(compound, reference)
        >>> aligned
        'Leu-----Pro---DPro'  # Gaps at pos 1 (LA03=5 chars) and 3 (Leu=3 chars)

        References
        ----------
        THEORY.md Section 2.3.4: Sequence alignment for visualization
        THEORY.md Section 1.5.6: Subsequence matching
        """
        # Use domain service to get subsequence alignment mapping
        mapping = self._lineage_finder.get_monomer_alignment_mapping(compound, reference)

        # Get reference monomers for gap sizing
        reference_monomers = []
        for bb in reversed(reference.building_blocks):
            reference_monomers.extend(bb.decompose_to_monomers())

        # Get candidate monomers
        candidate_monomers = []
        for bb in reversed(compound.building_blocks):
            candidate_monomers.extend(bb.decompose_to_monomers())

        # Handle all-null case
        if not mapping:
            # All gaps (all positions in reference)
            return "-".join("-" * len(m) for m in reference_monomers)

        # Build aligned sequence using mapping
        aligned = []
        cand_idx = 0

        for ref_idx, ref_monomer in enumerate(reference_monomers):
            if cand_idx < len(mapping) and mapping[cand_idx] == ref_idx:
                # Candidate has monomer at this reference position
                aligned.append(candidate_monomers[cand_idx])
                cand_idx += 1
            else:
                # Gap: candidate missing monomer at this position
                # Use dashes matching reference monomer length
                aligned.append("-" * len(ref_monomer))

        return "-".join(aligned)
