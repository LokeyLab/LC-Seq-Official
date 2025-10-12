"""Lineage finder domain service for hierarchy graph traversal.

This module provides domain services for finding compound lineages
(reference + descendants) within a hierarchy graph.

Terminology (THEORY.md Section 3.1):
    - Reference compound: The compound currently being analyzed (not "parent")
    - Lineage: All ancestors + descendants + self (not "family")
    - Descendant: Compound with fewer building blocks
    - Ancestor: Compound with more building blocks

Domain Responsibility:
    - Graph traversal to find descendants
    - Lineage extraction from hierarchy
    - Pure algorithmic logic using hierarchy operations

Not Responsible For:
    - Building the hierarchy (that's HierarchyBuilder)
    - Sorting compounds for display (that's presentation)
    - Visualization or plotting

References:
    THEORY.md Section 3.1: Terminology - Ancestry and Lineage
"""

from typing import List, Optional
from ..entities import Compound, BuildingBlock
from ..models import CompoundHierarchy, HierarchyMode


class LineageFinderService:
    """Domain service for finding compound lineages in hierarchy.

    This service provides algorithms for extracting lineages (reference +
    descendants) from a compound hierarchy. All methods are deterministic
    and side-effect free.

    Terminology aligned with THEORY.md Section 3.1:
    - Uses "reference" instead of "parent"
    - Uses "lineage" instead of "family"

    Examples:
        >>> finder = LineageFinderService()
        >>> lineage = finder.find_lineage(reference, hierarchy)
        >>> # Returns reference + all descendants
    """

    def find_principal_ideal(
        self,
        reference: Compound,
        compounds: List[Compound],
        mode: HierarchyMode
    ) -> List[Compound]:
        """Find principal ideal ↓X: reference + all descendants from compound library.

        This finds the lineage directly from the compound library WITHOUT building
        a full hierarchy first. This is the performance-optimized method for lineage
        analysis.

        Principal Ideal ↓X (THEORY.md Section 3.2):
        - Reference compound X
        - All compounds that are descendants of X (truncations)

        Building Block Mode vs Monomer Mode:
        - Building Block: Uses block support sequence for descendant checking
          (all positional variants included)
        - Monomer: Uses full monomer sequence for descendant checking
          (all positional variants included)

        Performance:
        - O(N * M) where N = library size, M = sequence length
        - Much faster than building full 64k-compound hierarchy first
        - Returns only lineage members for minimal hierarchy construction

        Parameters
        ----------
        reference : Compound
            Reference compound (root of lineage)
        compounds : List[Compound]
            Full compound library to search
        mode : HierarchyMode
            Building-block or monomer mode (determines descendant logic)

        Returns
        -------
        List[Compound]
            Lineage members: [reference] + descendants (principal ideal ↓X)
            All positional variants are included in both modes

        Examples
        --------
        >>> finder = LineageFinderService()
        >>> lineage = finder.find_principal_ideal(reference, compounds, HierarchyMode.MONOMER)
        >>> # Returns ~50 compounds from 64k library
        >>> len(lineage)
        47

        References
        ----------
        THEORY.md Section 3.2: Principal Ideal ↓X - All descendants of compound X
        """
        lineage = [reference]

        if mode == HierarchyMode.BUILDING_BLOCK:
            # Building-block mode: block support sequence descendant checking
            for compound in compounds:
                if compound == reference:
                    continue
                if self._is_building_block_descendant(compound, reference):
                    lineage.append(compound)
        else:  # MONOMER mode
            # Monomer mode: chemical identity (residue sequence) comparison
            for compound in compounds:
                if compound == reference:
                    continue
                if self._is_monomer_descendant(compound, reference):
                    lineage.append(compound)

        return lineage

    def _is_building_block_descendant(self, candidate: Compound, reference: Compound) -> bool:
        """Check if candidate is a building-block descendant of reference.

        Building-block descendant rules (THEORY.md Section 3.3):
        - Lower or equal level (candidate.level <= reference.level)
        - Candidate's block_support_sequence is a subsequence of reference's
        - Positions are IGNORED - only non-null blocks matter

        Examples:
            Reference: "Leu-Ala-Pro" (block_support_sequence: "Leu-Ala-Pro")
            Descendant: "Leu-Pro" (block_support_sequence: "Leu-Pro") ✓
            NOT descendant: "Leu-Null-Pro" (block_support_sequence: "Leu-Pro" at DIFFERENT positions) ✓
            NOT descendant: "Pro-Leu" (wrong order) ✗

        Returns True if candidate is descendant, False otherwise.
        """
        # Must have lower or equal level
        if candidate.level > reference.level:
            return False

        # Get block support sequences (non-null blocks only)
        candidate_blocks = [bb.code for bb in reversed(candidate.building_blocks) if not bb.is_null]
        reference_blocks = [bb.code for bb in reversed(reference.building_blocks) if not bb.is_null]

        # All-null compound is descendant of everything
        if not candidate_blocks:
            return True

        # Candidate must be subsequence of reference (same order, possible gaps)
        # Use greedy left-to-right matching
        ref_idx = 0
        for cand_block in candidate_blocks:
            # Find next occurrence of cand_block in remaining reference
            found = False
            while ref_idx < len(reference_blocks):
                if reference_blocks[ref_idx] == cand_block:
                    ref_idx += 1
                    found = True
                    break
                ref_idx += 1

            if not found:
                # Block not found in remaining reference sequence
                return False

        # All candidate blocks matched in order
        return True

    def _is_monomer_descendant(self, candidate: Compound, reference: Compound, debug: bool = False) -> bool:
        """Check if candidate is a monomer-level descendant of reference.

        Monomer descendant rules (THEORY.md Section 1.5.6, lines 461-473):
        - Descendants formed by REMOVING monomers from positions (not arbitrary subsets)
        - Result is a SUBSEQUENCE: same order, possible gaps
        - Example: "Leu-Ala-Val-Pro" → "Leu-Val-Pro" (removed Ala) ✓
        - Example: "Leu-Ala-Val-Pro" → "Val-Leu-Pro" (wrong order) ✗

        Returns True if candidate is descendant, False otherwise.
        """
        # Use alignment mapping to check descendant relationship
        mapping = self.get_monomer_alignment_mapping(candidate, reference)
        return mapping is not None

    def get_monomer_alignment_mapping(
        self,
        candidate: Compound,
        reference: Compound
    ) -> Optional[List[int]]:
        """Get monomer alignment mapping using greedy subsequence matching.

        Uses the same greedy left-to-right algorithm as descendant checking,
        but returns position mappings for visualization alignment.

        Parameters
        ----------
        candidate : Compound
            Candidate compound (potential descendant)
        reference : Compound
            Reference compound (alignment template)

        Returns
        -------
        Optional[List[int]]
            List of reference positions that matched each candidate monomer.
            None if candidate is not a valid subsequence of reference.
            Example: [0, 2, 4] means candidate monomers matched reference positions 0, 2, 4

        Examples
        --------
        >>> # Reference: Leu-LA03-Pro-Leu-DPro (5 monomers)
        >>> # Candidate: Leu-Pro-DPro (3 monomers)
        >>> mapping = finder.get_monomer_alignment_mapping(candidate, reference)
        >>> mapping
        [0, 2, 4]  # Candidate monomers at ref positions 0, 2, 4 (gaps at 1, 3)

        References
        ----------
        THEORY.md Section 1.5.6: Subsequence matching
        THEORY.md Section 2.3.4: Sequence alignment for visualization
        """
        # Must have lower or equal monomer level
        if candidate.monomer_level > reference.monomer_level:
            return None

        # Decompose building blocks to individual monomer sequences (THEORY.md 2.2)
        # CRITICAL: Reverse to maintain N→C order (building_blocks are stored C→N)
        candidate_monomers = []
        reference_monomers = []

        for bb in reversed(candidate.building_blocks):
            # Decompose composite blocks: "Leu-Ala-Val" → ["Leu", "Ala", "Val"]
            candidate_monomers.extend(bb.decompose_to_monomers())

        for bb in reversed(reference.building_blocks):
            reference_monomers.extend(bb.decompose_to_monomers())

        # All-null is descendant of everything (maps to no positions)
        if not candidate_monomers:
            return []

        # Greedy RIGHT-TO-LEFT subsequence matching (THEORY.md 1.5.6)
        # Match from C-terminus (right) to N-terminus (left) since C-terminus is anchor
        # Track which reference position each candidate monomer matched
        mapping = []
        ref_idx = len(reference_monomers) - 1

        # Process candidate monomers in reverse (C→N direction)
        for cand_monomer in reversed(candidate_monomers):
            # Find previous occurrence of cand_monomer in remaining reference sequence
            found = False
            while ref_idx >= 0:
                if reference_monomers[ref_idx] == cand_monomer:
                    mapping.insert(0, ref_idx)  # Insert at start to maintain order
                    ref_idx -= 1  # Move to previous monomer
                    found = True
                    break
                ref_idx -= 1

            if not found:
                # Monomer not found in remaining reference sequence
                return None

        # All candidate monomers matched in order
        return mapping

    def find_lineage(
        self,
        reference: Compound,
        hierarchy: CompoundHierarchy
    ) -> List[Compound]:
        """Find lineage: reference compound + all descendants.

        A lineage consists of:
        - The reference compound (root of subgraph)
        - All descendants (compounds reachable from reference)

        This uses the hierarchy's graph traversal to find all reachable
        descendants.

        Parameters
        ----------
        reference : Compound
            Reference compound (root of lineage) - THEORY.md Section 3.1:
            "The compound currently being analyzed"
        hierarchy : CompoundHierarchy
            Hierarchy containing the reference

        Returns
        -------
        List[Compound]
            Lineage members: [reference] + descendants

        Examples
        --------
        >>> finder = LineageFinderService()
        >>> lineage = finder.find_lineage(reference, hierarchy)
        >>> len(lineage)  # reference + descendants
        8

        References
        ----------
        THEORY.md Section 3.1: "Lineage: All ancestors + descendants + self"
        """
        lineage = [reference]
        descendants = hierarchy.get_descendants(reference)
        lineage.extend(descendants)
        return lineage

    def count_lineage_by_level(
        self,
        lineage: List[Compound],
        use_monomer_level: bool = False
    ) -> dict:
        """Count lineage members at each truncation level.

        Parameters
        ----------
        lineage : List[Compound]
            Lineage members to analyze
        use_monomer_level : bool, optional
            Use monomer level instead of building block level. Default is False.

        Returns
        -------
        dict
            Mapping of level → count

        Examples
        --------
        >>> finder = LineageFinderService()
        >>> counts = finder.count_lineage_by_level(lineage)
        >>> counts
        {3: 1, 2: 3, 1: 3, 0: 1}
        """
        counts = {}

        for compound in lineage:
            level = compound.monomer_level if use_monomer_level else compound.level
            counts[level] = counts.get(level, 0) + 1

        return counts

    def group_lineage_by_level(
        self,
        lineage: List[Compound],
        use_monomer_level: bool = False
    ) -> dict:
        """Group lineage members by truncation level.

        Parameters
        ----------
        lineage : List[Compound]
            Lineage members to group
        use_monomer_level : bool, optional
            Use monomer level instead of building block level. Default is False.

        Returns
        -------
        dict
            Mapping of level → List[Compound]

        Examples
        --------
        >>> finder = LineageFinderService()
        >>> by_level = finder.group_lineage_by_level(lineage)
        >>> by_level[3]  # All level-3 compounds
        [<Compound: Leu-Pro-Ala>]
        """
        groups = {}

        for compound in lineage:
            level = compound.monomer_level if use_monomer_level else compound.level
            if level not in groups:
                groups[level] = []
            groups[level].append(compound)

        return groups
