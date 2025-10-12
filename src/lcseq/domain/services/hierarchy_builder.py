"""
HierarchyBuilder - Constructs CompoundHierarchy from list of compounds.

Implementation based on THEORY.md Section 4.2, 3.3, 1.5.
"""

from typing import List
from ..entities.compound import Compound
from ..entities.building_block import BuildingBlock
from ..models.compound_hierarchy import CompoundHierarchy, HierarchyMode


class HierarchyBuilder:
    """
    Builds CompoundHierarchy from list of compounds.

    Automatically detects truncation relationships between compounds
    and constructs the DAG structure. Supports both building-block
    and monomer modes.

    Notes
    -----
    - Stateless service (no instance state)
    - Pure domain logic (no I/O operations)
    - Detects edges by comparing sequences
    - Validates DAG properties

    Examples
    --------
    >>> from lcseq.domain.entities.building_block import BuildingBlock
    >>> from lcseq.domain.entities.chromatogram import Chromatogram
    >>> import numpy as np
    >>>
    >>> # Create compounds
    >>> chromatogram = Chromatogram(
    ...     time_points=np.array([1.0, 2.0, 3.0]),
    ...     counts=np.array([100.0, 200.0, 150.0])
    ... )
    >>> bb0 = BuildingBlock.from_code(0, "Pro")
    >>> bb1 = BuildingBlock.from_code(1, "Leu")
    >>> bb_null = BuildingBlock.from_code(0, "Null")
    >>>
    >>> maximal = Compound([bb_null, bb1, bb0], chromatogram)
    >>> truncation = Compound([bb_null, bb_null, bb0], chromatogram)
    >>>
    >>> # Build hierarchy
    >>> builder = HierarchyBuilder()
    >>> hierarchy = builder.build([maximal, truncation], HierarchyMode.BUILDING_BLOCK)
    >>> hierarchy.get_descendants(maximal)
    [truncation]

    References
    ----------
    THEORY.md Section 4.2: Hierarchy Construction
    THEORY.md Section 3.3: Hierarchy Properties
    THEORY.md Section 1.5.2: Building-Block Level Truncation
    THEORY.md Section 1.5.3: Monomer-Level Truncation
    """

    def build(
        self,
        compounds: List[Compound],
        mode: HierarchyMode = HierarchyMode.BUILDING_BLOCK
    ) -> CompoundHierarchy:
        """
        Build hierarchy from list of compounds.

        Parameters
        ----------
        compounds : List[Compound]
            All compounds to include in hierarchy
        mode : HierarchyMode, optional
            Building-block or monomer mode (default: BUILDING_BLOCK)

        Returns
        -------
        CompoundHierarchy
            Constructed DAG with all truncation edges

        Notes
        -----
        - Detects all direct truncation edges (one building block/monomer removed)
        - Does not add transitive edges (follows transitive reduction principle)
        - Validates acyclic property and level ordering

        Algorithm
        ---------
        1. Create empty hierarchy with specified mode
        2. Add all compounds to hierarchy
        3. For each pair of compounds (ancestor, descendant):
           - Check if descendant is direct truncation of ancestor
           - If yes, add edge ancestor → descendant
        4. Return completed hierarchy

        Complexity: O(n² × m) where n = compounds, m = sequence length
        """
        # Create hierarchy
        hierarchy = CompoundHierarchy(mode=mode)

        # Add all compounds
        for compound in compounds:
            hierarchy.add_compound(compound)

        # Detect and add edges
        # Use transitive reduction: only add edges to NEAREST descendants
        for ancestor in compounds:
            # Find all potential descendants (any level below ancestor)
            potential_descendants = []
            for descendant in compounds:
                if ancestor == descendant:
                    continue

                # Check level ordering (descendant must be lower level)
                if mode == HierarchyMode.BUILDING_BLOCK:
                    if descendant.level >= ancestor.level:
                        continue
                else:  # MONOMER
                    if descendant.monomer_level >= ancestor.monomer_level:
                        continue

                # Check if this is a valid truncation (possibly multiple levels)
                if self._is_valid_descendant(ancestor, descendant, mode):
                    potential_descendants.append(descendant)

            # Add edges only to NEAREST descendants (direct or closest available)
            if potential_descendants:
                # Group by level
                if mode == HierarchyMode.BUILDING_BLOCK:
                    by_level = {}
                    for desc in potential_descendants:
                        level = desc.level
                        if level not in by_level:
                            by_level[level] = []
                        by_level[level].append(desc)
                else:  # MONOMER
                    by_level = {}
                    for desc in potential_descendants:
                        level = desc.monomer_level
                        if level not in by_level:
                            by_level[level] = []
                        by_level[level].append(desc)

                # Find nearest level with descendants
                nearest_level = max(by_level.keys())

                # Add edges to all descendants at nearest level
                for desc in by_level[nearest_level]:
                    hierarchy.add_edge(ancestor, desc)

        return hierarchy

    def _is_valid_descendant(
        self,
        ancestor: Compound,
        descendant: Compound,
        mode: HierarchyMode
    ) -> bool:
        """
        Check if descendant is a valid truncation of ancestor (any level).

        Unlike _detect_truncation_edge which only checks DIRECT truncations,
        this allows multiple levels of truncation (handles gaps in dataset).

        Parameters
        ----------
        ancestor : Compound
            Potential ancestor compound
        descendant : Compound
            Potential descendant compound
        mode : HierarchyMode
            Building-block or monomer mode

        Returns
        -------
        bool
            True if descendant is a truncation of ancestor (direct or transitive)

        Notes
        -----
        - Building-block mode: Descendant has nulls at positions where ancestor has blocks
        - Monomer mode: Descendant's monomer sequence is subsequence of ancestor's
        - Allows gaps in levels (handles missing intermediate compounds)
        """
        if mode == HierarchyMode.BUILDING_BLOCK:
            return self._is_building_block_descendant(ancestor, descendant)
        else:  # MONOMER
            return self._is_monomer_descendant(ancestor, descendant)

    def _detect_truncation_edge(
        self,
        ancestor: Compound,
        descendant: Compound,
        mode: HierarchyMode
    ) -> bool:
        """
        Detect if descendant is direct truncation of ancestor.

        Parameters
        ----------
        ancestor : Compound
            Potential ancestor compound
        descendant : Compound
            Potential descendant compound
        mode : HierarchyMode
            Building-block or monomer mode

        Returns
        -------
        bool
            True if descendant is direct truncation of ancestor

        Notes
        -----
        - Building-block mode: Exactly one building block replaced with null
        - Monomer mode: Exactly one monomer removed from sequence
        - Direct truncation only (not transitive)

        References
        ----------
        THEORY.md Section 1.5.2: Building-Block Level Truncation
        THEORY.md Section 1.5.3: Monomer-Level Truncation
        """
        if mode == HierarchyMode.BUILDING_BLOCK:
            return self._is_building_block_truncation(ancestor, descendant)
        else:  # MONOMER
            return self._is_monomer_truncation(ancestor, descendant)

    def _is_building_block_truncation(
        self,
        ancestor: Compound,
        descendant: Compound
    ) -> bool:
        """
        Check if descendant is direct building-block truncation of ancestor.

        A direct truncation means exactly one non-null building block
        in ancestor is replaced with null in descendant, all other positions match.

        Parameters
        ----------
        ancestor : Compound
            Potential ancestor
        descendant : Compound
            Potential descendant

        Returns
        -------
        bool
            True if descendant differs from ancestor by exactly one null substitution

        Examples
        --------
        >>> # Leu-Pro-Val → Leu-Null-Val (direct truncation)
        >>> ancestor = Compound([bb_leu, bb_pro, bb_val], chrom)
        >>> descendant = Compound([bb_leu, bb_null, bb_val], chrom)
        >>> builder._is_building_block_truncation(ancestor, descendant)
        True

        >>> # Leu-Pro-Val → Null-Null-Val (not direct - 2 removals)
        >>> descendant2 = Compound([bb_null, bb_null, bb_val], chrom)
        >>> builder._is_building_block_truncation(ancestor, descendant2)
        False
        """
        ancestor_blocks = ancestor.building_blocks
        descendant_blocks = descendant.building_blocks

        # Must have same number of positions
        if len(ancestor_blocks) != len(descendant_blocks):
            return False

        # Count differences
        differences = 0
        for anc_block, desc_block in zip(ancestor_blocks, descendant_blocks):
            # Must be same cycle
            if anc_block.cycle != desc_block.cycle:
                return False

            # Check if this position differs
            if anc_block.is_null and desc_block.is_null:
                # Both null - no difference
                continue
            elif anc_block.is_null and not desc_block.is_null:
                # Ancestor is null but descendant is not - invalid truncation
                return False
            elif not anc_block.is_null and desc_block.is_null:
                # Ancestor has block, descendant is null - this is a truncation
                differences += 1
            else:
                # Both non-null - must be same block
                if anc_block.code != desc_block.code:
                    return False

        # Direct truncation = exactly 1 difference
        return differences == 1

    def _is_monomer_truncation(
        self,
        ancestor: Compound,
        descendant: Compound
    ) -> bool:
        """
        Check if descendant is direct monomer-level truncation of ancestor.

        A direct monomer truncation means exactly one monomer is removed
        from ancestor's monomer sequence to yield descendant's monomer sequence.

        Parameters
        ----------
        ancestor : Compound
            Potential ancestor
        descendant : Compound
            Potential descendant

        Returns
        -------
        bool
            True if descendant's monomer sequence = ancestor's with one monomer removed

        Notes
        -----
        - Decomposes building blocks to monomer sequences
        - Checks if descendant is subsequence with exactly 1 monomer missing
        - Order-preserving subsequence check
        - Handles convergence (multiple ancestors → same descendant)

        Examples
        --------
        >>> # Leu-Ala-Val → Leu-Val (removed Ala)
        >>> ancestor_mono = "Leu-Ala-Val"
        >>> descendant_mono = "Leu-Val"
        >>> # This is a direct monomer truncation

        References
        ----------
        THEORY.md Section 1.5.3: Monomer-Level Truncation
        THEORY.md Section 1.3: Monomer Mode (DAG with Convergence)
        """
        # Get monomer sequences as strings
        ancestor_mono_str = ancestor.monomer_sequence
        descendant_mono_str = descendant.monomer_sequence

        # Handle empty sequences
        if not ancestor_mono_str:
            return False
        if not descendant_mono_str and ancestor_mono_str.count("-") == 0:
            # Ancestor has 1 monomer, descendant has 0
            return True
        if not descendant_mono_str:
            return False

        # Split into monomer lists
        ancestor_monomers = ancestor_mono_str.split("-")
        descendant_monomers = descendant_mono_str.split("-")

        # Descendant must have exactly 1 fewer monomer
        if len(ancestor_monomers) != len(descendant_monomers) + 1:
            return False

        # Check if descendant is ancestor with exactly one monomer removed
        # Try removing each monomer from ancestor and see if it matches descendant
        for i in range(len(ancestor_monomers)):
            # Create ancestor sequence with monomer at position i removed
            candidate = ancestor_monomers[:i] + ancestor_monomers[i+1:]

            # Check if this matches descendant
            if candidate == descendant_monomers:
                return True

        return False

    def _is_building_block_descendant(
        self,
        ancestor: Compound,
        descendant: Compound
    ) -> bool:
        """
        Check if descendant is a valid building-block truncation of ancestor (any level).

        Allows multiple nullification steps (handles gaps in dataset).
        """
        ancestor_blocks = ancestor.building_blocks
        descendant_blocks = descendant.building_blocks

        # Must have same number of positions
        if len(ancestor_blocks) != len(descendant_blocks):
            return False

        # Check each position
        for anc_block, desc_block in zip(ancestor_blocks, descendant_blocks):
            # Must be same cycle
            if anc_block.cycle != desc_block.cycle:
                return False

            # Check truncation validity
            if anc_block.is_null and not desc_block.is_null:
                # Ancestor is null but descendant is not - invalid
                return False

            if not anc_block.is_null and not desc_block.is_null:
                # Both non-null - must be same block
                if anc_block.code != desc_block.code:
                    return False

        # Valid if descendant has at least one more null than ancestor
        anc_nulls = sum(1 for b in ancestor_blocks if b.is_null)
        desc_nulls = sum(1 for b in descendant_blocks if b.is_null)
        return desc_nulls > anc_nulls

    def _is_monomer_descendant(
        self,
        ancestor: Compound,
        descendant: Compound
    ) -> bool:
        """
        Check if descendant is a valid monomer truncation of ancestor (any level).

        Checks if descendant's monomer sequence is a subsequence of ancestor's
        (allows multiple monomer removals for gaps in dataset).
        """
        # Get monomer sequences
        ancestor_mono_str = ancestor.monomer_sequence
        descendant_mono_str = descendant.monomer_sequence

        # Handle empty sequences
        if not ancestor_mono_str:
            return False
        if not descendant_mono_str:
            # Descendant has no monomers, ancestor has some - valid truncation
            return True

        # Split into monomer lists
        ancestor_monomers = ancestor_mono_str.split("-")
        descendant_monomers = descendant_mono_str.split("-")

        # Descendant must have fewer monomers
        if len(descendant_monomers) >= len(ancestor_monomers):
            return False

        # Check if descendant is a subsequence of ancestor (order-preserving)
        anc_idx = 0
        desc_idx = 0

        while anc_idx < len(ancestor_monomers) and desc_idx < len(descendant_monomers):
            if ancestor_monomers[anc_idx] == descendant_monomers[desc_idx]:
                desc_idx += 1
            anc_idx += 1

        # Valid if all descendant monomers were matched in order
        return desc_idx == len(descendant_monomers)
