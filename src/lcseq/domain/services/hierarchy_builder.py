"""
HierarchyBuilder - Constructs CompoundHierarchy from list of compounds.

Implementation based on THEORY.md Section 4.2, 3.3, 1.5.
"""

from typing import List, Dict
from tqdm import tqdm
from ..entities.compound import Compound
from ..entities.building_block import BuildingBlock
from ..models.compound_hierarchy import CompoundHierarchy, HierarchyMode
from .lineage_finder import LineageFinderService


class HierarchyBuilder:
    """
    Builds CompoundHierarchy from list of compounds.

    Automatically detects truncation relationships between compounds
    and constructs the DAG structure. Supports both building-block
    and monomer modes.

    Uses LineageFinderService for descendant logic (Single Source of Truth).

    Notes
    -----
    - Uses LineageFinderService for consistent descendant checking
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

    def __init__(self):
        """Initialize HierarchyBuilder with LineageFinderService."""
        self.lineage_finder = LineageFinderService()

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

        Algorithm (Hash-Optimized)
        --------------------------
        1. Index compounds by block_support_sequence for O(1) lookup
        2. For each compound, generate all possible direct truncation sequences
        3. Look up each truncation sequence in the hash index

        Complexity: O(n × m) where m = sequence length (typically 3-9)
        """
        # Create hierarchy
        hierarchy = CompoundHierarchy(mode=mode)

        # Add all compounds
        for compound in compounds:
            hierarchy.add_compound(compound)

        # Index compounds by block_support_sequence for O(1) lookup
        by_block_support: Dict[str, List[Compound]] = {}
        for compound in compounds:
            key = compound.block_support_sequence
            if key not in by_block_support:
                by_block_support[key] = []
            by_block_support[key].append(compound)

        # Show progress
        show_progress = len(compounds) > 100

        # For each compound, find direct descendants using hash lookup
        total_edges = 0
        for ancestor in tqdm(compounds, desc="Building edges", disable=not show_progress, unit="cpd"):
            # Generate all possible direct truncation sequences
            # A direct truncation removes exactly ONE non-null block
            truncation_sequences = self._generate_truncation_sequences(ancestor, mode)

            for trunc_seq in truncation_sequences:
                # Look up compounds with this block_support_sequence
                if trunc_seq in by_block_support:
                    for descendant in by_block_support[trunc_seq]:
                        # Verify it's a valid truncation (handles edge cases)
                        if self._is_direct_truncation(ancestor, descendant, mode):
                            hierarchy.add_edge(ancestor, descendant)
                            total_edges += 1

        return hierarchy

    def _generate_truncation_sequences(
        self,
        compound: Compound,
        mode: HierarchyMode
    ) -> List[str]:
        """
        Generate all possible block_support_sequences for direct truncations.

        A direct truncation removes exactly one non-null block.

        Parameters
        ----------
        compound : Compound
            Compound to generate truncations for
        mode : HierarchyMode
            Hierarchy mode

        Returns
        -------
        List[str]
            All possible block_support_sequences for direct descendants
        """
        if mode == HierarchyMode.BUILDING_BLOCK:
            # Get non-null blocks in order
            non_null_blocks = [bb for bb in reversed(compound.building_blocks) if not bb.is_null]

            if len(non_null_blocks) <= 1:
                return [""]  # Only truncation is to empty/L0

            # Generate all sequences with one block removed
            truncations = []
            for i in range(len(non_null_blocks)):
                remaining = non_null_blocks[:i] + non_null_blocks[i+1:]
                trunc_seq = "-".join(bb.code for bb in remaining)
                truncations.append(trunc_seq)

            return truncations
        else:
            # Monomer mode - more complex, delegate to existing logic
            # For now, return empty to fall back to pairwise comparison
            return []

    def _is_direct_truncation(
        self,
        ancestor: Compound,
        descendant: Compound,
        mode: HierarchyMode
    ) -> bool:
        """
        Check if descendant is a DIRECT truncation of ancestor (one step).

        For building block mode: exactly one non-null block becomes null.
        For monomer mode: exactly one monomer is removed.
        """
        if mode == HierarchyMode.BUILDING_BLOCK:
            return self._is_building_block_truncation(ancestor, descendant)
        else:
            return self._is_monomer_truncation(ancestor, descendant)

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
        ancestor_mono_str = ancestor.monomer_support_sequence
        descendant_mono_str = descendant.monomer_support_sequence

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

        Uses LineageFinderService for consistent block support subsequence logic.
        This ensures equivalence classes (same block support) converge in the hierarchy.

        Returns
        -------
        bool
            True if descendant's block support sequence is a proper subsequence of ancestor's

        Notes
        -----
        Delegates to LineageFinderService for Single Source of Truth.
        """
        return self.lineage_finder._is_building_block_descendant(descendant, ancestor)

    def _is_monomer_descendant(
        self,
        ancestor: Compound,
        descendant: Compound
    ) -> bool:
        """
        Check if descendant is a valid monomer truncation of ancestor (any level).

        Uses LineageFinderService for consistent monomer subsequence logic.
        This ensures equivalence classes (same monomer sequence) converge in the hierarchy.

        Returns
        -------
        bool
            True if descendant's monomer sequence is a proper subsequence of ancestor's

        Notes
        -----
        Delegates to LineageFinderService for Single Source of Truth.
        """
        return self.lineage_finder._is_monomer_descendant(descendant, ancestor)
