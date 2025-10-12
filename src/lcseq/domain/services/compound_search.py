"""Compound search domain service.

This module provides domain services for searching and filtering compounds
from a library based on various criteria (sequence, ID, building blocks, etc.).

Domain Responsibility:
    - Search for compounds by sequence (positional, residue, monomer)
    - Filter compounds by properties (level, building blocks, etc.)
    - Pure algorithmic logic with no I/O or presentation concerns

Not Responsible For:
    - Loading data from files (that's infrastructure)
    - Displaying search results (that's presentation)
    - Caching or indexing (that's infrastructure optimization)
"""

from typing import List, Optional
from ..entities import Compound


class CompoundSearchService:
    """Domain service for searching and filtering compounds.

    This service provides algorithms for finding compounds in a library
    based on various search criteria. All methods are deterministic and
    side-effect free.

    Examples:
        >>> searcher = CompoundSearchService()
        >>> reference = searcher.find_by_sequence(compounds, "Leu-Pro-Ala")
        >>> # Returns compound matching the sequence
    """

    def find_by_sequence(
        self,
        compounds: List[Compound],
        sequence: str,
        normalize_null: bool = True
    ) -> Optional[Compound]:
        """Find compound matching a given sequence.

        Searches for exact matches in:
        1. Positional block sequence (with building block positions)
        2. Block support sequence (without nulls)
        3. Normalized sequences (AgxNull → Null)

        Parameters
        ----------
        compounds : List[Compound]
            Library of compounds to search
        sequence : str
            Target sequence (e.g., "Leu-Pro-Ala" or "Leu-Null-Ala")
        normalize_null : bool, optional
            Normalize null representations (AgxNull → Null). Default is True.

        Returns
        -------
        Optional[Compound]
            First matching compound, or None if not found

        Examples
        --------
        >>> searcher = CompoundSearchService()
        >>> compound = searcher.find_by_sequence(compounds, "Leu-Pro-Ala")
        >>> compound = searcher.find_by_sequence(compounds, "Leu-Null-Ala")
        """
        sequence_parts = sequence.split("-")

        for compound in compounds:
            # Check positional block sequence (exact match)
            if compound.positional_block_sequence == sequence:
                return compound

            # Check block support sequence (match without nulls)
            if compound.block_support_sequence == sequence:
                return compound

            # Check with null normalization
            if normalize_null:
                cpd_parts = compound.positional_block_sequence.split("-")
                if len(cpd_parts) == len(sequence_parts):
                    match = all(
                        p1 == p2
                        or (p1.lower().replace("agxnull", "null") == p2.lower().replace("agxnull", "null"))
                        for p1, p2 in zip(cpd_parts, sequence_parts)
                    )
                    if match:
                        return compound

        return None

    def filter_by_level_range(
        self,
        compounds: List[Compound],
        min_level: Optional[int] = None,
        max_level: Optional[int] = None,
        use_monomer_level: bool = False
    ) -> List[Compound]:
        """Filter compounds by truncation level range.

        Parameters
        ----------
        compounds : List[Compound]
            Library of compounds to filter
        min_level : int, optional
            Minimum level (inclusive)
        max_level : int, optional
            Maximum level (inclusive)
        use_monomer_level : bool, optional
            Use monomer level instead of building block level. Default is False.

        Returns
        -------
        List[Compound]
            Filtered list of compounds

        Examples
        --------
        >>> searcher = CompoundSearchService()
        >>> truncated = searcher.filter_by_level_range(compounds, max_level=2)
        >>> # Returns compounds with level <= 2
        """
        result = []

        for compound in compounds:
            level = compound.monomer_level if use_monomer_level else compound.level

            if min_level is not None and level < min_level:
                continue
            if max_level is not None and level > max_level:
                continue

            result.append(compound)

        return result

    def filter_potential_descendants(
        self,
        compounds: List[Compound],
        reference: Compound
    ) -> List[Compound]:
        """Filter compounds to find potential descendants of a reference compound.

        In building-block mode, a descendant has:
        - Same or fewer building blocks (level <= reference.level)
        - At each position: same building block OR Null

        This is an optimization for building-block hierarchy construction.

        Parameters
        ----------
        compounds : List[Compound]
            Library of compounds to filter
        reference : Compound
            Reference compound to find descendants for - THEORY.md Section 3.1:
            "The compound currently being analyzed"

        Returns
        -------
        List[Compound]
            List containing reference + potential descendants

        Notes
        -----
        This optimization only works for building-block mode. In monomer mode,
        you need the full library due to convergence patterns.

        References
        ----------
        THEORY.md Section 3.1: Terminology - Ancestry and Lineage

        Examples
        --------
        >>> searcher = CompoundSearchService()
        >>> lineage_candidates = searcher.filter_potential_descendants(compounds, reference)
        >>> # Returns reference + compounds that could be descendants
        """
        potential_descendants = [reference]
        reference_bbs = reference.building_blocks

        for compound in compounds:
            if compound == reference:
                continue

            # Descendants have same or more truncation (lower or equal level)
            if compound.level > reference.level:
                continue

            # Check building block compatibility
            cpd_bbs = compound.building_blocks

            if len(cpd_bbs) != len(reference_bbs):
                continue

            # Each position must match reference OR be Null
            is_potential_descendant = True
            for i in range(len(reference_bbs)):
                reference_bb = reference_bbs[i]
                cpd_bb = cpd_bbs[i]

                # If reference has a building block at this position
                if not reference_bb.is_null:
                    # Compound must have same building block OR Null
                    if not cpd_bb.is_null:
                        # Both have building blocks - must match exactly
                        if cpd_bb.code != reference_bb.code:
                            is_potential_descendant = False
                            break

            if is_potential_descendant:
                potential_descendants.append(compound)

        return potential_descendants
