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
