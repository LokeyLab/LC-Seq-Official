"""
LevelAnalyzer - Analyzes hierarchy by truncation levels.

Implementation based on THEORY.md Section 3.3.
"""

from typing import Dict, List
from ..entities.compound import Compound
from ..models.compound_hierarchy import CompoundHierarchy, HierarchyMode


class LevelAnalyzer:
    """
    Analyzes hierarchy by truncation levels.

    Provides utilities for querying and analyzing compounds by their
    truncation level (number of non-null building blocks or monomers).

    Notes
    -----
    - Stateless service (no instance state)
    - Pure domain logic (no I/O operations)
    - Works with both building-block and monomer modes
    - Level = number of non-null building blocks (or monomers in monomer mode)

    Examples
    --------
    >>> from lcseq.domain.services.hierarchy_builder import HierarchyBuilder
    >>> from lcseq.domain.models.compound_hierarchy import HierarchyMode
    >>>
    >>> # Build hierarchy
    >>> builder = HierarchyBuilder()
    >>> hierarchy = builder.build(compounds, HierarchyMode.BUILDING_BLOCK)
    >>>
    >>> # Analyze levels
    >>> analyzer = LevelAnalyzer()
    >>> distribution = analyzer.get_level_distribution(hierarchy)
    >>> distribution
    {0: 1, 1: 3, 2: 3, 3: 1}
    >>> # 1 compound at level 0, 3 at level 1, etc.
    >>>
    >>> max_level = analyzer.get_max_level(hierarchy)
    >>> max_level
    3

    References
    ----------
    THEORY.md Section 3.3: Hierarchy Properties (Truncation Level)
    """

    def get_level_distribution(
        self,
        hierarchy: CompoundHierarchy
    ) -> Dict[int, int]:
        """
        Get distribution of compounds by level.

        Parameters
        ----------
        hierarchy : CompoundHierarchy
            The hierarchy to analyze

        Returns
        -------
        Dict[int, int]
            Mapping from level to count of compounds at that level

        Notes
        -----
        - Level 0 = all null (minimal compound)
        - Higher levels = more building blocks
        - Empty dict if hierarchy is empty

        Examples
        --------
        >>> distribution = analyzer.get_level_distribution(hierarchy)
        >>> distribution[2]
        5
        >>> # 5 compounds have exactly 2 non-null building blocks
        """
        distribution: Dict[int, int] = {}

        for compound in hierarchy.compounds:
            level = self._get_compound_level(compound, hierarchy.mode)

            if level not in distribution:
                distribution[level] = 0

            distribution[level] += 1

        return distribution

    def get_max_level(self, hierarchy: CompoundHierarchy) -> int:
        """
        Get maximum truncation level in hierarchy.

        Parameters
        ----------
        hierarchy : CompoundHierarchy
            The hierarchy to analyze

        Returns
        -------
        int
            Maximum level (most building blocks/monomers)

        Raises
        ------
        ValueError
            If hierarchy is empty

        Notes
        -----
        - Max level = maximal compounds in hierarchy
        - In building-block mode: most non-null building blocks
        - In monomer mode: most monomers

        Examples
        --------
        >>> max_level = analyzer.get_max_level(hierarchy)
        >>> max_level
        5
        >>> # Longest sequence has 5 building blocks
        """
        if not hierarchy.compounds:
            raise ValueError("Cannot get max level from empty hierarchy")

        levels = [
            self._get_compound_level(c, hierarchy.mode)
            for c in hierarchy.compounds
        ]

        return max(levels)

    def get_min_level(self, hierarchy: CompoundHierarchy) -> int:
        """
        Get minimum truncation level in hierarchy.

        Parameters
        ----------
        hierarchy : CompoundHierarchy
            The hierarchy to analyze

        Returns
        -------
        int
            Minimum level (fewest building blocks/monomers)

        Raises
        ------
        ValueError
            If hierarchy is empty

        Notes
        -----
        - Min level = minimal compounds in hierarchy
        - Often 0 (all-null compound L₀)

        Examples
        --------
        >>> min_level = analyzer.get_min_level(hierarchy)
        >>> min_level
        0
        >>> # L₀ compound present (all null)
        """
        if not hierarchy.compounds:
            raise ValueError("Cannot get min level from empty hierarchy")

        levels = [
            self._get_compound_level(c, hierarchy.mode)
            for c in hierarchy.compounds
        ]

        return min(levels)

    def get_compounds_at_level(
        self,
        hierarchy: CompoundHierarchy,
        level: int
    ) -> List[Compound]:
        """
        Get all compounds at a specific level.

        Parameters
        ----------
        hierarchy : CompoundHierarchy
            The hierarchy to query
        level : int
            Target level

        Returns
        -------
        List[Compound]
            All compounds at the specified level

        Notes
        -----
        - Returns empty list if no compounds at level
        - Level must be non-negative

        Examples
        --------
        >>> # Get all compounds with 2 building blocks
        >>> compounds = analyzer.get_compounds_at_level(hierarchy, 2)
        >>> len(compounds)
        3
        """
        if level < 0:
            return []

        return [
            c for c in hierarchy.compounds
            if self._get_compound_level(c, hierarchy.mode) == level
        ]

    def get_compounds_by_level_range(
        self,
        hierarchy: CompoundHierarchy,
        min_level: int,
        max_level: int
    ) -> List[Compound]:
        """
        Get compounds within a level range.

        Parameters
        ----------
        hierarchy : CompoundHierarchy
            The hierarchy to query
        min_level : int
            Minimum level (inclusive)
        max_level : int
            Maximum level (inclusive)

        Returns
        -------
        List[Compound]
            All compounds with level in [min_level, max_level]

        Notes
        -----
        - Range is inclusive on both ends
        - Returns empty list if no compounds in range

        Examples
        --------
        >>> # Get compounds with 2-4 building blocks
        >>> compounds = analyzer.get_compounds_by_level_range(hierarchy, 2, 4)
        >>> all(2 <= c.level <= 4 for c in compounds)
        True
        """
        if min_level < 0:
            min_level = 0

        return [
            c for c in hierarchy.compounds
            if min_level <= self._get_compound_level(c, hierarchy.mode) <= max_level
        ]

    def get_level_statistics(
        self,
        hierarchy: CompoundHierarchy
    ) -> Dict[str, float]:
        """
        Get statistical summary of levels.

        Parameters
        ----------
        hierarchy : CompoundHierarchy
            The hierarchy to analyze

        Returns
        -------
        Dict[str, float]
            Statistics including min, max, mean, median level

        Raises
        ------
        ValueError
            If hierarchy is empty

        Examples
        --------
        >>> stats = analyzer.get_level_statistics(hierarchy)
        >>> stats
        {
            'min': 0,
            'max': 5,
            'mean': 2.5,
            'median': 2.0,
            'count': 100
        }
        """
        if not hierarchy.compounds:
            raise ValueError("Cannot compute statistics from empty hierarchy")

        levels = [
            self._get_compound_level(c, hierarchy.mode)
            for c in hierarchy.compounds
        ]

        levels_sorted = sorted(levels)
        n = len(levels_sorted)

        return {
            'min': float(min(levels)),
            'max': float(max(levels)),
            'mean': sum(levels) / n,
            'median': levels_sorted[n // 2] if n % 2 == 1 else
                     (levels_sorted[n // 2 - 1] + levels_sorted[n // 2]) / 2.0,
            'count': float(n)
        }

    def _get_compound_level(
        self,
        compound: Compound,
        mode: HierarchyMode
    ) -> int:
        """
        Get level of a compound based on mode.

        Parameters
        ----------
        compound : Compound
            The compound
        mode : HierarchyMode
            Building-block or monomer mode

        Returns
        -------
        int
            Truncation level

        Notes
        -----
        - Building-block mode: uses compound.level
        - Monomer mode: uses compound.monomer_level
        """
        if mode == HierarchyMode.BUILDING_BLOCK:
            return compound.level
        else:  # MONOMER
            return compound.monomer_level
