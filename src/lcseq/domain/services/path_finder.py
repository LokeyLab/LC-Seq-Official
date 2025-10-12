"""
PathFinder - Finds paths in hierarchy DAG.

Implementation based on THEORY.md Section 3.1, 3.2.
"""

from typing import List, Optional
from ..entities.compound import Compound
from ..models.compound_hierarchy import CompoundHierarchy


class PathFinder:
    """
    Finds paths in hierarchy DAG.

    Provides algorithms for finding paths between compounds in the
    truncation hierarchy. Supports single path, all paths, and
    shortest path queries.

    Notes
    -----
    - Stateless service (no instance state)
    - Pure domain logic (no I/O operations)
    - Uses depth-first search for path finding
    - Handles both building-block and monomer modes

    Examples
    --------
    >>> from lcseq.domain.entities.building_block import BuildingBlock
    >>> from lcseq.domain.entities.chromatogram import Chromatogram
    >>> from lcseq.domain.models.compound_hierarchy import CompoundHierarchy, HierarchyMode
    >>> import numpy as np
    >>>
    >>> # Create hierarchy
    >>> hierarchy = CompoundHierarchy(mode=HierarchyMode.BUILDING_BLOCK)
    >>> chromatogram = Chromatogram(
    ...     time_points=np.array([1.0, 2.0, 3.0]),
    ...     counts=np.array([100.0, 200.0, 150.0])
    ... )
    >>>
    >>> # Add compounds and edges
    >>> c1 = Compound([bb1, bb2, bb3], chromatogram)
    >>> c2 = Compound([bb1, bb2, Null], chromatogram)
    >>> c3 = Compound([bb1, Null, Null], chromatogram)
    >>> hierarchy.add_compound(c1)
    >>> hierarchy.add_compound(c2)
    >>> hierarchy.add_compound(c3)
    >>> hierarchy.add_edge(c1, c2)
    >>> hierarchy.add_edge(c2, c3)
    >>>
    >>> # Find path
    >>> finder = PathFinder()
    >>> path = finder.find_path(hierarchy, c1, c3)
    >>> path
    [c1, c2, c3]

    References
    ----------
    THEORY.md Section 3.1: Terminology - Ancestry and Lineage
    THEORY.md Section 3.2: Mathematical Poset Terminology
    """

    def find_path(
        self,
        hierarchy: CompoundHierarchy,
        start: Compound,
        end: Compound
    ) -> Optional[List[Compound]]:
        """
        Find a path from start to end compound.

        Parameters
        ----------
        hierarchy : CompoundHierarchy
            The hierarchy to search
        start : Compound
            Starting compound
        end : Compound
            Target compound

        Returns
        -------
        Optional[List[Compound]]
            Path from start to end (inclusive), or None if no path exists

        Notes
        -----
        - Returns first path found (may not be shortest)
        - Path follows directed edges (start must be ancestor of end)
        - Returns None if no path exists
        - Includes both start and end in returned path

        Examples
        --------
        >>> path = finder.find_path(hierarchy, maximal, minimal)
        >>> path[0] == maximal
        True
        >>> path[-1] == minimal
        True
        """
        if start not in hierarchy.compounds or end not in hierarchy.compounds:
            return None

        # Special case: start == end
        if start == end:
            return [start]

        # DFS to find path
        return self._dfs_path(hierarchy, start, end, visited=[])

    def _dfs_path(
        self,
        hierarchy: CompoundHierarchy,
        current: Compound,
        target: Compound,
        visited: List[Compound]
    ) -> Optional[List[Compound]]:
        """
        Depth-first search to find path.

        Parameters
        ----------
        hierarchy : CompoundHierarchy
            The hierarchy
        current : Compound
            Current position in search
        target : Compound
            Target compound
        visited : List[Compound]
            Already visited compounds (cycle detection)

        Returns
        -------
        Optional[List[Compound]]
            Path if found, None otherwise
        """
        # Mark current as visited
        visited_copy = visited + [current]

        # Check if we reached target
        if current == target:
            return visited_copy

        # Try each descendant
        descendants = hierarchy.get_direct_descendants(current)
        for descendant in descendants:
            if descendant not in visited_copy:
                path = self._dfs_path(hierarchy, descendant, target, visited_copy)
                if path is not None:
                    return path

        return None

    def find_all_paths(
        self,
        hierarchy: CompoundHierarchy,
        start: Compound,
        end: Compound
    ) -> List[List[Compound]]:
        """
        Find all paths from start to end compound.

        Parameters
        ----------
        hierarchy : CompoundHierarchy
            The hierarchy to search
        start : Compound
            Starting compound
        end : Compound
            Target compound

        Returns
        -------
        List[List[Compound]]
            All paths from start to end, empty list if none exist

        Notes
        -----
        - Returns all possible paths (not just one)
        - Useful for analyzing convergence patterns
        - In building-block mode: typically one path per pair
        - In monomer mode: may have multiple paths (diamond patterns)

        Examples
        --------
        >>> # In monomer mode, multiple paths may exist
        >>> paths = finder.find_all_paths(hierarchy, maximal, minimal)
        >>> len(paths)
        3
        >>> # Three different truncation sequences lead to same minimal

        References
        ----------
        THEORY.md Section 1.3: Monomer Mode (DAG with Convergence)
        """
        if start not in hierarchy.compounds or end not in hierarchy.compounds:
            return []

        # Special case: start == end
        if start == end:
            return [[start]]

        # DFS to find all paths
        all_paths = []
        self._dfs_all_paths(hierarchy, start, end, [], all_paths)
        return all_paths

    def _dfs_all_paths(
        self,
        hierarchy: CompoundHierarchy,
        current: Compound,
        target: Compound,
        visited: List[Compound],
        all_paths: List[List[Compound]]
    ) -> None:
        """
        Depth-first search to find all paths.

        Parameters
        ----------
        hierarchy : CompoundHierarchy
            The hierarchy
        current : Compound
            Current position in search
        target : Compound
            Target compound
        visited : List[Compound]
            Current path being explored
        all_paths : List[List[Compound]]
            Accumulator for all found paths (modified in place)
        """
        # Add current to path
        visited_copy = visited + [current]

        # Check if we reached target
        if current == target:
            all_paths.append(visited_copy)
            return

        # Try each descendant
        descendants = hierarchy.get_direct_descendants(current)
        for descendant in descendants:
            if descendant not in visited_copy:
                self._dfs_all_paths(hierarchy, descendant, target, visited_copy, all_paths)

    def shortest_path(
        self,
        hierarchy: CompoundHierarchy,
        start: Compound,
        end: Compound
    ) -> Optional[List[Compound]]:
        """
        Find shortest path from start to end compound.

        Parameters
        ----------
        hierarchy : CompoundHierarchy
            The hierarchy to search
        start : Compound
            Starting compound
        end : Compound
            Target compound

        Returns
        -------
        Optional[List[Compound]]
            Shortest path from start to end, or None if no path exists

        Notes
        -----
        - Uses breadth-first search for shortest path
        - In DAG, shortest path = fewest edges
        - Useful for finding minimum truncation steps

        Examples
        --------
        >>> # Find minimum truncation steps
        >>> path = finder.shortest_path(hierarchy, maximal, minimal)
        >>> truncation_steps = len(path) - 1
        >>> truncation_steps
        3
        """
        if start not in hierarchy.compounds or end not in hierarchy.compounds:
            return None

        # Special case: start == end
        if start == end:
            return [start]

        # BFS for shortest path
        queue = [[start]]
        visited = {start}

        while queue:
            path = queue.pop(0)
            current = path[-1]

            # Get descendants
            descendants = hierarchy.get_direct_descendants(current)
            for descendant in descendants:
                if descendant == end:
                    return path + [descendant]

                if descendant not in visited:
                    visited.add(descendant)
                    queue.append(path + [descendant])

        return None
