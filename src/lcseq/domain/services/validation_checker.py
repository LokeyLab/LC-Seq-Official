"""
ValidationChecker - Validates hierarchy structural properties.

Implementation based on THEORY.md Section 1.1, 3.3.
"""

from typing import Set, List, Dict
from ..entities.compound import Compound
from ..models.compound_hierarchy import CompoundHierarchy, HierarchyMode


class ValidationChecker:
    """
    Validates hierarchy structural properties.

    Provides validation methods to ensure hierarchy maintains required
    mathematical properties: acyclic, proper level ordering, valid DAG.

    Notes
    -----
    - Stateless service (no instance state)
    - Pure domain logic (no I/O operations)
    - Validates DAG invariants
    - Checks structural consistency

    Examples
    --------
    >>> from lcseq.domain.services.hierarchy_builder import HierarchyBuilder
    >>> from lcseq.domain.models.compound_hierarchy import HierarchyMode
    >>>
    >>> # Build hierarchy
    >>> builder = HierarchyBuilder()
    >>> hierarchy = builder.build(compounds, HierarchyMode.BUILDING_BLOCK)
    >>>
    >>> # Validate structure
    >>> checker = ValidationChecker()
    >>> is_valid = checker.is_valid_dag(hierarchy)
    >>> is_valid
    True
    >>>
    >>> has_cycles = checker.has_cycles(hierarchy)
    >>> has_cycles
    False

    References
    ----------
    THEORY.md Section 1.1: Core Mathematical Model (DAG properties)
    THEORY.md Section 1.4: Why DAG Structure Holds (No Cycles)
    THEORY.md Section 3.3: Hierarchy Properties (Level ordering)
    """

    def is_valid_dag(self, hierarchy: CompoundHierarchy) -> bool:
        """
        Check if hierarchy is a valid DAG.

        Parameters
        ----------
        hierarchy : CompoundHierarchy
            The hierarchy to validate

        Returns
        -------
        bool
            True if hierarchy is valid DAG, False otherwise

        Notes
        -----
        Validates:
        - Acyclic (no cycles)
        - Proper level ordering (descendant level < ancestor level)
        - All edges connect compounds in hierarchy

        Examples
        --------
        >>> is_valid = checker.is_valid_dag(hierarchy)
        >>> is_valid
        True
        """
        # Check for cycles
        if self.has_cycles(hierarchy):
            return False

        # Check level ordering
        if not self.validate_level_ordering(hierarchy):
            return False

        # Check edge validity
        if not self._validate_edges(hierarchy):
            return False

        return True

    def has_cycles(self, hierarchy: CompoundHierarchy) -> bool:
        """
        Check if hierarchy contains cycles.

        Parameters
        ----------
        hierarchy : CompoundHierarchy
            The hierarchy to check

        Returns
        -------
        bool
            True if cycles exist, False otherwise

        Notes
        -----
        Uses depth-first search with three-color marking:
        - White: not visited
        - Gray: currently visiting (in recursion stack)
        - Black: completely visited

        A cycle exists if we encounter a gray node during DFS.

        References
        ----------
        THEORY.md Section 1.4: Why DAG Structure Holds (No Cycles)

        Examples
        --------
        >>> has_cycles = checker.has_cycles(hierarchy)
        >>> has_cycles
        False
        >>> # Valid hierarchy has no cycles
        """
        # Three-color DFS
        WHITE = 0
        GRAY = 1
        BLACK = 2

        colors: Dict[Compound, int] = {c: WHITE for c in hierarchy.compounds}

        def dfs(compound: Compound) -> bool:
            """Returns True if cycle detected."""
            colors[compound] = GRAY

            # Visit descendants
            descendants = hierarchy.get_direct_descendants(compound)
            for descendant in descendants:
                if colors[descendant] == GRAY:
                    # Found back edge - cycle detected
                    return True
                if colors[descendant] == WHITE:
                    if dfs(descendant):
                        return True

            colors[compound] = BLACK
            return False

        # Check all components
        for compound in hierarchy.compounds:
            if colors[compound] == WHITE:
                if dfs(compound):
                    return True

        return False

    def validate_level_ordering(self, hierarchy: CompoundHierarchy) -> bool:
        """
        Validate that descendant level < ancestor level for all edges.

        Parameters
        ----------
        hierarchy : CompoundHierarchy
            The hierarchy to validate

        Returns
        -------
        bool
            True if all edges respect level ordering, False otherwise

        Notes
        -----
        - In building-block mode: checks compound.level
        - In monomer mode: checks compound.monomer_level
        - Descendant must have strictly lower level than ancestor

        References
        ----------
        THEORY.md Section 3.3: Hierarchy Properties (Level ordering)

        Examples
        --------
        >>> valid = checker.validate_level_ordering(hierarchy)
        >>> valid
        True
        >>> # All descendants have lower levels than ancestors
        """
        level_attr = "level" if hierarchy.mode == HierarchyMode.BUILDING_BLOCK else "monomer_level"

        for ancestor in hierarchy.compounds:
            ancestor_level = getattr(ancestor, level_attr)
            descendants = hierarchy.get_direct_descendants(ancestor)

            for descendant in descendants:
                descendant_level = getattr(descendant, level_attr)

                if descendant_level >= ancestor_level:
                    return False

        return True

    def find_cycles(self, hierarchy: CompoundHierarchy) -> List[List[Compound]]:
        """
        Find all cycles in hierarchy.

        Parameters
        ----------
        hierarchy : CompoundHierarchy
            The hierarchy to check

        Returns
        -------
        List[List[Compound]]
            List of cycles found (empty if DAG is valid)

        Notes
        -----
        - Returns list of compound sequences forming cycles
        - Empty list indicates valid DAG
        - Used for debugging invalid hierarchies

        Examples
        --------
        >>> cycles = checker.find_cycles(hierarchy)
        >>> cycles
        []
        >>> # No cycles in valid hierarchy
        """
        cycles: List[List[Compound]] = []
        visited: Set[Compound] = set()
        rec_stack: List[Compound] = []

        def dfs(compound: Compound) -> None:
            """DFS to find cycles."""
            visited.add(compound)
            rec_stack.append(compound)

            descendants = hierarchy.get_direct_descendants(compound)
            for descendant in descendants:
                if descendant in rec_stack:
                    # Found cycle - extract it
                    cycle_start = rec_stack.index(descendant)
                    cycle = rec_stack[cycle_start:] + [descendant]
                    cycles.append(cycle)
                elif descendant not in visited:
                    dfs(descendant)

            rec_stack.pop()

        # Check all components
        for compound in hierarchy.compounds:
            if compound not in visited:
                dfs(compound)

        return cycles

    def validate_connectivity(self, hierarchy: CompoundHierarchy) -> bool:
        """
        Check if all compounds are in connected components.

        Parameters
        ----------
        hierarchy : CompoundHierarchy
            The hierarchy to validate

        Returns
        -------
        bool
            True if no isolated compounds exist, False otherwise

        Notes
        -----
        - Isolated compound = no ancestors and no descendants
        - Valid hierarchies may have multiple connected components (forest)
        - This checks if each compound participates in hierarchy

        Examples
        --------
        >>> connected = checker.validate_connectivity(hierarchy)
        >>> connected
        True
        """
        for compound in hierarchy.compounds:
            ancestors = hierarchy.get_ancestors(compound)
            descendants = hierarchy.get_descendants(compound)

            # Compound must have at least one relationship
            # (ancestor, descendant, or be maximal/minimal with edges)
            has_ancestors = len(hierarchy.get_direct_ancestors(compound)) > 0
            has_descendants = len(hierarchy.get_direct_descendants(compound)) > 0

            # Isolated if no edges at all
            if not has_ancestors and not has_descendants:
                # Check if it's a singleton (only compound)
                if len(hierarchy.compounds) > 1:
                    return False

        return True

    def _validate_edges(self, hierarchy: CompoundHierarchy) -> bool:
        """
        Validate that all edges connect compounds in hierarchy.

        Parameters
        ----------
        hierarchy : CompoundHierarchy
            The hierarchy to validate

        Returns
        -------
        bool
            True if all edges are valid, False otherwise

        Notes
        -----
        Checks that:
        - All edge endpoints are in hierarchy.compounds
        - No self-loops
        """
        for ancestor, descendants in hierarchy.edges.items():
            # Ancestor must be in hierarchy
            if ancestor not in hierarchy.compounds:
                return False

            # Check all descendants
            for descendant in descendants:
                # Descendant must be in hierarchy
                if descendant not in hierarchy.compounds:
                    return False

                # No self-loops
                if descendant == ancestor:
                    return False

        return True

    def get_validation_report(self, hierarchy: CompoundHierarchy) -> Dict[str, bool]:
        """
        Get comprehensive validation report.

        Parameters
        ----------
        hierarchy : CompoundHierarchy
            The hierarchy to validate

        Returns
        -------
        Dict[str, bool]
            Validation results for all checks

        Examples
        --------
        >>> report = checker.get_validation_report(hierarchy)
        >>> report
        {
            'is_valid_dag': True,
            'has_no_cycles': True,
            'level_ordering_valid': True,
            'edges_valid': True,
            'connectivity_valid': True
        }
        """
        return {
            'is_valid_dag': self.is_valid_dag(hierarchy),
            'has_no_cycles': not self.has_cycles(hierarchy),
            'level_ordering_valid': self.validate_level_ordering(hierarchy),
            'edges_valid': self._validate_edges(hierarchy),
            'connectivity_valid': self.validate_connectivity(hierarchy)
        }
