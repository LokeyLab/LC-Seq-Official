"""
CompoundHierarchy - DAG/Poset structure for truncation relationships.

Implementation based on THEORY.md Section 3.3, 4.2.2, 4.2.3.

Uses rustworkx for high-performance graph operations:
- O(1) node/edge lookup
- O(V+E) topological sort
- Efficient ancestor/descendant queries
"""

from dataclasses import dataclass, field
from typing import List, Dict, Set, Optional
from enum import Enum
import rustworkx as rx
from ..entities.compound import Compound


class HierarchyMode(Enum):
    """
    Hierarchy construction mode.

    Both modes use equivalence class logic based on atomic units:

    Attributes
    ----------
    BUILDING_BLOCK : str
        Building-block mode: DAG with convergence at block granularity
        (positional variants with same blocks → same equivalence class)
    MONOMER : str
        Monomer mode: DAG with convergence at monomer granularity
        (positional variants with same monomers → same equivalence class)
    """

    BUILDING_BLOCK = "building_block"
    MONOMER = "monomer"


class CompoundHierarchy:
    """
    Directed Acyclic Graph (DAG) representing truncation relationships.

    A poset structure where compounds are ordered by truncation relationships.
    Both modes use equivalence class logic (atomic unit principle):
    - Building-block mode: Convergence at block granularity (THEORY.md 4.2.2)
    - Monomer mode: Convergence at monomer granularity (THEORY.md 4.2.3)

    Uses rustworkx for high-performance graph operations:
    - O(1) node/edge lookup via hash maps
    - O(V+E) topological sort
    - Efficient ancestor/descendant queries

    Attributes
    ----------
    mode : HierarchyMode
        Building-block or monomer mode

    Notes
    -----
    - Graph is directed acyclic (DAG)
    - Edges flow from longer → shorter sequences (maximal → minimal)
    - Building-block mode: Positional variants with same blocks converge
    - Monomer mode: Positional variants with same monomers converge
    - Store only direct descendants (transitive reduction)
    - Ancestors computed via reverse traversal

    References
    ----------
    THEORY.md Section 1.1: Core Mathematical Model (DAG/Poset)
    THEORY.md Section 3.3: Hierarchy Properties
    THEORY.md Section 4.2.2: Building-Block Mode (Convergence at Block Granularity)
    THEORY.md Section 4.2.3: Monomer Mode (Convergence at Monomer Granularity)
    """

    def __init__(self, mode: HierarchyMode = HierarchyMode.BUILDING_BLOCK):
        """
        Initialize hierarchy with specified mode.

        Parameters
        ----------
        mode : HierarchyMode
            Building-block or monomer mode (default: BUILDING_BLOCK)
        """
        self.mode = mode
        self._graph = rx.PyDiGraph()
        self._compound_to_idx: Dict[Compound, int] = {}
        self._idx_to_compound: Dict[int, Compound] = {}
        self._level_index: Dict[int, List[Compound]] = {}  # Cache for O(1) level lookups

    @property
    def compounds(self) -> List[Compound]:
        """Get all compounds in hierarchy."""
        return list(self._compound_to_idx.keys())

    def add_compound(self, compound: Compound) -> None:
        """
        Add a compound to the hierarchy.

        Parameters
        ----------
        compound : Compound
            Compound to add

        Notes
        -----
        - Idempotent: adding same compound multiple times is safe
        - O(1) operation
        - Updates level index cache for O(1) level lookups
        """
        if compound not in self._compound_to_idx:
            idx = self._graph.add_node(compound)
            self._compound_to_idx[compound] = idx
            self._idx_to_compound[idx] = compound

            # Update level index cache
            level = compound.monomer_level if self.mode == HierarchyMode.MONOMER else compound.level
            if level not in self._level_index:
                self._level_index[level] = []
            self._level_index[level].append(compound)

    def add_edge(self, ancestor: Compound, descendant: Compound) -> None:
        """
        Add a truncation edge (ancestor → descendant).

        Parameters
        ----------
        ancestor : Compound
            Ancestor compound (longer sequence)
        descendant : Compound
            Descendant compound (shorter sequence)

        Raises
        ------
        ValueError
            If ancestor or descendant not in hierarchy
        ValueError
            If descendant has higher level than ancestor

        Notes
        -----
        - Validates level ordering (descendant.level < ancestor.level)
        - O(1) operation (cycle check skipped for performance - validated by level ordering)

        References
        ----------
        THEORY.md Section 1.1: DAG properties (acyclic, directed)
        THEORY.md Section 3.1: Terminology - Ancestry and Lineage
        THEORY.md Section 3.3: Level ordering
        """
        # Validate both compounds in hierarchy
        if ancestor not in self._compound_to_idx:
            raise ValueError(f"Ancestor compound not in hierarchy: {ancestor}")
        if descendant not in self._compound_to_idx:
            raise ValueError(f"Descendant compound not in hierarchy: {descendant}")

        # Validate level ordering (this guarantees no cycles)
        level_attr = "monomer_level" if self.mode == HierarchyMode.MONOMER else "level"
        ancestor_level = getattr(ancestor, level_attr)
        descendant_level = getattr(descendant, level_attr)

        if descendant_level >= ancestor_level:
            raise ValueError(
                f"Descendant level ({descendant_level}) must be < ancestor level ({ancestor_level})"
            )

        # Add edge
        ancestor_idx = self._compound_to_idx[ancestor]
        descendant_idx = self._compound_to_idx[descendant]
        self._graph.add_edge(ancestor_idx, descendant_idx, None)

    def get_descendants(self, compound: Compound) -> List[Compound]:
        """
        Get all descendants of a compound.

        Parameters
        ----------
        compound : Compound
            Query compound

        Returns
        -------
        List[Compound]
            All descendants (transitive closure)

        Notes
        -----
        Computes principal ideal ↓X in poset terminology (THEORY.md 3.2).
        Uses rustworkx's efficient BFS/DFS traversal.

        References
        ----------
        THEORY.md Section 3.2: Principal Ideal ↓X
        """
        if compound not in self._compound_to_idx:
            return []

        idx = self._compound_to_idx[compound]
        # rx.descendants returns set of node indices reachable from idx
        descendant_indices = rx.descendants(self._graph, idx)
        return [self._idx_to_compound[i] for i in descendant_indices]

    def get_ancestors(self, compound: Compound) -> List[Compound]:
        """
        Get all ancestors of a compound.

        Parameters
        ----------
        compound : Compound
            Query compound

        Returns
        -------
        List[Compound]
            All ancestors (transitive closure)

        Notes
        -----
        Computes principal filter ↑X in poset terminology (THEORY.md 3.2).
        Uses rustworkx's efficient reverse traversal.

        References
        ----------
        THEORY.md Section 3.2: Principal Filter ↑X
        """
        if compound not in self._compound_to_idx:
            return []

        idx = self._compound_to_idx[compound]
        # rx.ancestors returns set of node indices that can reach idx
        ancestor_indices = rx.ancestors(self._graph, idx)
        return [self._idx_to_compound[i] for i in ancestor_indices]

    def get_maximal_compounds(self) -> List[Compound]:
        """
        Get all maximal compounds (no ancestors in hierarchy).

        Returns
        -------
        List[Compound]
            Compounds with no ancestors (in-degree 0)

        Notes
        -----
        Typically the reference compound(s) at the root of the lineage.

        References
        ----------
        THEORY.md Section 3.2: Maximal Element
        """
        # Maximal = nodes with in-degree 0 (no incoming edges = no ancestors)
        return [
            self._idx_to_compound[idx]
            for idx in self._graph.node_indices()
            if self._graph.in_degree(idx) == 0
        ]

    def get_minimal_compounds(self) -> List[Compound]:
        """
        Get all minimal compounds (no descendants in hierarchy).

        Returns
        -------
        List[Compound]
            Compounds with no descendants (out-degree 0)

        Notes
        -----
        Often includes L₀ (all-null compound).

        References
        ----------
        THEORY.md Section 3.2: Minimal Element
        """
        # Minimal = nodes with out-degree 0 (no outgoing edges = no descendants)
        return [
            self._idx_to_compound[idx]
            for idx in self._graph.node_indices()
            if self._graph.out_degree(idx) == 0
        ]

    def get_level(self, level: int) -> List[Compound]:
        """
        Get all compounds at a specific truncation level.

        Parameters
        ----------
        level : int
            Truncation level (number of non-null building blocks or monomers)

        Returns
        -------
        List[Compound]
            All compounds at the specified level

        Notes
        -----
        O(1) lookup using cached level index, populated during add_compound().

        References
        ----------
        THEORY.md Section 3.3: Truncation Level
        """
        return self._level_index.get(level, [])

    def get_direct_descendants(self, compound: Compound) -> List[Compound]:
        """
        Get direct descendants (immediate truncations) only.

        Parameters
        ----------
        compound : Compound
            Query compound

        Returns
        -------
        List[Compound]
            Immediate descendants (one edge away)

        References
        ----------
        THEORY.md Section 3.1: "Descendant - Compound with fewer building blocks"
        """
        if compound not in self._compound_to_idx:
            return []

        idx = self._compound_to_idx[compound]
        # Get successor indices (direct outgoing edges)
        successor_indices = self._graph.successor_indices(idx)
        return [self._idx_to_compound[i] for i in successor_indices]

    def get_direct_ancestors(self, compound: Compound) -> List[Compound]:
        """
        Get direct ancestors (immediate extensions) only.

        Parameters
        ----------
        compound : Compound
            Query compound

        Returns
        -------
        List[Compound]
            Immediate ancestors (one edge away)

        References
        ----------
        THEORY.md Section 3.1: "Ancestor - Compound with more building blocks"
        """
        if compound not in self._compound_to_idx:
            return []

        idx = self._compound_to_idx[compound]
        # Get predecessor indices (direct incoming edges)
        predecessor_indices = self._graph.predecessor_indices(idx)
        return [self._idx_to_compound[i] for i in predecessor_indices]

    def size(self) -> int:
        """
        Get number of compounds in hierarchy.

        Returns
        -------
        int
            Number of vertices in DAG
        """
        return self._graph.num_nodes()

    def edge_count(self) -> int:
        """
        Get number of edges in hierarchy.

        Returns
        -------
        int
            Number of directed edges
        """
        return self._graph.num_edges()

    def topological_sort(self) -> List[Compound]:
        """
        Return compounds in topological order (L₀ first, bottom-up).

        Returns
        -------
        List[Compound]
            Compounds in topological order (bottom-up: minimal → maximal)

        Raises
        ------
        RuntimeError
            If cycle detected (should never happen with valid DAG)

        Notes
        -----
        Processing order guarantees:
        - L₀ (all-null) processed first
        - Every compound processed after all its descendants
        - O(V + E) complexity via rustworkx

        References
        ----------
        THEORY.md Section 7.1: Topological Sort (Kahn's Algorithm)
        """
        try:
            # rx.topological_sort returns indices in topological order
            # By default this is maximal → minimal, so we reverse for bottom-up
            sorted_indices = rx.topological_sort(self._graph)
            # Reverse to get minimal → maximal (bottom-up) order
            return [self._idx_to_compound[idx] for idx in reversed(sorted_indices)]
        except rx.DAGHasCycle:
            raise RuntimeError("Cycle detected in hierarchy!")

    def __contains__(self, compound: Compound) -> bool:
        """Check if compound is in hierarchy."""
        return compound in self._compound_to_idx

    def __repr__(self) -> str:
        """Detailed representation for debugging."""
        return (
            f"CompoundHierarchy(mode={self.mode.value}, "
            f"compounds={self.size()}, "
            f"edges={self.edge_count()})"
        )
