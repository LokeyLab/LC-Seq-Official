"""
CompoundHierarchy - DAG/Poset structure for truncation relationships.

Implementation based on THEORY.md Section 3.3, 4.2.2, 4.2.3.
"""

from dataclasses import dataclass, field
from typing import List, Dict, Set, Optional
from enum import Enum
from ..entities.compound import Compound


class HierarchyMode(Enum):
    """
    Hierarchy construction mode.

    Attributes
    ----------
    BUILDING_BLOCK : str
        Building-block mode: forest structure, no convergence
    MONOMER : str
        Monomer mode: DAG with convergence (diamond patterns)
    """

    BUILDING_BLOCK = "building_block"
    MONOMER = "monomer"


@dataclass
class CompoundHierarchy:
    """
    Directed Acyclic Graph (DAG) representing truncation relationships.

    A poset structure where compounds are ordered by truncation relationships.
    Supports two modes:
    - Building-block mode: Forest structure (THEORY.md 4.2.2)
    - Monomer mode: DAG with convergence (THEORY.md 4.2.3)

    Attributes
    ----------
    mode : HierarchyMode
        Building-block or monomer mode
    compounds : List[Compound]
        All compounds in the hierarchy
    edges : Dict[Compound, Set[Compound]]
        Direct descendant relationships: ancestor → {descendants}

    Notes
    -----
    - Graph is directed acyclic (DAG)
    - Edges flow from longer → shorter sequences (maximal → minimal)
    - Building-block mode: No convergence (forest of trees)
    - Monomer mode: Convergence allowed (multiple paths to same compound)
    - Store only direct descendants (transitive reduction)
    - Ancestors computed via reverse traversal

    Examples
    --------
    >>> from lcseq.domain.entities.building_block import BuildingBlock
    >>> from lcseq.domain.entities.chromatogram import Chromatogram
    >>> import numpy as np
    >>>
    >>> # Create simple hierarchy
    >>> hierarchy = CompoundHierarchy(mode=HierarchyMode.BUILDING_BLOCK)
    >>>
    >>> # Create compounds
    >>> chromatogram = Chromatogram(
    ...     time_points=np.array([1.0, 2.0, 3.0]),
    ...     counts=np.array([100.0, 200.0, 150.0])
    ... )
    >>> bb0 = BuildingBlock.from_code(0, "Pro")
    >>> bb1 = BuildingBlock.from_code(1, "Leu")
    >>> bb2 = BuildingBlock.from_code(2, "Val")
    >>> bb_null = BuildingBlock.from_code(0, "Null")
    >>>
    >>> # Maximal compound
    >>> maximal = Compound([bb0, bb1, bb2], chromatogram)
    >>> hierarchy.add_compound(maximal)
    >>>
    >>> # Add truncation
    >>> truncation = Compound([bb_null, bb1, bb2], chromatogram)
    >>> hierarchy.add_compound(truncation)
    >>> hierarchy.add_edge(maximal, truncation)
    >>>
    >>> # Query hierarchy
    >>> hierarchy.get_descendants(maximal)
    [truncation]

    References
    ----------
    THEORY.md Section 1.1: Core Mathematical Model (DAG/Poset)
    THEORY.md Section 3.3: Hierarchy Properties
    THEORY.md Section 4.2.2: Building-Block Mode (Forest)
    THEORY.md Section 4.2.3: Monomer Mode (DAG with Convergence)
    """

    mode: HierarchyMode
    compounds: List[Compound] = field(default_factory=list)
    edges: Dict[Compound, Set[Compound]] = field(default_factory=dict)

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
        """
        if compound not in self.compounds:
            self.compounds.append(compound)
        if compound not in self.edges:
            self.edges[compound] = set()

    def add_edge(self, ancestor: Compound, descendant: Compound) -> None:
        """
        Add a truncation edge (ancestor → descendant).

        Parameters
        ----------
        ancestor : Compound
            Ancestor compound (longer sequence) - THEORY.md Section 3.1:
            "Compound with more building blocks"
        descendant : Compound
            Descendant compound (shorter sequence) - THEORY.md Section 3.1:
            "Compound with fewer building blocks"

        Raises
        ------
        ValueError
            If ancestor or descendant not in hierarchy
        ValueError
            If edge would create cycle
        ValueError
            If descendant has higher level than ancestor

        Notes
        -----
        - Validates acyclic property (no cycles)
        - Validates level ordering (descendant.level < ancestor.level)
        - Only stores direct edges (transitive reduction)

        References
        ----------
        THEORY.md Section 1.1: DAG properties (acyclic, directed)
        THEORY.md Section 3.1: Terminology - Ancestry and Lineage
        THEORY.md Section 3.3: Level ordering
        """
        # Validate both compounds in hierarchy
        if ancestor not in self.compounds:
            raise ValueError(f"Ancestor compound not in hierarchy: {ancestor}")
        if descendant not in self.compounds:
            raise ValueError(f"Descendant compound not in hierarchy: {descendant}")

        # Validate level ordering
        if self.mode == HierarchyMode.BUILDING_BLOCK:
            level_attr = "level"
        else:  # MONOMER
            level_attr = "monomer_level"

        ancestor_level = getattr(ancestor, level_attr)
        descendant_level = getattr(descendant, level_attr)

        if descendant_level >= ancestor_level:
            raise ValueError(
                f"Descendant level ({descendant_level}) must be < ancestor level ({ancestor_level})"
            )

        # Check for cycles (would descendant → ancestor path exist?)
        if self._creates_cycle(ancestor, descendant):
            raise ValueError(
                f"Adding edge {ancestor} → {descendant} would create cycle"
            )

        # Add edge (initialize ancestor's descendants set if needed)
        if ancestor not in self.edges:
            self.edges[ancestor] = set()
        self.edges[ancestor].add(descendant)

    def _creates_cycle(self, ancestor: Compound, descendant: Compound) -> bool:
        """
        Check if adding ancestor → descendant edge would create cycle.

        Returns True if there exists path from descendant back to ancestor.
        """
        # If descendant can reach ancestor, then ancestor → descendant creates cycle
        return ancestor in self._get_all_descendants(descendant)

    def _get_all_descendants(self, compound: Compound) -> List[Compound]:
        """
        Get all descendants (transitive closure) via DFS.

        Parameters
        ----------
        compound : Compound
            Starting compound

        Returns
        -------
        List[Compound]
            All reachable descendants (not including compound itself)
        """
        descendants = []
        stack = [compound]
        visited = []

        while stack:
            current = stack.pop()
            if current in visited:
                continue
            visited.append(current)

            # Find descendants of current
            if current in self.edges:
                for descendant in self.edges[current]:
                    if descendant not in descendants:
                        descendants.append(descendant)
                    stack.append(descendant)

        return descendants

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

        References
        ----------
        THEORY.md Section 3.2: Principal Ideal ↓X
        """
        if compound not in self.compounds:
            return []
        return self._get_all_descendants(compound)

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

        References
        ----------
        THEORY.md Section 3.2: Principal Filter ↑X
        """
        if compound not in self.compounds:
            return []

        # DFS using reverse traversal of edges
        ancestors = []
        stack = [compound]
        visited = []

        while stack:
            current = stack.pop()
            if current in visited:
                continue
            visited.append(current)

            # Find ancestors of current (scan all edges)
            for ancestor, descendants_set in self.edges.items():
                if current in descendants_set:
                    if ancestor not in ancestors:
                        ancestors.append(ancestor)
                    stack.append(ancestor)

        return ancestors

    def get_maximal_compounds(self) -> List[Compound]:
        """
        Get all maximal compounds (no ancestors in hierarchy).

        Returns
        -------
        List[Compound]
            Compounds with no ancestors

        Notes
        -----
        In building-block mode: roots of trees in forest.
        In monomer mode: longest sequences in dataset.

        References
        ----------
        THEORY.md Section 3.2: Maximal Element
        """
        # Build set of all compounds that have ancestors
        has_ancestor = set()
        for ancestor, descendants_set in self.edges.items():
            has_ancestor.update(descendants_set)

        # Maximal = in hierarchy but not in has_ancestor
        return [c for c in self.compounds if c not in has_ancestor]

    def get_minimal_compounds(self) -> List[Compound]:
        """
        Get all minimal compounds (no descendants in hierarchy).

        Returns
        -------
        List[Compound]
            Compounds with no descendants

        Notes
        -----
        Often includes L₀ (all-null compound).

        References
        ----------
        THEORY.md Section 3.2: Minimal Element
        """
        # Minimal = compounds with no descendants
        has_descendants = set(ancestor for ancestor, descendants_set in self.edges.items() if descendants_set)
        return [c for c in self.compounds if c not in has_descendants]

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
        - In building-block mode: uses Compound.level
        - In monomer mode: uses Compound.monomer_level

        Examples
        --------
        >>> # Get all compounds with 2 non-null building blocks
        >>> hierarchy.get_level(2)
        [Leu-Val, Pro-Leu, ...]

        References
        ----------
        THEORY.md Section 3.3: Truncation Level
        """
        if self.mode == HierarchyMode.BUILDING_BLOCK:
            level_attr = "level"
        else:  # MONOMER
            level_attr = "monomer_level"

        return [c for c in self.compounds if getattr(c, level_attr) == level]

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

        Notes
        -----
        This returns stored edges (transitive reduction).

        References
        ----------
        THEORY.md Section 3.1: "Descendant - Compound with fewer building blocks"
        """
        if compound not in self.compounds:
            return []
        return list(self.edges.get(compound, set()))

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
        if compound not in self.compounds:
            return []

        return [ancestor for ancestor, descendants_set in self.edges.items() if compound in descendants_set]

    def size(self) -> int:
        """
        Get number of compounds in hierarchy.

        Returns
        -------
        int
            Number of vertices in DAG
        """
        return len(self.compounds)

    def edge_count(self) -> int:
        """
        Get number of edges in hierarchy.

        Returns
        -------
        int
            Number of directed edges
        """
        return sum(len(descendants_set) for descendants_set in self.edges.values())

    def topological_sort(self) -> List[Compound]:
        """
        Return compounds in topological order (L₀ first, bottom-up).

        Uses Kahn's algorithm (THEORY.md Section 7.1) to process compounds
        in guaranteed correct order for bottom-up traversal. Minimal elements
        (compounds with no descendants, including L₀) are processed first,
        followed by their ancestors in dependency order.

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
        - Automatically handles gaps in levels
        - O(V + E) complexity
        - Stable tie-breaking: within same level, sorted by canonical sequence

        Used for:
        - Peak classification constraint propagation (THEORY.md 2.4.6, 6.8)
        - Bottom-up analysis workflows
        - Ensuring correct dependency order

        Examples
        --------
        >>> sorted_compounds = hierarchy.topological_sort()
        >>> # Process in order: L₀ first, then ancestors
        >>> for compound in sorted_compounds:
        ...     process(compound)  # All descendants already processed

        References
        ----------
        THEORY.md Section 7.1: Topological Sort (Kahn's Algorithm)
        THEORY.md Section 2.4.6: L₀ as Anchor Point
        THEORY.md Section 6.8: DAG Constraint Propagation
        """
        # Determine level attribute based on hierarchy mode
        level_attr = "monomer_level" if self.mode == HierarchyMode.MONOMER else "level"

        # Count outgoing edges (out-degree) for each compound
        # Out-degree = number of descendants (compounds this one points to)
        # Minimal elements (L₀) have out-degree 0 (no descendants)
        out_degree = {compound: 0 for compound in self.compounds}

        for ancestor, descendants_set in self.edges.items():
            out_degree[ancestor] = len(descendants_set)

        # Start with minimal elements (out-degree 0 = no descendants, includes L₀)
        # Sort by (level, canonical_sequence) for stable, semantically meaningful order
        minimal = [c for c in self.compounds if out_degree[c] == 0]
        queue = sorted(minimal, key=lambda c: (getattr(c, level_attr), c.residue_sequence))
        result = []

        # Process compounds in dependency order (bottom-up)
        while queue:
            compound = queue.pop(0)
            result.append(compound)

            # Find ancestors of this compound (compounds that point to it)
            # Reduce their out-degree (they've had one descendant processed)
            # When ancestor's out-degree reaches 0, all its descendants processed
            for ancestor, descendants_set in self.edges.items():
                if compound in descendants_set:
                    out_degree[ancestor] -= 1
                    if out_degree[ancestor] == 0:
                        queue.append(ancestor)
                        # Maintain sorted order for stable tie-breaking
                        queue.sort(key=lambda c: (getattr(c, level_attr), c.residue_sequence))

        # Verify all compounds processed (detect cycles)
        if len(result) != len(self.compounds):
            raise RuntimeError(
                f"Cycle detected in hierarchy! "
                f"Processed {len(result)} of {len(self.compounds)} compounds"
            )

        return result

    def __repr__(self) -> str:
        """Detailed representation for debugging."""
        return (
            f"CompoundHierarchy(mode={self.mode.value}, "
            f"compounds={self.size()}, "
            f"edges={self.edge_count()})"
        )
