"""
Comprehensive tests for CompoundHierarchy model.

Tests DAG/Poset structure, both building-block and monomer modes.
"""

import pytest
import numpy as np
from lcseq.domain.models.compound_hierarchy import CompoundHierarchy, HierarchyMode
from lcseq.domain.entities.compound import Compound
from lcseq.domain.entities.building_block import BuildingBlock
from lcseq.domain.entities.chromatogram import Chromatogram


@pytest.fixture
def sample_chromatogram():
    """Create sample chromatogram for testing."""
    return Chromatogram(
        time_points=np.array([1.0, 2.0, 3.0, 4.0, 5.0]),
        counts=np.array([100.0, 200.0, 300.0, 200.0, 100.0]),
    )


@pytest.fixture
def building_blocks():
    """Create standard building blocks for testing."""
    return {
        "Pro": BuildingBlock.from_code(0, "Pro"),
        "Leu": BuildingBlock.from_code(1, "Leu"),
        "Val": BuildingBlock.from_code(2, "Val"),
        "Null0": BuildingBlock.from_code(0, "Null"),
        "Null1": BuildingBlock.from_code(1, "Null"),
        "Null2": BuildingBlock.from_code(2, "Null"),
    }


class TestHierarchyCreation:
    """Test hierarchy creation and initialization."""

    def test_create_empty_hierarchy_building_block_mode(self):
        """Test creating empty hierarchy in building-block mode."""
        hierarchy = CompoundHierarchy(mode=HierarchyMode.BUILDING_BLOCK)

        assert hierarchy.mode == HierarchyMode.BUILDING_BLOCK
        assert len(hierarchy.compounds) == 0
        assert len(hierarchy.edges) == 0
        assert hierarchy.size() == 0
        assert hierarchy.edge_count() == 0

    def test_create_empty_hierarchy_monomer_mode(self):
        """Test creating empty hierarchy in monomer mode."""
        hierarchy = CompoundHierarchy(mode=HierarchyMode.MONOMER)

        assert hierarchy.mode == HierarchyMode.MONOMER
        assert len(hierarchy.compounds) == 0
        assert len(hierarchy.edges) == 0

    def test_hierarchy_repr(self):
        """Test string representation."""
        hierarchy = CompoundHierarchy(mode=HierarchyMode.BUILDING_BLOCK)

        repr_str = repr(hierarchy)
        assert "CompoundHierarchy" in repr_str
        assert "building_block" in repr_str
        assert "compounds=0" in repr_str
        assert "edges=0" in repr_str


class TestAddCompound:
    """Test adding compounds to hierarchy."""

    def test_add_single_compound(self, sample_chromatogram, building_blocks):
        """Test adding a single compound."""
        hierarchy = CompoundHierarchy(mode=HierarchyMode.BUILDING_BLOCK)

        compound = Compound(
            [building_blocks["Pro"], building_blocks["Leu"]], sample_chromatogram
        )
        hierarchy.add_compound(compound)

        assert compound in hierarchy.compounds
        assert hierarchy.size() == 1
        assert compound in hierarchy.edges
        assert len(hierarchy.edges[compound]) == 0

    def test_add_multiple_compounds(self, sample_chromatogram, building_blocks):
        """Test adding multiple compounds."""
        hierarchy = CompoundHierarchy(mode=HierarchyMode.BUILDING_BLOCK)

        c1 = Compound([building_blocks["Pro"], building_blocks["Leu"]], sample_chromatogram)
        c2 = Compound([building_blocks["Null0"], building_blocks["Leu"]], sample_chromatogram)
        c3 = Compound([building_blocks["Pro"], building_blocks["Null1"]], sample_chromatogram)

        hierarchy.add_compound(c1)
        hierarchy.add_compound(c2)
        hierarchy.add_compound(c3)

        assert hierarchy.size() == 3
        assert c1 in hierarchy.compounds
        assert c2 in hierarchy.compounds
        assert c3 in hierarchy.compounds

    def test_add_compound_idempotent(self, sample_chromatogram, building_blocks):
        """Test that adding same compound multiple times is safe."""
        hierarchy = CompoundHierarchy(mode=HierarchyMode.BUILDING_BLOCK)

        compound = Compound([building_blocks["Pro"]], sample_chromatogram)
        hierarchy.add_compound(compound)
        hierarchy.add_compound(compound)
        hierarchy.add_compound(compound)

        assert hierarchy.size() == 1


class TestAddEdge:
    """Test adding truncation edges."""

    def test_add_valid_edge(self, sample_chromatogram, building_blocks):
        """Test adding valid truncation edge."""
        hierarchy = CompoundHierarchy(mode=HierarchyMode.BUILDING_BLOCK)

        ancestor = Compound([building_blocks["Pro"], building_blocks["Leu"]], sample_chromatogram)
        descendant = Compound([building_blocks["Null0"], building_blocks["Leu"]], sample_chromatogram)

        hierarchy.add_compound(ancestor)
        hierarchy.add_compound(descendant)
        hierarchy.add_edge(ancestor, descendant)

        assert descendant in hierarchy.edges[ancestor]
        assert hierarchy.edge_count() == 1

    def test_add_edge_parent_not_in_hierarchy(self, sample_chromatogram, building_blocks):
        """Test adding edge with ancestor not in hierarchy raises error."""
        hierarchy = CompoundHierarchy(mode=HierarchyMode.BUILDING_BLOCK)

        ancestor = Compound([building_blocks["Pro"], building_blocks["Leu"]], sample_chromatogram)
        descendant = Compound([building_blocks["Null0"], building_blocks["Leu"]], sample_chromatogram)
        hierarchy.add_compound(descendant)

        with pytest.raises(ValueError, match="Ancestor compound not in hierarchy"):
            hierarchy.add_edge(ancestor, descendant)

    def test_add_edge_child_not_in_hierarchy(self, sample_chromatogram, building_blocks):
        """Test adding edge with descendant not in hierarchy raises error."""
        hierarchy = CompoundHierarchy(mode=HierarchyMode.BUILDING_BLOCK)

        ancestor = Compound([building_blocks["Pro"], building_blocks["Leu"]], sample_chromatogram)
        descendant = Compound([building_blocks["Null0"], building_blocks["Leu"]], sample_chromatogram)
        hierarchy.add_compound(ancestor)

        with pytest.raises(ValueError, match="Descendant compound not in hierarchy"):
            hierarchy.add_edge(ancestor, descendant)

    def test_add_edge_invalid_level_ordering(self, sample_chromatogram, building_blocks):
        """Test adding edge with invalid level ordering raises error."""
        hierarchy = CompoundHierarchy(mode=HierarchyMode.BUILDING_BLOCK)

        # Descendant has same level as ancestor (invalid)
        ancestor = Compound([building_blocks["Pro"], building_blocks["Leu"]], sample_chromatogram)
        descendant = Compound([building_blocks["Pro"], building_blocks["Val"]], sample_chromatogram)

        hierarchy.add_compound(ancestor)
        hierarchy.add_compound(descendant)

        with pytest.raises(ValueError, match="Descendant level .* must be < ancestor level"):
            hierarchy.add_edge(ancestor, descendant)

    def test_add_edge_creates_cycle(self, sample_chromatogram, building_blocks):
        """Test adding edge that would create cycle raises error."""
        hierarchy = CompoundHierarchy(mode=HierarchyMode.BUILDING_BLOCK)

        c1 = Compound([building_blocks["Pro"], building_blocks["Leu"]], sample_chromatogram)
        c2 = Compound([building_blocks["Null0"], building_blocks["Leu"]], sample_chromatogram)
        c3 = Compound([building_blocks["Null0"], building_blocks["Null1"]], sample_chromatogram)

        hierarchy.add_compound(c1)
        hierarchy.add_compound(c2)
        hierarchy.add_compound(c3)

        # Add valid edges
        hierarchy.add_edge(c1, c2)
        hierarchy.add_edge(c2, c3)

        # Try to add edge that creates cycle: c3 -> c1
        # This will fail level check first (c3.level=0 < c1.level=2)
        with pytest.raises(ValueError, match="Descendant level .* must be < ancestor level"):
            hierarchy.add_edge(c3, c1)


class TestGetDescendants:
    """Test getting descendants (transitive closure)."""

    def test_get_descendants_no_children(self, sample_chromatogram, building_blocks):
        """Test getting descendants of leaf node."""
        hierarchy = CompoundHierarchy(mode=HierarchyMode.BUILDING_BLOCK)

        compound = Compound([building_blocks["Pro"]], sample_chromatogram)
        hierarchy.add_compound(compound)

        descendants = hierarchy.get_descendants(compound)
        assert len(descendants) == 0

    def test_get_descendants_direct_children(self, sample_chromatogram, building_blocks):
        """Test getting direct descendants."""
        hierarchy = CompoundHierarchy(mode=HierarchyMode.BUILDING_BLOCK)

        ancestor = Compound([building_blocks["Pro"], building_blocks["Leu"]], sample_chromatogram)
        descendant1 = Compound([building_blocks["Null0"], building_blocks["Leu"]], sample_chromatogram)
        descendant2 = Compound([building_blocks["Pro"], building_blocks["Null1"]], sample_chromatogram)

        hierarchy.add_compound(ancestor)
        hierarchy.add_compound(descendant1)
        hierarchy.add_compound(descendant2)
        hierarchy.add_edge(ancestor, descendant1)
        hierarchy.add_edge(ancestor, descendant2)

        descendants = hierarchy.get_descendants(ancestor)
        assert len(descendants) == 2
        assert descendant1 in descendants
        assert descendant2 in descendants

    def test_get_descendants_transitive(self, sample_chromatogram, building_blocks):
        """Test getting transitive descendants."""
        hierarchy = CompoundHierarchy(mode=HierarchyMode.BUILDING_BLOCK)

        c1 = Compound([building_blocks["Pro"], building_blocks["Leu"], building_blocks["Val"]], sample_chromatogram)
        c2 = Compound([building_blocks["Null0"], building_blocks["Leu"], building_blocks["Val"]], sample_chromatogram)
        c3 = Compound([building_blocks["Null0"], building_blocks["Null1"], building_blocks["Val"]], sample_chromatogram)
        c4 = Compound([building_blocks["Null0"], building_blocks["Null1"], building_blocks["Null2"]], sample_chromatogram)

        hierarchy.add_compound(c1)
        hierarchy.add_compound(c2)
        hierarchy.add_compound(c3)
        hierarchy.add_compound(c4)
        hierarchy.add_edge(c1, c2)
        hierarchy.add_edge(c2, c3)
        hierarchy.add_edge(c3, c4)

        # c1 should have all three as descendants
        descendants = hierarchy.get_descendants(c1)
        assert len(descendants) == 3
        assert c2 in descendants
        assert c3 in descendants
        assert c4 in descendants

    def test_get_descendants_compound_not_in_hierarchy(self, sample_chromatogram, building_blocks):
        """Test getting descendants of compound not in hierarchy."""
        hierarchy = CompoundHierarchy(mode=HierarchyMode.BUILDING_BLOCK)

        compound = Compound([building_blocks["Pro"]], sample_chromatogram)
        descendants = hierarchy.get_descendants(compound)

        assert len(descendants) == 0


class TestGetAncestors:
    """Test getting ancestors (reverse transitive closure)."""

    def test_get_ancestors_no_parents(self, sample_chromatogram, building_blocks):
        """Test getting ancestors of root node."""
        hierarchy = CompoundHierarchy(mode=HierarchyMode.BUILDING_BLOCK)

        compound = Compound([building_blocks["Pro"], building_blocks["Leu"]], sample_chromatogram)
        hierarchy.add_compound(compound)

        ancestors = hierarchy.get_ancestors(compound)
        assert len(ancestors) == 0

    def test_get_ancestors_direct_parents(self, sample_chromatogram, building_blocks):
        """Test getting direct ancestors."""
        hierarchy = CompoundHierarchy(mode=HierarchyMode.BUILDING_BLOCK)

        ancestor1 = Compound([building_blocks["Pro"], building_blocks["Leu"]], sample_chromatogram)
        ancestor2 = Compound([building_blocks["Pro"], building_blocks["Val"]], sample_chromatogram)
        descendant = Compound([building_blocks["Pro"], building_blocks["Null1"]], sample_chromatogram)

        hierarchy.add_compound(ancestor1)
        hierarchy.add_compound(ancestor2)
        hierarchy.add_compound(descendant)
        hierarchy.add_edge(ancestor1, descendant)
        hierarchy.add_edge(ancestor2, descendant)

        ancestors = hierarchy.get_ancestors(descendant)
        assert len(ancestors) == 2
        assert ancestor1 in ancestors
        assert ancestor2 in ancestors

    def test_get_ancestors_transitive(self, sample_chromatogram, building_blocks):
        """Test getting transitive ancestors."""
        hierarchy = CompoundHierarchy(mode=HierarchyMode.BUILDING_BLOCK)

        c1 = Compound([building_blocks["Pro"], building_blocks["Leu"], building_blocks["Val"]], sample_chromatogram)
        c2 = Compound([building_blocks["Null0"], building_blocks["Leu"], building_blocks["Val"]], sample_chromatogram)
        c3 = Compound([building_blocks["Null0"], building_blocks["Null1"], building_blocks["Val"]], sample_chromatogram)
        c4 = Compound([building_blocks["Null0"], building_blocks["Null1"], building_blocks["Null2"]], sample_chromatogram)

        hierarchy.add_compound(c1)
        hierarchy.add_compound(c2)
        hierarchy.add_compound(c3)
        hierarchy.add_compound(c4)
        hierarchy.add_edge(c1, c2)
        hierarchy.add_edge(c2, c3)
        hierarchy.add_edge(c3, c4)

        # c4 should have all three as ancestors
        ancestors = hierarchy.get_ancestors(c4)
        assert len(ancestors) == 3
        assert c1 in ancestors
        assert c2 in ancestors
        assert c3 in ancestors


class TestMaximalMinimalCompounds:
    """Test getting maximal and minimal compounds."""

    def test_get_maximal_compounds_single_root(self, sample_chromatogram, building_blocks):
        """Test getting maximal compounds with single root."""
        hierarchy = CompoundHierarchy(mode=HierarchyMode.BUILDING_BLOCK)

        root = Compound([building_blocks["Pro"], building_blocks["Leu"]], sample_chromatogram)
        descendant = Compound([building_blocks["Null0"], building_blocks["Leu"]], sample_chromatogram)

        hierarchy.add_compound(root)
        hierarchy.add_compound(descendant)
        hierarchy.add_edge(root, descendant)

        maximal = hierarchy.get_maximal_compounds()
        assert len(maximal) == 1
        assert root in maximal

    def test_get_maximal_compounds_multiple_roots(self, sample_chromatogram, building_blocks):
        """Test getting maximal compounds with forest structure."""
        hierarchy = CompoundHierarchy(mode=HierarchyMode.BUILDING_BLOCK)

        root1 = Compound([building_blocks["Pro"], building_blocks["Leu"]], sample_chromatogram)
        root2 = Compound([building_blocks["Pro"], building_blocks["Val"]], sample_chromatogram)
        descendant = Compound([building_blocks["Pro"], building_blocks["Null1"]], sample_chromatogram)

        hierarchy.add_compound(root1)
        hierarchy.add_compound(root2)
        hierarchy.add_compound(descendant)
        hierarchy.add_edge(root1, descendant)
        hierarchy.add_edge(root2, descendant)

        maximal = hierarchy.get_maximal_compounds()
        assert len(maximal) == 2
        assert root1 in maximal
        assert root2 in maximal

    def test_get_minimal_compounds_single_leaf(self, sample_chromatogram, building_blocks):
        """Test getting minimal compounds with single leaf."""
        hierarchy = CompoundHierarchy(mode=HierarchyMode.BUILDING_BLOCK)

        root = Compound([building_blocks["Pro"], building_blocks["Leu"]], sample_chromatogram)
        leaf = Compound([building_blocks["Null0"], building_blocks["Null1"]], sample_chromatogram)

        hierarchy.add_compound(root)
        hierarchy.add_compound(leaf)
        hierarchy.add_edge(root, leaf)

        minimal = hierarchy.get_minimal_compounds()
        assert len(minimal) == 1
        assert leaf in minimal

    def test_get_minimal_compounds_multiple_leaves(self, sample_chromatogram, building_blocks):
        """Test getting minimal compounds with multiple leaves."""
        hierarchy = CompoundHierarchy(mode=HierarchyMode.BUILDING_BLOCK)

        root = Compound([building_blocks["Pro"], building_blocks["Leu"]], sample_chromatogram)
        leaf1 = Compound([building_blocks["Null0"], building_blocks["Leu"]], sample_chromatogram)
        leaf2 = Compound([building_blocks["Pro"], building_blocks["Null1"]], sample_chromatogram)

        hierarchy.add_compound(root)
        hierarchy.add_compound(leaf1)
        hierarchy.add_compound(leaf2)
        hierarchy.add_edge(root, leaf1)
        hierarchy.add_edge(root, leaf2)

        minimal = hierarchy.get_minimal_compounds()
        assert len(minimal) == 2
        assert leaf1 in minimal
        assert leaf2 in minimal


class TestGetLevel:
    """Test getting compounds at specific level."""

    def test_get_level_building_block_mode(self, sample_chromatogram, building_blocks):
        """Test getting level in building-block mode."""
        hierarchy = CompoundHierarchy(mode=HierarchyMode.BUILDING_BLOCK)

        c0 = Compound([building_blocks["Null0"], building_blocks["Null1"], building_blocks["Null2"]], sample_chromatogram)
        c1 = Compound([building_blocks["Pro"], building_blocks["Null1"], building_blocks["Null2"]], sample_chromatogram)
        c2a = Compound([building_blocks["Pro"], building_blocks["Leu"], building_blocks["Null2"]], sample_chromatogram)
        c2b = Compound([building_blocks["Pro"], building_blocks["Null1"], building_blocks["Val"]], sample_chromatogram)
        c3 = Compound([building_blocks["Pro"], building_blocks["Leu"], building_blocks["Val"]], sample_chromatogram)

        for c in [c0, c1, c2a, c2b, c3]:
            hierarchy.add_compound(c)

        level0 = hierarchy.get_level(0)
        level1 = hierarchy.get_level(1)
        level2 = hierarchy.get_level(2)
        level3 = hierarchy.get_level(3)

        assert len(level0) == 1 and c0 in level0
        assert len(level1) == 1 and c1 in level1
        assert len(level2) == 2 and c2a in level2 and c2b in level2
        assert len(level3) == 1 and c3 in level3

    def test_get_level_empty(self, sample_chromatogram, building_blocks):
        """Test getting level with no compounds at that level."""
        hierarchy = CompoundHierarchy(mode=HierarchyMode.BUILDING_BLOCK)

        compound = Compound([building_blocks["Pro"], building_blocks["Leu"]], sample_chromatogram)
        hierarchy.add_compound(compound)

        level0 = hierarchy.get_level(0)
        level3 = hierarchy.get_level(3)

        assert len(level0) == 0
        assert len(level3) == 0


class TestDirectRelationships:
    """Test getting direct ancestors and descendants."""

    def test_get_direct_descendants(self, sample_chromatogram, building_blocks):
        """Test getting direct descendants only."""
        hierarchy = CompoundHierarchy(mode=HierarchyMode.BUILDING_BLOCK)

        c1 = Compound([building_blocks["Pro"], building_blocks["Leu"], building_blocks["Val"]], sample_chromatogram)
        c2 = Compound([building_blocks["Null0"], building_blocks["Leu"], building_blocks["Val"]], sample_chromatogram)
        c3 = Compound([building_blocks["Null0"], building_blocks["Null1"], building_blocks["Val"]], sample_chromatogram)

        hierarchy.add_compound(c1)
        hierarchy.add_compound(c2)
        hierarchy.add_compound(c3)
        hierarchy.add_edge(c1, c2)
        hierarchy.add_edge(c2, c3)

        # Direct descendants of c1
        direct = hierarchy.get_direct_descendants(c1)
        assert len(direct) == 1
        assert c2 in direct
        assert c3 not in direct  # c3 is transitive, not direct

    def test_get_direct_ancestors(self, sample_chromatogram, building_blocks):
        """Test getting direct ancestors only."""
        hierarchy = CompoundHierarchy(mode=HierarchyMode.BUILDING_BLOCK)

        c1 = Compound([building_blocks["Pro"], building_blocks["Leu"], building_blocks["Val"]], sample_chromatogram)
        c2 = Compound([building_blocks["Null0"], building_blocks["Leu"], building_blocks["Val"]], sample_chromatogram)
        c3 = Compound([building_blocks["Null0"], building_blocks["Null1"], building_blocks["Val"]], sample_chromatogram)

        hierarchy.add_compound(c1)
        hierarchy.add_compound(c2)
        hierarchy.add_compound(c3)
        hierarchy.add_edge(c1, c2)
        hierarchy.add_edge(c2, c3)

        # Direct ancestors of c3
        direct = hierarchy.get_direct_ancestors(c3)
        assert len(direct) == 1
        assert c2 in direct
        assert c1 not in direct  # c1 is transitive, not direct


class TestMonomerMode:
    """Test hierarchy in monomer mode."""

    def test_monomer_mode_level_calculation(self, sample_chromatogram):
        """Test level calculation uses monomer_level in monomer mode."""
        hierarchy = CompoundHierarchy(mode=HierarchyMode.MONOMER)

        # Create compound with composite building block
        bb_composite = BuildingBlock.from_code(1, "Leu-Ala-Val")  # 3 monomers
        bb_single = BuildingBlock.from_code(0, "Pro")  # 1 monomer

        compound = Compound([bb_single, bb_composite], sample_chromatogram)
        hierarchy.add_compound(compound)

        # Should use monomer_level (4) not level (2)
        level4 = hierarchy.get_level(4)
        assert len(level4) == 1
        assert compound in level4

    def test_monomer_mode_edge_validation(self, sample_chromatogram):
        """Test edge validation uses monomer_level in monomer mode."""
        hierarchy = CompoundHierarchy(mode=HierarchyMode.MONOMER)

        # Ancestor: 5 monomers
        bb1 = BuildingBlock.from_code(0, "Pro")
        bb2 = BuildingBlock.from_code(1, "Leu-Ala-Val")  # 3 monomers
        bb3 = BuildingBlock.from_code(2, "Phe")
        ancestor = Compound([bb1, bb2, bb3], sample_chromatogram)

        # Descendant: 4 monomers (removed Phe)
        bb3_null = BuildingBlock.from_code(2, "Null")
        descendant = Compound([bb1, bb2, bb3_null], sample_chromatogram)

        hierarchy.add_compound(ancestor)
        hierarchy.add_compound(descendant)
        hierarchy.add_edge(ancestor, descendant)  # Should succeed (5 > 4)

        assert descendant in hierarchy.edges[ancestor]
