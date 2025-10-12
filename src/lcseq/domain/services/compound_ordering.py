"""Compound ordering domain service.

This module provides domain services for ordering compounds based on
similarity metrics, as specified in THEORY.md Section 8.6.

Domain Responsibility:
    - Order compounds using hierarchical clustering on distance matrices
    - Apply average linkage with optimal leaf ordering
    - Pure algorithmic logic for determining compound order

Not Responsible For:
    - Computing distance metrics (use CompoundSimilarityAnalyzer)
    - Actual visualization or plotting
    - Color assignment or styling

References:
    THEORY.md Section 8.6: Similarity and Ordering
"""

from typing import List
import numpy as np
from scipy.spatial.distance import squareform
from scipy.cluster import hierarchy as scipy_hierarchy

from ..entities import Compound
from .compound_similarity import CompoundSimilarityAnalyzer


class CompoundOrderingService:
    """Domain service for ordering compounds based on similarity.

    As specified in THEORY.md Section 8.6:
    - Uses hierarchical clustering on distance matrix
    - Algorithm: Average linkage with optimal leaf ordering
    - THIS IS DOMAIN LOGIC (not presentation!)
    - Output: Ordered list for visualization

    Examples:
        >>> ordering = CompoundOrderingService()
        >>> ordered = ordering.order_by_similarity(compounds)
    """

    def __init__(self):
        """Initialize ordering service with similarity analyzer."""
        self.similarity_analyzer = CompoundSimilarityAnalyzer()

    def order_by_similarity(
        self,
        compounds: List[Compound],
        metric: str = "wasserstein",
        signal_variant: str = "raw"
    ) -> List[Compound]:
        """Order compounds by chromatographic similarity.

        Uses hierarchical clustering (average linkage with optimal leaf ordering)
        to find an optimal ordering of compounds based on their chromatographic
        similarity.

        As specified in THEORY.md Section 8.6, this is DOMAIN LOGIC because it
        determines the logical ordering of compounds based on their chemical
        properties, not presentation preferences.

        Parameters
        ----------
        compounds : List[Compound]
            Compounds to order
        metric : str, optional
            Distance metric to use. Default is "wasserstein".
        signal_variant : str, optional
            Signal variant to use. Default is "raw" (per THEORY.md - no baseline correction).

        Returns
        -------
        List[Compound]
            Ordered list of compounds

        Notes
        -----
        Algorithm (THEORY.md Section 8.6):
        1. Compute pairwise distance matrix using CompoundSimilarityAnalyzer
        2. Apply hierarchical clustering (average linkage)
        3. Use optimal leaf ordering to minimize adjacent distances
        4. Return compounds in clustered order

        Examples
        --------
        >>> ordering = CompoundOrderingService()
        >>> ordered_compounds = ordering.order_by_similarity(lineage)
        """
        if len(compounds) <= 1:
            return compounds

        # Compute similarity matrix using domain service
        distance_matrix = self.similarity_analyzer.compute_similarity_matrix_robust(
            compounds, metric=metric, signal_variant=signal_variant
        )

        # Apply hierarchical clustering (THEORY.md Section 8.6)
        # Algorithm: Average linkage with optimal leaf ordering
        condensed_dist = squareform(distance_matrix, checks=False)
        linkage_matrix = scipy_hierarchy.linkage(
            condensed_dist, method="average", optimal_ordering=True
        )
        dendrogram = scipy_hierarchy.dendrogram(linkage_matrix, no_plot=True)
        ordered_indices = dendrogram["leaves"]

        # Return ordered compounds
        return [compounds[i] for i in ordered_indices]
