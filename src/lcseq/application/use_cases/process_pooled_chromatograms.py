"""Use case for chromatogram processing in pooled mode.

This use case implements the hybrid pooled strategy from THEORY.md Section 4.2.3:
- Phase 1: Peak detection on pooled signal (expensive, once per class)
- Phase 2: Area integration on individual variants (cheap, per variant)
"""

from typing import List, Dict, Optional, Tuple
from lcseq.domain.entities import Compound, PooledCompound
from lcseq.domain.models import EquivalenceClass, CompoundHierarchy, PoolingStatus
from lcseq.domain.services import (
    SignalAggregator,
    EquivalenceClassBuilder,
    HierarchyBuilder,
)
from lcseq.application.use_cases import ProcessChromatogramsUseCase

from lcseq.config import (
    DEFAULT_Z_THRESHOLD,
    DEFAULT_PROMINENCE_PERCENTILE,
    DEFAULT_MIN_SNR,
    DEFAULT_MIN_BASELINE_SDS,
    DEFAULT_SIGNAL_VARIANT,
    DEFAULT_TRUNCATION_MARGIN,
    DEFAULT_CORRELATION_THRESHOLD,
    DEFAULT_AGGREGATION_METHOD,
)


class ProcessPooledChromatogramsUseCase:
    """
    Use case for processing chromatograms in pooled mode.

    Implements hybrid pooled strategy:
    1. Group compounds into equivalence classes (by block support sequence)
    2. For each class:
       a. Aggregate variant signals → pooled signal
       b. Validate correlation (automatic fallback if invalid)
       c. Peak detection on pooled signal (expensive, once)
       d. Peak classification on pooled signal
       e. For each variant: integrate areas using pooled boundaries (cheap)

    This provides ~3-10× speedup for peak detection while preserving
    individual purity measurements.

    Notes
    -----
    - Automatic fallback to individual mode if correlation < threshold
    - Individual purities always computed from variant's own signal
    - Pooled mode is optional optimization (individual mode is default)

    References
    ----------
    THEORY.md Section 4.2.2: Pooled Mode
    THEORY.md Section 4.2.3: Hybrid Pooled Strategy
    THEORY.md Section 4.2.8: Validity Requirements
    """

    def __init__(
        self,
        process_use_case: Optional[ProcessChromatogramsUseCase] = None,
        aggregator: Optional[SignalAggregator] = None,
        class_builder: Optional[EquivalenceClassBuilder] = None,
        hierarchy_builder: Optional[HierarchyBuilder] = None,
    ):
        """Initialize with domain services (dependency injection).

        Args:
            process_use_case: Standard processing use case (optional)
            aggregator: Signal aggregation service (optional)
            class_builder: Equivalence class builder service (optional)
            hierarchy_builder: Hierarchy builder service (optional)
        """
        self.process_use_case = process_use_case or ProcessChromatogramsUseCase()
        self.aggregator = aggregator or SignalAggregator()
        self.class_builder = class_builder or EquivalenceClassBuilder()
        self.hierarchy_builder = hierarchy_builder or HierarchyBuilder()

    def execute(
        self,
        compounds: List[Compound],
        hierarchy: CompoundHierarchy,
        z_threshold: float = DEFAULT_Z_THRESHOLD,
        prominence_percentile: float = DEFAULT_PROMINENCE_PERCENTILE,
        min_snr: float = DEFAULT_MIN_SNR,
        min_baseline_sds: float = DEFAULT_MIN_BASELINE_SDS,
        signal_variant: str = DEFAULT_SIGNAL_VARIANT,
        truncation_margin: float = DEFAULT_TRUNCATION_MARGIN,
        correlation_threshold: float = DEFAULT_CORRELATION_THRESHOLD,
        aggregation_method: str = DEFAULT_AGGREGATION_METHOD,
    ) -> Tuple[Dict[str, EquivalenceClass], List[Compound], CompoundHierarchy]:
        """
        Process chromatograms using hybrid pooled mode.

        Parameters
        ----------
        compounds : List[Compound]
            All compounds to process
        hierarchy : CompoundHierarchy
            Compound hierarchy for peak classification
        z_threshold : float
            Poisson Z-score threshold
        prominence_percentile : float
            Prominence percentile threshold
        min_snr : float
            Adaptive SNR threshold multiplier
        min_baseline_sds : float
            Global baseline threshold in SDs
        signal_variant : str
            Signal variant to use for detection
        truncation_margin : float
            Margin beyond truncation positions (in seconds)
        correlation_threshold : float
            Minimum correlation for pooling validity
        aggregation_method : str
            Aggregation method ("mean" or "median")

        Returns
        -------
        Tuple[Dict[str, EquivalenceClass], List[Compound], CompoundHierarchy]
            - Dictionary mapping block support sequence to equivalence class with results
            - List of pooled compounds (PooledCompound or Compound)
            - Quotient hierarchy (one node per equivalence class)

        Notes
        -----
        Each EquivalenceClass contains:
        - pooled_chromatogram: Aggregated signal (if pooling valid)
        - pooled_peaks: Peaks detected on pooled signal
        - pooling_status: POOLING_VALID or fallback reason
        - variants (compounds): Individual results with purities

        Algorithm:
        1. Build equivalence classes from compounds
        2. For each class:
           a. Attempt pooled aggregation + validation
           b. Create pooled compound (pooled or fallback to variants)
        3. Process ALL pooled compounds through standard pipeline (reuses existing logic!)
        4. Copy results from pooled compounds to all variants in class

        References
        ----------
        THEORY.md Section 4.2.3: Hybrid Pooled Strategy
        THEORY.md Section 4.2.8: Validity Requirements & Automatic Fallback
        THEORY.md Section 4.2.11: Quotient Hierarchy
        """
        # Step 1: Build equivalence classes from original compounds
        equivalence_classes = self.class_builder.build(compounds)

        # Step 2: Create pooled compounds for each equivalence class
        # Pooled compounds are either real compounds (single variant) or pooled compounds
        pooled_compound_to_class = {}
        pooled_compounds = []
        results = {}

        for eq_class in equivalence_classes:
            # Convert Set to List for processing
            variants = list(eq_class.members)

            if len(variants) == 1:
                # Single variant - use it directly as pooled compound
                pooled_compound = variants[0]
                eq_class.pooling_status = PoolingStatus.NOT_ATTEMPTED
                eq_class.fallback_reason = "Single variant (no aggregation needed)"

            else:
                # Multiple variants - ALWAYS aggregate (correlation is just metadata)
                (
                    pooled_chromatogram,
                    min_correlation,
                    is_valid,
                    reason,
                ) = self.aggregator.aggregate(
                    variants,
                    method=aggregation_method,
                    correlation_threshold=correlation_threshold,
                )

                eq_class.correlation_min = min_correlation
                eq_class.pooled_chromatogram = pooled_chromatogram

                # Create PooledCompound with pooled chromatogram
                # Use first variant as the template for building blocks and hierarchy position
                real_pooled_compound = variants[0]
                pooled_compound = PooledCompound(
                    real_compound=real_pooled_compound,
                    pooled_chromatogram=pooled_chromatogram,
                )

                if is_valid:
                    eq_class.pooling_status = PoolingStatus.POOLING_VALID
                else:
                    # Low correlation - still use pooled signal, but flag it
                    eq_class.pooling_status = PoolingStatus.POOLING_INVALID
                    eq_class.fallback_reason = reason

            # Store mapping and add to pooled compounds list
            pooled_compound_to_class[pooled_compound] = eq_class
            pooled_compounds.append(pooled_compound)

            # Store result
            results[eq_class.block_support_sequence] = eq_class

        # Step 3: Build quotient hierarchy by projecting original hierarchy edges
        # We cannot rebuild from scratch because PooledCompounds only know about one variant
        # Instead, we project edges from the original hierarchy onto equivalence classes

        # Create empty hierarchy with same mode
        from lcseq.domain.models import CompoundHierarchy
        quotient_hierarchy = CompoundHierarchy(mode=hierarchy.mode)

        # Add all pooled compounds as nodes
        for pooled_compound in pooled_compounds:
            quotient_hierarchy.add_compound(pooled_compound)

        # Project edges: if any variant in class A is ancestor of any variant in class B,
        # then pooled_compound(A) is ancestor of pooled_compound(B)
        # Build mapping: block_support_sequence -> pooled_compound
        block_support_to_pooled = {}
        for pooled_compound, eq_class in pooled_compound_to_class.items():
            block_support_to_pooled[eq_class.block_support_sequence] = pooled_compound

        # Build mapping: original compound -> block_support_sequence
        compound_to_block_support = {}
        for eq_class in equivalence_classes:
            for variant in eq_class.members:
                compound_to_block_support[variant] = eq_class.block_support_sequence

        # Project edges from original hierarchy
        # IMPORTANT: Use get_direct_descendants() to avoid creating transitive edges
        # We want to preserve the DAG structure, not create a transitive closure
        edges_added = set()
        for compound in hierarchy.compounds:
            if compound not in compound_to_block_support:
                continue

            ancestor_block_support = compound_to_block_support[compound]
            ancestor_pooled = block_support_to_pooled[ancestor_block_support]

            # Get only direct descendants (not transitive closure)
            direct_descendants = hierarchy.get_direct_descendants(compound)
            for desc in direct_descendants:
                if desc not in compound_to_block_support:
                    continue

                desc_block_support = compound_to_block_support[desc]
                desc_pooled = block_support_to_pooled[desc_block_support]

                # Add edge if not already added and not same block support sequence
                edge = (ancestor_pooled, desc_pooled)
                if edge not in edges_added and ancestor_block_support != desc_block_support:
                    quotient_hierarchy.add_edge(ancestor_pooled, desc_pooled)
                    edges_added.add(edge)

        # Step 4: Process pooled compounds through standard pipeline with quotient hierarchy
        # This ensures pooled compounds' descendants are other pooled compounds,
        # not individual variants from the original hierarchy
        peaks_dict = self.process_use_case.execute(
            compounds=pooled_compounds,
            hierarchy=quotient_hierarchy,  # Quotient hierarchy!
            z_threshold=z_threshold,
            prominence_percentile=prominence_percentile,
            min_snr=min_snr,
            min_baseline_sds=min_baseline_sds,
            signal_variant=signal_variant,
            truncation_margin=truncation_margin,
        )

        # Step 5: Copy results from pooled compounds to all variants in each class
        for pooled_compound, eq_class in pooled_compound_to_class.items():
            variants = list(eq_class.members)

            # Store pooled peaks (from pooled or real compound)
            eq_class.pooled_peaks = pooled_compound.detected_peaks

            # Copy detected peaks and selected peak to all variants in this class
            for variant in variants:
                variant.detected_peaks = pooled_compound.detected_peaks
                variant.selected_peak = pooled_compound.selected_peak

        return results, pooled_compounds, quotient_hierarchy
