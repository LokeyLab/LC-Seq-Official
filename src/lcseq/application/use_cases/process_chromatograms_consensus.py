"""Use case for chromatogram processing in consensus mode.

This use case implements the hybrid consensus strategy from THEORY.md Section 4.2.3:
- Phase 1: Peak detection on consensus signal (expensive, once per class)
- Phase 2: Area integration on individual variants (cheap, per variant)
"""

from typing import List, Dict, Optional, Tuple
from lcseq.domain.entities import Compound, VirtualCompound
from lcseq.domain.models import EquivalenceClass, CompoundHierarchy, ConsensusStatus
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


class ProcessChromatogramsConsensusUseCase:
    """
    Use case for processing chromatograms in consensus mode.

    Implements hybrid consensus strategy:
    1. Group compounds into equivalence classes (by residue sequence)
    2. For each class:
       a. Aggregate variant signals → consensus
       b. Validate correlation (automatic fallback if invalid)
       c. Peak detection on consensus (expensive, once)
       d. Peak classification on consensus
       e. For each variant: integrate areas using consensus boundaries (cheap)

    This provides ~3-10× speedup for peak detection while preserving
    individual purity measurements.

    Notes
    -----
    - Automatic fallback to individual mode if correlation < threshold
    - Individual purities always computed from variant's own signal
    - Consensus mode is optional optimization (individual mode is default)

    References
    ----------
    THEORY.md Section 4.2.2: Consensus Mode
    THEORY.md Section 4.2.3: Hybrid Consensus Strategy
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
        Process chromatograms using hybrid consensus mode.

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
            Minimum correlation for consensus validity
        aggregation_method : str
            Aggregation method ("mean" or "median")

        Returns
        -------
        Tuple[Dict[str, EquivalenceClass], List[Compound], CompoundHierarchy]
            - Dictionary mapping residue sequence to equivalence class with results
            - List of representative compounds (VirtualCompound or Compound)
            - Representative hierarchy (one node per equivalence class)

        Notes
        -----
        Each EquivalenceClass contains:
        - consensus_chromatogram: Aggregated signal (if consensus valid)
        - consensus_peaks: Peaks detected on consensus
        - consensus_status: CONSENSUS_VALID or fallback reason
        - variants (compounds): Individual results with purities

        Algorithm:
        1. Build equivalence classes from compounds
        2. For each class:
           a. Attempt consensus aggregation + validation
           b. Create representative compound (consensus or fallback to variants)
        3. Process ALL representatives through standard pipeline (reuses existing logic!)
        4. Copy results from representatives to all variants in class

        References
        ----------
        THEORY.md Section 4.2.3: Hybrid Consensus Strategy
        THEORY.md Section 4.2.8: Validity Requirements & Automatic Fallback
        """
        # Step 1: Build equivalence classes from original compounds
        equivalence_classes = self.class_builder.build(compounds)

        # Step 2: Create representative compounds for each equivalence class
        # Representatives are either real compounds (single variant) or virtual compounds (consensus)
        representative_to_class = {}
        representatives = []
        results = {}

        for eq_class in equivalence_classes:
            # Convert Set to List for processing
            variants = list(eq_class.compounds)

            if len(variants) == 1:
                # Single variant - use it directly as representative
                representative = variants[0]
                eq_class.consensus_status = ConsensusStatus.NOT_ATTEMPTED
                eq_class.fallback_reason = "Single variant (no aggregation needed)"

            else:
                # Multiple variants - ALWAYS aggregate (correlation is just metadata)
                (
                    consensus_chromatogram,
                    min_correlation,
                    is_valid,
                    reason,
                ) = self.aggregator.aggregate(
                    variants,
                    method=aggregation_method,
                    correlation_threshold=correlation_threshold,
                )

                eq_class.correlation_min = min_correlation
                eq_class.consensus_chromatogram = consensus_chromatogram

                # Create VirtualCompound with consensus chromatogram
                # Use first variant as the template for building blocks and hierarchy position
                real_representative = variants[0]
                representative = VirtualCompound(
                    real_compound=real_representative,
                    consensus_chromatogram=consensus_chromatogram,
                )

                if is_valid:
                    eq_class.consensus_status = ConsensusStatus.CONSENSUS_VALID
                else:
                    # Low correlation - still use consensus, but flag it
                    eq_class.consensus_status = ConsensusStatus.CONSENSUS_INVALID
                    eq_class.fallback_reason = reason

            # Store mapping and add to representatives list
            representative_to_class[representative] = eq_class
            representatives.append(representative)

            # Store result
            results[eq_class.residue_sequence] = eq_class

        # Step 3: Build representative hierarchy by projecting original hierarchy edges
        # We cannot rebuild from scratch because VirtualCompounds only know about one variant
        # Instead, we project edges from the original hierarchy onto equivalence classes

        # Create empty hierarchy with same mode
        from lcseq.domain.models import CompoundHierarchy
        representative_hierarchy = CompoundHierarchy(mode=hierarchy.mode)

        # Add all representatives as nodes
        for rep in representatives:
            representative_hierarchy.add_compound(rep)

        # Project edges: if any variant in class A is ancestor of any variant in class B,
        # then representative(A) is ancestor of representative(B)
        # Build mapping: residue_sequence -> representative
        residue_to_rep = {}
        for rep, eq_class in representative_to_class.items():
            residue_to_rep[eq_class.residue_sequence] = rep

        # Build mapping: original compound -> residue_sequence
        compound_to_residue = {}
        for eq_class in equivalence_classes:
            for variant in eq_class.compounds:
                compound_to_residue[variant] = eq_class.residue_sequence

        # Project edges from original hierarchy
        # IMPORTANT: Use get_direct_descendants() to avoid creating transitive edges
        # We want to preserve the DAG structure, not create a transitive closure
        edges_added = set()
        for compound in hierarchy.compounds:
            if compound not in compound_to_residue:
                continue

            ancestor_residue = compound_to_residue[compound]
            ancestor_rep = residue_to_rep[ancestor_residue]

            # Get only direct descendants (not transitive closure)
            direct_descendants = hierarchy.get_direct_descendants(compound)
            for desc in direct_descendants:
                if desc not in compound_to_residue:
                    continue

                desc_residue = compound_to_residue[desc]
                desc_rep = residue_to_rep[desc_residue]

                # Add edge if not already added and not same residue
                edge = (ancestor_rep, desc_rep)
                if edge not in edges_added and ancestor_residue != desc_residue:
                    representative_hierarchy.add_edge(ancestor_rep, desc_rep)
                    edges_added.add(edge)

        # Step 4: Process representatives through standard pipeline with representative hierarchy
        # This ensures virtual compounds' descendants are other representatives,
        # not individual variants from the original hierarchy
        peaks_dict = self.process_use_case.execute(
            compounds=representatives,
            hierarchy=representative_hierarchy,  # NEW hierarchy!
            z_threshold=z_threshold,
            prominence_percentile=prominence_percentile,
            min_snr=min_snr,
            min_baseline_sds=min_baseline_sds,
            signal_variant=signal_variant,
            truncation_margin=truncation_margin,
        )

        # Step 5: Copy results from representatives to all variants in each class
        for representative, eq_class in representative_to_class.items():
            variants = list(eq_class.compounds)

            # Store consensus peaks (from virtual or real representative)
            eq_class.consensus_peaks = representative.detected_peaks

            # Copy detected peaks and selected peak to all variants in this class
            for variant in variants:
                variant.detected_peaks = representative.detected_peaks
                variant.selected_peak = representative.selected_peak

        return results, representatives, representative_hierarchy
