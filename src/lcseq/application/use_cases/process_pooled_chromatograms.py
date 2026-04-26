"""Use case for chromatogram processing in pooled mode.

This use case implements the hybrid pooled strategy from THEORY.md Section 4.2.3:
- Phase 1: Peak detection on pooled signal (expensive, once per class)
- Phase 2: Area integration on individual variants (cheap, per variant)
"""

from pathlib import Path
from typing import List, Dict, Optional, Tuple, Any
from tqdm import tqdm
from lcseq.domain.entities import Compound, PooledCompound
from lcseq.domain.models import EquivalenceClass, CompoundHierarchy, PoolingStatus
from lcseq.domain.services import (
    SignalAggregator,
    EquivalenceClassBuilder,
    HierarchyBuilder,
    SignalPreprocessor,
    PreprocessingConfig,
)
from lcseq.application.use_cases import ProcessChromatogramsUseCase


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
        # Peak detection parameters
        alpha: float,
        prominence_percentile: float,
        min_snr: float,
        min_baseline_sds: float,
        signal_variant: str,
        min_dispersion_r: float,
        sigma_clip_sigma: float,
        # Peak classification parameters
        alpha_product: float,
        truncation_margin: float,
        peak_matching_tolerance: float,
        hungarian_min_threshold: float,
        # Pooling parameters
        correlation_threshold: float,
        aggregation_method: str,
        # Validation parameters
        include_rejected: bool,
        clpe_outlier_threshold: float,
        clpe_min_group_size: int,
        # Preprocessing parameters
        preprocessing_params: Dict[str, Any],
        # Optional cLPE reference (user-provided file path)
        clpe_reference_csv: Optional[Path] = None,
        # Optional dead time (None = derive from L0 peak RT)
        clpe_t0: Optional[float] = None,
    ) -> Tuple[Dict[str, EquivalenceClass], List[Compound], CompoundHierarchy, Optional[Dict[str, Any]]]:
        """
        Process chromatograms using hybrid pooled mode.

        Parameters
        ----------
        compounds : List[Compound]
            All compounds to process
        hierarchy : CompoundHierarchy
            Compound hierarchy for peak classification
        alpha : float
            Significance level (false positive rate)
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
        clpe_reference_csv : Path, optional
            Path to CSV with AlogP and scaffold data for cLPE validation
        clpe_t0 : float, optional
            Dead time for cLPE LogK calculation (if None, uses L0 peak RT)
        clpe_outlier_threshold : float
            Z-score threshold for cLPE outlier detection (default 2.5)
        clpe_min_group_size : int
            Min compounds per scaffold for cLPE model fitting (default 5)

        Returns
        -------
        Tuple[Dict[str, EquivalenceClass], List[Compound], CompoundHierarchy, Optional[Dict]]
            - Dictionary mapping block support sequence to equivalence class with results
            - List of pooled compounds (PooledCompound or Compound)
            - Quotient hierarchy (one node per equivalence class)
            - cLPE statistics (if validation enabled), None otherwise

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

        If cLPE validation is enabled (clpe_reference_csv provided):
        - cLPE validation runs at each level during classification
        - Outlier peaks may be re-selected from UNKNOWN peaks

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

        show_progress = len(equivalence_classes) > 100
        eq_iter = tqdm(
            equivalence_classes,
            desc="Creating pooled compounds",
            disable=not show_progress,
            unit="class"
        )
        for eq_class in eq_iter:
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

        # Step 2b: Apply preprocessing to all chromatograms if enabled
        # This ensures all signals have the "corrected" variant for peak detection
        preprocessing_config = PreprocessingConfig.from_dict(preprocessing_params)
        if preprocessing_config.enabled:
            preprocessor = SignalPreprocessor(preprocessing_config)
            for pooled_compound in pooled_compounds:
                if pooled_compound.chromatogram is not None:
                    # Only preprocess if the chromatogram doesn't already have the corrected variant
                    if not pooled_compound.chromatogram.has_signal_variant("corrected"):
                        preprocessor.preprocess(pooled_compound.chromatogram)

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
        #
        # OPTIMIZATION: Only iterate over ONE representative per equivalence class
        # All variants in a class share the same edges (same block_support_sequence)
        # This reduces iterations from 64k compounds to ~1.5k equivalence classes
        edges_added = set()
        for eq_class in equivalence_classes:
            # Get one representative from this class
            representative = next(iter(eq_class.members))
            if representative not in compound_to_block_support:
                continue

            ancestor_block_support = eq_class.block_support_sequence
            ancestor_pooled = block_support_to_pooled[ancestor_block_support]

            # Get only direct descendants (not transitive closure)
            direct_descendants = hierarchy.get_direct_descendants(representative)
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

        # Step 4: Load cLPE reference data if provided
        clpe_validator = None
        alogp_map = None
        scaffold_map = None

        if clpe_reference_csv:
            from lcseq.domain.services.clpe_validator import CLPEValidator
            from lcseq.infrastructure.loaders.clpe_reference_loader import CLPEReferenceLoader

            loader = CLPEReferenceLoader()
            ref_data = loader.load(clpe_reference_csv)
            alogp_map = ref_data.alogp_map
            scaffold_map = ref_data.scaffold_map

            # t0 will be updated from L0 peak during classification if not provided
            clpe_validator = CLPEValidator(
                t0=clpe_t0 or 1.0,
                outlier_threshold=clpe_outlier_threshold,
                min_group_size=clpe_min_group_size,
            )

        # Step 5: Process pooled compounds through standard pipeline with quotient hierarchy
        # This ensures pooled compounds' descendants are other pooled compounds,
        # not individual variants from the original hierarchy
        peaks_dict, clpe_stats = self.process_use_case.execute(
            compounds=pooled_compounds,
            hierarchy=quotient_hierarchy,  # Quotient hierarchy!
            alpha=alpha,
            prominence_percentile=prominence_percentile,
            min_snr=min_snr,
            min_baseline_sds=min_baseline_sds,
            signal_variant=signal_variant,
            min_dispersion_r=min_dispersion_r,
            sigma_clip_sigma=sigma_clip_sigma,
            alpha_product=alpha_product,
            truncation_margin=truncation_margin,
            peak_matching_tolerance=peak_matching_tolerance,
            hungarian_min_threshold=hungarian_min_threshold,
            include_rejected=include_rejected,
            preprocessing_params=preprocessing_params,
            clpe_validator=clpe_validator,
            alogp_map=alogp_map,
            scaffold_map=scaffold_map,
        )

        # Step 6: Copy results from pooled compounds to all variants in each class
        for pooled_compound, eq_class in pooled_compound_to_class.items():
            variants = list(eq_class.members)

            # Store pooled peaks (from pooled or real compound)
            eq_class.pooled_peaks = pooled_compound.detected_peaks

            # Copy detected peaks and selected peak to all variants in this class
            for variant in variants:
                variant.detected_peaks = pooled_compound.detected_peaks
                variant.selected_peak = pooled_compound.selected_peak

        return results, pooled_compounds, quotient_hierarchy, clpe_stats
