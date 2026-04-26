"""
Peak classification service using DAG constraint propagation.

Implementation based on THEORY.md Section 5.3-5.6.
"""

import numpy as np
from typing import List, Optional, Tuple, Dict, Any, TYPE_CHECKING
from scipy.optimize import linear_sum_assignment
from tqdm import tqdm
from ..entities.compound import Compound
from ..entities.peak import Peak, PeakType, RejectionReason
from ..models.compound_hierarchy import CompoundHierarchy

if TYPE_CHECKING:
    from .clpe_validator import CLPEValidator


class PeakClassifier:
    """
    Classifies peaks as NULL/TRUNCATION/PUTATIVE_PRODUCT/UNKNOWN.

    Peak classification is based on retention time position relative to
    expected elution order from DAG constraints. This is a POSITIONAL
    classification, NOT chemical validation.

    This is a stateless service - all methods are operations on input data with
    no instance state.

    Notes
    -----
    Classification types:
    - NULL: Peak at L₀ (full-null) position
    - TRUNCATION: Peak matching ancestor product or null position
    - PUTATIVE_PRODUCT: First persistent peak after truncations
    - UNKNOWN: Peak not matching any expected position

    CRITICAL: Putative product is a positional hypothesis, NOT chemically
    validated. Synthesis validation requires additional evidence (Part 6).

    Classification respects hierarchical constraints:
    - position(truncation) < position(product)
    - Ancestor products become descendant truncations
    - Global consistency across entire lineage DAG

    References
    ----------
    THEORY.md Section 5.3: Peak Type Classification
    THEORY.md Section 5.4: Global Classification via Constraint Propagation
    THEORY.md Section 5.6: Classification Scope and Limitations
    THEORY.md Section 6.13: Classification ≠ Validation

    Examples
    --------
    >>> classifier = PeakClassifier()
    >>> hierarchy = CompoundHierarchy(compounds, mode=HierarchyMode.BUILDING_BLOCK)
    >>> # Classify all compounds in hierarchy (bottom-up)
    >>> classifier.classify_hierarchy(hierarchy)
    >>> # Access classification results
    >>> compound.selected_peak.peak_type
    <PeakType.PUTATIVE_PRODUCT: 'PUTATIVE_PRODUCT'>
    """

    def _optimal_assignment(
        self,
        detected_peaks: List[Peak],
        expected_positions: List[float],
        signal_length: float,
        hungarian_min_threshold: float,
    ) -> Dict[int, int]:
        """
        Find optimal assignment of detected peaks to expected positions using Hungarian algorithm.

        Implements THEORY.md Section 5.5: Optimal Assignment and Matching.

        Parameters
        ----------
        detected_peaks : List[Peak]
            Detected peaks D = {d₁, d₂, ...}
        expected_positions : List[float]
            Expected positions E = {e₁, e₂, ...} (L₀ + truncation products)
        signal_length : float
            Total signal length for scale-invariant normalization

        Returns
        -------
        Dict[int, int]
            Mapping from detected_peak_index -> expected_position_index
            Only includes accepted assignments (cost < threshold)

        Notes
        -----
        Algorithm (THEORY.md Section 5.5):

        1. **Cost Matrix**: Cost(i,j) = |peak_i.position - expected_j| / signal_length
           - Scale-invariant normalization
           - Measures relative positional distance

        2. **Hungarian Algorithm**: Finds optimal assignment minimizing total cost
           - O(n³) time complexity
           - Globally optimal (not greedy)

        3. **Adaptive Threshold**: threshold = median(peak_spacings) / signal_length
           - Data-driven (no magic numbers)
           - Only accepts assignments where Cost(i,j) < threshold
           - Ensures peaks are genuinely close to expected positions

        References
        ----------
        THEORY.md Section 5.5: Optimal Assignment and Matching
        """
        if not detected_peaks or not expected_positions:
            return {}

        if signal_length <= 0:
            raise ValueError(f"Signal length must be positive, got {signal_length}")

        # Build cost matrix (THEORY.md Section 5.5)
        # OPTIMIZATION: Vectorized cost matrix using broadcasting (50-100x faster)
        # Cost(i,j) = |peak_i.position - expected_j| / signal_length
        peak_positions = np.array([p.position for p in detected_peaks])
        expected_array = np.array(expected_positions)
        # Broadcasting: (n_peaks, 1) - (n_expected,) -> (n_peaks, n_expected)
        cost_matrix = np.abs(peak_positions[:, np.newaxis] - expected_array) / signal_length

        # Solve optimal assignment using Hungarian algorithm
        row_ind, col_ind = linear_sum_assignment(cost_matrix)

        # Compute adaptive threshold (THEORY.md Section 5.5)
        # OPTIMIZATION: Vectorized peak spacing calculation
        sorted_expected = np.sort(expected_array)
        peak_spacings = np.diff(sorted_expected)  # Vectorized difference

        if len(peak_spacings) > 0:
            # Adaptive threshold based on median spacing between expected positions
            median_spacing = np.median(peak_spacings)
            threshold = max(median_spacing / signal_length, hungarian_min_threshold)
        else:
            # Single expected position: use minimum threshold from config
            threshold = hungarian_min_threshold

        # Accept only assignments below or equal to threshold
        # Use <= to ensure perfect matches (cost=0.0) are accepted when threshold=0
        assignments = {}
        for i, j in zip(row_ind, col_ind):
            if cost_matrix[i, j] <= threshold:
                assignments[i] = j

        return assignments

    def _find_l0_retention_time(self, hierarchy: CompoundHierarchy) -> Optional[float]:
        """
        Find L₀ (all-null) compound retention time from hierarchy.

        Classifies L0 compound if not already classified, then returns
        the selected peak position.

        Parameters
        ----------
        hierarchy : CompoundHierarchy
            Compound hierarchy

        Returns
        -------
        float or None
            L₀ peak position if found, None otherwise

        Notes
        -----
        Algorithm (THEORY.md Section 2.4.4):
        1. Get L0 compound from hierarchy
        2. Classify L0 if not already classified (selects earliest peak)
        3. Return selected peak position

        References
        ----------
        THEORY.md Section 2.4.4: L₀ Peak Detection Algorithm
        """
        # Find L₀ compound (level 0)
        l0_compounds = hierarchy.get_level(0)
        if not l0_compounds:
            return None

        l0 = l0_compounds[0]

        # Classify L0 if not already classified
        if l0.selected_peak is None and l0.detected_peaks:
            self._classify_l0_compound(l0)

        # Return selected peak position
        if l0.selected_peak is not None:
            return l0.selected_peak.position

        return None

    def _find_l0_retention_time_from_compounds(
        self,
        l0_compounds: List[Compound]
    ) -> Optional[float]:
        """
        Find L₀ retention time from pre-classified L0 compounds.

        Note: L0 compounds must be classified first via _classify_l0_compound()
        before calling this method.

        Parameters
        ----------
        l0_compounds : List[Compound]
            Pre-fetched and pre-classified L0 compounds

        Returns
        -------
        float or None
            L₀ peak position if found, None otherwise
        """
        if not l0_compounds:
            return None

        # Should be exactly one L₀
        l0 = l0_compounds[0]

        # Get selected peak (set by _classify_l0_compound)
        if l0.selected_peak is not None:
            return l0.selected_peak.position

        return None

    def _get_product_positions(self, compounds: List[Compound]) -> List[float]:
        """
        Extract product peak positions from list of compounds.

        Parameters
        ----------
        compounds : List[Compound]
            Compounds to extract positions from

        Returns
        -------
        List[float]
            Product peak positions (non-null)

        Notes
        -----
        Returns positions of selected_peak (putative product) for each compound.
        """
        positions = []
        for compound in compounds:
            if compound.selected_peak is not None:
                positions.append(compound.selected_peak.position)

        return positions

    def _get_all_descendant_peaks(
        self,
        compounds: List[Compound]
    ) -> List[Tuple[Compound, Peak]]:
        """
        Get all ACCEPTED detected peaks from compounds with their source compound.

        Parameters
        ----------
        compounds : List[Compound]
            Compounds to extract peaks from

        Returns
        -------
        List[Tuple[Compound, Peak]]
            List of (compound, peak) tuples for matching

        Notes
        -----
        Returns detected peaks for matching unknown peaks to descendant peaks.

        IMPORTANT: Only returns ACCEPTED peaks (not rejected). Rejected peaks
        represent noise/false positives and should not be used for truncation
        matching in ancestor compounds.

        Special handling for L0 (null compound):
        - Only includes the NULL peak, not UNKNOWN peaks
        - L0's UNKNOWN peaks are noise/contaminants, not meaningful species
        - L0 has no product to form oligomers from

        For all other levels (L1+):
        - Includes all accepted peaks (products, truncations, unknowns)
        - UNKNOWN peaks at L1+ could be oligomers worth propagating
        """
        result = []
        for compound in compounds:
            is_l0 = compound.level == 0
            for peak in compound.detected_peaks:
                # Skip rejected peaks - they're noise, not real peaks
                if peak.is_rejected:
                    continue
                # Skip L0's UNKNOWN peaks - they're noise, not meaningful species
                if is_l0 and peak.peak_type == PeakType.UNKNOWN:
                    continue
                result.append((compound, peak))
        return result

    def _find_matching_descendant_peak(
        self,
        peak: Peak,
        descendant_peaks: List[Tuple[Compound, Peak]],
        tolerance: float
    ) -> Optional[Tuple[Compound, Peak]]:
        """
        Find best matching peak among descendants.

        Parameters
        ----------
        peak : Peak
            Peak to match
        descendant_peaks : List[Tuple[Compound, Peak]]
            List of (compound, peak) tuples from descendants
        tolerance : float
            Relative tolerance for position matching

        Returns
        -------
        Optional[Tuple[Compound, Peak]]
            (compound, peak) tuple if match found, None otherwise

        Notes
        -----
        Finds the closest matching peak within tolerance.
        Uses relative distance: |peak.position - desc.position| / desc.position
        """
        best_match = None
        best_distance = float('inf')

        for compound, desc_peak in descendant_peaks:
            if desc_peak.position <= 0:
                continue
            rel_distance = abs(peak.position - desc_peak.position) / desc_peak.position
            if rel_distance <= tolerance and rel_distance < best_distance:
                best_distance = rel_distance
                best_match = (compound, desc_peak)

        return best_match

    def classify_hierarchy(
        self,
        hierarchy: CompoundHierarchy,
        tolerance: float,
        truncation_margin: float,
        alpha_product: float,
        hungarian_min_threshold: float,
        # Optional cLPE validation parameters
        clpe_validator: Optional["CLPEValidator"] = None,
        alogp_map: Optional[Dict[str, float]] = None,
        scaffold_map: Optional[Dict[str, str]] = None,
    ) -> Optional[Dict[str, Any]]:
        """
        Classify all peaks in hierarchy, processing bottom-up.

        Ensures descendants are classified before ancestors,
        enabling match propagation up the hierarchy.

        Optionally applies cLPE (chromatographic Linear Peptide Equation)
        validation at each level to validate and potentially re-select
        product peaks based on LogK ~ AlogP correlation.

        Parameters
        ----------
        hierarchy : CompoundHierarchy
            Compound hierarchy to classify
        tolerance : float, optional
            Relative tolerance for position matching
        truncation_margin : float, optional
            Margin beyond truncation positions (in seconds)
        alpha_product : float, optional
            Stricter significance threshold for product selection.
            Product candidates must have p_value < alpha_product.
        clpe_validator : CLPEValidator, optional
            If provided, enables cLPE validation at each level
        alogp_map : Dict[str, float], optional
            Mapping from compound identifier to AlogP value
        scaffold_map : Dict[str, str], optional
            Mapping from compound identifier to scaffold group

        Returns
        -------
        Optional[Dict[str, Any]]
            cLPE statistics if validation was enabled, None otherwise

        Notes
        -----
        Processing order: L0 -> L1 -> L2 -> ... -> Lmax

        This ensures that when classifying a compound at level N,
        all its descendants (levels 0 to N-1) have already been
        classified, enabling:
        1. Product positions are known for truncation matching
        2. Unknown peaks are identified for TRUNCATION_UNKNOWN matching
        3. Match chains are built progressively up the hierarchy

        If cLPE validation is enabled:
        - After classifying L0, t0 is extracted from the NULL peak
        - After classifying each level > 0, cLPE models are fit per scaffold
        - Product peaks are validated and potentially re-selected
        - Re-selection only considers UNKNOWN peaks (not TRUNCATION)

        References
        ----------
        THEORY.md Section 5.4: Global Classification via Constraint Propagation
        """
        # Find max level in hierarchy
        if hierarchy.mode.value == "building_block":
            max_level = max((c.level for c in hierarchy.compounds), default=0)
        else:
            max_level = max((c.monomer_level for c in hierarchy.compounds), default=0)

        # OPTIMIZATION: Cache L0 compounds ONCE (not per compound!)
        l0_compounds = hierarchy.get_level(0)
        l0_compounds_set = set(l0_compounds)  # O(1) lookup

        # Classify L0 compounds FIRST - this sets their selected_peak
        # L0 must be classified before we can get l0_retention_time
        for compound in l0_compounds:
            self._classify_l0_compound(compound)

        # Now get L0 retention time from classified L0 compound
        l0_retention_time = None
        if l0_compounds and l0_compounds[0].selected_peak is not None:
            l0_retention_time = l0_compounds[0].selected_peak.position

        # Initialize cLPE stats if validation enabled
        clpe_stats = {"levels": {}, "t0": None} if clpe_validator else None

        # Extract t0 from L0 for cLPE (already classified above)
        if clpe_validator and l0_retention_time is not None:
            t0 = l0_retention_time / 60.0  # sec to min
            clpe_validator.t0 = t0
            if clpe_stats:
                clpe_stats["t0"] = t0

        # Process bottom-up: L1, L2, ... (L0 already done)
        total_compounds = hierarchy.size()
        show_progress = total_compounds > 100

        with tqdm(total=total_compounds, desc="Classifying peaks", disable=not show_progress, unit="cpd") as pbar:
            # Update progress for already-classified L0 compounds
            pbar.update(len(l0_compounds))

            for level in range(1, max_level + 1):  # Start at 1, L0 already classified
                compounds_at_level = hierarchy.get_level(level)
                for compound in compounds_at_level:
                    self._classify_compound(
                        compound,
                        hierarchy,
                        l0_retention_time,
                        l0_compounds_set,  # Pass cached L0 set
                        tolerance,
                        truncation_margin,
                        alpha_product,
                        hungarian_min_threshold,
                    )
                    pbar.update(1)

                # cLPE validation after classifying this level (L1+)
                if clpe_validator and alogp_map and scaffold_map:
                    # Apply cLPE validation at this level
                    level_stats = self._apply_clpe_at_level(
                        compounds_at_level,
                        clpe_validator,
                        alogp_map,
                        scaffold_map,
                    )
                    clpe_stats["levels"][level] = level_stats

        return clpe_stats

    def _classify_compound(
        self,
        compound: Compound,
        hierarchy: CompoundHierarchy,
        l0_retention_time: Optional[float],
        l0_compounds_set: set,
        tolerance: float,
        truncation_margin: float,
        alpha_product: float,
        hungarian_min_threshold: float,
    ) -> None:
        """
        Classify peaks for a single compound with descendant matching.

        Parameters
        ----------
        compound : Compound
            Compound to classify
        hierarchy : CompoundHierarchy
            Compound hierarchy
        l0_retention_time : float or None
            L0 retention time (cached, computed once)
        l0_compounds_set : set
            Set of L0 compounds for O(1) lookup (cached, computed once)
        tolerance : float
            Relative tolerance for position matching
        truncation_margin : float
            Margin beyond truncation positions (in seconds)
        alpha_product : float, optional
            Stricter significance threshold for product selection.
            Product candidates must have p_value < alpha_product.

        Notes
        -----
        Classification flow:
        1. Handle L0 specially (NULL peak identification)
        2. Hungarian match to descendant product positions -> TRUNCATION
        3. Match remaining UNKNOWN to ANY descendant peak -> TRUNCATION_UNKNOWN
        4. Select product from unmatched peaks after boundary (filtered by alpha_product)
        """
        if not compound.detected_peaks:
            return

        # Special case: L0 (all-null) compound - O(1) lookup using cached set
        is_l0 = compound in l0_compounds_set

        if is_l0:
            self._classify_l0_compound(compound)
            return

        # Get only accepted peaks for main classification
        accepted_peaks = [p for p in compound.detected_peaks if p.is_accepted]
        if not accepted_peaks:
            # No accepted peaks - still classify rejected peaks
            descendants = hierarchy.get_descendants(compound)
            truncation_product_positions = self._get_product_positions(descendants)
            max_truncation_pos = 0.0
            if l0_retention_time is not None:
                max_truncation_pos = l0_retention_time
            if truncation_product_positions:
                max_truncation_pos = max(max_truncation_pos, max(truncation_product_positions))
            max_truncation_boundary = max_truncation_pos + truncation_margin
            self._classify_rejected_peaks(
                compound, descendants, l0_retention_time,
                max_truncation_boundary, tolerance
            )
            return

        # Get descendants
        descendants = hierarchy.get_descendants(compound)
        truncation_product_positions = self._get_product_positions(descendants)

        # Build expected truncation positions for Hungarian matching (L0 handled separately)
        truncation_positions = []
        truncation_compounds = []  # Track which compound each position came from

        for desc_compound in descendants:
            if desc_compound.selected_peak is not None:
                pos = desc_compound.selected_peak.position
                if pos > 0:
                    truncation_positions.append(pos)
                    truncation_compounds.append(desc_compound)

        # Calculate truncation boundary
        max_truncation_pos = 0.0
        if l0_retention_time is not None:
            max_truncation_pos = l0_retention_time
        if truncation_product_positions:
            max_truncation_pos = max(max_truncation_pos, max(truncation_product_positions))
        max_truncation_boundary = max_truncation_pos + truncation_margin

        # Get signal length for threshold calculations
        time_points = compound.chromatogram.time_points
        signal_length = time_points[-1] - time_points[0]

        assigned_peak_indices = set()

        # Step 1a: Match L0 FIRST with priority (before Hungarian matching)
        # This prevents the Hungarian algorithm from swapping L0 with truncation peaks
        if l0_retention_time is not None and l0_retention_time > 0:
            # Find the peak closest to L0 retention time
            best_l0_idx = None
            best_l0_cost = float('inf')

            for i, peak in enumerate(accepted_peaks):
                cost = abs(peak.position - l0_retention_time) / signal_length
                if cost < best_l0_cost:
                    best_l0_cost = cost
                    best_l0_idx = i

            # Accept if within threshold (use hungarian_min_threshold as baseline)
            if best_l0_idx is not None and best_l0_cost <= hungarian_min_threshold:
                accepted_peaks[best_l0_idx].peak_type = PeakType.NULL
                assigned_peak_indices.add(best_l0_idx)

        # Step 1b: Hungarian match truncation positions (excluding L0 peak)
        if truncation_positions:
            # Build list of remaining peaks (not assigned to L0)
            remaining_peaks = [
                (i, peak) for i, peak in enumerate(accepted_peaks)
                if i not in assigned_peak_indices
            ]

            if remaining_peaks:
                remaining_indices, remaining_peak_list = zip(*remaining_peaks)
                remaining_peak_list = list(remaining_peak_list)

                assignments = self._optimal_assignment(
                    remaining_peak_list,
                    truncation_positions,
                    signal_length,
                    hungarian_min_threshold,
                )

                for local_idx, trunc_idx in assignments.items():
                    original_idx = remaining_indices[local_idx]
                    peak = accepted_peaks[original_idx]
                    matched_compound = truncation_compounds[trunc_idx]

                    # Validate: TRUNCATION peaks must be before boundary
                    if peak.position > max_truncation_boundary:
                        continue  # Don't assign, let it be reconsidered

                    peak.peak_type = PeakType.TRUNCATION
                    assigned_peak_indices.add(original_idx)

                    # Record match tracking info
                    peak.matched_compound_sequence = matched_compound.block_support_sequence
                    peak.matched_peak_position = matched_compound.selected_peak.position
                    peak.matched_peak_type = PeakType.PUTATIVE_PRODUCT  # Matched to product

        # Step 2: Match remaining UNKNOWN accepted peaks to ANY descendant peak
        all_descendant_peaks = self._get_all_descendant_peaks(descendants)

        for i, peak in enumerate(accepted_peaks):
            if i in assigned_peak_indices:
                continue

            # Try to match to any descendant peak
            match = self._find_matching_descendant_peak(peak, all_descendant_peaks, tolerance)

            if match:
                matched_compound, matched_peak = match
                peak.peak_type = PeakType.TRUNCATION_UNKNOWN

                # Trace back to original source if matched peak is itself a propagated peak
                # This ensures we track where the peak ORIGINATED, not just the immediate match
                if matched_peak.matched_compound_sequence:
                    # The matched peak was itself matched from somewhere deeper - use that origin
                    peak.matched_compound_sequence = matched_peak.matched_compound_sequence
                    peak.matched_peak_position = matched_peak.matched_peak_position
                    peak.matched_peak_type = matched_peak.matched_peak_type
                else:
                    # The matched peak is the original source
                    peak.matched_compound_sequence = matched_compound.block_support_sequence
                    peak.matched_peak_position = matched_peak.position
                    peak.matched_peak_type = matched_peak.peak_type

                assigned_peak_indices.add(i)

        # Step 3: Select product from remaining unmatched accepted peaks after boundary
        product_candidates = []
        for i, peak in enumerate(accepted_peaks):
            if i in assigned_peak_indices:
                continue

            if peak.position > max_truncation_boundary:
                # Apply stricter alpha_product filter for product selection
                # Peaks for purity use permissive alpha, but product needs high confidence
                if peak.p_value is not None and peak.p_value < alpha_product:
                    product_candidates.append(peak)
                else:
                    # Doesn't meet product threshold - mark as UNKNOWN but keep for purity
                    peak.peak_type = PeakType.UNKNOWN
            else:
                peak.peak_type = PeakType.UNKNOWN

        # Select product by earliest position - product elutes just after truncations
        if product_candidates:
            product_peak = min(product_candidates, key=lambda p: p.position)
            product_peak.peak_type = PeakType.PUTATIVE_PRODUCT
            compound.selected_peak = product_peak

            # Mark remaining candidates as UNKNOWN
            for peak in product_candidates:
                if peak != product_peak:
                    peak.peak_type = PeakType.UNKNOWN
        else:
            compound.selected_peak = None

        # Step 4: Tentatively classify rejected peaks
        self._classify_rejected_peaks(
            compound, descendants, l0_retention_time,
            max_truncation_boundary, tolerance
        )

    def _classify_l0_compound(self, compound: Compound) -> None:
        """
        Classify L0 (null) compound - find the NULL peak directly from signal.

        L0 is special: the DNA tag is GUARANTEED to exist. We find the peak
        directly from the raw signal as the global maximum, bypassing all
        significance testing. L0 has no "unknown" peaks - only the NULL peak.

        Parameters
        ----------
        compound : Compound
            L0 compound to classify

        Notes
        -----
        L0 bypasses significance filtering because:
        1. The DNA tag always elutes (guaranteed to exist)
        2. We're not asking "is there a peak?" but "where is THE peak?"
        3. Global maximum is the L0 peak (largest signal)

        L0 does NOT have "unknown" peaks - any other signal is noise/artifacts.

        References
        ----------
        THEORY.md Section 2.4.4: L0 Peak Selection
        """
        chromatogram = compound.chromatogram
        time_points = chromatogram.time_points
        signal = chromatogram.get_signal("corrected")

        if len(signal) < 3:
            return

        # Find GLOBAL MAXIMUM - L0 peak is the largest peak
        peak_idx = int(np.argmax(signal))
        peak_time = float(time_points[peak_idx])
        peak_height = float(signal[peak_idx])

        # Find simple boundaries (extend to valleys or region edges)
        left_idx = peak_idx
        while left_idx > 0 and signal[left_idx - 1] < signal[left_idx]:
            left_idx -= 1

        right_idx = peak_idx
        while right_idx < len(signal) - 1 and signal[right_idx + 1] < signal[right_idx]:
            right_idx += 1

        left_base = time_points[left_idx]
        right_base = time_points[right_idx]

        # Compute area (simple sum)
        area = float(np.sum(signal[left_idx:right_idx + 1]))

        # Create the NULL peak
        null_peak = Peak(
            position=float(peak_time),
            left_base=float(left_base),
            right_base=float(right_base),
            height=float(peak_height),
            area=area,
            peak_type=PeakType.NULL,
            rejection_reason=RejectionReason.NONE,
        )

        # Set as selected peak and add to detected_peaks if not already there
        compound.selected_peak = null_peak
        if null_peak not in compound.detected_peaks:
            compound.detected_peaks.append(null_peak)

        # Clear detected_peaks of any other peaks - L0 only has the NULL peak
        compound.detected_peaks = [null_peak]

    def _classify_rejected_peaks(
        self,
        compound: Compound,
        descendants: List[Compound],
        l0_retention_time: Optional[float],
        max_truncation_boundary: float,
        tolerance: float
    ) -> None:
        """
        Tentatively classify rejected peaks based on position.

        Rejected peaks are classified to show what type they would have been,
        but they don't participate in optimal assignment or product selection.

        Parameters
        ----------
        compound : Compound
            Compound with rejected peaks to classify
        descendants : List[Compound]
            Descendant compounds for matching
        l0_retention_time : float or None
            L0 retention time
        max_truncation_boundary : float
            Truncation boundary (max descendant position + margin)
        tolerance : float
            Relative tolerance for position matching
        """
        # Get all descendant peaks for matching
        all_descendant_peaks = self._get_all_descendant_peaks(descendants)

        for peak in compound.detected_peaks:
            if not peak.is_rejected:
                continue

            # Check if matches L0 position
            if l0_retention_time is not None and l0_retention_time > 0:
                rel_diff = abs(peak.position - l0_retention_time) / l0_retention_time
                if rel_diff <= tolerance:
                    peak.peak_type = PeakType.NULL
                    continue

            # Check if matches any descendant product position
            matched_product = False
            for desc in descendants:
                if desc.selected_peak is not None:
                    pos = desc.selected_peak.position
                    if pos > 0:
                        rel_diff = abs(peak.position - pos) / pos
                        if rel_diff <= tolerance:
                            peak.peak_type = PeakType.TRUNCATION
                            peak.matched_compound_sequence = desc.block_support_sequence
                            peak.matched_peak_position = pos
                            peak.matched_peak_type = PeakType.PUTATIVE_PRODUCT
                            matched_product = True
                            break

            if matched_product:
                continue

            # Check if matches any descendant non-product peak
            match = self._find_matching_descendant_peak(peak, all_descendant_peaks, tolerance)
            if match:
                matched_compound, matched_peak = match
                peak.peak_type = PeakType.TRUNCATION_UNKNOWN

                # Trace to original source
                if matched_peak.matched_compound_sequence:
                    peak.matched_compound_sequence = matched_peak.matched_compound_sequence
                    peak.matched_peak_position = matched_peak.matched_peak_position
                    peak.matched_peak_type = matched_peak.matched_peak_type
                else:
                    peak.matched_compound_sequence = matched_compound.block_support_sequence
                    peak.matched_peak_position = matched_peak.position
                    peak.matched_peak_type = matched_peak.peak_type
                continue

            # Check if after truncation boundary (would be product candidate)
            if peak.position > max_truncation_boundary:
                peak.peak_type = PeakType.PUTATIVE_PRODUCT
            else:
                peak.peak_type = PeakType.UNKNOWN

    def _apply_clpe_at_level(
        self,
        compounds_at_level: List[Compound],
        clpe_validator: "CLPEValidator",
        alogp_map: Dict[str, float],
        scaffold_map: Dict[str, str],
    ) -> Dict[str, Any]:
        """
        Apply cLPE validation to all compounds at a single level.

        Groups compounds by scaffold, fits LogK ~ AlogP models per scaffold,
        validates product peaks against the model, and re-selects peaks if
        they are outliers.

        Called AFTER standard classification, BEFORE proceeding to next level.
        This ensures cLPE-corrected peaks are used for truncation matching
        in higher levels.

        Parameters
        ----------
        compounds_at_level : List[Compound]
            All compounds at the current level
        clpe_validator : CLPEValidator
            Validator with t0 already set from L0 peak
        alogp_map : Dict[str, float]
            Mapping from compound identifier to AlogP value
        scaffold_map : Dict[str, str]
            Mapping from compound identifier to scaffold group

        Returns
        -------
        Dict[str, Any]
            Statistics from cLPE validation at this level

        Notes
        -----
        Re-selection only considers UNKNOWN or PUTATIVE_PRODUCT peaks.
        TRUNCATION peaks are never re-selected (they're claimed by ancestors).
        """
        # Filter compounds with AlogP/scaffold data and selected peaks
        # Use local dicts instead of compound attributes (PooledCompound is immutable)
        compounds_with_data = []
        compound_alogp = {}  # Dict[Compound, float]
        compound_scaffold = {}  # Dict[Compound, str]

        for compound in compounds_at_level:
            if not compound.selected_peak:
                continue

            # Try to find AlogP and scaffold for this compound
            keys_to_try = [
                compound.positional_block_sequence,
                compound.block_support_sequence,
            ]
            if compound.compound_id:
                keys_to_try.insert(0, compound.compound_id)

            for key in keys_to_try:
                alogp = alogp_map.get(key)
                scaffold = scaffold_map.get(key)
                if alogp is not None and scaffold:
                    compound_alogp[compound] = alogp
                    compound_scaffold[compound] = scaffold
                    compounds_with_data.append(compound)
                    break

        if not compounds_with_data:
            return {
                "matched_compounds": 0,
                "validated": 0,
                "outliers": 0,
                "reselected": 0,
                "models_fitted": 0,
            }

        # Fit cLPE models for this level
        clpe_validator.fit_models(
            compounds_with_data,
            alogp_map,
            scaffold_map,
        )

        # Validate and optionally re-select peaks
        n_validated = 0
        n_outliers = 0
        n_reselected = 0

        for compound in compounds_with_data:
            alogp = compound_alogp.get(compound)
            scaffold = compound_scaffold.get(compound)
            if alogp is None or not scaffold:
                continue
            if not compound.selected_peak:
                continue

            result, new_peak = clpe_validator.validate_and_reselect(
                compound, alogp, scaffold
            )

            n_validated += 1

            # Store validation results on current peak
            if compound.selected_peak:
                compound.selected_peak.clpe_residual = result.residual
                compound.selected_peak.clpe_z_score = result.z_score
                compound.selected_peak.clpe_is_outlier = result.is_outlier

            if result.is_outlier:
                n_outliers += 1

            # Re-select peak if better alternative found
            if new_peak and new_peak != compound.selected_peak:
                old_peak = compound.selected_peak
                old_peak.peak_type = PeakType.UNKNOWN  # Demote old selection
                compound.selected_peak = new_peak
                new_peak.peak_type = PeakType.PUTATIVE_PRODUCT
                new_peak.clpe_reselected = True
                n_reselected += 1

        return {
            "matched_compounds": len(compounds_with_data),
            "validated": n_validated,
            "outliers": n_outliers,
            "reselected": n_reselected,
            "models_fitted": len(clpe_validator.models),
        }
