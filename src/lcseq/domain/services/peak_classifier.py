"""
Peak classification service using DAG constraint propagation.

Implementation based on THEORY.md Section 5.3-5.6.
"""

import numpy as np
from typing import List, Optional, Tuple, Dict
from scipy.optimize import linear_sum_assignment
from ..entities.compound import Compound
from ..entities.peak import Peak, PeakType
from ..models.compound_hierarchy import CompoundHierarchy
from ...config import (
    DEFAULT_PEAK_MATCHING_TOLERANCE,
    DEFAULT_TRUNCATION_MARGIN,
    DEFAULT_HUNGARIAN_MIN_THRESHOLD,
)


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
    >>> compound = Compound(building_blocks=[...], chromatogram=...)
    >>> hierarchy = CompoundHierarchy(mode=HierarchyMode.BUILDING_BLOCK)
    >>> # Detect peaks first
    >>> compound.detected_peaks = [...]
    >>> # Classify peaks
    >>> classified_peak = classifier.classify_peak(
    ...     compound, compound.detected_peaks[0], hierarchy
    ... )
    >>> classified_peak.peak_type
    <PeakType.PUTATIVE_PRODUCT: 'PUTATIVE_PRODUCT'>
    """

    def classify_peak(
        self,
        compound: Compound,
        peak: Peak,
        hierarchy: CompoundHierarchy,
        l0_retention_time: Optional[float] = None,
        tolerance: float = DEFAULT_PEAK_MATCHING_TOLERANCE
    ) -> Peak:
        """
        Classify a single peak for a compound.

        Parameters
        ----------
        compound : Compound
            Compound containing the peak
        peak : Peak
            Peak to classify
        hierarchy : CompoundHierarchy
            DAG structure with ancestor/descendant relationships
        l0_retention_time : float or None, optional
            Retention time of L₀ (null) peak. If None, search for L₀ in hierarchy.
        tolerance : float, optional
            Relative tolerance for retention time matching (fraction).
            Relative tolerance: |measured - expected| / expected

        Returns
        -------
        Peak
            Peak with updated peak_type classification

        Notes
        -----
        Classification priority:
        1. NULL: Matches L₀ position (±tolerance)
        2. TRUNCATION: Matches ancestor product or null position
        3. PUTATIVE_PRODUCT: First peak after truncations (positional hypothesis)
        4. UNKNOWN: No expected position match

        Tolerance is relative to position value: |measured - expected| / expected

        References
        ----------
        THEORY.md Section 5.3: Classification Logic
        THEORY.md Section 5.5: Optimal Assignment and Matching

        Examples
        --------
        >>> classifier = PeakClassifier()
        >>> peak = Peak(position=120.5, left_base=115, right_base=125, ...)
        >>> classified = classifier.classify_peak(compound, peak, hierarchy)
        """
        # Get L₀ retention time
        if l0_retention_time is None:
            l0_retention_time = self._find_l0_retention_time(hierarchy)

        # Check if NULL peak (matches L₀ position)
        if l0_retention_time is not None:
            if self._is_null_peak(peak, l0_retention_time, tolerance):
                peak.peak_type = PeakType.NULL
                return peak

        # Get truncation (descendant) positions for matching
        # Truncations are shorter sequences - descendants in the DAG
        truncations = hierarchy.get_descendants(compound)
        truncation_product_positions = self._get_product_positions(truncations)

        # Check if TRUNCATION peak (matches truncation product or null)
        if self._is_truncation_peak(peak, truncation_product_positions, l0_retention_time, tolerance):
            peak.peak_type = PeakType.TRUNCATION
            return peak

        # Check if PUTATIVE_PRODUCT (first peak after truncations)
        # This requires context of all peaks in compound
        if self._is_putative_product_peak(peak, compound, truncation_product_positions, l0_retention_time, tolerance):
            peak.peak_type = PeakType.PUTATIVE_PRODUCT
            return peak

        # Default: UNKNOWN
        peak.peak_type = PeakType.UNKNOWN
        return peak

    def classify_all_peaks(
        self,
        compound: Compound,
        hierarchy: CompoundHierarchy,
        l0_retention_time: Optional[float] = None,
        tolerance: float = DEFAULT_PEAK_MATCHING_TOLERANCE,
        truncation_margin: float = DEFAULT_TRUNCATION_MARGIN,
    ) -> None:
        """
        Classify all peaks for a compound in-place.

        Applies classification logic to all detected peaks and selects
        the putative product peak.

        Parameters
        ----------
        compound : Compound
            Compound with detected_peaks to classify
        hierarchy : CompoundHierarchy
            DAG structure
        l0_retention_time : float or None, optional
            L₀ retention time. If None, auto-detect.
        tolerance : float, optional
            Relative tolerance for position matching (fraction).
            Relative tolerance: |measured - expected| / expected
        truncation_margin : float, optional
            Margin beyond truncation positions (in seconds).
            Product peaks must be this far beyond max(truncation_positions).
            Accounts for retention time variability.

        Notes
        -----
        Updates compound.detected_peaks in-place with classifications.
        Sets compound.selected_peak to the PUTATIVE_PRODUCT peak if found.

        Special case for L₀ (THEORY.md Section 2.4.4):
        - L₀ has no truncations (minimal compound)
        - Selects earliest peak in early region (first 20% of chromatogram)
        - This peak is both NULL and the "product" for L₀

        Classification ensures:
        - At most one PUTATIVE_PRODUCT per compound
        - TRUNCATION peaks elute before PUTATIVE_PRODUCT
        - NULL peak identified if present
        - Remaining peaks marked UNKNOWN

        References
        ----------
        THEORY.md Section 2.4.4: L₀ Peak Detection Algorithm
        THEORY.md Section 5.3: Peak Type Classification
        THEORY.md Section 5.4: Global Classification
        """
        if not compound.detected_peaks:
            return

        # Special case: L₀ (all-null) compound (THEORY.md 2.4.4)
        # L₀ has no truncations - select earliest peak in early region
        # Check if compound is L₀ by checking if it's at level 0 in the hierarchy
        l0_compounds = hierarchy.get_level(0)
        is_l0 = compound in l0_compounds

        if is_l0:
            # Define early region (first 20% of chromatogram)
            time_points = compound.chromatogram.time_points
            t_min = time_points[0]
            t_max = time_points[-1]
            t_early_max = t_min + 0.2 * (t_max - t_min)

            # Filter peaks in early region
            early_peaks = [p for p in compound.detected_peaks if p.position <= t_early_max]

            if early_peaks:
                # Select earliest peak (minimum time) - THEORY.md 2.4.4
                earliest_peak = min(early_peaks, key=lambda p: p.position)
                earliest_peak.peak_type = PeakType.NULL
                compound.selected_peak = earliest_peak

                # Mark all other peaks as UNKNOWN
                for peak in compound.detected_peaks:
                    if peak != earliest_peak:
                        peak.peak_type = PeakType.UNKNOWN
            else:
                # Fallback: use earliest detected peak
                earliest_peak = min(compound.detected_peaks, key=lambda p: p.position)
                earliest_peak.peak_type = PeakType.NULL
                compound.selected_peak = earliest_peak

                # Mark all other peaks as UNKNOWN
                for peak in compound.detected_peaks:
                    if peak != earliest_peak:
                        peak.peak_type = PeakType.UNKNOWN

            return  # L₀ handling complete

        # Get L₀ retention time
        if l0_retention_time is None:
            l0_retention_time = self._find_l0_retention_time(hierarchy)

        # Get truncation (descendant) positions
        truncations = hierarchy.get_descendants(compound)
        truncation_product_positions = self._get_product_positions(truncations)

        # Build expected positions list (THEORY.md Section 5.5)
        # E = {L₀, truncation_products...}
        expected_positions = []
        expected_types = []  # Track which type each position represents

        if l0_retention_time is not None and l0_retention_time > 0:
            expected_positions.append(l0_retention_time)
            expected_types.append(PeakType.NULL)

        for trunc_pos in truncation_product_positions:
            if trunc_pos > 0:
                expected_positions.append(trunc_pos)
                expected_types.append(PeakType.TRUNCATION)

        # Optimal assignment using Hungarian algorithm (THEORY.md Section 5.5)
        if expected_positions:
            # Get signal length for normalization
            time_points = compound.chromatogram.time_points
            signal_length = time_points[-1] - time_points[0]

            # Find optimal assignment
            assignments = self._optimal_assignment(
                compound.detected_peaks,
                expected_positions,
                signal_length
            )

            # Classify peaks based on optimal assignment
            assigned_peak_indices = set()
            null_peak = None
            truncation_peaks = []

            for peak_idx, expected_idx in assignments.items():
                peak = compound.detected_peaks[peak_idx]
                peak_type = expected_types[expected_idx]

                peak.peak_type = peak_type
                assigned_peak_indices.add(peak_idx)

                if peak_type == PeakType.NULL:
                    null_peak = peak
                elif peak_type == PeakType.TRUNCATION:
                    truncation_peaks.append(peak)

            # CRITICAL: Validate TRUNCATION assignments against truncation boundary
            # Peaks assigned as TRUNCATION must actually elute BEFORE the boundary
            # This fixes cases where Hungarian algorithm matches peaks to distant expected positions
            max_truncation_pos = 0.0
            if l0_retention_time is not None:
                max_truncation_pos = l0_retention_time
            if truncation_product_positions:
                max_truncation_pos = max(max_truncation_pos, max(truncation_product_positions))
            max_truncation_boundary = max_truncation_pos + truncation_margin

            # Reclassify invalid TRUNCATION peaks
            invalid_truncations = []
            for peak in truncation_peaks:
                if peak.position > max_truncation_boundary:
                    # Peak elutes AFTER truncation boundary → cannot be truncation
                    # Mark as unassigned so it can be reconsidered for product
                    peak_idx = compound.detected_peaks.index(peak)
                    assigned_peak_indices.remove(peak_idx)
                    invalid_truncations.append(peak)

            # Remove invalid truncations from truncation_peaks list
            for peak in invalid_truncations:
                truncation_peaks.remove(peak)

            # Unassigned peaks could be PRODUCT or UNKNOWN
            unassigned_peaks = [
                compound.detected_peaks[i]
                for i in range(len(compound.detected_peaks))
                if i not in assigned_peak_indices
            ]
        else:
            # No expected positions - all peaks unassigned
            unassigned_peaks = list(compound.detected_peaks)
            null_peak = None
            truncation_peaks = []

        # Determine product candidates from unassigned peaks (THEORY.md 5.3)
        # Constraint: position > max(ALL possible descendant positions) + margin
        product_candidates = []

        # Find max truncation position (reuse calculation from above if available)
        if 'max_truncation_boundary' in locals():
            # Already computed during TRUNCATION validation
            max_truncation_pos_with_margin = max_truncation_boundary
        else:
            # First calculation (when no expected positions)
            max_truncation_pos = 0.0
            if l0_retention_time is not None:
                max_truncation_pos = l0_retention_time
            if truncation_product_positions:
                max_truncation_pos = max(max_truncation_pos, max(truncation_product_positions))

            # Apply truncation margin to account for retention time variability
            # Margin is absolute time in seconds
            max_truncation_pos_with_margin = max_truncation_pos + truncation_margin

        for peak in unassigned_peaks:
            # Must elute after all possible truncations (descendant products) + margin
            # Margin accounts for retention time variability from peak matching
            if peak.position > max_truncation_pos_with_margin:
                product_candidates.append(peak)
            else:
                # Peak elutes before expected product position (including margin)
                peak.peak_type = PeakType.UNKNOWN

        # Select PUTATIVE_PRODUCT: first by position among candidates (THEORY.md 5.3)
        product_peak = None
        if product_candidates:
            # Select first peak by position (earliest elution after truncations)
            product_peak = min(product_candidates, key=lambda p: p.position)
            product_peak.peak_type = PeakType.PUTATIVE_PRODUCT

            # Mark non-selected candidates as UNKNOWN
            for peak in product_candidates:
                if peak != product_peak:
                    peak.peak_type = PeakType.UNKNOWN

        # Set selected peak
        compound.selected_peak = product_peak

    def _optimal_assignment(
        self,
        detected_peaks: List[Peak],
        expected_positions: List[float],
        signal_length: float
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
        # Cost(i,j) = |peak_i.position - expected_j| / signal_length
        n_peaks = len(detected_peaks)
        n_expected = len(expected_positions)
        cost_matrix = np.zeros((n_peaks, n_expected))

        for i, peak in enumerate(detected_peaks):
            for j, expected_pos in enumerate(expected_positions):
                cost_matrix[i, j] = abs(peak.position - expected_pos) / signal_length

        # Solve optimal assignment using Hungarian algorithm
        row_ind, col_ind = linear_sum_assignment(cost_matrix)

        # Compute adaptive threshold (THEORY.md Section 5.5)
        # threshold = median(peak_spacings) / signal_length
        sorted_expected = sorted(expected_positions)
        peak_spacings = [sorted_expected[i+1] - sorted_expected[i]
                        for i in range(len(sorted_expected)-1)]

        if peak_spacings:
            median_spacing = np.median(peak_spacings)
            threshold = median_spacing / signal_length
        else:
            # Fallback: use adaptive threshold if only one expected position
            threshold = DEFAULT_PEAK_MATCHING_TOLERANCE

        # Minimum threshold to handle LC-Seq's discrete fractionation
        # Accounts for retention time variability in discrete fraction collection
        threshold = max(threshold, DEFAULT_HUNGARIAN_MIN_THRESHOLD)

        # Accept only assignments below or equal to threshold
        # Use <= to ensure perfect matches (cost=0.0) are accepted when threshold=0
        assignments = {}
        for i, j in zip(row_ind, col_ind):
            if cost_matrix[i, j] <= threshold:
                assignments[i] = j

        return assignments

    def _is_null_peak(
        self,
        peak: Peak,
        l0_retention_time: float,
        tolerance: float
    ) -> bool:
        """
        Check if peak matches L₀ (null) position.

        Parameters
        ----------
        peak : Peak
            Peak to check
        l0_retention_time : float
            L₀ retention time
        tolerance : float
            Relative tolerance (fraction)

        Returns
        -------
        bool
            True if peak matches L₀ position within tolerance
        """
        if l0_retention_time <= 0:
            return False

        relative_diff = abs(peak.position - l0_retention_time) / l0_retention_time
        return relative_diff <= tolerance

    def _is_truncation_peak(
        self,
        peak: Peak,
        truncation_positions: List[float],
        l0_retention_time: Optional[float],
        tolerance: float
    ) -> bool:
        """
        Check if peak matches any truncation product position or null position.

        Parameters
        ----------
        peak : Peak
            Peak to check
        truncation_positions : List[float]
            Product positions from truncations (descendants)
        l0_retention_time : float or None
            L₀ retention time
        tolerance : float
            Relative tolerance

        Returns
        -------
        bool
            True if peak matches any expected truncation position
        """
        # Check against truncation product positions
        for pos in truncation_positions:
            if pos <= 0:
                continue
            relative_diff = abs(peak.position - pos) / pos
            if relative_diff <= tolerance:
                return True

        # Check against null position (if not already classified as NULL)
        if l0_retention_time is not None and l0_retention_time > 0:
            relative_diff = abs(peak.position - l0_retention_time) / l0_retention_time
            if relative_diff <= tolerance:
                return True

        return False

    def _is_putative_product_peak(
        self,
        peak: Peak,
        compound: Compound,
        truncation_positions: List[float],
        l0_retention_time: Optional[float],
        tolerance: float
    ) -> bool:
        """
        Check if peak is putative product (first peak after truncations).

        Parameters
        ----------
        peak : Peak
            Peak to check
        compound : Compound
            Compound containing peak
        truncation_positions : List[float]
            Truncation product positions (from descendants)
        l0_retention_time : float or None
            L₀ retention time
        tolerance : float
            Relative tolerance

        Returns
        -------
        bool
            True if peak is putative product candidate

        Notes
        -----
        Putative product criteria:
        - Not NULL or TRUNCATION
        - Elutes after all truncation positions
        - First such peak (if multiple candidates)
        """
        # Must not match null or truncation positions
        if l0_retention_time is not None and self._is_null_peak(peak, l0_retention_time, tolerance):
            return False

        if self._is_truncation_peak(peak, truncation_positions, l0_retention_time, tolerance):
            return False

        # Must elute after truncations
        # Find maximum truncation position
        max_truncation_pos = 0.0
        if l0_retention_time is not None:
            max_truncation_pos = l0_retention_time
        if truncation_positions:
            max_truncation_pos = max(max_truncation_pos, max(truncation_positions))

        # Product should elute after truncations
        if peak.position <= max_truncation_pos:
            return False

        return True

    def _find_l0_retention_time(self, hierarchy: CompoundHierarchy) -> Optional[float]:
        """
        Find L₀ (all-null) compound retention time from hierarchy.

        Uses earliest significant peak in early region (first 20% of chromatogram)
        per THEORY.md Section 2.4.4.

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
        1. Define early region: first 20% of chromatogram
        2. Find peaks in early region
        3. Select earliest peak (minimum retention time)

        Rationale: Earlier peak = more hydrophilic → more likely pure DNA tag

        References
        ----------
        THEORY.md Section 2.4.4: L₀ Peak Detection Algorithm
        """
        # Find L₀ compound (level 0)
        l0_compounds = hierarchy.get_level(0)

        if not l0_compounds:
            return None

        # Should be exactly one L₀
        l0 = l0_compounds[0]

        # Get selected peak (putative product for L₀ is the null peak)
        if l0.selected_peak is not None:
            return l0.selected_peak.position

        # Fallback: find earliest significant peak in early region (THEORY.md 2.4.4)
        if l0.detected_peaks:
            # Define early region (first 20% of chromatogram)
            time_points = l0.chromatogram.time_points
            t_min = time_points[0]
            t_max = time_points[-1]
            t_early_max = t_min + 0.2 * (t_max - t_min)

            # Filter peaks in early region
            early_peaks = [p for p in l0.detected_peaks if p.position <= t_early_max]

            if early_peaks:
                # Select earliest peak (minimum time) - THEORY.md 2.4.4
                earliest_peak = min(early_peaks, key=lambda p: p.position)
                return earliest_peak.position

            # Fallback: if no peaks in early region, use earliest detected peak
            earliest_peak = min(l0.detected_peaks, key=lambda p: p.position)
            return earliest_peak.position

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
