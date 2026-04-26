"""
Coupling efficiency calculation service.

Computes per-cycle coupling efficiency from classified peak data.

Implementation based on docs/SYNTHETIC_EFFICIENCY.md.
"""

from typing import Dict, List, Tuple, Optional
import math

from ..entities.compound import Compound
from ..entities.peak import Peak, PeakType
from ..models.compound_hierarchy import CompoundHierarchy
from ..models.coupling_efficiency import CycleEfficiency, CompoundEfficiency


class CouplingEfficiencyCalculator:
    """
    Calculates coupling efficiency from classified peak data.

    For a compound with truncation peaks, computes efficiency of each
    synthesis cycle based on area partitioning by truncation level.

    This is a stateless service - all methods are operations on input
    data with no instance state.

    Notes
    -----
    Efficiency calculation:
    - efficiency[k] = area(passed cycle k) / area(reached cycle k)
    - For cycle 0: efficiency = (total - null_area) / total
    - For cycle k: efficiency = sum(areas at levels > k) / sum(areas at levels >= k)

    Requires:
    - Compound with classified peaks (peak_type set)
    - Hierarchy to identify descendant levels

    References
    ----------
    docs/SYNTHETIC_EFFICIENCY.md: Per-Cycle Efficiency
    THEORY.md Section 5.3: Peak Type Classification

    Examples
    --------
    >>> calculator = CouplingEfficiencyCalculator()
    >>> result = calculator.calculate(compound, hierarchy)
    >>> result.overall_efficiency
    0.823
    >>> result.cycle_efficiencies[0].efficiency
    0.95
    """

    def calculate(
        self,
        compound: Compound,
        hierarchy: CompoundHierarchy,
        peak_matching_tolerance: float,
    ) -> CompoundEfficiency:
        """
        Calculate per-cycle efficiency for a compound.

        Main entry point for efficiency calculation.

        Parameters
        ----------
        compound : Compound
            Compound with classified peaks (peak_type and area set)
        hierarchy : CompoundHierarchy
            Hierarchy containing this compound and its descendants

        Returns
        -------
        CompoundEfficiency
            Complete efficiency analysis for the compound

        Notes
        -----
        Requires peaks to be classified with:
        - PeakType.NULL for L0 peak
        - PeakType.TRUNCATION for truncation peaks
        - PeakType.PUTATIVE_PRODUCT for product peak

        Examples
        --------
        >>> result = calculator.calculate(compound, hierarchy)
        >>> print(f"Overall efficiency: {result.overall_efficiency:.1%}")
        Overall efficiency: 82.3%
        """
        # Step 1: Get areas by level
        areas_by_level = self.get_areas_by_level(compound, hierarchy, peak_matching_tolerance)

        # Step 2: Extract transitions from building blocks
        transitions = self.extract_transitions(compound)

        # Step 3: Compute per-cycle efficiency
        cycle_efficiencies = []
        for prev_block, next_block, cycle in transitions:
            efficiency, area_reached, area_passed = self._compute_cycle_efficiency(
                cycle, compound.level, areas_by_level
            )
            cycle_efficiencies.append(CycleEfficiency(
                cycle=cycle,
                prev_block=prev_block,
                next_block=next_block,
                efficiency=efficiency,
                area_reached=area_reached,
                area_passed=area_passed
            ))

        return CompoundEfficiency(
            compound_id=compound.compound_id or str(compound),
            sequence=compound.positional_block_sequence,
            level=compound.level,
            cycle_efficiencies=cycle_efficiencies,
            areas_by_level=areas_by_level
        )

    def get_areas_by_level(
        self,
        compound: Compound,
        hierarchy: CompoundHierarchy,
        peak_matching_tolerance: float,
    ) -> Dict[int, float]:
        """
        Map classified peaks to their truncation levels.

        Parameters
        ----------
        compound : Compound
            Compound with classified peaks
        hierarchy : CompoundHierarchy
            Hierarchy for descendant lookup

        Returns
        -------
        Dict[int, float]
            Mapping from level to total peak area at that level:
            - level 0 = NULL peak area (L0)
            - level k = truncation at level k
            - level n = PUTATIVE_PRODUCT area (compound's own level)

        Notes
        -----
        For TRUNCATION peaks, we match to descendants by retention time
        to determine which level each truncation represents.
        """
        areas_by_level: Dict[int, float] = {}

        # Get descendants for truncation matching
        descendants = hierarchy.get_descendants(compound)
        descendant_positions = self._get_descendant_position_map(descendants)

        for peak in compound.detected_peaks:
            if peak.peak_type == PeakType.NULL:
                # NULL peak = level 0
                level = 0
            elif peak.peak_type == PeakType.PUTATIVE_PRODUCT:
                # Product peak = compound's own level
                level = compound.level
            elif peak.peak_type == PeakType.TRUNCATION:
                # Match truncation to descendant by retention time
                level = self._match_truncation_to_level(
                    peak, descendant_positions, peak_matching_tolerance
                )
            else:
                # UNKNOWN peaks are excluded from efficiency calculation
                continue

            # Accumulate area at this level
            if level not in areas_by_level:
                areas_by_level[level] = 0.0
            areas_by_level[level] += peak.area

        return areas_by_level

    def extract_transitions(
        self,
        compound: Compound
    ) -> List[Tuple[str, str, int]]:
        """
        Extract (prev_block, next_block, cycle) transitions from building blocks.

        Parameters
        ----------
        compound : Compound
            Compound with building blocks

        Returns
        -------
        List[Tuple[str, str, int]]
            List of (prev_block_code, next_block_code, cycle) tuples.
            prev_block is "N" for the first transition (null → first block).

        Notes
        -----
        Building blocks are ordered by cycle (0, 1, 2, ...).
        Null blocks are skipped in the transition chain.

        Examples
        --------
        For compound Leu-Pro-Ala (cycles 0, 1, 2):
        - ("N", "Ala", 0)      # null → first block
        - ("Ala", "Pro", 1)    # first → second
        - ("Pro", "Leu", 2)    # second → third
        """
        transitions = []
        blocks = compound.building_blocks  # Ordered by cycle

        # Track previous non-null block for transitions
        prev_block_code = "N"  # Start with null
        prev_cycle = -1  # Before first cycle

        for block in blocks:
            if block.is_null:
                # Skip null blocks in transition chain
                continue

            # Record transition from previous to current
            # The cycle is the current block's cycle (where coupling happened)
            transitions.append((prev_block_code, block.code, block.cycle))

            # Update previous for next iteration
            prev_block_code = block.code
            prev_cycle = block.cycle

        return transitions

    def _compute_cycle_efficiency(
        self,
        cycle: int,
        compound_level: int,
        areas_by_level: Dict[int, float]
    ) -> Tuple[float, float, float]:
        """
        Compute efficiency for a specific cycle.

        Parameters
        ----------
        cycle : int
            Cycle number (0-indexed)
        compound_level : int
            Level of the compound being analyzed
        areas_by_level : Dict[int, float]
            Peak areas by truncation level

        Returns
        -------
        Tuple[float, float, float]
            (efficiency, area_reached, area_passed)

        Notes
        -----
        efficiency[k] = area(passed cycle k) / area(reached cycle k)

        For cycle k:
        - area_reached = sum of areas at levels >= k
        - area_passed = sum of areas at levels > k

        The level corresponds to how many blocks a compound has:
        - Level 0 = no blocks (L0/null)
        - Level 1 = 1 block (truncation at cycle 0 means failed before adding block 1)
        - Level k = k blocks

        So truncation at cycle k produces a compound at level k (has blocks 0..k-1).
        Material that "reached" cycle k has at least k blocks (level >= k).
        Material that "passed" cycle k has at least k+1 blocks (level > k).
        """
        # area_reached = sum of areas at levels >= cycle
        # (material that had at least 'cycle' blocks before this coupling)
        area_reached = sum(
            area for lvl, area in areas_by_level.items()
            if lvl >= cycle
        )

        # area_passed = sum of areas at levels > cycle
        # (material that successfully added the block at this cycle)
        area_passed = sum(
            area for lvl, area in areas_by_level.items()
            if lvl > cycle
        )

        if area_reached == 0:
            return (float('nan'), area_reached, area_passed)

        efficiency = area_passed / area_reached
        return (efficiency, area_reached, area_passed)

    def _get_descendant_position_map(
        self,
        descendants: List[Compound]
    ) -> Dict[float, int]:
        """
        Build mapping from retention time to level for descendants.

        Parameters
        ----------
        descendants : List[Compound]
            Descendant compounds (truncations)

        Returns
        -------
        Dict[float, int]
            Mapping from selected_peak.position to compound level
        """
        position_to_level = {}
        for desc in descendants:
            if desc.selected_peak is not None:
                position_to_level[desc.selected_peak.position] = desc.level
        return position_to_level

    def _match_truncation_to_level(
        self,
        peak: Peak,
        descendant_positions: Dict[float, int],
        tolerance: float,
    ) -> int:
        """
        Match a truncation peak to its corresponding level.

        Parameters
        ----------
        peak : Peak
            Truncation peak to match
        descendant_positions : Dict[float, int]
            Mapping from descendant retention times to levels
        tolerance : float
            Relative tolerance for position matching

        Returns
        -------
        int
            Matched level

        Raises
        ------
        ValueError
            If no descendant positions provided or no match within tolerance

        Notes
        -----
        Uses closest position match within tolerance.
        """
        if not descendant_positions:
            raise ValueError(
                f"Cannot match truncation peak at {peak.position}: no descendant positions"
            )

        # Find closest descendant position
        best_match_level = None
        best_distance = float('inf')

        for desc_pos, level in descendant_positions.items():
            distance = abs(peak.position - desc_pos)
            rel_distance = distance / max(desc_pos, 1.0)

            if rel_distance < best_distance:
                best_distance = rel_distance
                best_match_level = level

        # Accept match if within tolerance
        if best_match_level is not None and best_distance <= tolerance:
            return best_match_level

        raise ValueError(
            f"Cannot match truncation peak at {peak.position}: "
            f"closest match at distance {best_distance:.4f} exceeds tolerance {tolerance}"
        )
