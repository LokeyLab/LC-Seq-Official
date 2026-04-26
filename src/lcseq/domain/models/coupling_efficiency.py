"""
Coupling efficiency models for synthetic efficiency analysis.

Implementation based on docs/SYNTHETIC_EFFICIENCY.md.
"""

from dataclasses import dataclass, field
from typing import Dict, List, Optional
import math


@dataclass
class CycleEfficiency:
    """
    Efficiency measurement for a single synthesis cycle.

    Represents the efficiency of coupling one building block after another
    during synthesis. Computed as the ratio of material that passed a cycle
    to material that reached that cycle.

    Attributes
    ----------
    cycle : int
        Synthesis cycle number (0-indexed)
    prev_block : str
        Previous block code ("N" for null/start, or building block code)
    next_block : str
        Block code that was coupled at this cycle
    efficiency : float
        Coupling efficiency (area_passed / area_reached), in [0, 1]
    area_reached : float
        Denominator: total area of material that reached this cycle
    area_passed : float
        Numerator: total area of material that successfully passed this cycle

    Notes
    -----
    - efficiency = area_passed / area_reached
    - For cycle 0: prev_block = "N" (null/starting material)
    - efficiency of NaN indicates no material reached this cycle (divide by zero)

    Examples
    --------
    >>> ce = CycleEfficiency(
    ...     cycle=1,
    ...     prev_block="Leu",
    ...     next_block="Pro",
    ...     efficiency=0.92,
    ...     area_reached=100000.0,
    ...     area_passed=92000.0
    ... )
    >>> ce.efficiency
    0.92

    References
    ----------
    docs/SYNTHETIC_EFFICIENCY.md: Per-Cycle Efficiency
    """

    cycle: int
    prev_block: str
    next_block: str
    efficiency: float
    area_reached: float
    area_passed: float

    def __post_init__(self) -> None:
        """Validate cycle efficiency properties."""
        if self.cycle < 0:
            raise ValueError(f"Cycle must be non-negative, got {self.cycle}")
        if not self.prev_block:
            raise ValueError("prev_block cannot be empty")
        if not self.next_block:
            raise ValueError("next_block cannot be empty")

    @property
    def is_valid(self) -> bool:
        """
        Check if efficiency is a valid number (not NaN).

        Returns
        -------
        bool
            True if efficiency is a valid number
        """
        return not math.isnan(self.efficiency)

    @property
    def transition(self) -> str:
        """
        Human-readable transition string.

        Returns
        -------
        str
            Transition in format "prev -> next"
        """
        return f"{self.prev_block} -> {self.next_block}"

    def __repr__(self) -> str:
        """Detailed representation for debugging."""
        return (
            f"CycleEfficiency(cycle={self.cycle}, "
            f"transition='{self.transition}', "
            f"efficiency={self.efficiency:.3f})"
        )


@dataclass
class CompoundEfficiency:
    """
    Complete efficiency analysis for a single compound.

    Contains per-cycle efficiency measurements and supporting data
    for computing overall synthetic efficiency.

    Attributes
    ----------
    compound_id : str
        Unique identifier for the compound
    sequence : str
        Positional block sequence (e.g., "Leu-Pro-Ala")
    level : int
        Compound level (number of non-null building blocks)
    cycle_efficiencies : List[CycleEfficiency]
        Per-cycle efficiency measurements
    areas_by_level : Dict[int, float]
        Peak areas grouped by truncation level (for debugging/verification)

    Notes
    -----
    - overall_efficiency = product of all cycle efficiencies
    - areas_by_level[0] = L0 (null) peak area
    - areas_by_level[level] = product peak area
    - Missing levels in areas_by_level have area 0 (no truncation detected)

    Examples
    --------
    >>> ce = CompoundEfficiency(
    ...     compound_id="cpd_001",
    ...     sequence="Leu-Pro-Ala",
    ...     level=3,
    ...     cycle_efficiencies=[...],
    ...     areas_by_level={0: 1000, 1: 500, 2: 200, 3: 8000}
    ... )
    >>> ce.overall_efficiency
    0.823

    References
    ----------
    docs/SYNTHETIC_EFFICIENCY.md: Per-Cycle Efficiency
    """

    compound_id: str
    sequence: str
    level: int
    cycle_efficiencies: List[CycleEfficiency] = field(default_factory=list)
    areas_by_level: Dict[int, float] = field(default_factory=dict)

    def __post_init__(self) -> None:
        """Validate compound efficiency properties."""
        if self.level < 0:
            raise ValueError(f"Level must be non-negative, got {self.level}")

    @property
    def overall_efficiency(self) -> float:
        """
        Compute overall synthetic efficiency as product of cycle efficiencies.

        Returns
        -------
        float
            Product of all cycle efficiencies, in [0, 1]
            Returns NaN if any cycle efficiency is NaN

        Notes
        -----
        overall_efficiency = efficiency[0] * efficiency[1] * ... * efficiency[n-1]

        This equals: area_product / area_total (when all cycles have valid data)
        """
        if not self.cycle_efficiencies:
            return float('nan')

        result = 1.0
        for ce in self.cycle_efficiencies:
            if not ce.is_valid:
                return float('nan')
            result *= ce.efficiency
        return result

    @property
    def synthetic_efficiency(self) -> float:
        """
        Compute synthetic efficiency directly from areas.

        Returns
        -------
        float
            product_area / (product_area + truncation_areas)
            Excludes unknown peaks, only considers synthesis-related signal.

        Notes
        -----
        This should equal overall_efficiency when all cycle efficiencies are valid.
        Useful for verification and when some cycle data is missing.
        """
        if not self.areas_by_level:
            return float('nan')

        product_area = self.areas_by_level.get(self.level, 0.0)
        # Truncation areas = all levels except product (level 0 to level-1)
        truncation_area = sum(
            area for lvl, area in self.areas_by_level.items()
            if lvl < self.level
        )
        total = product_area + truncation_area

        if total == 0:
            return float('nan')

        return product_area / total

    @property
    def total_area(self) -> float:
        """
        Get total area across all levels.

        Returns
        -------
        float
            Sum of all peak areas
        """
        return sum(self.areas_by_level.values())

    @property
    def product_area(self) -> float:
        """
        Get product peak area.

        Returns
        -------
        float
            Area at compound's level (product), or 0 if not found
        """
        return self.areas_by_level.get(self.level, 0.0)

    @property
    def null_area(self) -> float:
        """
        Get L0 (null) peak area.

        Returns
        -------
        float
            Area at level 0, or 0 if not found
        """
        return self.areas_by_level.get(0, 0.0)

    def get_efficiency_for_cycle(self, cycle: int) -> Optional[CycleEfficiency]:
        """
        Get efficiency measurement for a specific cycle.

        Parameters
        ----------
        cycle : int
            Cycle number to query

        Returns
        -------
        Optional[CycleEfficiency]
            Efficiency for that cycle, or None if not found
        """
        for ce in self.cycle_efficiencies:
            if ce.cycle == cycle:
                return ce
        return None

    def __repr__(self) -> str:
        """Detailed representation for debugging."""
        eff = self.overall_efficiency
        eff_str = f"{eff:.3f}" if not math.isnan(eff) else "NaN"
        return (
            f"CompoundEfficiency("
            f"id='{self.compound_id}', "
            f"sequence='{self.sequence}', "
            f"level={self.level}, "
            f"overall_efficiency={eff_str}, "
            f"n_cycles={len(self.cycle_efficiencies)})"
        )
