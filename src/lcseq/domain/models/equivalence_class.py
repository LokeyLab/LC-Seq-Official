"""
EquivalenceClass - groups positional variants by chemical identity.

Implementation based on THEORY.md Section 4.2.
"""

from dataclasses import dataclass, field
from typing import List, Set, Optional, Dict, Any
from enum import Enum
from ..entities.compound import Compound
from ..entities.peak import Peak
from ..entities.chromatogram import Chromatogram


class PoolingStatus(Enum):
    """
    Status of pooled mode processing for an equivalence class.

    References
    ----------
    THEORY.md Section 4.2.8: Validity Requirements
    """

    NOT_ATTEMPTED = "not_attempted"
    POOLING_VALID = "pooling_valid"
    POOLING_INVALID = "pooling_invalid"
    HETEROGENEOUS = "heterogeneous"
    POOLING_INVALID_BUT_SIMILAR = "pooling_invalid_but_similar"


@dataclass
class EquivalenceClass:
    """
    Groups positional variants representing the same chemical peptide.

    An equivalence class under the relation R: "has same block support sequence".
    Multiple positional block sequences (different synthesis paths) that produce
    the same chemical molecule at block granularity.

    Supports both individual mode (default) and optional pooled mode
    for performance optimization.

    Attributes
    ----------
    block_support_sequence : str
        Equivalence class identifier (non-null building blocks only)
    members : Set[Compound]
        All positional variants (members of this equivalence class)
    pooled_chromatogram : Optional[Chromatogram]
        Aggregated pooled chromatogram (if pooled mode used)
    pooled_peaks : List[Peak]
        Peaks detected on pooled signal
    pooling_status : PoolingStatus
        Status of pooled processing
    correlation_min : Optional[float]
        Minimum pairwise correlation between variants
    fallback_reason : Optional[str]
        Reason for fallback to individual mode (if applicable)
    metadata : Dict[str, Any]
        Additional metadata for this class

    Notes
    -----
    - Equivalence relation properties (THEORY.md 4.2.1):
      - Reflexive: C R C
      - Symmetric: A R B ⟹ B R A
      - Transitive: A R B ∧ B R C ⟹ A R C
    - Same block support sequence = same chemical molecule at block granularity
    - Different positional block sequences = different synthesis paths
    - Pooled mode is optional optimization (THEORY.md 4.2.2)
    - Automatic fallback to individual mode if correlation < threshold

    Examples
    --------
    >>> from lcseq.domain.entities.building_block import BuildingBlock
    >>> from lcseq.domain.entities.chromatogram import Chromatogram
    >>> import numpy as np
    >>>
    >>> # Create equivalence class for "Val"
    >>> eq_class = EquivalenceClass(block_support_sequence="Val")
    >>>
    >>> # Create positional variants
    >>> chromatogram = Chromatogram(
    ...     time_points=np.array([1.0, 2.0, 3.0]),
    ...     counts=np.array([100.0, 200.0, 150.0])
    ... )
    >>>
    >>> # Variant 1: [Val, Null, Null]
    >>> bb0 = BuildingBlock.from_code(0, "Val")
    >>> bb1 = BuildingBlock.from_code(1, "Null")
    >>> bb2 = BuildingBlock.from_code(2, "Null")
    >>> variant1 = Compound([bb0, bb1, bb2], chromatogram)
    >>> eq_class.add_compound(variant1)
    >>>
    >>> # Variant 2: [Null, Val, Null]
    >>> bb0_null = BuildingBlock.from_code(0, "Null")
    >>> bb1_val = BuildingBlock.from_code(1, "Val")
    >>> variant2 = Compound([bb0_null, bb1_val, bb2], chromatogram)
    >>> eq_class.add_compound(variant2)
    >>>
    >>> # Both have same block_support_sequence = "Val"
    >>> len(eq_class.members)
    2

    References
    ----------
    THEORY.md Section 4.2.1: EquivalenceClass Definition
    THEORY.md Section 1.2: Chemical Identity vs Positional Identity
    THEORY.md Section 4.2.2: Pooled Mode
    THEORY.md Section 4.2.3: Hybrid Pooled Strategy
    THEORY.md Section 4.2.8: Validity Requirements
    """

    block_support_sequence: str
    members: Set[Compound] = field(default_factory=set)
    pooled_chromatogram: Optional[Chromatogram] = None
    pooled_peaks: List[Peak] = field(default_factory=list)
    pooling_status: PoolingStatus = PoolingStatus.NOT_ATTEMPTED
    correlation_min: Optional[float] = None
    fallback_reason: Optional[str] = None
    metadata: Dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        """Validate equivalence class properties."""
        # Validate all members have same block support sequence
        for member in self.members:
            if member.block_support_sequence != self.block_support_sequence:
                raise ValueError(
                    f"Compound block support sequence '{member.block_support_sequence}' "
                    f"does not match class block support sequence '{self.block_support_sequence}'"
                )

    def add_compound(self, compound: Compound) -> None:
        """
        Add a positional variant (member) to the equivalence class.

        Parameters
        ----------
        compound : Compound
            Positional variant to add as member

        Raises
        ------
        ValueError
            If compound's block support sequence doesn't match class

        Notes
        -----
        - Idempotent: adding same compound multiple times is safe
        - Validates block support sequence matches class identifier
        """
        if compound.block_support_sequence != self.block_support_sequence:
            raise ValueError(
                f"Cannot add compound with block support sequence "
                f"'{compound.block_support_sequence}' to class '{self.block_support_sequence}'"
            )
        self.members.add(compound)

    def size(self) -> int:
        """
        Get number of members in equivalence class.

        Returns
        -------
        int
            Number of members in equivalence class

        Examples
        --------
        >>> eq_class = EquivalenceClass(block_support_sequence="Val")
        >>> # Add 3 positional variants...
        >>> eq_class.size()
        3
        """
        return len(self.members)

    def is_empty(self) -> bool:
        """
        Check if equivalence class is empty.

        Returns
        -------
        bool
            True if no members in class
        """
        return len(self.members) == 0

    def get_positional_block_sequences(self) -> List[str]:
        """
        Get all distinct positional block sequences in class.

        Returns
        -------
        List[str]
            All positional block sequences

        Examples
        --------
        >>> eq_class = EquivalenceClass(block_support_sequence="Leu-Pro")
        >>> # Add variants: [Leu-Pro-Null], [Leu-Null-Pro], [Null-Leu-Pro]
        >>> eq_class.get_positional_block_sequences()
        ['Leu-Pro-Null', 'Leu-Null-Pro', 'Null-Leu-Pro']
        """
        return [c.positional_block_sequence for c in self.members]

    @property
    def is_pooling_valid(self) -> bool:
        """
        Check if pooled mode was successfully applied.

        Returns
        -------
        bool
            True if pooled mode succeeded

        References
        ----------
        THEORY.md Section 4.2.8: Validity Requirements
        """
        return self.pooling_status == PoolingStatus.POOLING_VALID

    @property
    def mean_purity(self) -> Optional[float]:
        """
        Mean purity across all members.

        Returns None if no members have computed purity.

        Returns
        -------
        Optional[float]
            Mean purity or None

        References
        ----------
        THEORY.md Section 4.2.7: Aggregate Statistics
        """
        purities = []
        for member in self.members:
            if member.selected_peak is not None:
                purity = member.selected_peak.metadata.get("purity")
                if purity is not None:
                    purities.append(purity)

        if not purities:
            return None
        return sum(purities) / len(purities)

    @property
    def purity_std(self) -> Optional[float]:
        """
        Standard deviation of purity across members.

        Returns None if fewer than 2 members have computed purity.

        Returns
        -------
        Optional[float]
            Purity standard deviation or None

        References
        ----------
        THEORY.md Section 4.2.7: Aggregate Statistics
        """
        purities = []
        for member in self.members:
            if member.selected_peak is not None:
                purity = member.selected_peak.metadata.get("purity")
                if purity is not None:
                    purities.append(purity)

        if len(purities) < 2:
            return None

        mean = sum(purities) / len(purities)
        variance = sum((p - mean) ** 2 for p in purities) / len(purities)
        return variance ** 0.5

    def __len__(self) -> int:
        """Return number of members (for len() operator)."""
        return len(self.members)

    def __repr__(self) -> str:
        """Detailed representation for debugging."""
        return (
            f"EquivalenceClass(block_support_sequence='{self.block_support_sequence}', "
            f"members={len(self.members)})"
        )

    def __str__(self) -> str:
        """String representation using block support sequence."""
        return f"[{self.block_support_sequence}] ({len(self.members)} members)"
