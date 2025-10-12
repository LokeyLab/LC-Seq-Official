"""
PeakBoundaries value object - immutable peak integration region.

Implementation based on THEORY.md Section 5.2.4.
"""

from dataclasses import dataclass
from typing import Optional


@dataclass(frozen=True)
class PeakBoundaries:
    """
    Represents the integration boundaries for a chromatographic peak.

    Peak boundaries define the region over which peak area is integrated.
    Boundaries are determined by valley detection or intensity thresholds
    (5% of peak height) as described in THEORY.md Section 5.2.4.

    Attributes
    ----------
    left_base : float
        Left boundary time point (start of integration)
    right_base : float
        Right boundary time point (end of integration)
    left_valley : Optional[float]
        Time of left valley (local minimum), if detected
    right_valley : Optional[float]
        Time of right valley (local minimum), if detected

    Notes
    -----
    - Immutable value object (frozen dataclass)
    - left_base < right_base (validated)
    - Valleys are optional (may be threshold-based instead)
    - Width = right_base - left_base
    - Boundaries used for peak area integration

    Examples
    --------
    >>> # Peak with detected valleys
    >>> bounds = PeakBoundaries(
    ...     left_base=44.0,
    ...     right_base=46.0,
    ...     left_valley=44.2,
    ...     right_valley=45.8
    ... )
    >>> bounds.width()
    2.0

    >>> # Peak with threshold-based boundaries (no valleys)
    >>> bounds = PeakBoundaries(
    ...     left_base=44.0,
    ...     right_base=46.0,
    ...     left_valley=None,
    ...     right_valley=None
    ... )
    >>> bounds.has_valleys()
    False

    References
    ----------
    THEORY.md Section 5.2.4: Peak Boundary Determination
    """

    left_base: float
    right_base: float
    left_valley: Optional[float] = None
    right_valley: Optional[float] = None

    def __post_init__(self) -> None:
        """Validate peak boundary properties."""
        # Validate left < right
        if self.left_base >= self.right_base:
            raise ValueError(
                f"left_base ({self.left_base}) must be < right_base ({self.right_base})"
            )

        # Validate valleys are within boundaries if present
        if self.left_valley is not None:
            if self.left_valley < self.left_base or self.left_valley > self.right_base:
                raise ValueError(
                    f"left_valley ({self.left_valley}) must be within "
                    f"[{self.left_base}, {self.right_base}]"
                )

        if self.right_valley is not None:
            if self.right_valley < self.left_base or self.right_valley > self.right_base:
                raise ValueError(
                    f"right_valley ({self.right_valley}) must be within "
                    f"[{self.left_base}, {self.right_base}]"
                )

        # Validate valley ordering if both present
        if self.left_valley is not None and self.right_valley is not None:
            if self.left_valley >= self.right_valley:
                raise ValueError(
                    f"left_valley ({self.left_valley}) must be < "
                    f"right_valley ({self.right_valley})"
                )

    @classmethod
    def from_valleys(
        cls,
        left_valley: float,
        right_valley: float,
    ) -> "PeakBoundaries":
        """
        Create boundaries from detected valleys.

        Uses valleys as both the base boundaries and valley positions.

        Parameters
        ----------
        left_valley : float
            Left valley position (becomes left_base)
        right_valley : float
            Right valley position (becomes right_base)

        Returns
        -------
        PeakBoundaries
            New boundaries with valleys at base positions

        Examples
        --------
        >>> bounds = PeakBoundaries.from_valleys(44.2, 45.8)
        >>> bounds.left_base
        44.2
        >>> bounds.has_valleys()
        True
        """
        return cls(
            left_base=left_valley,
            right_base=right_valley,
            left_valley=left_valley,
            right_valley=right_valley,
        )

    @classmethod
    def from_threshold(
        cls,
        left_base: float,
        right_base: float,
    ) -> "PeakBoundaries":
        """
        Create boundaries from threshold-based detection.

        Used when valleys are not detected and boundaries are determined
        by intensity threshold (5% of peak height).

        Parameters
        ----------
        left_base : float
            Left boundary from threshold
        right_base : float
            Right boundary from threshold

        Returns
        -------
        PeakBoundaries
            New boundaries without valley information

        Examples
        --------
        >>> bounds = PeakBoundaries.from_threshold(44.0, 46.0)
        >>> bounds.has_valleys()
        False

        References
        ----------
        THEORY.md Section 5.2.4: "threshold_fraction = 0.05 (5% of peak height)"
        """
        return cls(
            left_base=left_base,
            right_base=right_base,
            left_valley=None,
            right_valley=None,
        )

    def width(self) -> float:
        """
        Calculate peak width (right_base - left_base).

        Returns
        -------
        float
            Peak width (time units)

        Examples
        --------
        >>> bounds = PeakBoundaries(left_base=44.0, right_base=46.0)
        >>> bounds.width()
        2.0
        """
        return self.right_base - self.left_base

    def contains(self, time: float) -> bool:
        """
        Check if time point is within boundaries.

        Parameters
        ----------
        time : float
            Time point to check

        Returns
        -------
        bool
            True if left_base <= time <= right_base

        Examples
        --------
        >>> bounds = PeakBoundaries(left_base=44.0, right_base=46.0)
        >>> bounds.contains(45.0)
        True
        >>> bounds.contains(47.0)
        False
        """
        return self.left_base <= time <= self.right_base

    def has_valleys(self) -> bool:
        """
        Check if valleys were detected.

        Returns
        -------
        bool
            True if both valleys are present

        Examples
        --------
        >>> bounds = PeakBoundaries.from_valleys(44.0, 46.0)
        >>> bounds.has_valleys()
        True

        >>> bounds = PeakBoundaries.from_threshold(44.0, 46.0)
        >>> bounds.has_valleys()
        False
        """
        return self.left_valley is not None and self.right_valley is not None

    def has_partial_valleys(self) -> bool:
        """
        Check if only one valley was detected.

        This can occur at signal boundaries where only one side has a valley.

        Returns
        -------
        bool
            True if exactly one valley is present

        Examples
        --------
        >>> bounds = PeakBoundaries(
        ...     left_base=44.0,
        ...     right_base=46.0,
        ...     left_valley=44.2,
        ...     right_valley=None
        ... )
        >>> bounds.has_partial_valleys()
        True
        """
        return (self.left_valley is None) != (self.right_valley is None)

    def valley_width(self) -> Optional[float]:
        """
        Calculate width between valleys.

        Returns
        -------
        Optional[float]
            Width between valleys, or None if valleys not detected

        Examples
        --------
        >>> bounds = PeakBoundaries.from_valleys(44.2, 45.8)
        >>> bounds.valley_width()
        1.6

        >>> bounds = PeakBoundaries.from_threshold(44.0, 46.0)
        >>> bounds.valley_width() is None
        True
        """
        if self.has_valleys():
            return self.right_valley - self.left_valley
        return None

    def __str__(self) -> str:
        """String representation shows boundaries."""
        if self.has_valleys():
            return f"[{self.left_base:.2f}, {self.right_base:.2f}] (valleys: {self.left_valley:.2f}, {self.right_valley:.2f})"
        else:
            return f"[{self.left_base:.2f}, {self.right_base:.2f}]"

    def __repr__(self) -> str:
        """Detailed representation for debugging."""
        return (
            f"PeakBoundaries(left_base={self.left_base}, right_base={self.right_base}, "
            f"left_valley={self.left_valley}, right_valley={self.right_valley})"
        )
