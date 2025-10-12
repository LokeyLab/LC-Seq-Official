"""
RetentionTime value object - scalar time value with unit handling.

Implementation based on THEORY.md Section 2.3.1.
"""

from dataclasses import dataclass
from enum import Enum


class TimeUnit(Enum):
    """Time units for retention time."""

    SECONDS = "s"
    MINUTES = "min"


@dataclass(frozen=True)
class RetentionTime:
    """
    Represents an absolute retention time with unit handling.

    Retention times are scalar values representing the time at which a
    compound elutes from the chromatography column. All comparisons and
    operations are performed on absolute time values, independent of
    signal boundaries or sampling rates.

    Attributes
    ----------
    value : float
        Time value (must be non-negative)
    unit : TimeUnit
        Time unit (SECONDS or MINUTES)

    Notes
    -----
    - Immutable value object (frozen dataclass)
    - Non-negative values only (time from injection)
    - Unit conversion methods for interoperability
    - Comparison operators based on absolute time
    - Physically meaningful values (not array indices)

    Examples
    --------
    >>> # Create retention time in seconds
    >>> rt = RetentionTime(value=45.2, unit=TimeUnit.SECONDS)
    >>> rt.in_seconds()
    45.2

    >>> # Create in minutes and convert
    >>> rt_min = RetentionTime(value=0.75, unit=TimeUnit.MINUTES)
    >>> rt_min.in_seconds()
    45.0

    >>> # Comparison (automatic unit handling)
    >>> rt1 = RetentionTime(value=45.0, unit=TimeUnit.SECONDS)
    >>> rt2 = RetentionTime(value=0.75, unit=TimeUnit.MINUTES)
    >>> rt1 == rt2
    True

    References
    ----------
    THEORY.md Section 2.3.1: Absolute Time Representation
    """

    value: float
    unit: TimeUnit

    def __post_init__(self) -> None:
        """Validate retention time properties."""
        if self.value < 0:
            raise ValueError(
                f"Retention time must be non-negative, got {self.value}"
            )

        # Reasonable upper bound check (10 hours = 600 minutes = 36000 seconds)
        max_seconds = 36000.0
        if self.in_seconds() > max_seconds:
            raise ValueError(
                f"Retention time {self.value} {self.unit.value} exceeds "
                f"reasonable limit (10 hours)"
            )

    @classmethod
    def from_seconds(cls, seconds: float) -> "RetentionTime":
        """
        Create retention time from seconds.

        Parameters
        ----------
        seconds : float
            Time in seconds

        Returns
        -------
        RetentionTime
            New retention time in seconds

        Examples
        --------
        >>> rt = RetentionTime.from_seconds(45.2)
        >>> rt.value
        45.2
        >>> rt.unit
        <TimeUnit.SECONDS: 's'>
        """
        return cls(value=seconds, unit=TimeUnit.SECONDS)

    @classmethod
    def from_minutes(cls, minutes: float) -> "RetentionTime":
        """
        Create retention time from minutes.

        Parameters
        ----------
        minutes : float
            Time in minutes

        Returns
        -------
        RetentionTime
            New retention time in minutes

        Examples
        --------
        >>> rt = RetentionTime.from_minutes(0.75)
        >>> rt.value
        0.75
        >>> rt.unit
        <TimeUnit.MINUTES: 'min'>
        """
        return cls(value=minutes, unit=TimeUnit.MINUTES)

    def in_seconds(self) -> float:
        """
        Get time value in seconds.

        Returns
        -------
        float
            Time in seconds (converts if needed)

        Examples
        --------
        >>> rt = RetentionTime.from_minutes(1.5)
        >>> rt.in_seconds()
        90.0
        """
        if self.unit == TimeUnit.SECONDS:
            return self.value
        else:  # MINUTES
            return self.value * 60.0

    def in_minutes(self) -> float:
        """
        Get time value in minutes.

        Returns
        -------
        float
            Time in minutes (converts if needed)

        Examples
        --------
        >>> rt = RetentionTime.from_seconds(90.0)
        >>> rt.in_minutes()
        1.5
        """
        if self.unit == TimeUnit.MINUTES:
            return self.value
        else:  # SECONDS
            return self.value / 60.0

    def to_unit(self, target_unit: TimeUnit) -> "RetentionTime":
        """
        Convert to different time unit.

        Parameters
        ----------
        target_unit : TimeUnit
            Target unit for conversion

        Returns
        -------
        RetentionTime
            New retention time in target unit

        Examples
        --------
        >>> rt = RetentionTime.from_seconds(90.0)
        >>> rt_min = rt.to_unit(TimeUnit.MINUTES)
        >>> rt_min.value
        1.5
        >>> rt_min.unit
        <TimeUnit.MINUTES: 'min'>
        """
        if self.unit == target_unit:
            return self

        if target_unit == TimeUnit.SECONDS:
            return RetentionTime.from_seconds(self.in_seconds())
        else:  # MINUTES
            return RetentionTime.from_minutes(self.in_minutes())

    def matches(self, other: "RetentionTime", tolerance: float) -> bool:
        """
        Check if retention time matches another within tolerance.

        Parameters
        ----------
        other : RetentionTime
            Other retention time to compare
        tolerance : float
            Maximum allowed difference (in seconds)

        Returns
        -------
        bool
            True if |self - other| <= tolerance

        Examples
        --------
        >>> rt1 = RetentionTime.from_seconds(45.0)
        >>> rt2 = RetentionTime.from_seconds(45.2)
        >>> rt1.matches(rt2, tolerance=0.5)
        True
        >>> rt1.matches(rt2, tolerance=0.1)
        False

        Notes
        -----
        Comparison is always in seconds regardless of original units.

        References
        ----------
        THEORY.md Section 2.3.1: "Does peak at 45.2s match expected 45.0s?
        → |45.2 - 45.0| = 0.2s"
        """
        diff = abs(self.in_seconds() - other.in_seconds())
        return diff <= tolerance

    def __eq__(self, other: object) -> bool:
        """
        Check equality (same absolute time).

        Compares absolute time in seconds, regardless of original units.
        """
        if not isinstance(other, RetentionTime):
            return NotImplemented

        # Compare in seconds (handles unit conversion)
        return abs(self.in_seconds() - other.in_seconds()) < 1e-9

    def __lt__(self, other: "RetentionTime") -> bool:
        """Less than comparison."""
        return self.in_seconds() < other.in_seconds()

    def __le__(self, other: "RetentionTime") -> bool:
        """Less than or equal comparison."""
        return self.in_seconds() <= other.in_seconds()

    def __gt__(self, other: "RetentionTime") -> bool:
        """Greater than comparison."""
        return self.in_seconds() > other.in_seconds()

    def __ge__(self, other: "RetentionTime") -> bool:
        """Greater than or equal comparison."""
        return self.in_seconds() >= other.in_seconds()

    def __hash__(self) -> int:
        """Hash based on absolute time in seconds."""
        # Round to avoid floating point issues
        return hash(round(self.in_seconds(), 9))

    def __str__(self) -> str:
        """String representation with units."""
        return f"{self.value}{self.unit.value}"

    def __repr__(self) -> str:
        """Detailed representation for debugging."""
        return f"RetentionTime({self.value}, {self.unit.name})"
