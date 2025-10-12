"""
Chromatogram entity - represents elution profile across fractions/time.

Implementation based on THEORY.md Section 2.1, 5.0.1.
"""

from dataclasses import dataclass, field
from typing import Dict, Optional
import numpy as np
from numpy.typing import NDArray


@dataclass
class Chromatogram:
    """
    Represents elution profile across fractions/time for a DNA-encoded library member.

    A chromatogram contains raw signal data and multiple signal variants derived
    from preprocessing operations (baseline correction, derivatives, etc.).

    Attributes
    ----------
    time_points : NDArray[np.float64]
        Time points in seconds or minutes (absolute time from injection)
    counts : NDArray[np.float64]
        Raw signal intensities at each time point
    signal_variants : Dict[str, NDArray[np.float64]]
        Additional signal representations:
        - "corrected": Baseline-corrected signal
        - "derivative": First derivative
        - "derivative_2": Second derivative
        - Custom processing can add more variants

    Notes
    -----
    - Time units must be consistent across dataset (seconds or minutes)
    - Time points are absolute (not array indices)
    - All signals must have same length as time_points
    - Raw counts stored separately from processed variants
    - Signals may have different start/end times across experimental runs

    References
    ----------
    THEORY.md Section 2.1: Core Entities
    THEORY.md Section 5.0.1: DEL Signal Characteristics
    THEORY.md Section 2.3.1: Absolute Time Representation
    THEORY.md Section 2.3.2: Variable Signal Boundaries
    """

    time_points: NDArray[np.float64]
    counts: NDArray[np.float64]
    signal_variants: Dict[str, NDArray[np.float64]] = field(default_factory=dict)

    def __post_init__(self) -> None:
        """Validate chromatogram data."""
        # Convert to numpy arrays with proper dtype
        if not isinstance(self.time_points, np.ndarray):
            object.__setattr__(self, 'time_points', np.array(self.time_points, dtype=np.float64))
        else:
            # Ensure existing arrays are float64
            if self.time_points.dtype != np.float64:
                object.__setattr__(self, 'time_points', self.time_points.astype(np.float64))

        if not isinstance(self.counts, np.ndarray):
            object.__setattr__(self, 'counts', np.array(self.counts, dtype=np.float64))
        else:
            # Ensure existing arrays are float64
            if self.counts.dtype != np.float64:
                object.__setattr__(self, 'counts', self.counts.astype(np.float64))

        # Validate arrays
        if len(self.time_points) == 0:
            raise ValueError("Chromatogram must have at least one time point")

        if len(self.time_points) != len(self.counts):
            raise ValueError(
                f"Time points ({len(self.time_points)}) and counts ({len(self.counts)}) "
                f"must have same length"
            )

        # Validate time points are strictly increasing
        if len(self.time_points) > 1 and not np.all(np.diff(self.time_points) > 0):
            raise ValueError("Time points must be strictly increasing")

        # Validate all signal variants have correct length
        for variant_name, signal in self.signal_variants.items():
            if len(signal) != len(self.time_points):
                raise ValueError(
                    f"Signal variant '{variant_name}' has length {len(signal)}, "
                    f"expected {len(self.time_points)}"
                )

    @property
    def duration(self) -> float:
        """
        Total duration of chromatogram in time units.

        Returns
        -------
        float
            Duration from first to last time point

        Examples
        --------
        >>> chrom = Chromatogram(time_points=[0.0, 30.0, 60.0], counts=[10, 50, 20])
        >>> chrom.duration
        60.0
        """
        return float(self.time_points[-1] - self.time_points[0])

    @property
    def time_range(self) -> tuple[float, float]:
        """
        Time range of chromatogram.

        Returns
        -------
        tuple[float, float]
            (start_time, end_time)

        Examples
        --------
        >>> chrom = Chromatogram(time_points=[600.0, 1200.0, 1800.0], counts=[10, 50, 20])
        >>> chrom.time_range
        (600.0, 1800.0)
        """
        return (float(self.time_points[0]), float(self.time_points[-1]))

    def get_signal(self, variant: str = "raw") -> NDArray[np.float64]:
        """
        Get signal by variant name.

        Parameters
        ----------
        variant : str, optional
            Signal variant name. Use "raw" for original counts.
            Default is "raw".

        Returns
        -------
        NDArray[np.float64]
            Signal array

        Raises
        ------
        KeyError
            If variant not found

        Examples
        --------
        >>> chrom = Chromatogram(time_points=[0, 1, 2], counts=[10, 20, 15])
        >>> chrom.get_signal("raw")
        array([10., 20., 15.])

        >>> chrom.signal_variants["corrected"] = np.array([5, 15, 10])
        >>> chrom.get_signal("corrected")
        array([ 5., 15., 10.])
        """
        if variant == "raw":
            return self.counts
        
        if variant not in self.signal_variants:
            raise KeyError(
                f"Signal variant '{variant}' not found. "
                f"Available: ['raw'] + {list(self.signal_variants.keys())}"
            )
        
        return self.signal_variants[variant]

    def add_signal_variant(self, name: str, signal: NDArray[np.float64]) -> None:
        """
        Add a new signal variant (e.g., baseline-corrected, derivative).

        Parameters
        ----------
        name : str
            Variant name (e.g., "corrected", "derivative", "derivative_2")
        signal : NDArray[np.float64]
            Signal array (must match length of time_points)

        Raises
        ------
        ValueError
            If signal length doesn't match time_points length

        Examples
        --------
        >>> chrom = Chromatogram(time_points=[0, 1, 2], counts=[10, 20, 15])
        >>> corrected = np.array([5, 15, 10])
        >>> chrom.add_signal_variant("corrected", corrected)
        >>> chrom.get_signal("corrected")
        array([ 5., 15., 10.])
        """
        signal_array = np.asarray(signal, dtype=np.float64)
        
        if len(signal_array) != len(self.time_points):
            raise ValueError(
                f"Signal length ({len(signal_array)}) must match "
                f"time_points length ({len(self.time_points)})"
            )
        
        self.signal_variants[name] = signal_array

    def has_signal_variant(self, name: str) -> bool:
        """
        Check if signal variant exists.

        Parameters
        ----------
        name : str
            Variant name to check

        Returns
        -------
        bool
            True if variant exists, False otherwise

        Examples
        --------
        >>> chrom = Chromatogram(time_points=[0, 1], counts=[10, 20])
        >>> chrom.has_signal_variant("corrected")
        False
        >>> chrom.add_signal_variant("corrected", np.array([5, 15]))
        >>> chrom.has_signal_variant("corrected")
        True
        """
        if name == "raw":
            return True
        return name in self.signal_variants

    def __len__(self) -> int:
        """Number of time points in chromatogram."""
        return len(self.time_points)

    def __eq__(self, other: object) -> bool:
        """
        Test equality based on array values.

        Parameters
        ----------
        other : object
            Object to compare with

        Returns
        -------
        bool
            True if chromatograms have equal data
        """
        if not isinstance(other, Chromatogram):
            return NotImplemented

        # Compare arrays using np.array_equal
        if not np.array_equal(self.time_points, other.time_points):
            return False
        if not np.array_equal(self.counts, other.counts):
            return False

        # Compare signal variants
        if set(self.signal_variants.keys()) != set(other.signal_variants.keys()):
            return False

        for key in self.signal_variants:
            if not np.array_equal(self.signal_variants[key], other.signal_variants[key]):
                return False

        return True

    def __repr__(self) -> str:
        """Detailed representation for debugging."""
        variants = ["raw"] + list(self.signal_variants.keys())
        return (
            f"Chromatogram(n_points={len(self)}, "
            f"time_range={self.time_range}, "
            f"variants={variants})"
        )
