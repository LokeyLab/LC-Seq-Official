"""
Virtual Compound - Proxy for consensus chromatogram processing.

This module provides VirtualCompound, a safe immutable proxy that allows
processing consensus chromatograms through the standard pipeline while
preserving the original compound data.
"""

from typing import List, Optional
from .compound import Compound
from .chromatogram import Chromatogram
from .peak import Peak


class VirtualCompound:
    """
    Virtual compound for consensus chromatogram processing.

    VirtualCompound acts as a proxy that:
    - Uses a consensus chromatogram for peak detection
    - Delegates all other properties to the real compound
    - Stores pipeline results separately (detected_peaks, selected_peak)
    - Maintains compatibility with hierarchy lookups via __eq__ and __hash__

    This allows consensus mode to process aggregated signals through the
    standard pipeline without mutating the original compound data.

    Safety Features
    ---------------
    - Immutable: original compound is never modified
    - Automatic delegation: __getattr__ delegates all non-overridden attributes
    - Restricted mutation: __setattr__ only allows setting detected_peaks/selected_peak
    - __slots__: prevents accidental attribute creation
    - Hierarchy compatible: __eq__ and __hash__ delegate to real compound

    Examples
    --------
    >>> # Create virtual compound with consensus chromatogram
    >>> virtual = VirtualCompound(real_compound, consensus_chrom)
    >>>
    >>> # Process through pipeline (uses consensus_chrom)
    >>> peaks = detector.detect_peaks(virtual.chromatogram)
    >>> virtual.detected_peaks = peaks
    >>>
    >>> # Virtual compound works with hierarchy
    >>> assert virtual in hierarchy.compounds  # True (via __hash__)
    >>> descendants = hierarchy.get_descendants(virtual)  # Works!
    >>>
    >>> # Transfer results to real compounds
    >>> for variant in equivalence_class.compounds:
    ...     variant.detected_peaks = virtual.detected_peaks

    Notes
    -----
    VirtualCompound is a temporary processing artifact. After pipeline
    execution, results should be transferred to the real compounds and
    the virtual compound discarded.

    References
    ----------
    THEORY.md Section 4.2: Consensus Mode
    """

    # Restrict attributes to prevent accidental creation
    __slots__ = ('_real', '_consensus_chromatogram', 'detected_peaks', 'selected_peak')

    def __init__(self, real_compound: Compound, consensus_chromatogram: Chromatogram):
        """
        Create virtual compound with consensus chromatogram.

        Parameters
        ----------
        real_compound : Compound
            Real compound to delegate properties to
        consensus_chromatogram : Chromatogram
            Consensus chromatogram to use for peak detection
        """
        # Use object.__setattr__ to bypass our custom __setattr__
        object.__setattr__(self, '_real', real_compound)
        object.__setattr__(self, '_consensus_chromatogram', consensus_chromatogram)
        object.__setattr__(self, 'detected_peaks', [])
        object.__setattr__(self, 'selected_peak', None)

    @property
    def chromatogram(self) -> Chromatogram:
        """
        Return consensus chromatogram instead of real compound's chromatogram.

        This is the key override that allows consensus signal processing.
        """
        return self._consensus_chromatogram

    def __getattr__(self, name: str):
        """
        Delegate all other attributes to real compound.

        This provides automatic delegation for:
        - building_blocks
        - level
        - monomer_level
        - positional_block_sequence
        - block_support_sequence
        - monomer_support_sequence
        - Any future properties added to Compound

        Safety: raises AttributeError for private attributes to prevent
        accessing internal state incorrectly.
        """
        if name.startswith('_'):
            raise AttributeError(f"VirtualCompound has no attribute '{name}'")
        return getattr(self._real, name)

    def __setattr__(self, name: str, value):
        """
        Only allow setting detected_peaks and selected_peak.

        This prevents accidental mutation of delegated attributes while
        allowing the pipeline to store results.

        Raises
        ------
        AttributeError
            If attempting to set any attribute except detected_peaks or selected_peak
        """
        if name in ('detected_peaks', 'selected_peak'):
            object.__setattr__(self, name, value)
        else:
            raise AttributeError(
                f"Cannot set attribute '{name}' on VirtualCompound. "
                f"Only 'detected_peaks' and 'selected_peak' can be modified."
            )

    def __eq__(self, other):
        """
        Equality based on real compound for hierarchy lookups.

        This ensures:
        - virtual_compound == real_compound returns True
        - virtual_compound in hierarchy.compounds returns True
        - Hierarchy edge lookups work correctly

        Parameters
        ----------
        other : object
            Object to compare with

        Returns
        -------
        bool
            True if other is the same real compound or a virtual compound
            wrapping the same real compound
        """
        if isinstance(other, VirtualCompound):
            return self._real == other._real
        return self._real == other

    def __hash__(self):
        """
        Hash based on real compound for hierarchy lookups.

        This ensures:
        - hash(virtual_compound) == hash(real_compound)
        - Virtual compound can be used as dict key
        - Hierarchy set operations work correctly

        Returns
        -------
        int
            Hash of the real compound
        """
        return hash(self._real)

    def __repr__(self):
        """
        String representation for debugging.

        Returns
        -------
        str
            Representation showing this is a virtual compound wrapping a real one
        """
        return f"VirtualCompound({self._real.positional_block_sequence})"
