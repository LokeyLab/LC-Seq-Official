"""
MonomerSequence value object - immutable decomposed chemical sequence.

Implementation based on THEORY.md Section 2.2, 1.5.3.
"""

from dataclasses import dataclass
from typing import List, Tuple


@dataclass(frozen=True)
class MonomerSequence:
    """
    Represents the fully decomposed chemical sequence as individual monomers.

    The monomer sequence is the chemical identity of the peptide, obtained by
    decomposing all building blocks (including composite blocks) into their
    constituent monomers. This represents what molecule actually exists
    chemically, independent of how it was synthesized.

    Attributes
    ----------
    monomers : Tuple[str, ...]
        Ordered tuple of monomer codes from N→C (left to right)
        Tuple for immutability and hashability

    Notes
    -----
    - N→C ordering: First monomer = N-terminus, Last = C-terminus
    - Immutable value object (frozen dataclass with tuple)
    - Composite blocks are fully decomposed
    - Null blocks contribute no monomers
    - Equality based on chemical sequence only
    - Multiple positional sequences can produce same monomer sequence

    Examples
    --------
    >>> # Simple sequence
    >>> seq = MonomerSequence.from_list(["Leu", "Ala", "Pro"])
    >>> seq.to_string()
    'Leu-Ala-Pro'
    >>> len(seq)
    3

    >>> # Empty sequence (all nulls)
    >>> empty = MonomerSequence.from_list([])
    >>> len(empty)
    0
    >>> empty.is_empty()
    True

    References
    ----------
    THEORY.md Section 2.2: Monomer Sequence definition
    THEORY.md Section 1.5.3: Monomer-Level Decomposition
    THEORY.md Section 1.2: Chemical Identity vs Positional Identity
    """

    monomers: Tuple[str, ...]

    def __post_init__(self) -> None:
        """Validate monomer sequence properties."""
        # Validate no empty monomer codes
        for i, monomer in enumerate(self.monomers):
            if not monomer or not monomer.strip():
                raise ValueError(f"Monomer at position {i} cannot be empty or whitespace")

    @classmethod
    def from_list(cls, monomers: List[str]) -> "MonomerSequence":
        """
        Create sequence from ordered list of monomers.

        Parameters
        ----------
        monomers : List[str]
            Monomer codes in N→C order (left to right)

        Returns
        -------
        MonomerSequence
            New immutable sequence

        Examples
        --------
        >>> seq = MonomerSequence.from_list(["Leu", "Ala", "Pro"])
        >>> len(seq)
        3
        """
        return cls(monomers=tuple(monomers))

    @classmethod
    def from_string(cls, sequence_str: str) -> "MonomerSequence":
        """
        Create sequence from hyphen-separated string.

        Parameters
        ----------
        sequence_str : str
            Hyphen-separated monomer codes (N→C order)
            e.g., "Leu-Ala-Pro"

        Returns
        -------
        MonomerSequence
            New sequence parsed from string

        Examples
        --------
        >>> seq = MonomerSequence.from_string("Leu-Ala-Pro")
        >>> seq.to_string()
        'Leu-Ala-Pro'

        >>> # Empty string produces empty sequence
        >>> empty = MonomerSequence.from_string("")
        >>> empty.is_empty()
        True
        """
        if not sequence_str or not sequence_str.strip():
            return cls(monomers=())

        monomers = [m.strip() for m in sequence_str.split("-")]
        return cls(monomers=tuple(monomers))

    def to_string(self) -> str:
        """
        Convert to hyphen-separated string (N→C order).

        Returns
        -------
        str
            Hyphen-separated sequence
            e.g., "Leu-Ala-Pro"

        Examples
        --------
        >>> seq = MonomerSequence.from_list(["Leu", "Ala", "Pro"])
        >>> seq.to_string()
        'Leu-Ala-Pro'
        """
        return "-".join(self.monomers)

    def is_empty(self) -> bool:
        """
        Check if sequence is empty (no monomers).

        Returns
        -------
        bool
            True if sequence has no monomers (e.g., all null blocks)

        Examples
        --------
        >>> empty = MonomerSequence.from_list([])
        >>> empty.is_empty()
        True

        >>> seq = MonomerSequence.from_list(["Leu"])
        >>> seq.is_empty()
        False
        """
        return len(self.monomers) == 0

    def get_n_terminus(self) -> str:
        """
        Get N-terminal monomer.

        Returns
        -------
        str
            First monomer in sequence (N-terminus)

        Raises
        ------
        IndexError
            If sequence is empty

        Examples
        --------
        >>> seq = MonomerSequence.from_list(["Leu", "Ala", "Pro"])
        >>> seq.get_n_terminus()
        'Leu'
        """
        return self.monomers[0]

    def get_c_terminus(self) -> str:
        """
        Get C-terminal monomer.

        Returns
        -------
        str
            Last monomer in sequence (C-terminus)

        Raises
        ------
        IndexError
            If sequence is empty

        Examples
        --------
        >>> seq = MonomerSequence.from_list(["Leu", "Ala", "Pro"])
        >>> seq.get_c_terminus()
        'Pro'
        """
        return self.monomers[-1]

    def get_monomer_at(self, position: int) -> str:
        """
        Get monomer at position (0-indexed from N-terminus).

        Parameters
        ----------
        position : int
            Position in sequence (0 = N-terminus)

        Returns
        -------
        str
            Monomer at position

        Raises
        ------
        IndexError
            If position out of range

        Examples
        --------
        >>> seq = MonomerSequence.from_list(["Leu", "Ala", "Pro"])
        >>> seq.get_monomer_at(1)
        'Ala'
        """
        return self.monomers[position]

    def reverse(self) -> "MonomerSequence":
        """
        Create reversed sequence (C→N instead of N→C).

        Returns
        -------
        MonomerSequence
            New sequence with monomers in reverse order

        Examples
        --------
        >>> seq = MonomerSequence.from_list(["Leu", "Ala", "Pro"])
        >>> rev = seq.reverse()
        >>> rev.to_string()
        'Pro-Ala-Leu'
        """
        return MonomerSequence(monomers=tuple(reversed(self.monomers)))

    def __len__(self) -> int:
        """Return number of monomers."""
        return len(self.monomers)

    def __iter__(self):
        """Iterate over monomers (N→C order)."""
        return iter(self.monomers)

    def __getitem__(self, index):
        """Access monomer by index or slice."""
        result = self.monomers[index]
        if isinstance(index, slice):
            return MonomerSequence(monomers=result)
        return result

    def __str__(self) -> str:
        """String representation shows monomer sequence."""
        return self.to_string()

    def __repr__(self) -> str:
        """Detailed representation for debugging."""
        return f"MonomerSequence('{self.to_string()}')"
