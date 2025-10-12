"""
BuildingBlock entity - represents a chemical building block in library synthesis.

Implementation based on THEORY.md Section 2.1, 1.5.3.
"""

from dataclasses import dataclass
from typing import List


@dataclass(frozen=True)
class BuildingBlock:
    """
    Represents a chemical building block used in DNA-encoded library synthesis.

    A building block may be:
    - A single monomer (e.g., "Leu")
    - A composite block containing multiple monomers (e.g., "Leu-Ala-Val")
    - A null block representing no synthesis at this position

    Composite blocks are written identically to regular sequences. Only library
    design metadata distinguishes them from individual monomers.

    Attributes
    ----------
    cycle : int
        Synthesis cycle/round number (0-indexed from C-terminus)
    code : str
        Building block identifier (e.g., "Leu", "Leu-Ala-Val", "Null")
    is_null : bool
        Whether this is a null block (no synthesis)

    Notes
    -----
    - N→C sequence convention: Position 0 = C-terminus (rightmost)
    - Null detection: Any code containing "null" (case-insensitive) → null block
    - Composite blocks: "Leu-Ala-Val" looks identical to tripeptide sequence
      (only library metadata distinguishes them)

    References
    ----------
    THEORY.md Section 2.1: Core Entities
    THEORY.md Section 1.5.3: Monomer-Level Decomposition
    """

    cycle: int
    code: str
    is_null: bool

    def __post_init__(self) -> None:
        """Validate building block properties."""
        if self.cycle < 0:
            raise ValueError(f"Cycle must be non-negative, got {self.cycle}")
        if not self.code:
            raise ValueError("Building block code cannot be empty")

    @classmethod
    def from_code(cls, cycle: int, code: str) -> "BuildingBlock":
        """
        Create BuildingBlock with automatic null detection.

        Parameters
        ----------
        cycle : int
            Synthesis cycle number
        code : str
            Building block identifier

        Returns
        -------
        BuildingBlock
            New building block instance with is_null automatically set

        Examples
        --------
        >>> BuildingBlock.from_code(0, "Leu")
        BuildingBlock(cycle=0, code='Leu', is_null=False)

        >>> BuildingBlock.from_code(1, "Null")
        BuildingBlock(cycle=1, code='Null', is_null=True)

        >>> BuildingBlock.from_code(2, "null_variant")
        BuildingBlock(cycle=2, code='null_variant', is_null=True)
        """
        is_null = "null" in code.lower()
        return cls(cycle=cycle, code=code, is_null=is_null)

    @classmethod
    def null(cls, cycle: int) -> "BuildingBlock":
        """
        Create a null building block (represents no synthesis at this position).

        Parameters
        ----------
        cycle : int
            Synthesis cycle number

        Returns
        -------
        BuildingBlock
            Null building block instance

        Examples
        --------
        >>> BuildingBlock.null(2)
        BuildingBlock(cycle=2, code='AgxNull', is_null=True)
        """
        return cls(cycle=cycle, code="AgxNull", is_null=True)

    def decompose_to_monomers(self) -> List[str]:
        """
        Decompose building block into constituent monomers.

        For composite blocks (e.g., "Leu-Ala-Val"), splits into individual monomers.
        For single monomers, returns single-element list.
        For null blocks, returns empty list.

        Returns
        -------
        List[str]
            List of monomer codes (empty for null blocks)

        Examples
        --------
        >>> bb = BuildingBlock.from_code(0, "Leu")
        >>> bb.decompose_to_monomers()
        ['Leu']

        >>> bb = BuildingBlock.from_code(1, "Leu-Ala-Val")
        >>> bb.decompose_to_monomers()
        ['Leu', 'Ala', 'Val']

        >>> bb = BuildingBlock.from_code(2, "Null")
        >>> bb.decompose_to_monomers()
        []

        Notes
        -----
        Composite blocks are split on hyphen ('-') character. This is the
        standard notation for peptide sequences (N→C convention).

        References
        ----------
        THEORY.md Section 1.5.3: Monomer-Level Decomposition
        """
        if self.is_null:
            return []

        # Split on hyphen to get individual monomers
        # "Leu-Ala-Val" → ["Leu", "Ala", "Val"]
        # "Leu" → ["Leu"]
        return [monomer.strip() for monomer in self.code.split("-")]

    def __str__(self) -> str:
        """String representation of building block."""
        return self.code

    def __repr__(self) -> str:
        """Detailed representation for debugging."""
        return f"BuildingBlock(cycle={self.cycle}, code='{self.code}', is_null={self.is_null})"
