"""
BuildingBlockSequence value object - immutable position-to-block mapping.

Implementation based on THEORY.md Section 2.2, 1.5.1.
"""

from dataclasses import dataclass
from typing import Dict, List, Tuple
from lcseq.domain.entities.building_block import BuildingBlock


@dataclass(frozen=True)
class BuildingBlockSequence:
    """
    Represents the mapping of synthesis positions to building blocks.

    A building block sequence defines what building block was used at each
    synthesis cycle/position. This is the positional representation that
    captures the synthesis path, independent of the final chemical identity.

    Attributes
    ----------
    blocks : Dict[int, BuildingBlock]
        Mapping from position (cycle number) to BuildingBlock
        Position 0 = C-terminus (synthesized first)
        Higher positions = toward N-terminus

    Notes
    -----
    - N→C sequence convention: Position 0 = C-terminus (rightmost)
    - Immutable value object (frozen dataclass)
    - Positions must be contiguous: 0, 1, 2, ..., n-1
    - Multiple sequences can produce same chemical peptide (positional variants)

    Examples
    --------
    >>> from lcseq.domain.entities.building_block import BuildingBlock
    >>> bb0 = BuildingBlock.from_code(0, "Pro")
    >>> bb1 = BuildingBlock.from_code(1, "Leu")
    >>> bb2 = BuildingBlock.from_code(2, "Null")
    >>> seq = BuildingBlockSequence.from_blocks([bb0, bb1, bb2])
    >>> len(seq)
    3
    >>> seq.get_block_at(0).code
    'Pro'

    References
    ----------
    THEORY.md Section 2.2: Sequence Representations
    THEORY.md Section 1.5.1: N→C Sequence Convention
    """

    blocks: Dict[int, BuildingBlock]

    def __post_init__(self) -> None:
        """Validate building block sequence properties."""
        if not self.blocks:
            raise ValueError("Building block sequence cannot be empty")

        # Validate positions are non-negative
        for position in self.blocks.keys():
            if position < 0:
                raise ValueError(f"Position must be non-negative, got {position}")

        # Validate positions are contiguous (0, 1, 2, ..., n-1)
        positions = sorted(self.blocks.keys())
        expected = list(range(len(positions)))
        if positions != expected:
            raise ValueError(
                f"Positions must be contiguous starting from 0. "
                f"Expected {expected}, got {positions}"
            )

        # Validate each block's cycle matches its position
        for position, block in self.blocks.items():
            if block.cycle != position:
                raise ValueError(
                    f"Block at position {position} has mismatched cycle {block.cycle}"
                )

    @classmethod
    def from_blocks(cls, blocks: List[BuildingBlock]) -> "BuildingBlockSequence":
        """
        Create sequence from ordered list of building blocks.

        Parameters
        ----------
        blocks : List[BuildingBlock]
            Building blocks in order from C-terminus (pos 0) to N-terminus

        Returns
        -------
        BuildingBlockSequence
            New sequence with blocks mapped to positions

        Examples
        --------
        >>> bb0 = BuildingBlock.from_code(0, "Pro")
        >>> bb1 = BuildingBlock.from_code(1, "Leu")
        >>> seq = BuildingBlockSequence.from_blocks([bb0, bb1])
        >>> len(seq)
        2
        """
        block_dict = {block.cycle: block for block in blocks}
        return cls(blocks=block_dict)

    @classmethod
    def from_codes(cls, codes: List[str]) -> "BuildingBlockSequence":
        """
        Create sequence from ordered list of building block codes.

        Automatically detects null blocks and assigns cycles.

        Parameters
        ----------
        codes : List[str]
            Building block codes in order from C-terminus to N-terminus

        Returns
        -------
        BuildingBlockSequence
            New sequence with auto-generated BuildingBlock instances

        Examples
        --------
        >>> seq = BuildingBlockSequence.from_codes(["Pro", "Leu", "Null"])
        >>> seq.get_block_at(0).code
        'Pro'
        >>> seq.get_block_at(2).is_null
        True
        """
        blocks = [BuildingBlock.from_code(i, code) for i, code in enumerate(codes)]
        return cls.from_blocks(blocks)

    def get_block_at(self, position: int) -> BuildingBlock:
        """
        Get building block at specified position.

        Parameters
        ----------
        position : int
            Position to retrieve (0 = C-terminus)

        Returns
        -------
        BuildingBlock
            Building block at the position

        Raises
        ------
        KeyError
            If position is not in sequence
        """
        return self.blocks[position]

    def get_non_null_blocks(self) -> List[BuildingBlock]:
        """
        Get all non-null building blocks in sequence order.

        Returns blocks in order from C-terminus (lowest position) to
        N-terminus (highest position), excluding null blocks.

        Returns
        -------
        List[BuildingBlock]
            Non-null blocks in position order

        Examples
        --------
        >>> seq = BuildingBlockSequence.from_codes(["Pro", "Null", "Leu"])
        >>> blocks = seq.get_non_null_blocks()
        >>> [b.code for b in blocks]
        ['Pro', 'Leu']
        """
        positions = sorted(self.blocks.keys())
        return [
            self.blocks[pos]
            for pos in positions
            if not self.blocks[pos].is_null
        ]

    def to_positional_string(self) -> str:
        """
        Convert to positional sequence string (includes nulls).

        Returns
        -------
        str
            Hyphen-separated sequence from N→C (left to right)
            e.g., "Leu-Null-Pro" (Leu at position 2, Pro at position 0)

        Notes
        -----
        String representation is N→C (left to right) for readability,
        but internal storage is C→N (position 0 = C-terminus).

        Examples
        --------
        >>> seq = BuildingBlockSequence.from_codes(["Pro", "Null", "Leu"])
        >>> seq.to_positional_string()
        'Leu-Null-Pro'
        """
        # Get positions in reverse order for N→C display
        positions = sorted(self.blocks.keys(), reverse=True)
        codes = [self.blocks[pos].code for pos in positions]
        return "-".join(codes)

    def to_residue_string(self) -> str:
        """
        Convert to block support sequence string (non-null blocks only).

        Returns
        -------
        str
            Hyphen-separated non-null sequence from N→C
            e.g., "Leu-Pro" (from "Leu-Null-Pro")

        Notes
        -----
        This is the block support sequence used for grouping
        positional variants. Multiple positional sequences can have
        the same block support sequence.

        Examples
        --------
        >>> seq = BuildingBlockSequence.from_codes(["Pro", "Null", "Leu"])
        >>> seq.to_residue_string()
        'Leu-Pro'

        References
        ----------
        THEORY.md Section 2.2: Block Support Sequence definition
        """
        non_null = self.get_non_null_blocks()
        # Reverse for N→C display
        codes = [block.code for block in reversed(non_null)]
        return "-".join(codes)

    def __len__(self) -> int:
        """Return number of positions (including nulls)."""
        return len(self.blocks)

    def __str__(self) -> str:
        """String representation shows positional sequence."""
        return self.to_positional_string()

    def __repr__(self) -> str:
        """Detailed representation for debugging."""
        positions = sorted(self.blocks.keys())
        blocks_str = ", ".join(
            f"{pos}: {self.blocks[pos].code}" for pos in positions
        )
        return f"BuildingBlockSequence({blocks_str})"
