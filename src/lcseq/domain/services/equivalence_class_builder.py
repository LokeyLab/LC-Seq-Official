"""
EquivalenceClassBuilder - Groups compounds by block support sequence.

Implementation based on THEORY.md Section 4.2.1.
"""

from typing import List, Dict
from ..entities.compound import Compound
from ..models.equivalence_class import EquivalenceClass


class EquivalenceClassBuilder:
    """
    Groups compounds into equivalence classes by block support sequence.

    Compounds with the same block support sequence (non-null building blocks)
    but different positional block sequences are grouped into the same
    equivalence class. This identifies positional variants of the
    same chemical peptide at block granularity.

    Notes
    -----
    - Stateless service (no instance state)
    - Pure domain logic (no I/O operations)
    - Groups by block_support_sequence
    - Each class contains positional variant members

    Examples
    --------
    >>> from lcseq.domain.entities.building_block import BuildingBlock
    >>> from lcseq.domain.entities.chromatogram import Chromatogram
    >>> import numpy as np
    >>>
    >>> # Create compounds with same block support sequence
    >>> chromatogram = Chromatogram(
    ...     time_points=np.array([1.0, 2.0, 3.0]),
    ...     counts=np.array([100.0, 200.0, 150.0])
    ... )
    >>> bb_leu = BuildingBlock.from_code(0, "Leu")
    >>> bb_pro = BuildingBlock.from_code(1, "Pro")
    >>> bb_null = BuildingBlock.from_code(2, "Null")
    >>>
    >>> # Leu-Pro-Null and Leu-Null-Pro have same block support sequence "Leu-Pro"
    >>> compound1 = Compound([bb_null, bb_pro, bb_leu], chromatogram)
    >>> compound2 = Compound([bb_null, bb_null, bb_leu], chromatogram)
    >>>
    >>> builder = EquivalenceClassBuilder()
    >>> classes = builder.build([compound1, compound2])
    >>> # Classes will group by unique block support sequences

    References
    ----------
    THEORY.md Section 4.2.1: Equivalence Classes (Positional Variants)
    THEORY.md Section 2.2: Sequence Representations
    THEORY.md Section 1.2: Chemical Identity vs Positional Identity
    """

    def build(self, compounds: List[Compound]) -> List[EquivalenceClass]:
        """
        Build equivalence classes from list of compounds.

        Parameters
        ----------
        compounds : List[Compound]
            All compounds to group

        Returns
        -------
        List[EquivalenceClass]
            Equivalence classes, each containing members with same block support sequence

        Notes
        -----
        - Groups by block_support_sequence (non-null building blocks only)
        - Each class has unique block support sequence
        - Order of classes is deterministic (sorted by block support sequence)

        Algorithm
        ---------
        1. Group compounds by block support sequence
        2. Create EquivalenceClass for each unique sequence
        3. Return sorted list of classes

        Examples
        --------
        >>> compounds = [
        ...     Compound([Leu, Pro, Null], chrom),    # block_support_sequence: "Leu-Pro"
        ...     Compound([Leu, Null, Pro], chrom),    # block_support_sequence: "Leu-Pro"
        ...     Compound([Val, Pro, Null], chrom),    # block_support_sequence: "Val-Pro"
        ... ]
        >>> classes = builder.build(compounds)
        >>> len(classes)
        2
        >>> # One class for "Leu-Pro" (2 members), one for "Val-Pro" (1 member)
        """
        # Group by block support sequence
        groups: Dict[str, List[Compound]] = {}

        for compound in compounds:
            block_support_seq = self._get_block_support_sequence(compound)

            if block_support_seq not in groups:
                groups[block_support_seq] = []

            groups[block_support_seq].append(compound)

        # Create equivalence classes
        classes = []
        for block_support_seq in sorted(groups.keys()):
            eq_class = EquivalenceClass(
                block_support_sequence=block_support_seq,
                members=set(groups[block_support_seq])
            )
            classes.append(eq_class)

        return classes

    def _get_block_support_sequence(self, compound: Compound) -> str:
        """
        Get block support sequence for a compound.

        Parameters
        ----------
        compound : Compound
            Compound to get block support sequence from

        Returns
        -------
        str
            Block support sequence (non-null building blocks only)

        Notes
        -----
        Uses the compound's block_support_sequence property which filters
        out null building blocks and joins remaining codes.

        Examples
        --------
        >>> compound = Compound([Leu, Null, Pro], chrom)
        >>> builder._get_block_support_sequence(compound)
        'Leu-Pro'

        >>> all_null = Compound([Null, Null, Null], chrom)
        >>> builder._get_block_support_sequence(all_null)
        ''
        """
        return compound.block_support_sequence
