"""
AnalysisResult - complete analysis output.

Implementation based on THEORY.md Part 5, 6.
"""

from dataclasses import dataclass, field
from typing import Dict, List, Optional
from .compound_hierarchy import CompoundHierarchy
from .equivalence_class import EquivalenceClass
from .peak_classification import PeakClassification
from ..entities.compound import Compound
from ..entities.peak import ValidationStatus


@dataclass
class AnalysisResult:
    """
    Complete output from LC-Seq analysis run.

    Aggregates all analysis results: hierarchy structure, equivalence classes,
    peak classifications, and validation results.

    Attributes
    ----------
    hierarchy : CompoundHierarchy
        DAG/Poset of truncation relationships
    equivalence_classes : List[EquivalenceClass]
        Grouped positional variants (used in pooled mode)
    peak_classifications : Dict[Compound, List[PeakClassification]]
        All peak classifications per compound
    validation_results : Dict[Compound, ValidationStatus]
        Synthesis validation status per compound
    metadata : Dict[str, Any]
        Optional metadata (timestamps, parameters, etc.)

    Notes
    -----
    - Immutable once created (use dataclass defaults)
    - Contains complete analysis state for reporting/export
    - Validation results independent of peak classifications
    - Equivalence classes may be empty in individual mode

    Examples
    --------
    >>> from lcseq.domain.models.compound_hierarchy import HierarchyMode
    >>> from lcseq.domain.entities.peak import ValidationStatus
    >>>
    >>> # Create analysis result
    >>> hierarchy = CompoundHierarchy(mode=HierarchyMode.BUILDING_BLOCK)
    >>> result = AnalysisResult(
    ...     hierarchy=hierarchy,
    ...     equivalence_classes=[],
    ...     peak_classifications={},
    ...     validation_results={}
    ... )
    >>>
    >>> # Add validation result
    >>> result.validation_results[compound] = ValidationStatus.VALIDATED

    References
    ----------
    THEORY.md Section 5.3: Peak Classification
    THEORY.md Section 6.10: Validation Classification
    THEORY.md Section 4.2.1: Equivalence Classes
    """

    hierarchy: CompoundHierarchy
    equivalence_classes: List[EquivalenceClass] = field(default_factory=list)
    peak_classifications: Dict[Compound, List[PeakClassification]] = field(
        default_factory=dict
    )
    validation_results: Dict[Compound, ValidationStatus] = field(default_factory=dict)
    metadata: Dict[str, any] = field(default_factory=dict)

    def add_peak_classification(
        self, compound: Compound, classification: PeakClassification
    ) -> None:
        """
        Add a peak classification for a compound.

        Parameters
        ----------
        compound : Compound
            The compound being analyzed
        classification : PeakClassification
            Peak classification result

        Notes
        -----
        Compound can have multiple peak classifications (multiple peaks detected).
        """
        if compound not in self.peak_classifications:
            self.peak_classifications[compound] = []
        self.peak_classifications[compound].append(classification)

    def set_validation_status(
        self, compound: Compound, status: ValidationStatus
    ) -> None:
        """
        Set validation status for a compound.

        Parameters
        ----------
        compound : Compound
            The compound being validated
        status : ValidationStatus
            Validation result

        Notes
        -----
        Overwrites previous validation status if already set.
        """
        self.validation_results[compound] = status

    def get_product_peak(self, compound: Compound) -> Optional[PeakClassification]:
        """
        Get the putative product peak classification for a compound.

        Parameters
        ----------
        compound : Compound
            Query compound

        Returns
        -------
        Optional[PeakClassification]
            Product peak classification, or None if not found

        Notes
        -----
        Returns first PUTATIVE_PRODUCT classification found.
        Most compounds should have at most one product peak.
        """
        classifications = self.peak_classifications.get(compound, [])
        for classification in classifications:
            if classification.is_product:
                return classification
        return None

    def get_truncation_peaks(
        self, compound: Compound
    ) -> List[PeakClassification]:
        """
        Get all truncation peak classifications for a compound.

        Parameters
        ----------
        compound : Compound
            Query compound

        Returns
        -------
        List[PeakClassification]
            All TRUNCATION classifications (may be empty)
        """
        classifications = self.peak_classifications.get(compound, [])
        return [c for c in classifications if c.is_truncation]

    def get_validated_compounds(self) -> List[Compound]:
        """
        Get all compounds with VALIDATED status.

        Returns
        -------
        List[Compound]
            Compounds successfully validated
        """
        return [
            compound
            for compound, status in self.validation_results.items()
            if status == ValidationStatus.VALIDATED
        ]

    def get_failed_compounds(self) -> List[Compound]:
        """
        Get all compounds with FAILED status.

        Returns
        -------
        List[Compound]
            Compounds that failed validation
        """
        return [
            compound
            for compound, status in self.validation_results.items()
            if status == ValidationStatus.FAILED
        ]

    def get_validation_summary(self) -> Dict[ValidationStatus, int]:
        """
        Get summary counts of validation statuses.

        Returns
        -------
        Dict[ValidationStatus, int]
            Count of compounds per validation status

        Examples
        --------
        >>> result.get_validation_summary()
        {
            <ValidationStatus.VALIDATED: 'VALIDATED'>: 45,
            <ValidationStatus.LIKELY_SUCCESS: 'LIKELY_SUCCESS'>: 12,
            <ValidationStatus.UNCERTAIN: 'UNCERTAIN'>: 5,
            <ValidationStatus.LIKELY_FAILURE: 'LIKELY_FAILURE'>: 3,
            <ValidationStatus.FAILED: 'FAILED'>: 7
        }
        """
        summary = {}
        for status in ValidationStatus:
            count = sum(1 for s in self.validation_results.values() if s == status)
            if count > 0:
                summary[status] = count
        return summary

    def get_equivalence_class(
        self, block_support_sequence: str
    ) -> Optional[EquivalenceClass]:
        """
        Get equivalence class by block support sequence.

        Parameters
        ----------
        block_support_sequence : str
            Block support sequence (non-null blocks)

        Returns
        -------
        Optional[EquivalenceClass]
            Equivalence class, or None if not found
        """
        for eq_class in self.equivalence_classes:
            if eq_class.block_support_sequence == block_support_sequence:
                return eq_class
        return None

    def total_compounds(self) -> int:
        """
        Get total number of compounds analyzed.

        Returns
        -------
        int
            Number of compounds in hierarchy
        """
        return self.hierarchy.size()

    def total_peaks_detected(self) -> int:
        """
        Get total number of peaks detected across all compounds.

        Returns
        -------
        int
            Total peak count
        """
        return sum(
            len(classifications) for classifications in self.peak_classifications.values()
        )

    def __repr__(self) -> str:
        """Detailed representation for debugging."""
        n_validated = sum(
            1
            for s in self.validation_results.values()
            if s == ValidationStatus.VALIDATED
        )
        return (
            f"AnalysisResult("
            f"compounds={self.total_compounds()}, "
            f"peaks={self.total_peaks_detected()}, "
            f"validated={n_validated}, "
            f"eq_classes={len(self.equivalence_classes)})"
        )
