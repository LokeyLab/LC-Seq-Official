"""
End-to-end validation pipeline orchestrating classification and validation.

Implementation based on THEORY.md Section 6.13.
"""

from typing import List, Dict, Optional
from ...entities.compound import Compound
from ...models.compound_hierarchy import CompoundHierarchy
from ..bayesian_validator import BayesianValidator
from ..purity_calculator import PurityCalculator
from ..snr_calculator import SNRCalculator
from .adaptive_validator import AdaptiveValidator
from .pooling_validator import PoolingValidator


class ValidationWorkflow:
    """
    Coordinates end-to-end validation pipeline.

    Orchestrates:
    1. Dataset statistics computation (adaptive thresholds)
    2. Pooled mode validation (if applicable)
    3. Bayesian validation with adaptive thresholds
    4. Result aggregation

    This is the main entry point for synthesis validation.

    Notes
    -----
    Workflow steps:
    1. Compute dataset-wide statistics
    2. Validate each compound using Bayesian framework
    3. Apply adaptive thresholds
    4. Generate comprehensive validation report

    References
    ----------
    THEORY.md Section 6.13: Validation Workflow

    Examples
    --------
    >>> workflow = ValidationWorkflow()
    >>> hierarchy = CompoundHierarchy(...)
    >>> compounds = [...]  # Analyzed compounds
    >>> results = workflow.validate_library(compounds, hierarchy)
    >>> results['dataset_stats']
    {'purity_p25': 0.5, 'purity_p50': 0.7, 'purity_p75': 0.85, ...}
    """

    def __init__(self):
        """Initialize validation workflow with required validators."""
        self.bayesian_validator = BayesianValidator()
        self.adaptive_validator = AdaptiveValidator()
        self.pooling_validator = PoolingValidator()

    def validate_library(
        self,
        compounds: List[Compound],
        hierarchy: CompoundHierarchy,
        retention_precision: float = 0.5
    ) -> Dict[str, any]:
        """
        Validate entire library with adaptive thresholds.

        Parameters
        ----------
        compounds : List[Compound]
            All compounds in library (must have detected_peaks)
        hierarchy : CompoundHierarchy
            DAG structure for constraint checking
        retention_precision : float, optional
            Minimum resolvable retention time difference

        Returns
        -------
        Dict[str, any]
            Comprehensive validation results with keys:
            - 'dataset_stats': Dataset-wide statistics
            - 'validation_results': Per-compound validation
            - 'summary': Aggregate metrics

        Notes
        -----
        This is the main validation entry point. Coordinates all validation
        services to produce complete validation assessment.

        Examples
        --------
        >>> results = workflow.validate_library(compounds, hierarchy)
        >>> results['summary']['validated_count']
        450
        >>> results['summary']['validation_rate']
        0.75
        """
        # Step 1: Compute dataset statistics
        dataset_stats = self.adaptive_validator.compute_dataset_statistics(compounds)

        # Step 2: Validate each compound
        validation_results = []
        for compound in compounds:
            result = self.validate_compound(
                compound=compound,
                hierarchy=hierarchy,
                dataset_stats=dataset_stats,
                retention_precision=retention_precision
            )
            validation_results.append(result)

        # Step 3: Generate summary
        summary = self._generate_summary(validation_results, dataset_stats)

        return {
            'dataset_stats': dataset_stats,
            'validation_results': validation_results,
            'summary': summary
        }

    def validate_compound(
        self,
        compound: Compound,
        hierarchy: CompoundHierarchy,
        dataset_stats: Dict[str, float],
        retention_precision: float = 0.5
    ) -> Dict[str, any]:
        """
        Validate single compound with full diagnostics.

        Parameters
        ----------
        compound : Compound
            Compound to validate
        hierarchy : CompoundHierarchy
            DAG structure
        dataset_stats : Dict[str, float]
            Dataset statistics from adaptive_validator
        retention_precision : float, optional
            Retention time precision

        Returns
        -------
        Dict[str, any]
            Detailed validation result with keys:
            - 'compound_id': Compound identifier
            - 'validation_status': ValidationStatus
            - 'purity': float
            - 'purity_category': str (adaptive category)
            - 'purity_ci': Tuple[float, float] (confidence interval)
            - 'snr': float
            - 'retention_order_valid': bool
            - 'descendant_fraction': float
            - 'likelihood_ratio': float

        Examples
        --------
        >>> result = workflow.validate_compound(compound, hierarchy, stats)
        >>> result
        {
            'compound_id': 'compound_123',
            'validation_status': <ValidationStatus.VALIDATED>,
            'purity': 0.85,
            'purity_category': 'high',
            'snr': 12.5,
            ...
        }
        """
        from ...entities.peak import ValidationStatus

        # Validate using Bayesian framework
        validation_status = self.bayesian_validator.validate(
            compound=compound,
            hierarchy=hierarchy,
            dataset_stats=dataset_stats,
            retention_precision=retention_precision
        )

        # Calculate purity using shared service
        purity = PurityCalculator.calculate(compound)

        # Get adaptive category
        purity_category = self.adaptive_validator.get_adaptive_category(
            purity, dataset_stats
        )

        # Calculate confidence interval
        total_counts = sum(peak.area for peak in compound.detected_peaks)
        purity_ci = self.adaptive_validator.compute_purity_confidence_interval(
            purity, total_counts
        )

        # Calculate SNR using shared service
        background = dataset_stats.get('background', 10.0)
        snr = SNRCalculator.calculate(compound, background)

        # Check retention order
        retention_order_valid = self.bayesian_validator._check_retention_order(
            compound, hierarchy, retention_precision
        )

        # Get descendant validation
        descendant_fraction = self.bayesian_validator._get_descendant_validation_fraction(
            compound, hierarchy
        )

        # Compute likelihood ratio
        likelihood_ratio = self.bayesian_validator.compute_likelihood_ratio(
            purity=purity,
            retention_order_valid=retention_order_valid,
            descendant_fraction=descendant_fraction
        )

        return {
            'compound_id': str(compound),
            'validation_status': validation_status,
            'purity': purity,
            'purity_category': purity_category,
            'purity_ci': purity_ci,
            'snr': snr,
            'retention_order_valid': retention_order_valid,
            'descendant_fraction': descendant_fraction,
            'likelihood_ratio': likelihood_ratio
        }

    def _generate_summary(
        self,
        validation_results: List[Dict[str, any]],
        dataset_stats: Dict[str, float]
    ) -> Dict[str, any]:
        """
        Generate summary statistics from validation results.

        Parameters
        ----------
        validation_results : List[Dict[str, any]]
            Per-compound validation results
        dataset_stats : Dict[str, float]
            Dataset statistics

        Returns
        -------
        Dict[str, any]
            Summary with keys:
            - 'total_compounds': int
            - 'validated_count': int
            - 'likely_success_count': int
            - 'uncertain_count': int
            - 'likely_failure_count': int
            - 'failed_count': int
            - 'validation_rate': float
            - 'stringency_level': str
        """
        from ...entities.peak import ValidationStatus

        total = len(validation_results)
        if total == 0:
            return {
                'total_compounds': 0,
                'validated_count': 0,
                'likely_success_count': 0,
                'uncertain_count': 0,
                'likely_failure_count': 0,
                'failed_count': 0,
                'validation_rate': 0.0,
                'stringency_level': 'unknown'
            }

        # Count by status
        status_counts = {
            ValidationStatus.VALIDATED: 0,
            ValidationStatus.LIKELY_SUCCESS: 0,
            ValidationStatus.UNCERTAIN: 0,
            ValidationStatus.LIKELY_FAILURE: 0,
            ValidationStatus.FAILED: 0
        }

        for result in validation_results:
            status = result['validation_status']
            if status in status_counts:
                status_counts[status] += 1

        # Get stringency level
        _, stringency = self.adaptive_validator.get_validation_stringency(dataset_stats)

        # Calculate validation rate (VALIDATED + LIKELY_SUCCESS)
        success_count = (
            status_counts[ValidationStatus.VALIDATED] +
            status_counts[ValidationStatus.LIKELY_SUCCESS]
        )
        validation_rate = success_count / total

        return {
            'total_compounds': total,
            'validated_count': status_counts[ValidationStatus.VALIDATED],
            'likely_success_count': status_counts[ValidationStatus.LIKELY_SUCCESS],
            'uncertain_count': status_counts[ValidationStatus.UNCERTAIN],
            'likely_failure_count': status_counts[ValidationStatus.LIKELY_FAILURE],
            'failed_count': status_counts[ValidationStatus.FAILED],
            'validation_rate': validation_rate,
            'stringency_level': stringency
        }

    def validate_with_pooling_check(
        self,
        compounds: List[Compound],
        hierarchy: CompoundHierarchy,
        retention_precision: float = 0.5
    ) -> Dict[str, any]:
        """
        Validate library with pooled mode checking.

        Parameters
        ----------
        compounds : List[Compound]
            Compounds (may include variants)
        hierarchy : CompoundHierarchy
            DAG structure
        retention_precision : float, optional
            Retention time precision

        Returns
        -------
        Dict[str, any]
            Validation results with pooling status

        Notes
        -----
        This method should be used when compounds include multiple variants
        that need pooling checking. For single variants, use validate_library().

        Examples
        --------
        >>> results = workflow.validate_with_pooling_check(compounds, hierarchy)
        >>> results['pooling_report']
        {'status': 'POOLING_VALID', 'min_correlation': 0.92, ...}
        """
        # Standard validation
        validation_results = self.validate_library(
            compounds, hierarchy, retention_precision
        )

        # Add pooling information (if applicable)
        # Note: This is a placeholder - full pooling integration would require
        # chromatogram access at this level
        validation_results['pooling_report'] = {
            'status': 'NOT_APPLICABLE',
            'note': 'Pooling checking requires chromatogram access'
        }

        return validation_results
