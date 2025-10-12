"""
Bayesian validation service for synthesis success assessment.

Implementation based on THEORY.md Section 6.7-6.10.
"""

import numpy as np
from typing import Dict, List, Optional
from ..entities.compound import Compound
from ..entities.peak import Peak, PeakType, ValidationStatus
from ..models.compound_hierarchy import CompoundHierarchy
from .purity_calculator import PurityCalculator
from .snr_calculator import SNRCalculator


class BayesianValidator:
    """
    Implements Bayesian synthesis validation framework.

    Validates synthesis success using multiple evidence sources:
    - Purity (product area / total area)
    - Signal-to-noise ratio (peak height / background)
    - Retention time ordering (chromatographic physics)
    - Descendant validation status (DAG constraints)

    This is a stateless service - all methods are operations on input data with
    no instance state.

    Notes
    -----
    The validation framework uses Bayesian inference:

        P(synthesis_succeeded | evidence) ∝
            P(purity | synthesis_succeeded) ×
            P(retention_order | synthesis_succeeded) ×
            P(descendants | synthesis_succeeded) ×
            P(synthesis_succeeded)

    Validation categories:
    - VALIDATED: Very high confidence (>95%)
    - LIKELY_SUCCESS: High confidence (80-95%)
    - UNCERTAIN: Ambiguous evidence
    - LIKELY_FAILURE: High confidence of failure
    - FAILED: Very high confidence of failure

    CRITICAL: This is synthesis validation, NOT peak classification.
    Classification (NULL/TRUNCATION/PUTATIVE_PRODUCT) is a separate step
    that identifies peak positions. Validation assesses synthesis success.

    References
    ----------
    THEORY.md Section 6.7: Bayesian Validation Framework
    THEORY.md Section 6.10: Validation Classification
    THEORY.md Section 6.13: Classification ≠ Validation

    Examples
    --------
    >>> validator = BayesianValidator()
    >>> compound = Compound(...)
    >>> hierarchy = CompoundHierarchy(...)
    >>> dataset_stats = {
    ...     'purity_p25': 0.5,
    ...     'purity_p50': 0.7,
    ...     'purity_p75': 0.85,
    ...     'background': 10.0
    ... }
    >>> status = validator.validate(compound, hierarchy, dataset_stats)
    >>> status
    <ValidationStatus.VALIDATED: 'VALIDATED'>
    """

    def validate(
        self,
        compound: Compound,
        hierarchy: CompoundHierarchy,
        dataset_stats: Dict[str, float],
        retention_precision: float = 0.5
    ) -> ValidationStatus:
        """
        Validate synthesis success for a compound.

        Parameters
        ----------
        compound : Compound
            Compound to validate
        hierarchy : CompoundHierarchy
            DAG structure for constraint checking
        dataset_stats : Dict[str, float]
            Dataset statistics for adaptive thresholds.
            Required keys:
            - 'purity_p25': 25th percentile purity
            - 'purity_p50': median purity
            - 'purity_p75': 75th percentile purity
            - 'background': Background signal level
        retention_precision : float, optional
            Minimum resolvable retention time difference.
            Default is 0.5 (seconds or minutes).

        Returns
        -------
        ValidationStatus
            Validation category (VALIDATED, LIKELY_SUCCESS, etc.)

        Notes
        -----
        Validation requires:
        1. Compound has selected_peak (putative product)
        2. Purity can be calculated from detected peaks
        3. Dataset statistics available for thresholds
        4. Descendant validation status known (if applicable)

        If selected_peak is None, returns NOT_VALIDATED.

        References
        ----------
        THEORY.md Section 6.10: Validation Classification

        Examples
        --------
        >>> validator = BayesianValidator()
        >>> status = validator.validate(compound, hierarchy, dataset_stats)
        """
        # Check if compound has been analyzed
        if compound.selected_peak is None:
            return ValidationStatus.NOT_VALIDATED

        # Calculate purity using shared service
        purity = PurityCalculator.calculate(compound)

        # Calculate SNR using shared service
        background = dataset_stats.get('background', 10.0)
        snr = SNRCalculator.calculate(compound, background)

        # Check retention time ordering
        retention_order_valid = self._check_retention_order(
            compound, hierarchy, retention_precision
        )

        # Get descendant validation
        descendant_fraction = self._get_descendant_validation_fraction(compound, hierarchy)

        # Get thresholds
        p25 = dataset_stats.get('purity_p25', 0.5)
        p50 = dataset_stats.get('purity_p50', 0.7)
        p75 = dataset_stats.get('purity_p75', 0.85)

        # Apply decision logic
        return self._classify_validation_status(
            purity=purity,
            snr=snr,
            retention_order_valid=retention_order_valid,
            descendant_fraction=descendant_fraction,
            p25=p25,
            p50=p50,
            p75=p75
        )

    def _check_retention_order(
        self,
        compound: Compound,
        hierarchy: CompoundHierarchy,
        precision: float
    ) -> bool:
        """
        Check retention time ordering constraint.

        Validates that compound elutes AFTER all its descendants
        (chromatographic physics: longer peptides elute later).

        Parameters
        ----------
        compound : Compound
            Compound to check
        hierarchy : CompoundHierarchy
            DAG structure
        precision : float
            Minimum resolvable time difference

        Returns
        -------
        bool
            True if retention order is valid, False if violated

        Notes
        -----
        Physical law: retention_time(compound) > retention_time(descendant)

        For confident ordering:
            t_compound - t_descendant > 2 × precision

        If difference < precision → ambiguous (return True, but flag uncertainty)

        Retention violations are HARD constraints - if violated, synthesis
        definitely failed OR peak assignment is incorrect.

        References
        ----------
        THEORY.md Section 6.6: Retention Time Constraints
        """
        if compound.selected_peak is None:
            return True  # Cannot check, assume valid

        compound_rt = compound.selected_peak.position

        # Get all descendants
        descendants = hierarchy.get_descendants(compound)

        for desc in descendants:
            if desc.selected_peak is not None:
                desc_rt = desc.selected_peak.position

                # Compound must elute AFTER descendant
                if compound_rt <= desc_rt + precision:
                    # Violation or ambiguous
                    return False

        return True

    def _get_descendant_validation_fraction(
        self,
        compound: Compound,
        hierarchy: CompoundHierarchy
    ) -> float:
        """
        Calculate fraction of descendants that are validated.

        Parameters
        ----------
        compound : Compound
            Compound to check
        hierarchy : CompoundHierarchy
            DAG structure

        Returns
        -------
        float
            Fraction of descendants with VALIDATED status [0, 1]

        Notes
        -----
        If compound has no descendants, returns 1.0 (no constraint).

        Used in Bayesian inference:
        - All descendants validated → increases confidence
        - Any descendant failed → decreases confidence

        References
        ----------
        THEORY.md Section 6.8: DAG Constraint Propagation
        """
        descendants = hierarchy.get_descendants(compound)

        if not descendants:
            return 1.0  # No descendants = no constraint

        validated_count = 0
        total_count = 0

        for desc in descendants:
            if desc.selected_peak is not None:
                total_count += 1
                if desc.selected_peak.validation_status == ValidationStatus.VALIDATED:
                    validated_count += 1

        if total_count == 0:
            return 1.0

        return validated_count / total_count

    def _classify_validation_status(
        self,
        purity: float,
        snr: float,
        retention_order_valid: bool,
        descendant_fraction: float,
        p25: float,
        p50: float,
        p75: float
    ) -> ValidationStatus:
        """
        Classify validation status based on evidence.

        Parameters
        ----------
        purity : float
            Purity value [0, 1]
        snr : float
            Signal-to-noise ratio
        retention_order_valid : bool
            Whether retention order is satisfied
        descendant_fraction : float
            Fraction of validated descendants
        p25, p50, p75 : float
            Dataset percentiles for adaptive thresholds

        Returns
        -------
        ValidationStatus
            Classification result

        Notes
        -----
        Decision framework (THEORY.md Section 6.10):

        VALIDATED: Very high confidence (>95%)
        - Retention order correct
        - Purity > P₇₅
        - SNR > 5
        - All descendants validated

        LIKELY_SUCCESS: High confidence (80-95%)
        - Retention order correct
        - Purity > P₅₀
        - SNR > 3
        - Majority descendants validated

        UNCERTAIN: Ambiguous
        - P₂₅ < Purity < P₇₅ OR
        - SNR ≈ 3 OR
        - Mixed descendant results

        LIKELY_FAILURE: High confidence of failure
        - Purity < P₂₅ OR
        - SNR < 3 OR
        - Retention order ambiguous

        FAILED: Very high confidence of failure
        - Retention order violated OR
        - Purity very low OR
        - SNR < 2

        References
        ----------
        THEORY.md Section 6.10: Validation Classification
        """
        # FAILED: Hard constraints violated
        if not retention_order_valid:
            return ValidationStatus.FAILED

        if snr < 2.0:
            return ValidationStatus.FAILED

        if purity < 0.1:
            return ValidationStatus.FAILED

        # VALIDATED: All evidence positive
        if (purity > p75 and
            snr > 5.0 and
            retention_order_valid and
            descendant_fraction > 0.8):
            return ValidationStatus.VALIDATED

        # LIKELY_SUCCESS: Most evidence positive
        if (purity > p50 and
            snr > 3.0 and
            retention_order_valid and
            descendant_fraction > 0.5):
            return ValidationStatus.LIKELY_SUCCESS

        # LIKELY_FAILURE: Most evidence negative
        if purity < p25 or snr < 3.0:
            return ValidationStatus.LIKELY_FAILURE

        # UNCERTAIN: Mixed or ambiguous evidence
        return ValidationStatus.UNCERTAIN

    def compute_likelihood_ratio(
        self,
        purity: float,
        retention_order_valid: bool,
        descendant_fraction: float
    ) -> float:
        """
        Compute Bayesian likelihood ratio for synthesis success.

        Parameters
        ----------
        purity : float
            Purity value [0, 1]
        retention_order_valid : bool
            Whether retention order is satisfied
        descendant_fraction : float
            Fraction of validated descendants

        Returns
        -------
        float
            Likelihood ratio P(evidence | succeeded) / P(evidence | failed)

        Notes
        -----
        Uses likelihood functions from THEORY.md Section 6.7:

        Purity likelihood:
        - P(purity | succeeded) ~ Beta(α=19, β=1) [mode ≈ 0.95]
        - P(purity | failed) ~ Beta(α=2, β=8) [mode ≈ 0.20]

        Retention order likelihood:
        - P(order_correct | succeeded) = 0.95
        - P(order_correct | failed) = 0.05

        Descendant evidence:
        - P(descendants_validated | succeeded) = 0.90^n
        - P(descendants_validated | failed) = 0.10^n

        Returns likelihood ratio (LR):
        - LR > 10: Strong evidence for success
        - 1 < LR < 10: Moderate evidence for success
        - LR ≈ 1: No evidence either way
        - 0.1 < LR < 1: Moderate evidence for failure
        - LR < 0.1: Strong evidence for failure

        References
        ----------
        THEORY.md Section 6.7: Bayesian Validation Framework
        """
        from scipy.stats import beta

        # Purity likelihood (Beta distributions)
        p_purity_given_success = beta.pdf(purity, a=19, b=1)
        p_purity_given_failure = beta.pdf(purity, a=2, b=8)

        purity_lr = p_purity_given_success / max(p_purity_given_failure, 1e-10)

        # Retention order likelihood
        if retention_order_valid:
            retention_lr = 0.95 / 0.05  # = 19.0
        else:
            retention_lr = 0.05 / 0.95  # = 0.053

        # Descendant evidence likelihood
        # Simplified: use fraction as proxy for number validated
        if descendant_fraction > 0:
            descendant_lr = (0.90 ** descendant_fraction) / (0.10 ** descendant_fraction)
        else:
            descendant_lr = 1.0  # Neutral if no descendants

        # Combined likelihood ratio (multiply independent evidence)
        likelihood_ratio = purity_lr * retention_lr * descendant_lr

        return float(likelihood_ratio)
