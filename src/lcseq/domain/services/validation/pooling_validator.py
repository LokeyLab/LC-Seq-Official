"""
Pooled mode validation with automatic fallback.

Implementation based on THEORY.md Section 4.2.8.1.
"""

import numpy as np
from typing import List, Dict, Optional, Tuple
from ...entities.compound import Compound
from ...entities.chromatogram import Chromatogram


class PoolingValidator:
    """
    Implements pooled mode validation with correlation checking.

    Validates that signal variants are sufficiently correlated (r > 0.8) for
    pooled mode analysis. Provides automatic fallback to individual mode
    if correlation check fails.

    Core Principle: Pooled mode is an optional optimization. If variants
    are heterogeneous (low correlation), automatically fall back to
    individual mode processing.

    Notes
    -----
    Workflow:
    1. Check correlation between all variant pairs
    2. If min(correlation) > 0.8: POOLING_VALID
    3. If min(correlation) ≤ 0.8: HETEROGENEOUS → fallback
    4. Process in individual mode
    5. Aggregate results with class-level status

    Status flags:
    - POOLING_VALID: Correlation check passed
    - HETEROGENEOUS: Correlation check failed, used individual mode
    - POOLING_INVALID_BUT_SIMILAR: Fallback successful

    References
    ----------
    THEORY.md Section 4.2.8.1: Operational Fallback Workflow

    Examples
    --------
    >>> validator = PoolingValidator()
    >>> chromatograms = [chrom1, chrom2, chrom3]  # 3 variants
    >>> is_valid, min_corr = validator.check_correlation(chromatograms)
    >>> is_valid
    True
    >>> min_corr
    0.92
    """

    def __init__(self, correlation_threshold: float = 0.8):
        """
        Initialize pooling validator.

        Parameters
        ----------
        correlation_threshold : float, optional
            Minimum correlation required for pooled mode.
            Default is 0.8 (THEORY.md Section 4.2.8.1).
        """
        self.correlation_threshold = correlation_threshold

    def check_correlation(
        self,
        chromatograms: List[Chromatogram],
        signal_key: str = "corrected"
    ) -> Tuple[bool, float]:
        """
        Check if chromatograms are sufficiently correlated for pooled mode.

        Parameters
        ----------
        chromatograms : List[Chromatogram]
            Signal variants to check
        signal_key : str, optional
            Which signal variant to use (default: "corrected")

        Returns
        -------
        Tuple[bool, float]
            (is_valid, min_correlation)
            - is_valid: True if min(correlation) > threshold
            - min_correlation: Minimum pairwise correlation

        Notes
        -----
        Computes Pearson correlation for all pairs of chromatograms.
        Valid for pooled mode if min(correlation) > 0.8.

        If correlation < 0.8, signals are heterogeneous and should be
        processed individually.

        References
        ----------
        THEORY.md Section 4.2.8.1: Operational Fallback Workflow

        Examples
        --------
        >>> is_valid, min_corr = validator.check_correlation([c1, c2, c3])
        >>> is_valid
        True
        >>> min_corr
        0.92
        """
        if len(chromatograms) < 2:
            # Single variant - always valid for pooling (trivial case)
            return True, 1.0

        # Extract signals
        signals = []
        for chrom in chromatograms:
            signal = chrom.signals.get(signal_key)
            if signal is None:
                raise ValueError(f"Signal '{signal_key}' not found in chromatogram")
            signals.append(signal)

        # Ensure all same length (interpolation should handle this upstream)
        lengths = [len(s) for s in signals]
        if len(set(lengths)) > 1:
            raise ValueError(
                f"Chromatograms have different lengths: {lengths}. "
                "Apply interpolation first."
            )

        # Compute all pairwise correlations
        correlations = []
        n = len(signals)
        for i in range(n):
            for j in range(i + 1, n):
                corr = self._compute_correlation(signals[i], signals[j])
                correlations.append(corr)

        if not correlations:
            return True, 1.0

        min_correlation = float(np.min(correlations))
        is_valid = min_correlation > self.correlation_threshold

        return is_valid, min_correlation

    def _compute_correlation(
        self,
        signal1: np.ndarray,
        signal2: np.ndarray
    ) -> float:
        """
        Compute Pearson correlation between two signals.

        Parameters
        ----------
        signal1, signal2 : np.ndarray
            Signals to correlate (same length)

        Returns
        -------
        float
            Pearson correlation coefficient [-1, 1]

        Notes
        -----
        Uses numpy corrcoef for numerical stability.
        Returns 0.0 if either signal has zero variance.
        """
        # Check for zero variance
        if np.std(signal1) == 0 or np.std(signal2) == 0:
            return 0.0

        # Compute correlation matrix
        corr_matrix = np.corrcoef(signal1, signal2)

        # Extract correlation coefficient
        correlation = corr_matrix[0, 1]

        # Handle NaN (shouldn't happen with variance check, but be safe)
        if np.isnan(correlation):
            return 0.0

        return float(correlation)

    def validate_pooling_eligibility(
        self,
        chromatograms: List[Chromatogram],
        signal_key: str = "corrected"
    ) -> Dict[str, any]:
        """
        Full validation check for pooled mode eligibility.

        Parameters
        ----------
        chromatograms : List[Chromatogram]
            Signal variants to validate
        signal_key : str, optional
            Which signal variant to use

        Returns
        -------
        Dict[str, any]
            Validation report with keys:
            - 'is_eligible': bool
            - 'min_correlation': float
            - 'status': str ('POOLING_VALID', 'HETEROGENEOUS', etc.)
            - 'recommendation': str (human-readable)
            - 'all_correlations': List[float] (all pairwise correlations)

        Examples
        --------
        >>> report = validator.validate_pooling_eligibility(chroms)
        >>> report
        {
            'is_eligible': True,
            'min_correlation': 0.92,
            'status': 'POOLING_VALID',
            'recommendation': 'Use pooled mode - signals are highly correlated',
            'all_correlations': [0.92, 0.95, 0.94]
        }
        """
        is_valid, min_corr = self.check_correlation(chromatograms, signal_key)

        # Extract all pairwise correlations for detailed report
        signals = [chrom.signals.get(signal_key) for chrom in chromatograms]
        all_correlations = []
        n = len(signals)
        for i in range(n):
            for j in range(i + 1, n):
                corr = self._compute_correlation(signals[i], signals[j])
                all_correlations.append(corr)

        # Determine status and recommendation
        if is_valid:
            status = 'POOLING_VALID'
            recommendation = (
                f'Use pooled mode - signals are highly correlated '
                f'(min r = {min_corr:.3f})'
            )
        else:
            status = 'HETEROGENEOUS'
            recommendation = (
                f'Fallback to individual mode - signals are heterogeneous '
                f'(min r = {min_corr:.3f} < {self.correlation_threshold})'
            )

        return {
            'is_eligible': is_valid,
            'min_correlation': min_corr,
            'status': status,
            'recommendation': recommendation,
            'all_correlations': all_correlations
        }

    def aggregate_individual_results(
        self,
        compounds: List[Compound]
    ) -> Dict[str, any]:
        """
        Aggregate validation results from individual mode processing.

        Parameters
        ----------
        compounds : List[Compound]
            Compounds from individual mode (one per variant)

        Returns
        -------
        Dict[str, any]
            Aggregated results with keys:
            - 'pooling_status': str
            - 'variant_count': int
            - 'validated_count': int
            - 'validation_fraction': float

        Notes
        -----
        After fallback to individual mode, aggregate results to provide
        class-level validation status.

        Examples
        --------
        >>> results = validator.aggregate_individual_results(compounds)
        >>> results
        {
            'pooling_status': 'POOLING_INVALID_BUT_SIMILAR',
            'variant_count': 3,
            'validated_count': 3,
            'validation_fraction': 1.0
        }
        """
        from ...entities.peak import ValidationStatus

        variant_count = len(compounds)
        validated_count = sum(
            1 for c in compounds
            if c.selected_peak and
            c.selected_peak.validation_status == ValidationStatus.VALIDATED
        )

        validation_fraction = (
            validated_count / variant_count if variant_count > 0 else 0.0
        )

        # Determine pooling status based on validation fraction
        if validation_fraction >= 0.8:
            pooling_status = 'POOLING_INVALID_BUT_SIMILAR'
        elif validation_fraction >= 0.5:
            pooling_status = 'HETEROGENEOUS_MIXED'
        else:
            pooling_status = 'HETEROGENEOUS_POOR'

        return {
            'pooling_status': pooling_status,
            'variant_count': variant_count,
            'validated_count': validated_count,
            'validation_fraction': validation_fraction
        }
