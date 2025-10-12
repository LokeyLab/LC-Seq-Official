"""
Adaptive validation with dataset-relative thresholds.

Implementation based on THEORY.md Section 6.2-6.4.
"""

import numpy as np
from typing import Dict, List, Tuple
from ...entities.compound import Compound
from ...entities.peak import ValidationStatus
from ..purity_calculator import PurityCalculator


class AdaptiveValidator:
    """
    Implements adaptive validation with dataset-relative thresholds.

    Computes dataset-wide statistics (P25, P50, P75 for purity) and applies
    likelihood ratio thresholds adaptively based on data distribution.

    Core Principle: All validation metrics must be dataset-relative and
    scale-invariant to accommodate different sequencing depths, scaling
    methods, and library quality levels.

    Notes
    -----
    The adaptive approach:
    1. Learn dataset distribution (bootstrap phase)
    2. Define adaptive thresholds (percentiles, MAD)
    3. Classify relative to dataset characteristics
    4. Adjust stringency based on data quality

    NO fixed thresholds - all parameters derived from data distribution.

    References
    ----------
    THEORY.md Section 6.2: Adaptive Validation Principle
    THEORY.md Section 6.4: Distribution-Based Thresholds

    Examples
    --------
    >>> validator = AdaptiveValidator()
    >>> compounds = [...]  # List of analyzed compounds
    >>> stats = validator.compute_dataset_statistics(compounds)
    >>> stats
    {
        'purity_p10': 0.3,
        'purity_p25': 0.5,
        'purity_p50': 0.7,
        'purity_p75': 0.85,
        'purity_p90': 0.95,
        'purity_mad': 0.15,
        'background': 10.0
    }
    """

    def compute_dataset_statistics(
        self,
        compounds: List[Compound]
    ) -> Dict[str, float]:
        """
        Compute dataset-wide statistics for adaptive thresholds.

        Parameters
        ----------
        compounds : List[Compound]
            All compounds in dataset (must have detected_peaks)

        Returns
        -------
        Dict[str, float]
            Dataset statistics with keys:
            - 'purity_p10', 'purity_p25', 'purity_p50', 'purity_p75', 'purity_p90'
            - 'purity_mad': Median absolute deviation
            - 'background': Background signal level

        Notes
        -----
        Computes:
        1. Purity for each compound
        2. Percentiles (P10, P25, P50, P75, P90)
        3. MAD (median absolute deviation)
        4. Background estimate from low-count tail

        References
        ----------
        THEORY.md Section 6.4: Distribution-Based Thresholds

        Examples
        --------
        >>> stats = validator.compute_dataset_statistics(compounds)
        >>> stats['purity_p50']  # Median purity
        0.7
        >>> stats['purity_mad']  # Spread measure
        0.15
        """
        # Compute purity for all compounds
        purities = []
        all_counts = []

        for compound in compounds:
            if not compound.detected_peaks:
                continue

            purity = PurityCalculator.calculate(compound)
            purities.append(purity)

            # Collect all peak heights for background estimation
            for peak in compound.detected_peaks:
                all_counts.append(peak.height)

        if not purities:
            # No valid compounds - return defaults
            return {
                'purity_p10': 0.0,
                'purity_p25': 0.0,
                'purity_p50': 0.0,
                'purity_p75': 0.0,
                'purity_p90': 0.0,
                'purity_mad': 0.0,
                'background': 1.0
            }

        purities_array = np.array(purities)

        # Compute percentiles
        p10 = float(np.percentile(purities_array, 10))
        p25 = float(np.percentile(purities_array, 25))
        p50 = float(np.percentile(purities_array, 50))
        p75 = float(np.percentile(purities_array, 75))
        p90 = float(np.percentile(purities_array, 90))

        # Compute MAD (Median Absolute Deviation)
        mad = float(np.median(np.abs(purities_array - p50)))

        # Estimate background from low-count tail
        if all_counts:
            counts_array = np.array(all_counts)
            background = float(np.percentile(counts_array, 10))
        else:
            background = 1.0

        return {
            'purity_p10': p10,
            'purity_p25': p25,
            'purity_p50': p50,
            'purity_p75': p75,
            'purity_p90': p90,
            'purity_mad': mad,
            'background': max(background, 1.0)  # Ensure > 0
        }

    def get_adaptive_category(
        self,
        purity: float,
        dataset_stats: Dict[str, float]
    ) -> str:
        """
        Classify purity into adaptive category based on dataset distribution.

        Parameters
        ----------
        purity : float
            Purity value [0, 1]
        dataset_stats : Dict[str, float]
            Dataset statistics from compute_dataset_statistics()

        Returns
        -------
        str
            Category: 'exceptional', 'high', 'moderate', 'low', 'very_low'

        Notes
        -----
        Categories (dataset-relative):
        - Exceptional: purity > P90 (top 10%)
        - High: P75 < purity ≤ P90
        - Moderate: P50 < purity ≤ P75
        - Low: P25 < purity ≤ P50
        - Very low: purity ≤ P25

        References
        ----------
        THEORY.md Section 6.4: Distribution-Based Thresholds

        Examples
        --------
        >>> category = validator.get_adaptive_category(0.85, stats)
        >>> category
        'high'
        """
        p25 = dataset_stats.get('purity_p25', 0.5)
        p50 = dataset_stats.get('purity_p50', 0.7)
        p75 = dataset_stats.get('purity_p75', 0.85)
        p90 = dataset_stats.get('purity_p90', 0.95)

        if purity > p90:
            return 'exceptional'
        elif purity > p75:
            return 'high'
        elif purity > p50:
            return 'moderate'
        elif purity > p25:
            return 'low'
        else:
            return 'very_low'

    def get_validation_stringency(
        self,
        dataset_stats: Dict[str, float]
    ) -> Tuple[float, str]:
        """
        Determine validation stringency based on dataset quality.

        Parameters
        ----------
        dataset_stats : Dict[str, float]
            Dataset statistics

        Returns
        -------
        Tuple[float, str]
            (threshold_percentile, stringency_level)
            - threshold_percentile: Which percentile to use (50 or 75)
            - stringency_level: 'strict' or 'lenient'

        Notes
        -----
        If MAD(purity) < 0.1: High-quality library → strict (P75)
        If MAD(purity) > 0.2: Variable library → lenient (P50)

        This ensures fair evaluation regardless of library-wide quality.

        References
        ----------
        THEORY.md Section 6.4: Distribution-Based Thresholds

        Examples
        --------
        >>> threshold, stringency = validator.get_validation_stringency(stats)
        >>> threshold
        75.0
        >>> stringency
        'strict'
        """
        mad = dataset_stats.get('purity_mad', 0.15)

        if mad < 0.1:
            # High-quality library - use strict thresholds
            return 75.0, 'strict'
        elif mad > 0.2:
            # Variable library - use lenient thresholds
            return 50.0, 'lenient'
        else:
            # Moderate quality - use P75 but with moderate stringency
            return 75.0, 'moderate'

    def compute_purity_confidence_interval(
        self,
        purity: float,
        total_counts: float
    ) -> Tuple[float, float]:
        """
        Compute 95% confidence interval for purity estimate.

        Parameters
        ----------
        purity : float
            Purity estimate [0, 1]
        total_counts : float
            Total scaled counts

        Returns
        -------
        Tuple[float, float]
            (lower_bound, upper_bound) for 95% CI

        Notes
        -----
        Standard error:
        SE(purity) = √[purity × (1-purity) / total_counts]

        95% CI:
        CI = purity ± 1.96 × SE(purity)

        Minimum count threshold: total_counts > 100 for CI width < 0.2

        References
        ----------
        THEORY.md Section 6.3: Purity Definition and Calculation

        Examples
        --------
        >>> lower, upper = validator.compute_purity_confidence_interval(0.7, 500)
        >>> (lower, upper)
        (0.66, 0.74)
        """
        if total_counts <= 0:
            return (0.0, 1.0)

        # Standard error
        se = np.sqrt(purity * (1 - purity) / total_counts)

        # 95% CI
        lower = max(0.0, purity - 1.96 * se)
        upper = min(1.0, purity + 1.96 * se)

        return (float(lower), float(upper))
