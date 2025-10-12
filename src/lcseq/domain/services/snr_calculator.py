"""
Signal-to-noise ratio (SNR) calculation service for synthesis validation.

This service provides a single implementation of SNR calculation used by
all validation services (BayesianValidator, ValidationWorkflow, etc.).

Implementation based on THEORY.md Section 6.5.
"""

from ..entities.compound import Compound


class SNRCalculator:
    """
    Domain service for calculating signal-to-noise ratio.

    SNR is a universal scale-invariant metric that quantifies whether a peak
    is significantly above the background noise floor.

    This is a stateless service with a single static method to ensure
    consistent SNR calculation across all validation workflows.

    Notes
    -----
    SNR = peak_height / background_level

    Interpretation (THEORY.md Section 6.5):
    - SNR > 10: High confidence signal
    - 3 ≤ SNR ≤ 10: Moderate confidence signal
    - SNR < 3: Near noise floor (low confidence)

    SNR is scale-invariant because both numerator and denominator
    are in the same units and scale proportionally.

    References
    ----------
    THEORY.md Section 6.5: Signal-to-Noise Ratio

    Examples
    --------
    >>> from lcseq.domain.services import SNRCalculator
    >>> snr = SNRCalculator.calculate(compound, background=10.0)
    >>> snr
    12.5
    """

    @staticmethod
    def calculate(compound: Compound, background: float) -> float:
        """
        Calculate signal-to-noise ratio for compound's selected peak.

        Parameters
        ----------
        compound : Compound
            Compound with selected_peak populated (putative product peak)
        background : float
            Background signal level (typically from dataset P10 or P25)

        Returns
        -------
        float
            SNR = peak_height / background

        Notes
        -----
        Calculation:
        1. Get height of selected_peak (putative product)
        2. Divide by background level
        3. Return ratio

        Edge cases:
        - No selected peak → SNR = 0.0
        - Background ≤ 0 → SNR = 0.0 (invalid background)
        - Valid peak and background → SNR = height / background

        The background level should be obtained from dataset-wide
        statistics (e.g., P10 of all peak heights) to ensure
        scale-invariance and comparability across experiments.

        References
        ----------
        THEORY.md Section 6.5: Signal-to-Noise Ratio
        THEORY.md Section 6.4: Distribution-Based Thresholds

        Examples
        --------
        >>> compound = Compound(...)
        >>> compound.selected_peak = Peak(height=125.0, ...)
        >>> SNRCalculator.calculate(compound, background=10.0)
        12.5

        >>> # Low SNR (near noise floor)
        >>> compound.selected_peak = Peak(height=25.0, ...)
        >>> SNRCalculator.calculate(compound, background=10.0)
        2.5
        """
        if compound.selected_peak is None or background <= 0:
            return 0.0

        snr = compound.selected_peak.height / background
        return float(snr)
