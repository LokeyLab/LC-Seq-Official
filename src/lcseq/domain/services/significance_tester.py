"""
Significance testing service using Negative Binomial distribution.

Provides statistical significance testing for peak detection using the
Negative Binomial distribution, which properly handles overdispersion
(Var > Mean) common in scaled/normalized count data.

The NB distribution naturally reduces to Poisson when no overdispersion
is present (dispersion parameter r → ∞).

References
----------
THEORY.md Section 5.2: Statistical Significance Testing

Robinson, M.D. and Smyth, G.K. (2007). Moderated statistical tests for
assessing differences in tag abundance. Bioinformatics 23(21):2881-2887.

Anders, S. and Huber, W. (2010). Differential expression analysis for
sequence count data. Genome Biology 11:R106.
"""

from typing import Tuple
import numpy as np
from scipy import stats


class SignificanceTesterService:
    """
    Test statistical significance of peaks above background using Negative Binomial.

    The Negative Binomial distribution is the most principled choice for count data
    because it:
    1. Properly handles overdispersion (Var > Mean)
    2. Reduces to Poisson when no overdispersion (r → ∞)
    3. Uses exact CDF (no Normal approximation)
    4. Is the standard in genomics (DESeq2, edgeR)

    The user specifies a false positive rate (alpha) rather than an
    arbitrary threshold, making the choice statistically principled.

    Notes
    -----
    For Negative Binomial with mean μ and dispersion r:
    - Variance = μ + μ²/r (quadratic variance function)
    - When r → ∞: Var → μ (Poisson limit)
    - When r is small: high overdispersion

    The NB arises naturally as a Poisson-Gamma mixture:
    if λ ~ Gamma(r, r/μ) and X | λ ~ Poisson(λ), then X ~ NB(r, r/(r+μ)).

    References
    ----------
    THEORY.md Section 5.2: Statistical Significance Testing
    """

    def test(
        self,
        observed: float,
        background: float,
        alpha: float,
        dispersion: float,
    ) -> Tuple[bool, float]:
        """
        Test if observed count is significantly above background.

        Parameters
        ----------
        observed : float
            Observed count (peak height)
        background : float
            Expected background level (baseline)
        alpha : float
            Significance level (false positive rate).
            Default 0.001 = 0.1% false positive rate per test.
        dispersion : float
            Negative Binomial dispersion parameter r.
            Larger values = closer to Poisson (less overdispersion).
            Default 1e6 effectively gives Poisson behavior.

        Returns
        -------
        Tuple[bool, float]
            (is_significant, p_value) where:
            - is_significant: True if p_value < alpha
            - p_value: Probability of observing this count or higher
              under the null hypothesis (just background)

        Notes
        -----
        This is a one-tailed test: we only care if the count is
        significantly ABOVE background (not below).

        The test uses the Negative Binomial distribution:
        P(X >= observed | NB(μ = background, r = dispersion))

        scipy.stats.nbinom uses (n, p) parameterization where:
        - n = r (dispersion parameter)
        - p = r / (r + μ) (success probability)

        This ensures:
        - Mean = μ
        - Var = μ + μ²/r
        """
        if background <= 0:
            # Can't compute significance with zero/negative background
            # Treat any positive observation as significant
            return observed > 0, 0.0 if observed > 0 else 1.0

        # For continuous observed values, use floor for discrete NB
        k = int(np.floor(observed))

        if k <= 0:
            # Zero or negative counts are never significant above background
            return False, 1.0

        # Convert to scipy's (n, p) parameterization
        r = dispersion
        mu = background
        p = r / (r + mu)  # Success probability

        # P(X >= k) = 1 - P(X <= k-1) = sf(k-1)
        # Using survival function for numerical stability at tails
        p_value = stats.nbinom.sf(k - 1, n=r, p=p)

        return p_value < alpha, float(p_value)

    def compute_p_value(
        self, observed: float, background: float, dispersion: float
    ) -> float:
        """
        Compute p-value for observed count given background.

        Parameters
        ----------
        observed : float
            Observed count (peak height)
        background : float
            Expected background level (baseline)
        dispersion : float
            NB dispersion parameter r

        Returns
        -------
        float
            P-value: probability of observing this count or higher
            under the null hypothesis
        """
        _, p_value = self.test(observed, background, alpha=1.0, dispersion=dispersion)
        return p_value

    def compute_threshold(
        self, background: float, alpha: float, dispersion: float
    ) -> float:
        """
        Compute the minimum count needed for significance at given alpha.

        Parameters
        ----------
        background : float
            Expected background level
        alpha : float
            Significance level
        dispersion : float
            NB dispersion parameter r

        Returns
        -------
        float
            Minimum count that would be significant at level alpha

        Notes
        -----
        This is the inverse of the p-value calculation:
        find X such that P(X >= x | NB(μ, r)) < alpha
        """
        if background <= 0:
            return 0.0

        # Find smallest k such that P(X >= k | NB(μ, r)) < alpha
        # scipy.stats.nbinom uses (n, p) parameterization:
        # n = r, p = r / (r + μ)
        r = dispersion
        p = r / (r + background)
        return float(stats.nbinom.ppf(1 - alpha, n=r, p=p) + 1)
