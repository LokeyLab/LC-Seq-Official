"""
ComputeGlobalScalesUseCase - Compute uniform scale parameters for persistent homology.

Implementation based on THEORY.md Section 5.7.1 (lines 2212-2225).

CRITICAL: All parameters must be uniform across the entire dataset for comparability.
"""

import numpy as np
from typing import List, Dict
from ...domain.entities.compound import Compound


class ComputeGlobalScalesUseCase:
    """
    Compute uniform scale parameters for entire dataset.

    Per THEORY.md line 2210: "All parameters must be uniform across the
    entire dataset for comparability."

    This use case computes the scale-space filtration parameters (σ_min, σ_max, num_scales)
    that will be used consistently for peak detection across all compounds in the dataset.

    Notes
    -----
    - σ_min is derived from actual sampling intervals (data-driven)
    - σ_max = σ_min × 50 (covers 2 orders of magnitude in smoothing)
    - num_scales = 20 (linearly spaced scales)
    - These parameters must be computed ONCE and used for ALL compounds

    References
    ----------
    THEORY.md Section 5.7.1: Uniformity Requirements
    THEORY.md Section 5.2: Persistent Homology for Peak Significance

    Examples
    --------
    >>> use_case = ComputeGlobalScalesUseCase()
    >>> scales = use_case.execute(all_compounds)
    >>> scales
    {'sigma_min': 0.5, 'sigma_max': 25.0, 'num_scales': 20}
    """

    def execute(self, compounds: List[Compound]) -> Dict[str, float]:
        """
        Compute global scale parameters from dataset.

        Parameters
        ----------
        compounds : List[Compound]
            All compounds in the dataset

        Returns
        -------
        Dict[str, float]
            Dictionary with keys:
            - 'sigma_min': Minimum smoothing scale (in array index units)
            - 'sigma_max': Maximum smoothing scale (in array index units)
            - 'num_scales': Number of scales in filtration (20)

        Notes
        -----
        Algorithm (THEORY.md lines 2215-2225):

        IMPORTANT: Gaussian filter operates on array INDICES, not time units.

        Step 1: Determine minimum smoothing scale
          σ_min = 0.5 (smooths over immediate neighbors)
          This is data-driven: small enough to preserve sharp peaks

        Step 2: Set maximum scale based on signal length
          σ_max = median_signal_length / 10
          Smooths over ~10% of signal, enough to identify broad features
          Typically ~10-20 for signals with ~100 points

        Step 3: Define scale sequence
          num_scales = 20
          scales = linspace(σ_min, σ_max, num_scales)
          (Linearly spaced smoothing scales in index units)

        Rationale:
        - ✅ Uniform: Same scales for every compound in dataset
        - ✅ Data-driven: σ_max based on actual signal lengths
        - ✅ No circular dependency: Doesn't require peak detection first
        - ✅ Simple: One computation, apply to all
        - ✅ Index-based: Works with scipy.ndimage.gaussian_filter1d
        """
        if not compounds:
            raise ValueError("Cannot compute scales for empty compound list")

        # Step 1: Compute median signal length across dataset
        signal_lengths = []
        for compound in compounds:
            n_points = len(compound.chromatogram.time_points)
            if n_points >= 3:
                signal_lengths.append(n_points)

        if not signal_lengths:
            raise ValueError("No valid signals found in dataset")

        median_length = float(np.median(signal_lengths))

        # Step 2: Set index-based smoothing scales
        # σ_min: Small enough to preserve sharp peaks (0.5 array indices)
        sigma_min = 0.5

        # σ_max: Smooths over ~10% of signal (data-driven based on signal length)
        # For ~96 point signals, this gives σ_max ≈ 9.6
        sigma_max = median_length / 10.0

        # Step 3: Define number of scales (THEORY.md line 2223)
        num_scales = 20

        return {
            'sigma_min': sigma_min,
            'sigma_max': sigma_max,
            'num_scales': num_scales
        }
