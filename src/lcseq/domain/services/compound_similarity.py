"""Compound chromatographic similarity analysis domain service.

This module provides domain services for computing chromatographic similarity
between compounds using various distance metrics on their chromatogram signals.

Implementation based on THEORY.md Section 4.4 (Similarity Analysis).

Domain Responsibility:
    - Compute chromatographic similarity/distance metrics between compounds
    - Provide various distance metrics (Wasserstein, Euclidean, etc.)
    - Signal preprocessing for similarity computation
    - Pure algorithmic logic with no presentation concerns

Not Responsible For:
    - Sorting compounds for display (that's presentation)
    - Hierarchical clustering for visualization
    - Visualization or plotting
    - File I/O or data loading

References:
    THEORY.md Section 4.4: Similarity Analysis
    Wasserstein Distance: https://en.wikipedia.org/wiki/Wasserstein_metric
    Also known as Earth Mover's Distance (EMD)
"""

import numpy as np
from numpy.typing import NDArray
from typing import List
from scipy.stats import wasserstein_distance

from ..entities import Compound


class CompoundSimilarityAnalyzer:
    """Domain service for computing chromatographic similarity between compounds.

    This service provides algorithms for computing distance/similarity metrics
    between compound chromatograms using various statistical measures. All methods
    are deterministic and side-effect free.

    This is a stateless service - all configuration is passed as method parameters.

    Examples:
        >>> analyzer = CompoundSimilarityAnalyzer()
        >>> distance = analyzer.compute_pairwise_distance(
        ...     compound1, compound2, metric="wasserstein"
        ... )
        >>> # distance is a float representing chromatographic dissimilarity
        >>>
        >>> matrix = analyzer.compute_similarity_matrix(
        ...     compounds, metric="wasserstein"
        ... )
        >>> # matrix[i][j] is distance between compounds[i] and compounds[j]
    """

    def compute_pairwise_distance(
        self,
        compound1: Compound,
        compound2: Compound,
        metric: str = "wasserstein",
        signal_variant: str = "raw"
    ) -> float:
        """Compute distance between two compounds' chromatograms.

        This computes the distance between two chromatograms using the specified
        metric. The Wasserstein metric measures the minimum "work" required to
        transform one probability distribution into another, capturing similarity
        in elution profiles.

        Parameters
        ----------
        compound1 : Compound
            First compound
        compound2 : Compound
            Second compound
        metric : str, optional
            Distance metric to use. Options:
            - 'wasserstein': Wasserstein distance (Earth Mover's Distance)
            - 'euclidean': Euclidean distance (requires same-length signals)
            Default is 'wasserstein'.
        signal_variant : str, optional
            Which chromatogram signal variant to use. Default is 'raw' (per THEORY.md - no baseline correction).

        Returns
        -------
        float
            Non-negative float representing distance (0 = identical, higher = more different)
            Returns 1.0 if either signal is invalid (all zeros, negative, etc.)

        Raises
        ------
        ValueError
            If metric is not recognized

        Notes
        -----
        - Wasserstein distance treats signals as probability distributions
        - Euclidean distance requires equal-length signals
        - Invalid signals return fallback distance of 1.0

        References
        ----------
        THEORY.md Section 4.4: Similarity Metrics

        Examples
        --------
        >>> analyzer = CompoundSimilarityAnalyzer()
        >>> dist = analyzer.compute_pairwise_distance(cpd1, cpd2, metric="wasserstein")
        >>> # dist is in range [0, +inf), typically [0, 1] for normalized signals
        """
        if metric not in ["wasserstein", "euclidean"]:
            raise ValueError(f"Unknown metric: {metric}. Options: 'wasserstein', 'euclidean'")

        if metric == "wasserstein":
            return self._compute_wasserstein_distance(compound1, compound2, signal_variant)
        elif metric == "euclidean":
            return self._compute_euclidean_distance(compound1, compound2, signal_variant)
        else:
            raise ValueError(f"Unknown metric: {metric}")

    def _compute_wasserstein_distance(
        self,
        compound1: Compound,
        compound2: Compound,
        signal_variant: str
    ) -> float:
        """Compute Wasserstein distance between two chromatograms.

        The Wasserstein distance treats the chromatogram signal as a probability
        distribution over retention time. It measures how much "mass" needs to be
        moved (and how far) to transform one distribution into the other.

        Algorithm:
            1. Get time points and signal intensities for both compounds
            2. Normalize signals to be non-negative
            3. Normalize signals to sum to 1 (probability distribution)
            4. Compute Wasserstein distance using scipy
            5. Handle edge cases (zero signals, invalid data)

        Args:
            compound1: First compound
            compound2: Second compound

        Returns:
            Wasserstein distance in range [0, +inf)
            Returns 1.0 for invalid signals as a fallback
        """
        # Get signals
        signal1 = compound1.chromatogram.get_signal(signal_variant)
        signal2 = compound2.chromatogram.get_signal(signal_variant)
        time1 = compound1.chromatogram.time_points
        time2 = compound2.chromatogram.time_points

        # Normalize signals to be non-negative
        # (subtract minimum to shift any negative values to zero)
        signal1_normalized = np.maximum(signal1 - signal1.min(), 0)
        signal2_normalized = np.maximum(signal2 - signal2.min(), 0)

        # Check if signals are valid (non-zero sum)
        sum1 = signal1_normalized.sum()
        sum2 = signal2_normalized.sum()

        if sum1 <= 0 or sum2 <= 0:
            # Invalid signals - return maximum distance
            return 1.0

        # Normalize to probability distributions (sum = 1)
        prob1 = signal1_normalized / sum1
        prob2 = signal2_normalized / sum2

        # Compute Wasserstein distance
        try:
            distance = wasserstein_distance(time1, time2, prob1, prob2)
            return float(distance)
        except (ValueError, RuntimeError):
            # Numerical issues - return fallback
            return 1.0

    def _compute_euclidean_distance(
        self,
        compound1: Compound,
        compound2: Compound,
        signal_variant: str
    ) -> float:
        """Compute Euclidean distance between two chromatograms.

        This requires both chromatograms to have the same length and time points.
        It measures the direct point-by-point difference between signals.

        Args:
            compound1: First compound
            compound2: Second compound

        Returns:
            Euclidean distance, or raises ValueError if lengths don't match

        Raises:
            ValueError: If chromatograms have different lengths
        """
        signal1 = compound1.chromatogram.get_signal(signal_variant)
        signal2 = compound2.chromatogram.get_signal(signal_variant)

        if len(signal1) != len(signal2):
            raise ValueError(
                f"Euclidean distance requires equal-length signals: "
                f"{len(signal1)} vs {len(signal2)}"
            )

        # Normalize signals
        signal1_norm = np.maximum(signal1 - signal1.min(), 0)
        signal2_norm = np.maximum(signal2 - signal2.min(), 0)

        if signal1_norm.sum() > 0:
            signal1_norm = signal1_norm / signal1_norm.sum()
        if signal2_norm.sum() > 0:
            signal2_norm = signal2_norm / signal2_norm.sum()

        # Compute Euclidean distance
        return float(np.linalg.norm(signal1_norm - signal2_norm))

    def compute_similarity_matrix(
        self,
        compounds: List[Compound],
        metric: str = "wasserstein",
        signal_variant: str = "raw"
    ) -> NDArray[np.float64]:
        """Compute pairwise distance matrix for all compounds.

        This creates an n×n symmetric matrix where matrix[i][j] is the distance
        between compounds[i] and compounds[j]. The diagonal is all zeros (compounds
        are identical to themselves).

        Parameters
        ----------
        compounds : List[Compound]
            List of compounds to compare
        metric : str, optional
            Distance metric to use. Default is 'wasserstein'.
        signal_variant : str, optional
            Signal variant to use. Default is 'raw' (per THEORY.md - no baseline correction).

        Returns
        -------
        NDArray[np.float64]
            Symmetric n×n numpy array of distances
            matrix[i][j] = matrix[j][i] = distance between compounds i and j

        Raises
        ------
        ValueError
            If compounds list is empty

        Notes
        -----
        - Matrix is symmetric: matrix[i][j] == matrix[j][i]
        - Diagonal is zero: matrix[i][i] == 0
        - Only upper triangle is computed for efficiency

        Examples
        --------
        >>> analyzer = CompoundSimilarityAnalyzer()
        >>> compounds = [cpd1, cpd2, cpd3]
        >>> matrix = analyzer.compute_similarity_matrix(compounds, metric="wasserstein")
        >>> matrix.shape
        (3, 3)
        >>> matrix[0][1] == matrix[1][0]  # Symmetric
        True
        >>> matrix[0][0]  # Diagonal is zero
        0.0
        """
        if not compounds:
            raise ValueError("Compounds list cannot be empty")

        n = len(compounds)
        distance_matrix = np.zeros((n, n), dtype=np.float64)

        # Compute upper triangle (matrix is symmetric)
        for i in range(n):
            for j in range(i + 1, n):
                dist = self.compute_pairwise_distance(
                    compounds[i], compounds[j], metric=metric, signal_variant=signal_variant
                )
                distance_matrix[i, j] = dist
                distance_matrix[j, i] = dist  # Symmetric

        return distance_matrix

    def compute_similarity_matrix_robust(
        self,
        compounds: List[Compound],
        metric: str = "wasserstein",
        signal_variant: str = "raw"
    ) -> NDArray[np.float64]:
        """Compute pairwise distance matrix with additional robustness checks.

        This is like compute_similarity_matrix() but with additional numerical
        stability guarantees:
        - Clips distances to [0, 1] range
        - Replaces NaN/inf values with 1.0 (maximum distance)
        - Ensures perfect symmetry

        This is useful when working with noisy or edge-case data.

        Parameters
        ----------
        compounds : List[Compound]
            List of compounds to compare
        metric : str, optional
            Distance metric to use. Default is 'wasserstein'.
        signal_variant : str, optional
            Signal variant to use. Default is 'raw' (per THEORY.md - no baseline correction).

        Returns
        -------
        NDArray[np.float64]
            Robust n×n distance matrix with guaranteed numerical properties

        Notes
        -----
        Additional robustness features:
        - All distances clipped to [0, 1] range
        - NaN values replaced with 1.0 (maximum distance)
        - Infinity values replaced with 1.0
        - Perfect symmetry enforced by averaging with transpose
        """
        # Compute base matrix
        distance_matrix = self.compute_similarity_matrix(
            compounds, metric=metric, signal_variant=signal_variant
        )

        # Ensure valid distance range [0, 1]
        distance_matrix = np.clip(distance_matrix, 0, 1)

        # Replace any NaN or inf values with maximum distance
        distance_matrix = np.nan_to_num(
            distance_matrix,
            nan=1.0,
            posinf=1.0,
            neginf=0.0
        )

        # Ensure perfect symmetry (average with transpose)
        distance_matrix = (distance_matrix + distance_matrix.T) / 2

        return distance_matrix
