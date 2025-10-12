"""Use cases for chromatogram processing workflows.

These use cases encapsulate the business logic of processing chromatograms:
- Peak detection
- Peak classification
- Peak integration

These workflows orchestrate multiple domain services to accomplish complete
processing tasks.
"""

from typing import List, Dict, Optional
from lcseq.domain.entities import Compound, Peak
from lcseq.domain.models import CompoundHierarchy
from lcseq.domain.services import (
    PeakDetector,
    PeakClassifier,
    PeakIntegrator,
)
from lcseq.config import (
    DEFAULT_Z_THRESHOLD,
    DEFAULT_PROMINENCE_PERCENTILE,
    DEFAULT_MIN_SNR,
    DEFAULT_MIN_BASELINE_SDS,
    DEFAULT_SIGNAL_VARIANT,
    DEFAULT_TRUNCATION_MARGIN,
)
class ProcessChromatogramsUseCase:
    """Use case for processing chromatograms with detection and classification.

    This use case orchestrates:
    1. Peak detection using Discrete Morse Theory + Poisson statistics
    2. Peak classification using hierarchy

    Used by: lineage mode
    """

    def __init__(
        self,
        detector: Optional[PeakDetector] = None,
        classifier: Optional[PeakClassifier] = None,
    ):
        """Initialize with domain services (dependency injection).

        Args:
            detector: Peak detection service (optional)
            classifier: Peak classification service (optional)
        """
        self.detector = detector or PeakDetector()
        self.classifier = classifier or PeakClassifier()

    def execute(
        self,
        compounds: List[Compound],
        hierarchy: CompoundHierarchy,
        z_threshold: float = DEFAULT_Z_THRESHOLD,
        prominence_percentile: float = DEFAULT_PROMINENCE_PERCENTILE,
        min_snr: float = DEFAULT_MIN_SNR,
        min_baseline_sds: float = DEFAULT_MIN_BASELINE_SDS,
        signal_variant: str = DEFAULT_SIGNAL_VARIANT,
        truncation_margin: float = DEFAULT_TRUNCATION_MARGIN,
    ) -> Dict[Compound, List[Peak]]:
        """Process all chromatograms: detect and classify peaks.

        Args:
            compounds: List of compounds to process
            hierarchy: Compound hierarchy for peak classification
            z_threshold: Poisson Z-score threshold
            prominence_percentile: Prominence percentile threshold
            min_snr: Adaptive SNR threshold multiplier
            min_baseline_sds: Global baseline threshold in SDs
            signal_variant: Signal variant to use for detection
            truncation_margin: Margin beyond truncation positions (in seconds)

        Returns:
            Dictionary mapping compounds to their detected peaks

        Notes
        -----
        Processes compounds in topological order (bottom-up, L₀ first) per
        THEORY.md Section 5.4 to ensure descendants are processed before ancestors.
        This allows peak classification to use descendant product positions as
        truncation constraints.

        References
        ----------
        THEORY.md Section 5.2: Statistical Significance Testing (Poisson + Prominence)
        THEORY.md Section 5.4: Global Classification via Constraint Propagation
        THEORY.md Section 7.1: Topological Sort (Kahn's Algorithm)
        """
        peaks_dict = {}

        # Process in topological order: L₀ first, then ancestors (THEORY.md 5.4)
        # This ensures descendants are processed before ancestors
        sorted_compounds = hierarchy.topological_sort()

        # Filter to only include compounds in the input list
        compounds_set = set(compounds)
        sorted_compounds = [c for c in sorted_compounds if c in compounds_set]

        for compound in sorted_compounds:
            # Step 1: Peak detection with Poisson statistics
            peaks = self.detector.detect_peaks(
                compound.chromatogram,
                z_threshold=z_threshold,
                prominence_percentile=prominence_percentile,
                min_snr=min_snr,
                min_baseline_sds=min_baseline_sds,
                signal_variant=signal_variant,
            )

            # Step 2: Peak classification
            compound.detected_peaks = peaks
            self.classifier.classify_all_peaks(
                compound,
                hierarchy,
                truncation_margin=truncation_margin,
            )

            # Step 3: Store results
            peaks_dict[compound] = compound.detected_peaks

        return peaks_dict


class ProcessChromatogramsWithIntegrationUseCase:
    """Use case for processing chromatograms with detection and integration.

    This use case orchestrates:
    1. Peak detection using Discrete Morse Theory + Poisson statistics
    2. Peak integration (area calculation)
    3. Results aggregation

    Used by: batch mode
    """

    def __init__(
        self,
        detector: Optional[PeakDetector] = None,
        integrator: Optional[PeakIntegrator] = None,
    ):
        """Initialize with domain services (dependency injection).

        Args:
            detector: Peak detection service (optional)
            integrator: Peak integration service (optional)
        """
        self.detector = detector or PeakDetector()
        self.integrator = integrator or PeakIntegrator()

    def execute(
        self,
        compounds: List[Compound],
        z_threshold: float = DEFAULT_Z_THRESHOLD,
        prominence_percentile: float = DEFAULT_PROMINENCE_PERCENTILE,
        min_snr: float = DEFAULT_MIN_SNR,
        min_baseline_sds: float = DEFAULT_MIN_BASELINE_SDS,
        signal_variant: str = DEFAULT_SIGNAL_VARIANT,
    ) -> List[Dict]:
        """Process all chromatograms: detect and integrate peaks.

        Args:
            compounds: List of compounds to process
            z_threshold: Poisson Z-score threshold
            prominence_percentile: Prominence percentile threshold
            min_snr: Adaptive SNR threshold multiplier
            min_baseline_sds: Global baseline threshold in SDs
            signal_variant: Signal variant to use for detection

        Returns:
            List of result dictionaries with keys:
            - sequence: Compound positional sequence
            - n_peaks: Number of peaks detected
            - total_area: Total integrated peak area

        References
        ----------
        THEORY.md Section 5.2: Statistical Significance Testing (Poisson + Prominence)
        """
        results = []

        for compound in compounds:
            # Step 1: Peak detection with Poisson statistics
            peaks = self.detector.detect_peaks(
                compound.chromatogram,
                z_threshold=z_threshold,
                prominence_percentile=prominence_percentile,
                min_snr=min_snr,
                min_baseline_sds=min_baseline_sds,
                signal_variant=signal_variant,
            )

            # Step 2: Peak integration
            total_area = 0.0
            for peak in peaks:
                left, right, area = self.integrator.integrate_peak(
                    compound.chromatogram,
                    peak.position,
                    signal_variant=signal_variant,
                )
                total_area += area

            # Step 3: Aggregate results
            results.append({
                "sequence": compound.positional_sequence,
                "n_peaks": len(peaks),
                "total_area": total_area,
            })

        return results
