"""Use cases for chromatogram processing workflows.

These use cases encapsulate the business logic of processing chromatograms:
- Peak detection
- Peak classification
- Peak integration

These workflows orchestrate multiple domain services to accomplish complete
processing tasks.
"""

import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import List, Dict, Optional, Any, TYPE_CHECKING, Tuple
from tqdm import tqdm
from lcseq.domain.entities import Compound, Peak
from lcseq.domain.models import CompoundHierarchy
from lcseq.domain.services import (
    PeakDetector,
    PeakClassifier,
    PeakIntegrator,
    SignalPreprocessor,
    PreprocessingConfig,
)


def _detect_peaks_for_compound(
    compound_data: Tuple[int, Any, Any],  # (index, chromatogram, detection_params)
) -> Tuple[int, List[Peak]]:
    """
    Worker function for parallel peak detection.

    Must be top-level function for pickling in multiprocessing.

    Parameters
    ----------
    compound_data : tuple
        (compound_index, chromatogram, detection_params_dict)

    Returns
    -------
    tuple
        (compound_index, detected_peaks)
    """
    idx, chromatogram, params = compound_data
    detector = PeakDetector(sigma_clip_sigma=params["sigma_clip_sigma"])
    peaks = detector.detect_peaks(
        chromatogram,
        alpha=params["alpha"],
        prominence_percentile=params["prominence_percentile"],
        min_snr=params["min_snr"],
        min_baseline_sds=params["min_baseline_sds"],
        signal_variant=params["signal_variant"],
        min_dispersion_r=params["min_dispersion_r"],
        include_rejected=params["include_rejected"],
    )
    return idx, peaks

if TYPE_CHECKING:
    from lcseq.domain.services.clpe_validator import CLPEValidator
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
        sigma_clip_sigma: float = 2.0,
    ):
        """Initialize with domain services (dependency injection).

        Args:
            detector: Peak detection service (optional, created if not provided)
            classifier: Peak classification service (optional, created if not provided)
            sigma_clip_sigma: Sigma for baseline estimation (default 2.0 = 95% CI)
        """
        self._sigma_clip_sigma = sigma_clip_sigma
        self.detector = detector or PeakDetector(sigma_clip_sigma=sigma_clip_sigma)
        self.classifier = classifier or PeakClassifier()

    def execute(
        self,
        compounds: List[Compound],
        hierarchy: CompoundHierarchy,
        # Peak detection parameters
        alpha: float,
        prominence_percentile: float,
        min_snr: float,
        min_baseline_sds: float,
        signal_variant: str,
        min_dispersion_r: float,
        sigma_clip_sigma: float,
        # Peak classification parameters
        alpha_product: float,
        truncation_margin: float,
        peak_matching_tolerance: float,
        hungarian_min_threshold: float,
        # Validation parameters
        include_rejected: bool,
        # Preprocessing parameters
        preprocessing_params: Dict[str, Any],
        # Performance parameters
        num_workers: Optional[int] = None,
        # Optional cLPE validator (dependency injection)
        clpe_validator: Optional["CLPEValidator"] = None,
        alogp_map: Optional[Dict[str, float]] = None,
        scaffold_map: Optional[Dict[str, str]] = None,
    ) -> Tuple[Dict[Compound, List[Peak]], Optional[Dict[str, Any]]]:
        """Process all chromatograms: detect and classify peaks.

        Args:
            compounds: List of compounds to process
            hierarchy: Compound hierarchy for peak classification
            alpha: Significance level for general peak acceptance (permissive)
            alpha_product: Significance level for product selection (strict)
            prominence_percentile: Prominence percentile threshold
            min_snr: Adaptive SNR threshold multiplier
            min_baseline_sds: Global baseline threshold in SDs
            signal_variant: Signal variant to use for detection
            truncation_margin: Margin beyond truncation positions (in seconds)
            peak_matching_tolerance: Relative tolerance for peak matching
            clpe_validator: Optional cLPE validator for peak validation
            alogp_map: Optional mapping from compound ID to AlogP
            scaffold_map: Optional mapping from compound ID to scaffold
            num_workers: Number of parallel workers for peak detection.
                None or 1 = sequential processing (default).
                >1 = parallel processing with specified workers.
                -1 = use all available CPU cores.

        Returns:
            Tuple of:
            - Dictionary mapping compounds to their detected peaks
            - cLPE statistics (if validation enabled), None otherwise

        Notes
        -----
        Processing is done in two phases:
        1. Peak detection: Detect peaks for all compounds (parallelizable)
        2. Classification: Classify peaks hierarchically (bottom-up) with
           descendant matching to propagate peak origins up the hierarchy

        If cLPE validation is enabled, it runs at each level during
        classification to validate and potentially re-select peaks.

        This two-phase approach ensures all descendant peaks are detected
        before classification, enabling TRUNCATION_UNKNOWN matching where
        peaks are matched to any descendant peak (not just products).

        Parallel processing (num_workers > 1) provides 4-8x speedup on
        multi-core systems for large datasets (1000+ compounds).

        References
        ----------
        THEORY.md Section 5.2: Statistical Significance Testing (Poisson + Prominence)
        THEORY.md Section 5.4: Global Classification via Constraint Propagation
        THEORY.md Section 7.1: Topological Sort (Kahn's Algorithm)
        """
        peaks_dict = {}

        # Phase 0: Apply preprocessing if enabled
        preprocessing_config = PreprocessingConfig.from_dict(preprocessing_params)
        if preprocessing_config.enabled:
            preprocessor = SignalPreprocessor(preprocessing_config)
            for compound in compounds:
                if compound.chromatogram is not None:
                    if not compound.chromatogram.has_signal_variant("corrected"):
                        preprocessor.preprocess(compound.chromatogram)

        # Phase 1: Detect peaks for all compounds
        show_progress = len(compounds) > 100

        # Determine number of workers
        if num_workers == -1:
            num_workers = os.cpu_count() or 1
        use_parallel = num_workers is not None and num_workers > 1

        if use_parallel and len(compounds) > num_workers:
            # Parallel peak detection
            detection_params = {
                "alpha": alpha,
                "prominence_percentile": prominence_percentile,
                "min_snr": min_snr,
                "min_baseline_sds": min_baseline_sds,
                "signal_variant": signal_variant,
                "min_dispersion_r": min_dispersion_r,
                "sigma_clip_sigma": sigma_clip_sigma,
                "include_rejected": include_rejected,
            }

            # Prepare work items: (index, chromatogram, params)
            work_items = [
                (i, compound.chromatogram, detection_params)
                for i, compound in enumerate(compounds)
            ]

            # Process in parallel with progress bar
            with ProcessPoolExecutor(max_workers=num_workers) as executor:
                futures = {
                    executor.submit(_detect_peaks_for_compound, item): item[0]
                    for item in work_items
                }

                with tqdm(
                    total=len(compounds),
                    desc=f"Detecting peaks ({num_workers} workers)",
                    disable=not show_progress,
                    unit="cpd"
                ) as pbar:
                    for future in as_completed(futures):
                        idx, peaks = future.result()
                        compounds[idx].detected_peaks = peaks
                        pbar.update(1)
        else:
            # Sequential peak detection (default)
            compounds_iter = tqdm(
                compounds,
                desc="Detecting peaks",
                disable=not show_progress,
                unit="cpd"
            )
            for compound in compounds_iter:
                peaks = self.detector.detect_peaks(
                    compound.chromatogram,
                    alpha=alpha,
                    prominence_percentile=prominence_percentile,
                    min_snr=min_snr,
                    min_baseline_sds=min_baseline_sds,
                    signal_variant=signal_variant,
                    min_dispersion_r=min_dispersion_r,
                    include_rejected=include_rejected,
                )
                compound.detected_peaks = peaks

        # Phase 2: Classify peaks hierarchically with descendant matching
        # This processes bottom-up (L0 -> L1 -> ...) to enable match propagation
        # If cLPE is enabled, validation happens at each level
        clpe_stats = self.classifier.classify_hierarchy(
            hierarchy,
            tolerance=peak_matching_tolerance,
            truncation_margin=truncation_margin,
            alpha_product=alpha_product,
            hungarian_min_threshold=hungarian_min_threshold,
            clpe_validator=clpe_validator,
            alogp_map=alogp_map,
            scaffold_map=scaffold_map,
        )

        # Collect results
        for compound in compounds:
            peaks_dict[compound] = compound.detected_peaks

        return peaks_dict, clpe_stats


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
        sigma_clip_sigma: float = 2.0,
    ):
        """Initialize with domain services (dependency injection).

        Args:
            detector: Peak detection service (optional)
            integrator: Peak integration service (optional)
            sigma_clip_sigma: Sigma for baseline estimation (default 2.0 = 95% CI)
        """
        self.detector = detector or PeakDetector(sigma_clip_sigma=sigma_clip_sigma)
        self.integrator = integrator or PeakIntegrator()

    def execute(
        self,
        compounds: List[Compound],
        alpha: float,
        prominence_percentile: float,
        min_snr: float,
        min_baseline_sds: float,
        signal_variant: str,
        min_dispersion_r: float,
        include_rejected: bool,
    ) -> List[Dict]:
        """Process all chromatograms: detect and integrate peaks.

        Args:
            compounds: List of compounds to process
            alpha: Significance level (false positive rate)
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
                alpha=alpha,
                prominence_percentile=prominence_percentile,
                min_snr=min_snr,
                min_baseline_sds=min_baseline_sds,
                signal_variant=signal_variant,
                min_dispersion_r=min_dispersion_r,
                include_rejected=include_rejected,
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
                "sequence": compound.positional_block_sequence,
                "n_peaks": len(peaks),
                "total_area": total_area,
            })

        return results
