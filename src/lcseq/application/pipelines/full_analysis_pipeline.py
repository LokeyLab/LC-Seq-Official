"""
Full end-to-end analysis pipeline.

Orchestrates: peak detection → integration → classification → validation.
"""

import os
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime
from typing import List, Dict, Any, Tuple, Optional
from uuid import uuid4

from tqdm import tqdm

from ...domain.entities.compound import Compound
from ...domain.entities.peak import PeakType
from ...domain.models.compound_hierarchy import CompoundHierarchy, HierarchyMode
from ...domain.services.hierarchy_builder import HierarchyBuilder
from ...domain.services.peak_detector import PeakDetector
from ...domain.services.peak_integrator import PeakIntegrator
from ...domain.services.peak_classifier import PeakClassifier
from ...domain.services.clpe_validator import CLPEValidator
from ...domain.services.signal_preprocessor import SignalPreprocessor, PreprocessingConfig
from ...domain.services.validation.validation_workflow import ValidationWorkflow
from ...infrastructure.loaders.clpe_reference_loader import CLPEReferenceLoader
from ..dtos.analysis_request import AnalysisRequest, CLPEParams
from ..dtos.analysis_response import (
    AnalysisResponse,
    CompoundResult,
    ValidationSummary
)


def _detect_and_integrate_peaks(
    work_item: Tuple[int, Any, Dict[str, Any], str],
) -> Tuple[int, List[Any]]:
    """
    Worker function for parallel peak detection and integration.

    Must be top-level function for pickling in multiprocessing.

    Parameters
    ----------
    work_item : tuple
        (compound_index, chromatogram, detection_params, signal_variant)

    Returns
    -------
    tuple
        (compound_index, detected_peaks_with_areas)
    """
    idx, chromatogram, detection_params, signal_variant = work_item

    # Create fresh instances (not shared across processes)
    detector = PeakDetector()
    integrator = PeakIntegrator()

    # Peak detection - all params required from config (no fallback defaults)
    peaks = detector.detect_peaks(
        chromatogram,
        alpha=detection_params['alpha'],
        prominence_percentile=detection_params['prominence_percentile'],
        min_snr=detection_params['min_snr'],
        min_baseline_sds=detection_params['min_baseline_sds'],
        signal_variant=signal_variant,
        min_dispersion_r=detection_params['min_dispersion_r'],
        include_rejected=detection_params['include_rejected'],
    )

    # Peak integration
    for peak in peaks:
        left, right, area = integrator.integrate_peak(
            chromatogram, peak.position,
            signal_variant=signal_variant
        )
        peak.area = area
        peak.left_boundary = left
        peak.right_boundary = right

    return idx, peaks


class FullAnalysisPipeline:
    """
    Complete end-to-end LC-Seq analysis pipeline.

    Orchestrates all processing stages from raw chromatograms to validated
    synthesis results.

    Processing stages:
    1. Hierarchy construction from library design
    2. Peak detection (Discrete Morse theory + Poisson statistics)
    3. Peak integration (area calculation)
    4. Peak classification (DAG constraint propagation)
    5. Synthesis validation (Bayesian framework with adaptive thresholds)

    Notes
    -----
    This pipeline implements the complete workflow from THEORY.md Parts 1-7.

    References
    ----------
    THEORY.md Part 5: Peak Detection Mathematical Foundations
    THEORY.md Part 6: Synthesis Validation Theory
    THEORY.md Section 5.7: Three-Stage Pipeline

    Examples
    --------
    >>> pipeline = FullAnalysisPipeline()
    >>> request = AnalysisRequest(...)
    >>> compounds = load_compounds(request.data_path)
    >>> response = pipeline.execute(compounds, request)
    >>> response.validation_summary.validation_rate
    0.75
    """

    def __init__(self):
        """Initialize pipeline with required services."""
        self.hierarchy_builder = HierarchyBuilder()
        self.signal_preprocessor = None  # Created per-request with params
        self.peak_detector = PeakDetector()
        self.peak_integrator = PeakIntegrator()
        self.peak_classifier = PeakClassifier()
        self.validation_workflow = ValidationWorkflow()

    def execute(
        self,
        compounds: List[Compound],
        request: AnalysisRequest
    ) -> AnalysisResponse:
        """
        Execute full analysis pipeline.

        Parameters
        ----------
        compounds : List[Compound]
            Compounds to analyze (with chromatograms)
        request : AnalysisRequest
            Analysis parameters

        Returns
        -------
        AnalysisResponse
            Complete analysis results

        Notes
        -----
        Pipeline stages:
        1. Build hierarchy
        2. Process chromatograms (detect peaks, integrate areas)
        3. Classify peaks
        4. Validate synthesis
        5. Aggregate results

        Examples
        --------
        >>> response = pipeline.execute(compounds, request)
        """
        start_time = time.time()
        request_id = f"analysis_{datetime.now().strftime('%Y%m%d_%H%M%S')}_{uuid4().hex[:8]}"

        errors = []
        warnings = []

        try:
            # Stage 1: Build hierarchy
            hierarchy_mode = (
                HierarchyMode.MONOMER if request.hierarchy_mode == 'monomer'
                else HierarchyMode.BUILDING_BLOCK
            )
            hierarchy = self.hierarchy_builder.build(compounds, hierarchy_mode)

            # Stage 1b: Preprocess signals (filtering + baseline correction)
            if request.preprocessing_params and request.preprocessing_params.enabled:
                self._preprocess_signals(compounds, request.preprocessing_params)

            # Stage 2: Process chromatograms
            processed_compounds = self._process_chromatograms(
                compounds, request.detection_params, request.preprocessing_params,
                num_workers=request.num_workers
            )

            # Stage 3: Classify peaks
            classified_compounds = self._classify_peaks(
                processed_compounds, hierarchy
            )

            # Stage 3b: cLPE validation (optional)
            clpe_stats = {}
            if request.clpe_params and request.clpe_params.enabled:
                classified_compounds, clpe_stats = self._run_clpe_validation(
                    classified_compounds, request.clpe_params
                )

            # Stage 4: Validate synthesis
            validation_results = self.validation_workflow.validate_library(
                classified_compounds,
                hierarchy,
                retention_precision=request.validation_params['retention_precision']
            )

            # Stage 5: Build response
            compound_results = self._build_compound_results(
                classified_compounds, validation_results
            )
            validation_summary = self._build_validation_summary(
                validation_results
            )

            end_time = time.time()
            preprocessing_info = None
            if request.preprocessing_params and request.preprocessing_params.enabled:
                preprocessing_info = {
                    'baseline_order': request.preprocessing_params.baseline_order,
                }
            processing_metadata = {
                'runtime_seconds': end_time - start_time,
                'compound_count': len(compounds),
                'hierarchy_mode': request.hierarchy_mode,
                'variant_mode': request.variant_mode,
                'detection_params': request.detection_params,
                'preprocessing': preprocessing_info,
                'clpe_validation': clpe_stats if clpe_stats else None,
                'num_workers': request.num_workers,
            }

            return AnalysisResponse(
                request_id=request_id,
                timestamp=datetime.now(),
                compound_results=compound_results,
                validation_summary=validation_summary,
                dataset_stats=validation_results['dataset_stats'],
                processing_metadata=processing_metadata,
                errors=errors,
                warnings=warnings
            )

        except Exception as e:
            errors.append(f"Pipeline execution failed: {str(e)}")
            # Return error response
            return AnalysisResponse(
                request_id=request_id,
                timestamp=datetime.now(),
                compound_results=[],
                validation_summary=ValidationSummary(
                    total_compounds=0,
                    validated_count=0,
                    likely_success_count=0,
                    uncertain_count=0,
                    likely_failure_count=0,
                    failed_count=0,
                    validation_rate=0.0,
                    median_purity=0.0,
                    dataset_quality='ERROR'
                ),
                dataset_stats={},
                processing_metadata={'runtime_seconds': time.time() - start_time},
                errors=errors,
                warnings=warnings
            )

    def _preprocess_signals(
        self,
        compounds: List[Compound],
        preprocessing_params: PreprocessingConfig
    ) -> None:
        """
        Apply signal preprocessing to all chromatograms.

        Applies baseline correction, storing corrected signal as a variant
        on each chromatogram.

        Parameters
        ----------
        compounds : List[Compound]
            Compounds with chromatograms (modified in-place)
        preprocessing_params : PreprocessingConfig
            Preprocessing configuration
        """
        preprocessor = SignalPreprocessor(preprocessing_params)

        for compound in compounds:
            if compound.chromatogram is not None:
                preprocessor.preprocess(compound.chromatogram)

    def _process_chromatograms(
        self,
        compounds: List[Compound],
        detection_params: Dict[str, Any],
        preprocessing_params: PreprocessingConfig = None,
        num_workers: Optional[int] = None
    ) -> List[Compound]:
        """
        Process chromatograms: peak detection and integration.

        Parameters
        ----------
        compounds : List[Compound]
            Compounds with chromatograms
        detection_params : Dict[str, Any]
            Peak detection parameters
        preprocessing_params : PreprocessingConfig, optional
            If provided and enabled, use corrected signal for peak detection
        num_workers : Optional[int]
            Number of parallel workers. None/1=sequential, >1=parallel, -1=all cores.

        Returns
        -------
        List[Compound]
            Compounds with detected and integrated peaks
        """
        # Determine which signal variant to use
        use_corrected = (
            preprocessing_params is not None
            and preprocessing_params.enabled
        )
        signal_variant = "corrected" if use_corrected else "raw"

        # Filter compounds with chromatograms
        compounds_with_chrom = [c for c in compounds if c.chromatogram is not None]

        # Determine number of workers
        if num_workers == -1:
            num_workers = os.cpu_count() or 1
        use_parallel = num_workers is not None and num_workers > 1

        show_progress = len(compounds_with_chrom) > 100

        if use_parallel and len(compounds_with_chrom) > num_workers:
            # Parallel processing
            work_items = [
                (i, c.chromatogram, detection_params, signal_variant)
                for i, c in enumerate(compounds_with_chrom)
            ]

            with ProcessPoolExecutor(max_workers=num_workers) as executor:
                futures = {
                    executor.submit(_detect_and_integrate_peaks, item): item[0]
                    for item in work_items
                }

                with tqdm(
                    total=len(compounds_with_chrom),
                    desc=f"Processing chromatograms ({num_workers} workers)",
                    disable=not show_progress,
                    unit="cpd"
                ) as pbar:
                    for future in as_completed(futures):
                        idx, peaks = future.result()
                        compounds_with_chrom[idx].detected_peaks = peaks
                        pbar.update(1)
        else:
            # Sequential processing
            compounds_iter = tqdm(
                compounds_with_chrom,
                desc="Processing chromatograms",
                disable=not show_progress,
                unit="cpd"
            )
            for compound in compounds_iter:
                # Peak detection - all params required from config (no fallback defaults)
                peaks = self.peak_detector.detect_peaks(
                    compound.chromatogram,
                    alpha=detection_params['alpha'],
                    prominence_percentile=detection_params['prominence_percentile'],
                    min_snr=detection_params['min_snr'],
                    min_baseline_sds=detection_params['min_baseline_sds'],
                    signal_variant=signal_variant,
                    min_dispersion_r=detection_params['min_dispersion_r'],
                    include_rejected=detection_params['include_rejected'],
                )

                # Peak integration
                for peak in peaks:
                    left, right, area = self.peak_integrator.integrate_peak(
                        compound.chromatogram, peak.position,
                        signal_variant=signal_variant
                    )
                    peak.area = area
                    peak.left_boundary = left
                    peak.right_boundary = right

                compound.detected_peaks = peaks

        return compounds_with_chrom

    def _classify_peaks(
        self,
        compounds: List[Compound],
        hierarchy: CompoundHierarchy
    ) -> List[Compound]:
        """
        Classify detected peaks using DAG constraints.

        Processes compounds bottom-up (L0 -> L1 -> ...) to enable
        descendant matching and peak origin tracking.

        Parameters
        ----------
        compounds : List[Compound]
            Compounds with detected peaks
        hierarchy : CompoundHierarchy
            DAG structure for constraint propagation

        Returns
        -------
        List[Compound]
            Compounds with classified peaks
        """
        # Use hierarchical classification with descendant matching
        self.peak_classifier.classify_hierarchy(hierarchy)
        return compounds

    def _run_clpe_validation(
        self,
        compounds: List[Compound],
        clpe_params: CLPEParams
    ) -> tuple:
        """
        Run cLPE validation to verify peak selection consistency.

        Uses the chromatographic Linear Peptide Equation to validate that
        selected peaks have retention times consistent with compound
        lipophilicity (AlogP) based on regression models per scaffold.

        Parameters
        ----------
        compounds : List[Compound]
            Compounds with classified peaks
        clpe_params : CLPEParams
            cLPE validation parameters

        Returns
        -------
        tuple
            (compounds, stats_dict) where compounds may have re-selected peaks
        """
        if not clpe_params.reference_csv_path:
            return compounds, {'error': 'No reference CSV path provided'}

        # Load reference data (AlogP and scaffold only - LogK computed from observed RT)
        loader = CLPEReferenceLoader()
        ref_data = loader.load(clpe_params.reference_csv_path)

        # Initialize validator
        validator = CLPEValidator(
            t0=clpe_params.t0,
            outlier_threshold=clpe_params.outlier_threshold,
            min_group_size=clpe_params.min_group_size
        )

        # Match compounds to reference data
        matched_compounds = []
        for compound in compounds:
            # Try to match by sequence or block_support_sequence
            keys_to_try = [
                compound.positional_block_sequence,
                compound.block_support_sequence
            ]
            if compound.compound_id:
                keys_to_try.insert(0, compound.compound_id)

            for key in keys_to_try:
                alogp = ref_data.alogp_map.get(key)
                scaffold = ref_data.scaffold_map.get(key)
                if alogp is not None and scaffold:
                    compound.alogp = alogp
                    compound.scaffold_group = scaffold
                    matched_compounds.append(compound)
                    break

        # Fit cLPE models
        validator.fit_models(
            matched_compounds,
            ref_data.alogp_map,
            ref_data.scaffold_map
        )

        # Validate and optionally re-select peaks
        n_validated = 0
        n_outliers = 0
        n_reselected = 0

        for compound in matched_compounds:
            if not compound.selected_peak or compound.alogp is None or not compound.scaffold_group:
                continue

            result, new_peak = validator.validate_and_reselect(
                compound, compound.alogp, compound.scaffold_group
            )

            n_validated += 1

            # Store validation results on peak
            if compound.selected_peak:
                compound.selected_peak.clpe_residual = result.residual
                compound.selected_peak.clpe_z_score = result.z_score
                compound.selected_peak.clpe_is_outlier = result.is_outlier

            if result.is_outlier:
                n_outliers += 1

            # Re-select peak if requested and better alternative found
            if clpe_params.reselect_peaks and new_peak:
                old_peak = compound.selected_peak
                old_peak.peak_type = PeakType.UNKNOWN  # Demote old selection
                compound.selected_peak = new_peak
                new_peak.peak_type = PeakType.PUTATIVE_PRODUCT
                new_peak.clpe_reselected = True
                n_reselected += 1

        # Build stats
        stats = {
            'matched_compounds': len(matched_compounds),
            'validated_compounds': n_validated,
            'outliers': n_outliers,
            'outlier_rate': n_outliers / n_validated if n_validated > 0 else 0.0,
            'reselected_peaks': n_reselected,
            'models_fitted': len(validator.models)
        }

        return compounds, stats

    def _build_compound_results(
        self,
        compounds: List[Compound],
        validation_results: Dict[str, Any]
    ) -> List[CompoundResult]:
        """
        Build compound result DTOs.

        Parameters
        ----------
        compounds : List[Compound]
            Analyzed compounds
        validation_results : Dict[str, Any]
            Validation results from workflow

        Returns
        -------
        List[CompoundResult]
            DTO representations
        """
        results = []
        validation_map = {
            vr['compound_id']: vr
            for vr in validation_results['validation_results']
        }

        for compound in compounds:
            compound_id = str(compound)
            vr = validation_map.get(compound_id, {})

            # Count peaks by type
            # TRUNCATION: matched to descendant product peak
            # TRUNCATION_UNKNOWN: matched to descendant non-product peak
            truncation_count = sum(
                1 for p in compound.detected_peaks
                if p.peak_type == PeakType.TRUNCATION
            )
            truncation_unknown_count = sum(
                1 for p in compound.detected_peaks
                if p.peak_type == PeakType.TRUNCATION_UNKNOWN
            )
            unknown_count = sum(
                1 for p in compound.detected_peaks
                if p.peak_type == PeakType.UNKNOWN
            )

            result = CompoundResult(
                compound_id=compound_id,
                sequence=str(compound.positional_block_sequence),
                level=compound.level,
                validation_status=str(vr.get('validation_status', 'NOT_VALIDATED')),
                purity=vr.get('purity', 0.0),
                purity_category=vr.get('purity_category', 'unknown'),
                snr=vr.get('snr', 0.0),
                retention_time=compound.selected_peak.position if compound.selected_peak else None,
                peak_count=len(compound.detected_peaks),
                truncation_count=truncation_count,
                unknown_count=unknown_count
            )
            results.append(result)

        return results

    def _build_validation_summary(
        self,
        validation_results: Dict[str, Any]
    ) -> ValidationSummary:
        """
        Build validation summary DTO.

        Parameters
        ----------
        validation_results : Dict[str, Any]
            Validation results from workflow

        Returns
        -------
        ValidationSummary
            Summary DTO
        """
        summary = validation_results['summary']
        stats = validation_results['dataset_stats']

        # Determine dataset quality
        median_purity = stats.get('purity_p50', 0.0)
        if median_purity > 0.8:
            quality = 'HIGH'
        elif median_purity > 0.6:
            quality = 'MODERATE'
        elif median_purity > 0.4:
            quality = 'LOW'
        else:
            quality = 'POOR'

        return ValidationSummary(
            total_compounds=summary['total_compounds'],
            validated_count=summary['validated_count'],
            likely_success_count=summary['likely_success_count'],
            uncertain_count=summary['uncertain_count'],
            likely_failure_count=summary['likely_failure_count'],
            failed_count=summary['failed_count'],
            validation_rate=summary['validation_rate'],
            median_purity=median_purity,
            dataset_quality=quality
        )
