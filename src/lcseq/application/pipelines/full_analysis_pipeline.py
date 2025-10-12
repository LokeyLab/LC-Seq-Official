"""
Full end-to-end analysis pipeline.

Orchestrates: peak detection → integration → classification → validation.
"""

import time
from datetime import datetime
from typing import List, Dict, Any
from uuid import uuid4

from ...domain.entities.compound import Compound
from ...domain.models.compound_hierarchy import CompoundHierarchy, HierarchyMode
from ...domain.services.hierarchy_builder import HierarchyBuilder
from ...domain.services.peak_detector import PeakDetector
from ...domain.services.peak_integrator import PeakIntegrator
from ...domain.services.peak_classifier import PeakClassifier
from ...domain.services.validation.validation_workflow import ValidationWorkflow
from ..dtos.analysis_request import AnalysisRequest
from ..dtos.analysis_response import (
    AnalysisResponse,
    CompoundResult,
    ValidationSummary
)


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

            # Stage 2: Process chromatograms
            processed_compounds = self._process_chromatograms(
                compounds, request.detection_params
            )

            # Stage 3: Classify peaks
            classified_compounds = self._classify_peaks(
                processed_compounds, hierarchy
            )

            # Stage 4: Validate synthesis
            validation_results = self.validation_workflow.validate_library(
                classified_compounds,
                hierarchy,
                retention_precision=request.validation_params.get('retention_precision', 0.5)
            )

            # Stage 5: Build response
            compound_results = self._build_compound_results(
                classified_compounds, validation_results
            )
            validation_summary = self._build_validation_summary(
                validation_results
            )

            end_time = time.time()
            processing_metadata = {
                'runtime_seconds': end_time - start_time,
                'compound_count': len(compounds),
                'hierarchy_mode': request.hierarchy_mode,
                'variant_mode': request.variant_mode,
                'detection_params': request.detection_params
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

    def _process_chromatograms(
        self,
        compounds: List[Compound],
        detection_params: Dict[str, Any]
    ) -> List[Compound]:
        """
        Process chromatograms: peak detection and integration.

        Parameters
        ----------
        compounds : List[Compound]
            Compounds with chromatograms
        detection_params : Dict[str, Any]
            Peak detection parameters

        Returns
        -------
        List[Compound]
            Compounds with detected and integrated peaks
        """
        processed = []

        for compound in compounds:
            if compound.chromatogram is None:
                continue

            # Peak detection
            peaks = self.peak_detector.detect_peaks(
                compound.chromatogram,
                z_threshold=detection_params.get('z_threshold', 3.0),
                prominence_percentile=detection_params.get('prominence_percentile', 0.2)
            )

            # Peak integration
            integrated_peaks = [
                self.peak_integrator.integrate_peak(compound.chromatogram, peak)
                for peak in peaks
            ]

            # Update compound
            compound.detected_peaks = integrated_peaks
            processed.append(compound)

        return processed

    def _classify_peaks(
        self,
        compounds: List[Compound],
        hierarchy: CompoundHierarchy
    ) -> List[Compound]:
        """
        Classify detected peaks using DAG constraints.

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
        return self.peak_classifier.classify_library(compounds, hierarchy)

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
        from ...entities.peak import PeakType

        results = []
        validation_map = {
            vr['compound_id']: vr
            for vr in validation_results['validation_results']
        }

        for compound in compounds:
            compound_id = str(compound)
            vr = validation_map.get(compound_id, {})

            # Count peaks by type
            truncation_count = sum(
                1 for p in compound.detected_peaks
                if p.peak_type == PeakType.TRUNCATION
            )
            unknown_count = sum(
                1 for p in compound.detected_peaks
                if p.peak_type == PeakType.UNKNOWN
            )

            result = CompoundResult(
                compound_id=compound_id,
                sequence=str(compound.positional_sequence),
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
