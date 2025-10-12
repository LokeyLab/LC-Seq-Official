"""
Analysis response DTO for application boundaries.

Defines the output structure from LC-Seq analysis workflows.
"""

from dataclasses import dataclass, field
from typing import List, Dict, Any, Optional
from datetime import datetime


@dataclass
class CompoundResult:
    """
    Analysis result for a single compound.

    Attributes
    ----------
    compound_id : str
        Compound identifier
    sequence : str
        Compound sequence (N→C order)
    level : int
        Compound level in hierarchy
    validation_status : str
        Validation status (VALIDATED, LIKELY_SUCCESS, etc.)
    purity : float
        Purity value [0, 1]
    purity_category : str
        Adaptive purity category
    snr : float
        Signal-to-noise ratio
    retention_time : Optional[float]
        Product peak retention time (if detected)
    peak_count : int
        Number of detected peaks
    truncation_count : int
        Number of truncation peaks
    unknown_count : int
        Number of unknown peaks
    """

    compound_id: str
    sequence: str
    level: int
    validation_status: str
    purity: float
    purity_category: str
    snr: float
    retention_time: Optional[float] = None
    peak_count: int = 0
    truncation_count: int = 0
    unknown_count: int = 0

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary."""
        return {
            'compound_id': self.compound_id,
            'sequence': self.sequence,
            'level': self.level,
            'validation_status': self.validation_status,
            'purity': self.purity,
            'purity_category': self.purity_category,
            'snr': self.snr,
            'retention_time': self.retention_time,
            'peak_count': self.peak_count,
            'truncation_count': self.truncation_count,
            'unknown_count': self.unknown_count
        }


@dataclass
class ValidationSummary:
    """
    Summary of validation results across dataset.

    Attributes
    ----------
    total_compounds : int
        Total number of compounds analyzed
    validated_count : int
        Number with VALIDATED status
    likely_success_count : int
        Number with LIKELY_SUCCESS status
    uncertain_count : int
        Number with UNCERTAIN status
    likely_failure_count : int
        Number with LIKELY_FAILURE status
    failed_count : int
        Number with FAILED status
    validation_rate : float
        Fraction validated or likely success
    median_purity : float
        Median purity across dataset
    dataset_quality : str
        Overall dataset quality assessment
    """

    total_compounds: int
    validated_count: int
    likely_success_count: int
    uncertain_count: int
    likely_failure_count: int
    failed_count: int
    validation_rate: float
    median_purity: float
    dataset_quality: str

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary."""
        return {
            'total_compounds': self.total_compounds,
            'validated_count': self.validated_count,
            'likely_success_count': self.likely_success_count,
            'uncertain_count': self.uncertain_count,
            'likely_failure_count': self.likely_failure_count,
            'failed_count': self.failed_count,
            'validation_rate': self.validation_rate,
            'median_purity': self.median_purity,
            'dataset_quality': self.dataset_quality
        }


@dataclass
class AnalysisResponse:
    """
    Complete response from LC-Seq library analysis.

    This DTO defines the output boundary for the analysis workflow.

    Attributes
    ----------
    request_id : str
        Unique identifier for this analysis
    timestamp : datetime
        When analysis completed
    compound_results : List[CompoundResult]
        Per-compound analysis results
    validation_summary : ValidationSummary
        Dataset-wide validation summary
    dataset_stats : Dict[str, float]
        Dataset statistics (percentiles, MAD, etc.)
    processing_metadata : Dict[str, Any]
        Processing parameters and metadata
    errors : List[str]
        Any errors encountered during analysis
    warnings : List[str]
        Any warnings generated during analysis

    Examples
    --------
    >>> response = AnalysisResponse(
    ...     request_id='analysis_20251008_123456',
    ...     timestamp=datetime.now(),
    ...     compound_results=[...],
    ...     validation_summary=ValidationSummary(...),
    ...     dataset_stats={'purity_p50': 0.7, ...},
    ...     processing_metadata={'runtime_seconds': 120.5},
    ...     errors=[],
    ...     warnings=[]
    ... )
    """

    request_id: str
    timestamp: datetime
    compound_results: List[CompoundResult]
    validation_summary: ValidationSummary
    dataset_stats: Dict[str, float]
    processing_metadata: Dict[str, Any] = field(default_factory=dict)
    errors: List[str] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)

    def to_dict(self) -> Dict[str, Any]:
        """
        Convert to dictionary representation.

        Returns
        -------
        Dict[str, Any]
            Dictionary with all response data
        """
        return {
            'request_id': self.request_id,
            'timestamp': self.timestamp.isoformat(),
            'compound_results': [r.to_dict() for r in self.compound_results],
            'validation_summary': self.validation_summary.to_dict(),
            'dataset_stats': self.dataset_stats,
            'processing_metadata': self.processing_metadata,
            'errors': self.errors,
            'warnings': self.warnings
        }

    def get_success_rate(self) -> float:
        """Get fraction of compounds successfully validated."""
        return self.validation_summary.validation_rate

    def get_failed_compounds(self) -> List[CompoundResult]:
        """Get all compounds with FAILED status."""
        return [
            r for r in self.compound_results
            if r.validation_status == 'FAILED'
        ]

    def get_high_purity_compounds(self, threshold: float = 0.8) -> List[CompoundResult]:
        """Get compounds with purity above threshold."""
        return [
            r for r in self.compound_results
            if r.purity >= threshold
        ]
