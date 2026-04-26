"""
Analysis request DTO for application boundaries.

Defines the input structure for LC-Seq analysis workflows.
"""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, Dict, Any

# Import PreprocessingConfig from domain - single source of truth
from lcseq.domain.services.signal_preprocessor import PreprocessingConfig


@dataclass(frozen=True)
class CLPEParams:
    """
    Parameters for cLPE (chromatographic Linear Peptide Equation) validation.

    cLPE validates peak selection by checking if observed retention times
    are consistent with compound lipophilicity (AlogP) based on a linear
    regression model fitted per scaffold group.

    LogK is computed from OBSERVED RT: LogK = log10((RT - t0) / t0)
    NOT loaded from reference data.

    Attributes
    ----------
    enabled : bool
        Whether to run cLPE validation
    outlier_threshold : float
        Z-score threshold for outlier detection (from config)
    min_group_size : int
        Minimum compounds per scaffold for model fitting (from config)
    reference_csv_path : Optional[Path]
        Path to CSV with pre-computed AlogP and scaffold data (user-provided)
    t0 : Optional[float]
        Column dead time in minutes (None = derive from L0 peak)
    reselect_peaks : bool
        Whether to re-select peaks for outliers
    """
    enabled: bool
    outlier_threshold: float
    min_group_size: int
    reselect_peaks: bool
    reference_csv_path: Optional[Path] = None
    t0: Optional[float] = None

    @classmethod
    def from_dict(cls, params: dict) -> "CLPEParams":
        """Create CLPEParams from a dictionary. All keys required except paths."""
        if not params:
            raise ValueError("CLPEParams requires a non-empty dictionary")
        ref_path = params.get("reference_csv_path")
        return cls(
            enabled=params["enabled"],
            outlier_threshold=params["outlier_threshold"],
            min_group_size=params["min_group_size"],
            reselect_peaks=params["reselect_peaks"],
            reference_csv_path=Path(ref_path) if ref_path else None,
            t0=params.get("t0"),
        )

    def to_dict(self) -> dict:
        """Convert to dictionary for serialization."""
        return {
            "enabled": self.enabled,
            "reference_csv_path": str(self.reference_csv_path) if self.reference_csv_path else None,
            "t0": self.t0,
            "outlier_threshold": self.outlier_threshold,
            "min_group_size": self.min_group_size,
            "reselect_peaks": self.reselect_peaks,
        }


@dataclass(frozen=True)
class AnalysisRequest:
    """
    Request for LC-Seq library analysis.

    This DTO defines the input boundary for the analysis workflow,
    encapsulating all parameters needed to perform analysis.

    All required parameters must be provided from config - no hardcoded defaults.

    Attributes
    ----------
    data_path : Path
        Path to chromatogram data (HDF5, CSV, or Excel)
    output_path : Path
        Where to save results
    detection_params : Dict[str, Any]
        Peak detection parameters (from config, required)
    validation_params : Dict[str, Any]
        Validation parameters (from config, required)
    preprocessing_params : PreprocessingConfig
        Preprocessing configuration (from config, required)
    library_path : Optional[Path]
        Path to library design file (user-provided)
    variant_mode : str
        'individual' or 'pooled' for variant handling
    hierarchy_mode : str
        'building_block' or 'monomer' for hierarchy construction
    export_formats : list
        Which formats to export ('csv', 'excel', 'json')
    clpe_params : Optional[CLPEParams]
        cLPE validation parameters (if validation enabled)
    num_workers : Optional[int]
        Number of parallel workers (None/1=sequential, >1=parallel, -1=all cores)

    Examples
    --------
    >>> request = AnalysisRequest(
    ...     data_path=Path('data/chromatograms.hdf5'),
    ...     output_path=Path('results/'),
    ...     detection_params=config.peak_detection_params,
    ...     validation_params=config.validation_params,
    ...     preprocessing_params=PreprocessingConfig.from_dict(config.preprocessing_params),
    ... )
    """

    # Required parameters (from config)
    data_path: Path
    output_path: Path
    detection_params: Dict[str, Any]
    validation_params: Dict[str, Any]
    preprocessing_params: PreprocessingConfig
    # Optional parameters
    library_path: Optional[Path] = None
    variant_mode: str = 'individual'
    hierarchy_mode: str = 'building_block'
    export_formats: list = None
    clpe_params: Optional[CLPEParams] = None
    num_workers: Optional[int] = None  # None/1=sequential, >1=parallel, -1=all cores

    def __post_init__(self):
        """Set default export formats if not provided."""
        # Use object.__setattr__ because dataclass is frozen
        if self.export_formats is None:
            object.__setattr__(self, 'export_formats', ['csv', 'excel'])

    def to_dict(self) -> Dict[str, Any]:
        """
        Convert to dictionary representation.

        Returns
        -------
        Dict[str, Any]
            Dictionary with all request parameters
        """
        return {
            'data_path': str(self.data_path),
            'library_path': str(self.library_path) if self.library_path else None,
            'output_path': str(self.output_path),
            'variant_mode': self.variant_mode,
            'hierarchy_mode': self.hierarchy_mode,
            'detection_params': self.detection_params,
            'validation_params': self.validation_params,
            'export_formats': self.export_formats,
            'preprocessing_params': self.preprocessing_params.to_dict() if self.preprocessing_params else None,
            'clpe_params': self.clpe_params.to_dict() if self.clpe_params else None,
            'num_workers': self.num_workers,
        }

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> 'AnalysisRequest':
        """
        Create from dictionary. Required fields must be present.

        Parameters
        ----------
        data : Dict[str, Any]
            Dictionary with request parameters

        Returns
        -------
        AnalysisRequest
            Constructed request object

        Raises
        ------
        KeyError
            If required fields are missing
        """
        # Parse optional clpe_params (may be None)
        clpe_data = data.get('clpe_params')
        clpe_params = CLPEParams.from_dict(clpe_data) if clpe_data else None

        return cls(
            data_path=Path(data['data_path']),
            output_path=Path(data['output_path']),
            detection_params=data['detection_params'],
            validation_params=data['validation_params'],
            preprocessing_params=PreprocessingConfig.from_dict(data['preprocessing_params']),
            library_path=Path(data['library_path']) if data.get('library_path') else None,
            variant_mode=data.get('variant_mode', 'individual'),
            hierarchy_mode=data.get('hierarchy_mode', 'building_block'),
            export_formats=data.get('export_formats'),
            clpe_params=clpe_params,
            num_workers=data.get('num_workers'),
        )
