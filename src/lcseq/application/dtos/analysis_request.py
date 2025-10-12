"""
Analysis request DTO for application boundaries.

Defines the input structure for LC-Seq analysis workflows.
"""

from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Dict, Any


@dataclass(frozen=True)
class AnalysisRequest:
    """
    Request for LC-Seq library analysis.

    This DTO defines the input boundary for the analysis workflow,
    encapsulating all parameters needed to perform analysis.

    Attributes
    ----------
    data_path : Path
        Path to chromatogram data (HDF5, CSV, or Excel)
    library_path : Optional[Path]
        Path to library design file
    output_path : Path
        Where to save results
    variant_mode : str
        'individual' or 'pooled' for variant handling
    hierarchy_mode : str
        'building_block' or 'monomer' for hierarchy construction
    detection_params : Dict[str, Any]
        Peak detection parameters (min_persistence, etc.)
    validation_params : Dict[str, Any]
        Validation parameters (thresholds, etc.)
    export_formats : list
        Which formats to export ('csv', 'excel', 'json')

    Examples
    --------
    >>> request = AnalysisRequest(
    ...     data_path=Path('data/chromatograms.hdf5'),
    ...     library_path=Path('data/library_design.csv'),
    ...     output_path=Path('results/'),
    ...     variant_mode='individual',
    ...     hierarchy_mode='building_block',
    ...     detection_params={'min_persistence': 0.05},
    ...     validation_params={'correlation_threshold': 0.8},
    ...     export_formats=['csv', 'excel']
    ... )
    """

    data_path: Path
    output_path: Path
    library_path: Optional[Path] = None
    variant_mode: str = 'individual'
    hierarchy_mode: str = 'building_block'
    detection_params: Dict[str, Any] = None
    validation_params: Dict[str, Any] = None
    export_formats: list = None

    def __post_init__(self):
        """Set default parameters if not provided."""
        # Use object.__setattr__ because dataclass is frozen
        if self.detection_params is None:
            object.__setattr__(self, 'detection_params', {
                'min_persistence': 0.05,
                'boundary_method': 'valley_or_5pct'
            })

        if self.validation_params is None:
            object.__setattr__(self, 'validation_params', {
                'correlation_threshold': 0.8,
                'retention_precision': 0.5
            })

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
            'export_formats': self.export_formats
        }

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> 'AnalysisRequest':
        """
        Create from dictionary.

        Parameters
        ----------
        data : Dict[str, Any]
            Dictionary with request parameters

        Returns
        -------
        AnalysisRequest
            Constructed request object
        """
        return cls(
            data_path=Path(data['data_path']),
            library_path=Path(data['library_path']) if data.get('library_path') else None,
            output_path=Path(data['output_path']),
            variant_mode=data.get('variant_mode', 'individual'),
            hierarchy_mode=data.get('hierarchy_mode', 'building_block'),
            detection_params=data.get('detection_params'),
            validation_params=data.get('validation_params'),
            export_formats=data.get('export_formats')
        )
