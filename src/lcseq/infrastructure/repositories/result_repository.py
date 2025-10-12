"""
Repository for saving analysis results.

Implements persistence operations for analysis results.
"""

import json
from pathlib import Path
from typing import Dict, Any
from ...application.dtos.analysis_response import AnalysisResponse


class ResultRepository:
    """
    Repository for analysis result persistence.

    Handles saving analysis results in various formats.

    Examples
    --------
    >>> repo = ResultRepository()
    >>> response = AnalysisResponse(...)
    >>> repo.save_json(response, Path('results.json'))
    """

    def save_json(
        self,
        response: AnalysisResponse,
        file_path: Path
    ) -> None:
        """
        Save analysis response as JSON.

        Parameters
        ----------
        response : AnalysisResponse
            Analysis results to save
        file_path : Path
            Output file path

        Examples
        --------
        >>> repo.save_json(response, Path('results/analysis.json'))
        """
        # Convert to dict
        data = response.to_dict()

        # Ensure parent directory exists
        file_path.parent.mkdir(parents=True, exist_ok=True)

        # Write JSON
        with open(file_path, 'w') as f:
            json.dump(data, f, indent=2, default=str)

    def load_json(self, file_path: Path) -> Dict[str, Any]:
        """
        Load analysis results from JSON.

        Parameters
        ----------
        file_path : Path
            JSON file path

        Returns
        -------
        Dict[str, Any]
            Loaded analysis data

        Examples
        --------
        >>> data = repo.load_json(Path('results.json'))
        """
        if not file_path.exists():
            raise FileNotFoundError(f"Results file not found: {file_path}")

        with open(file_path, 'r') as f:
            data = json.load(f)

        return data
