"""
CSV exporter for analysis results.

Exports analysis results to CSV format for easy viewing and processing.
"""

import pandas as pd
from pathlib import Path
from typing import Dict
from ...application.dtos.analysis_response import AnalysisResponse


class CSVExporter:
    """
    Export analysis results to CSV format.

    Creates CSV files with compound-level results and summary statistics.

    Examples
    --------
    >>> exporter = CSVExporter()
    >>> response = AnalysisResponse(...)
    >>> exporter.export(response, Path('results/'))
    """

    def export(
        self,
        response: AnalysisResponse,
        output_dir: Path
    ) -> Path:
        """
        Export analysis results to CSV.

        Parameters
        ----------
        response : AnalysisResponse
            Analysis results to export
        output_dir : Path
            Output directory

        Returns
        -------
        Path
            Path to created CSV file

        Examples
        --------
        >>> csv_path = exporter.export(response, Path('results/'))
        >>> csv_path
        Path('results/compound_results.csv')
        """
        # Ensure output directory exists
        output_dir.mkdir(parents=True, exist_ok=True)

        # Convert compound results to DataFrame
        data = [result.to_dict() for result in response.compound_results]
        df = pd.DataFrame(data)

        # Save to CSV
        output_path = output_dir / 'compound_results.csv'
        df.to_csv(output_path, index=False)

        return output_path

    def export_summary(
        self,
        response: AnalysisResponse,
        output_dir: Path
    ) -> Path:
        """
        Export validation summary to CSV.

        Parameters
        ----------
        response : AnalysisResponse
            Analysis results
        output_dir : Path
            Output directory

        Returns
        -------
        Path
            Path to summary CSV

        Examples
        --------
        >>> summary_path = exporter.export_summary(response, Path('results/'))
        """
        output_dir.mkdir(parents=True, exist_ok=True)

        # Create summary DataFrame
        summary_data = response.validation_summary.to_dict()
        df = pd.DataFrame([summary_data])

        # Save to CSV
        output_path = output_dir / 'validation_summary.csv'
        df.to_csv(output_path, index=False)

        return output_path

    def export_dataset_stats(
        self,
        response: AnalysisResponse,
        output_dir: Path
    ) -> Path:
        """
        Export dataset statistics to CSV.

        Parameters
        ----------
        response : AnalysisResponse
            Analysis results
        output_dir : Path
            Output directory

        Returns
        -------
        Path
            Path to stats CSV

        Examples
        --------
        >>> stats_path = exporter.export_dataset_stats(response, Path('results/'))
        """
        output_dir.mkdir(parents=True, exist_ok=True)

        # Create stats DataFrame
        stats_data = response.dataset_stats
        df = pd.DataFrame([stats_data])

        # Save to CSV
        output_path = output_dir / 'dataset_statistics.csv'
        df.to_csv(output_path, index=False)

        return output_path

    def export_all(
        self,
        response: AnalysisResponse,
        output_dir: Path
    ) -> Dict[str, Path]:
        """
        Export all results to CSV files.

        Parameters
        ----------
        response : AnalysisResponse
            Analysis results
        output_dir : Path
            Output directory

        Returns
        -------
        Dict[str, Path]
            Map of file type -> path

        Examples
        --------
        >>> paths = exporter.export_all(response, Path('results/'))
        >>> paths
        {
            'compound_results': Path('results/compound_results.csv'),
            'summary': Path('results/validation_summary.csv'),
            'stats': Path('results/dataset_statistics.csv')
        }
        """
        paths = {}
        paths['compound_results'] = self.export(response, output_dir)
        paths['summary'] = self.export_summary(response, output_dir)
        paths['stats'] = self.export_dataset_stats(response, output_dir)

        return paths
