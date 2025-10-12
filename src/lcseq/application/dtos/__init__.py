"""Data Transfer Objects for application boundaries."""

from .analysis_request import AnalysisRequest
from .analysis_response import AnalysisResponse, CompoundResult, ValidationSummary

__all__ = [
    'AnalysisRequest',
    'AnalysisResponse',
    'CompoundResult',
    'ValidationSummary',
]
