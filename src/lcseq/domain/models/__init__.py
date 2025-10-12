"""
Domain models for LC-Seq analysis.

Graph structures and analysis results.
"""

from .compound_hierarchy import CompoundHierarchy, HierarchyMode
from .equivalence_class import EquivalenceClass, ConsensusStatus
from .peak_classification import PeakClassification
from .analysis_configuration import (
    AnalysisConfiguration,
    AnalysisMode,
    HierarchyMode as ConfigHierarchyMode,
)
from .analysis_result import AnalysisResult

__all__ = [
    "CompoundHierarchy",
    "HierarchyMode",
    "EquivalenceClass",
    "ConsensusStatus",
    "PeakClassification",
    "AnalysisConfiguration",
    "AnalysisMode",
    "ConfigHierarchyMode",
    "AnalysisResult",
]
