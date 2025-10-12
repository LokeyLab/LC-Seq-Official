"""
LC-Seq: DNA-Encoded Library Chromatographic Data Analysis

A Clean Architecture implementation for analyzing chromatography data
from DNA-encoded library screens.

Based on theoretical foundations in docs/THEORY.md.
"""

__version__ = "0.1.0"
__author__ = "Adam"

# Expose key domain entities
from lcseq.domain.entities.building_block import BuildingBlock
from lcseq.domain.entities.chromatogram import Chromatogram
from lcseq.domain.entities.compound import Compound
from lcseq.domain.entities.peak import Peak, PeakType, ValidationStatus

# Expose main configuration
from lcseq.domain.models.analysis_configuration import (
    AnalysisConfiguration,
    AnalysisMode,
    HierarchyMode,
)

__all__ = [
    "__version__",
    "__author__",
    "BuildingBlock",
    "Chromatogram",
    "Compound",
    "Peak",
    "PeakType",
    "ValidationStatus",
    "AnalysisConfiguration",
    "AnalysisMode",
    "HierarchyMode",
]
