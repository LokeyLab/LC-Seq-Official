"""
Domain services for LC-Seq analysis.

This package contains stateless domain services that orchestrate domain entities
and models. Services provide graph operations, hierarchy analysis, algorithms,
and validation.

Graph Operations
----------------
- HierarchyBuilder: Constructs CompoundHierarchy from list of compounds
- EquivalenceClassBuilder: Groups compounds by residue sequence
- PathFinder: Finds paths in hierarchy DAG
- LevelAnalyzer: Analyzes hierarchy by truncation levels
- ValidationChecker: Validates hierarchy structural properties

Algorithm Services
------------------
- BaselineEstimatorService: Swappable baseline estimation methods (THEORY.md 5.1)
- PeakDetector: Morse theory + persistent homology on raw signals (THEORY.md 5.1-5.2)
- PeakIntegrator: Area integration and boundary detection on raw signals (THEORY.md 5.0.7)
- PeakClassifier: Peak type classification via DAG constraints (THEORY.md 5.3-5.6)
- BayesianValidator: Synthesis validation framework (THEORY.md 6.7-6.10)
- PurityCalculator: Single-source-of-truth purity calculation (THEORY.md 6.3)
- SNRCalculator: Single-source-of-truth SNR calculation (THEORY.md 6.5)
- SequenceAligner: Single-source-of-truth MSA-style sequence alignment (THEORY.md 2.3.4)
- SequenceSimilarityAnalyzer: Sequence distance metrics (Levenshtein, etc.)
- CompoundSimilarityAnalyzer: Chromatographic similarity (Wasserstein, etc.)
- SignalAggregator: Pooled signal aggregation for positional variants (THEORY.md 4.2.4)
- QualityAssessor: Signal quality metrics for library pre-filtering (THEORY.md 4.2.8)

References
----------
THEORY.md Section 4.2: Hierarchy Construction
THEORY.md Section 3.3: Hierarchy Properties
THEORY.md Section 5: Peak Detection Mathematical Foundations
THEORY.md Section 6: Synthesis Validation Theory
ARCHITECTURE.md: Domain Services
"""

from .hierarchy_builder import HierarchyBuilder
from .equivalence_class_builder import EquivalenceClassBuilder
from .path_finder import PathFinder
from .level_analyzer import LevelAnalyzer
from .validation_checker import ValidationChecker
from .baseline_estimator import BaselineEstimatorService
from .significance_tester import SignificanceTesterService
from .peak_detector import PeakDetector
from .peak_integrator import PeakIntegrator
from .peak_classifier import PeakClassifier
from .bayesian_validator import BayesianValidator
from .purity_calculator import PurityCalculator
from .snr_calculator import SNRCalculator
from .sequence_aligner import SequenceAligner
from .sequence_similarity import SequenceSimilarityAnalyzer
from .compound_similarity import CompoundSimilarityAnalyzer
from .compound_search import CompoundSearchService
from .lineage_finder import LineageFinderService
from .compound_ordering import CompoundOrderingService
from .pooling.signal_aggregator import SignalAggregator
from .coupling_efficiency_calculator import CouplingEfficiencyCalculator
from .signal_preprocessor import SignalPreprocessor, PreprocessingConfig
from .quality_assessor import QualityAssessor, SignalQualityMetrics, EquivalenceClassQuality

__all__ = [
    "HierarchyBuilder",
    "EquivalenceClassBuilder",
    "PathFinder",
    "LevelAnalyzer",
    "ValidationChecker",
    "BaselineEstimatorService",
    "SignificanceTesterService",
    "PeakDetector",
    "PeakIntegrator",
    "PeakClassifier",
    "BayesianValidator",
    "PurityCalculator",
    "SNRCalculator",
    "SequenceAligner",
    "SequenceSimilarityAnalyzer",
    "CompoundSimilarityAnalyzer",
    "CompoundSearchService",
    "LineageFinderService",
    "CompoundOrderingService",
    "SignalAggregator",
    "CouplingEfficiencyCalculator",
    "SignalPreprocessor",
    "PreprocessingConfig",
    "QualityAssessor",
    "SignalQualityMetrics",
    "EquivalenceClassQuality",
]
