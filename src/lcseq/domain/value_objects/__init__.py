"""
Value objects for LC-Seq domain.

Value objects are immutable objects that represent descriptive aspects of the
domain with no conceptual identity. They are compared based on their values,
not their identity.
"""

from lcseq.domain.value_objects.building_block_sequence import BuildingBlockSequence
from lcseq.domain.value_objects.monomer_sequence import MonomerSequence
from lcseq.domain.value_objects.peak_boundaries import PeakBoundaries
from lcseq.domain.value_objects.retention_time import RetentionTime, TimeUnit

__all__ = [
    "BuildingBlockSequence",
    "MonomerSequence",
    "PeakBoundaries",
    "RetentionTime",
    "TimeUnit",
]
