"""Infrastructure layer for LC-Seq.

This module provides infrastructure concerns like data loading.
"""

from .loaders import HDF5CompoundLoader

__all__ = ["HDF5CompoundLoader"]
