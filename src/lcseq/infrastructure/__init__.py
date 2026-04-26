"""Infrastructure layer for LC-Seq.

This module provides infrastructure concerns like data loading and result export.

Note: Exports are imported lazily to avoid circular dependencies.
The circular import chain was:
  config -> infrastructure -> exporters -> application -> config

To use exporters, import them directly:
  from lcseq.infrastructure.exporters import JSONExporter, CSVExporter
"""

from .loaders import HDF5CompoundLoader

# Lazy imports for __all__ - these are only imported when explicitly requested
# to avoid circular dependency with application layer
def __getattr__(name):
    if name == "JSONExporter":
        from .exporters import JSONExporter
        return JSONExporter
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")

__all__ = ["HDF5CompoundLoader", "JSONExporter"]
