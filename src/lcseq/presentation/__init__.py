"""Presentation Layer for LC-Seq.

This layer contains all external-facing adapters that present
the system to users and external systems.

Components
----------
- cli: Command-line interface (Typer-based)
- visualization: Matplotlib-based plotters

Architecture Notes
------------------
The presentation layer:
- Depends on the application layer
- Transforms external requests into use case calls
- Transforms domain results into external formats (CLI output, plots)
- Contains NO business logic
- Is the outermost layer in Clean Architecture

References
----------
ARCHITECTURE.md: Presentation Layer (formerly Interface Layer)
"""

__all__ = []
