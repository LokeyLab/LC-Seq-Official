#!/usr/bin/env python3
"""LC-Seq Analysis - Lineage analysis for reference compounds.

This script contains NO business logic - it only orchestrates domain services
and handles visualization (a presentation concern, not domain logic).

Analyzes a reference compound and its complete lineage (ancestors + descendants).

Hierarchy Modes:
  building_block - Poset structure based on positional sequences
  monomer        - DAG with convergence based on chemical identity (default, per THEORY.md)

Usage:
    # Analyze reference compound with its lineage
    python examples/analyze.py --reference "Phe-DNvl-DPhe"

    # Use building-block hierarchy mode
    python examples/analyze.py --reference "Phe-DNvl-DPhe" --hierarchy-mode building_block
"""

from pathlib import Path
import sys
import argparse

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from lcseq.domain.services import (
    HierarchyBuilder,
    CompoundSearchService,
    LineageFinderService,
)
from lcseq.domain.models import HierarchyMode
from lcseq.application.use_cases import ProcessChromatogramsUseCase
from lcseq.presentation.visualization.plotters import LineageOffsetPlotter
from lcseq.infrastructure import HDF5CompoundLoader


def analyze_lineage(args):
    """Analyze reference compound and its complete lineage (ancestors + descendants)."""
    print("=" * 80)
    print("LC-Seq Lineage Analysis")
    print("=" * 80)
    print(f"\nReference Compound: {args.reference}")
    print(f"Data: {args.data}")
    print(f"Output: {args.output}")

    args.output.mkdir(exist_ok=True, parents=True)
    plots_dir = args.output / "plots"
    plots_dir.mkdir(exist_ok=True, parents=True)

    # STEP 1: Load data
    print("\n" + "=" * 80)
    print("STEP 1: Load Data")
    print("=" * 80)
    loader = HDF5CompoundLoader()
    compounds = loader.load_all(args.data)

    # STEP 2: Find reference compound (using domain service)
    print("\n" + "=" * 80)
    print("STEP 2: Find Reference Compound")
    print("=" * 80)
    search_service = CompoundSearchService()
    reference = search_service.find_by_sequence(compounds, args.reference)

    if reference is None:
        print(f"\n❌ ERROR: Reference compound '{args.reference}' not found!")
        print("\nAvailable sequences (first 20):")
        for i, cpd in enumerate(compounds[:20]):
            print(f"  {i+1}. {cpd.positional_sequence}")
        sys.exit(1)

    print(f"✓ Found: {reference.positional_sequence}")
    print(f"  Canonical Sequence: {reference.residue_sequence}")

    # Display correct level based on hierarchy mode (THEORY.md Section 3.3)
    if args.hierarchy_mode == "building_block":
        print(f"  Level: {reference.level} (Block Mode)")
    else:  # monomer
        print(f"  Level: {reference.monomer_level} (Monomer Mode)")

    # STEP 3: Find lineage directly (Principal Ideal ↓X)
    print("\n" + "=" * 80)
    print("STEP 3: Find Lineage (Principal Ideal ↓X)")
    print("=" * 80)
    print(f"Scanning {len(compounds):,} compounds to find lineage members...")
    lineage_finder = LineageFinderService()
    lineage = lineage_finder.find_principal_ideal(reference, compounds, args.hierarchy_mode_enum)
    print(f"✓ Found: {len(lineage)} compounds (reference + {len(lineage)-1} descendants)")
    print(f"  Reduction: {len(compounds):,} → {len(lineage)} ({100*len(lineage)/len(compounds):.2f}%)")

    # STEP 4: Build minimal hierarchy (lineage only)
    print("\n" + "=" * 80)
    print("STEP 4: Build Minimal Hierarchy")
    print("=" * 80)
    print(f"Building hierarchy for lineage ({len(lineage)} compounds)...")
    print(f"  Mode: {args.hierarchy_mode}")
    builder = HierarchyBuilder()
    hierarchy = builder.build(lineage, args.hierarchy_mode_enum)
    print(f"✓ Built: {hierarchy.size():,} compounds, {hierarchy.edge_count():,} edges")

    # Show structure (using domain service to group by level)
    use_monomer_level = args.hierarchy_mode == "monomer"
    by_level = lineage_finder.group_lineage_by_level(lineage, use_monomer_level=use_monomer_level)
    mode_label = "Monomer Mode" if use_monomer_level else "Block Mode"
    print(f"\n  Structure ({mode_label} levels):")
    for level in sorted(by_level.keys(), reverse=True):
        print(f"    Level {level}: {len(by_level[level])} compounds")

    # STEP 5: Process chromatograms (using application layer use case)
    print("\n" + "=" * 80)
    print("STEP 5: Process Chromatograms")
    print("=" * 80)
    print("Detecting peaks using Discrete Morse Theory + Poisson statistics...")
    process_use_case = ProcessChromatogramsUseCase()
    peaks_dict = process_use_case.execute(
        compounds=lineage,
        hierarchy=hierarchy,
        # All parameters use defaults from config
    )

    total_peaks = sum(len(peaks) for peaks in peaks_dict.values())
    print(f"✓ Processed: {len(lineage)} compounds")
    print(f"  Peaks: {total_peaks} total, {total_peaks/len(lineage):.1f} avg/compound")

    # STEP 6: Generate plot (using presentation layer)
    print("\n" + "=" * 80)
    print("STEP 6: Generate Plot")
    print("=" * 80)
    plot_path = plots_dir / f"lineage_{reference.residue_sequence}.png"
    plotter = LineageOffsetPlotter()
    plotter.plot(
        lineage,
        peaks_dict,
        output_path=plot_path,
        reference=reference,
        hierarchy_mode=args.hierarchy_mode_enum,
        hierarchy=hierarchy,  # Enable dashed lines for truncation regions
        # min_baseline_sds uses default from config
    )

    # Export CSV (sorted to match visualization order)
    csv_path = args.output / f"lineage_{reference.residue_sequence}.csv"
    with open(csv_path, "w") as f:
        # Header with mode-specific level name (THEORY.md Section 3.3)
        level_attr = "monomer_level" if use_monomer_level else "level"
        level_label = "monomer_level" if use_monomer_level else "block_level"
        f.write(f"sequence,{level_label},n_peaks\n")

        # Sort lineage same way as visualization (maximal→minimal, grouped by canonical sequence)
        lineage_sorted = sorted(lineage, key=lambda c: (-getattr(c, level_attr), c.residue_sequence))

        for cpd in lineage_sorted:
            peaks = peaks_dict.get(cpd, [])
            level_value = getattr(cpd, level_attr)
            f.write(f"{cpd.positional_sequence},{level_value},{len(peaks)}\n")
    print(f"✓ CSV exported: {csv_path.name}")

    # Summary
    print("\n" + "=" * 80)
    print("✅ COMPLETE")
    print("=" * 80)
    print(f"\nReference: {reference.positional_sequence}")
    print(f"Lineage: {len(lineage)} compounds")
    print(f"Peaks: {total_peaks} detected")
    print(f"\nResults: {args.output}")
    print(f"  - {plot_path.name}")
    print(f"  - {csv_path.name}")


# =============================================================================
# Main Entry Point
# =============================================================================


def main():
    """Main entry point for lineage analysis."""
    parser = argparse.ArgumentParser(
        description="LC-Seq Lineage Analysis - Analyze reference compound and descendants",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Analyze reference compound with monomer hierarchy (default, per THEORY.md)
  python examples/analyze.py --reference "Phe-DNvl-DPhe"

  # Use building-block hierarchy (poset structure)
  python examples/analyze.py --reference "Phe-DNvl-DPhe" --hierarchy-mode building_block

  # Custom data file
  python examples/analyze.py --reference "Leu-Pro-Ala" --data my_data.h5

Hierarchy Modes:
  monomer        - DAG with convergence, chemical identity (default, per THEORY.md)
  building_block - Poset structure based on positional sequences

Terminology (per THEORY.md Section 3.1):
  - Reference Compound: The compound currently being analyzed
  - Lineage: All ancestors + descendants + self
  - Descendant: Compound with fewer building blocks
  - Ancestor: Compound with more building blocks

Architecture:
  This script contains NO business logic. It only:
  - Loads data (infrastructure)
  - Calls domain services (HierarchyBuilder, LineageFinderService, PeakDetector)
  - Handles visualization (presentation concern)
        """,
    )

    parser.add_argument(
        "--reference",
        type=str,
        default="Leu-LA03-Pro-Leu-DLeuMe-DPro-Leu-Leu-DPro",
        help='Reference compound sequence (e.g., "Phe-DNvl-DPhe")',
    )
    parser.add_argument(
        "--data",
        type=Path,
        default=Path("test_data/processed_data.h5"),
        help="HDF5 data file (default: test_data/processed_data.h5)",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("results/"),
        help="Output directory (default: examples/lineage_results)",
    )
    parser.add_argument(
        "--hierarchy-mode",
        choices=["building_block", "monomer"],
        default="monomer",
        help="Hierarchy mode: monomer (DAG, default per THEORY.md) or building_block (poset)",
    )

    args = parser.parse_args()

    # Convert hierarchy mode string to enum
    if args.hierarchy_mode == "building_block":
        args.hierarchy_mode_enum = HierarchyMode.BUILDING_BLOCK
    else:  # monomer
        args.hierarchy_mode_enum = HierarchyMode.MONOMER

    # Run lineage analysis
    analyze_lineage(args)


if __name__ == "__main__":
    main()
