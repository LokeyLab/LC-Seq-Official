#!/usr/bin/env python3
"""LC-Seq Consensus Mode Demo - demonstrates hybrid consensus strategy.

This script demonstrates the consensus mode workflow from THEORY.md Section 4.2:
1. Group compounds into equivalence classes (by residue sequence)
2. Aggregate variant signals → consensus
3. Detect peaks on consensus (expensive, once per class)
4. Integrate areas on individual variants (cheap, per variant)

This provides ~3-10× speedup for peak detection while preserving individual
purity measurements.

Usage:
    python examples/analyze_consensus.py --reference "Leu-Pro-Ala"
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
from lcseq.domain.models import HierarchyMode, ConsensusStatus
from lcseq.application.use_cases import ProcessChromatogramsConsensusUseCase
from lcseq.presentation.visualization.plotters import LineageOffsetPlotter
from lcseq.infrastructure import HDF5CompoundLoader


def analyze_consensus(args):
    """Demonstrate consensus mode on reference compound lineage."""
    print("=" * 80)
    print("LC-Seq Consensus Mode Demo")
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
    print(f"✓ Loaded {len(compounds):,} compounds")

    # STEP 2: Find reference compound
    print("\n" + "=" * 80)
    print("STEP 2: Find Reference Compound")
    print("=" * 80)
    search_service = CompoundSearchService()
    reference = search_service.find_by_sequence(compounds, args.reference)

    if reference is None:
        print(f"\n❌ ERROR: Reference compound '{args.reference}' not found!")
        sys.exit(1)

    print(f"✓ Found: {reference.positional_sequence}")
    print(f"  Canonical Sequence: {reference.residue_sequence}")
    print(f"  Level: {reference.monomer_level} (Monomer Mode)")

    # STEP 3: Find lineage
    print("\n" + "=" * 80)
    print("STEP 3: Find Lineage")
    print("=" * 80)
    lineage_finder = LineageFinderService()
    lineage = lineage_finder.find_principal_ideal(
        reference, compounds, HierarchyMode.MONOMER
    )
    print(f"✓ Found: {len(lineage)} compounds")

    # STEP 4: Build hierarchy
    print("\n" + "=" * 80)
    print("STEP 4: Build Hierarchy")
    print("=" * 80)
    builder = HierarchyBuilder()
    hierarchy = builder.build(lineage, HierarchyMode.MONOMER)
    print(f"✓ Built: {hierarchy.size():,} compounds, {hierarchy.edge_count():,} edges")

    # STEP 5: Process with consensus mode
    print("\n" + "=" * 80)
    print("STEP 5: Process Chromatograms (Consensus Mode)")
    print("=" * 80)
    print("Using hybrid consensus strategy (THEORY.md Section 4.2.3):")
    print("  Phase 1: Peak detection on consensus signal (expensive, once per class)")
    print("  Phase 2: Area integration on variants (cheap, per variant)")

    process_use_case = ProcessChromatogramsConsensusUseCase()
    equivalence_classes, representatives, representative_hierarchy = process_use_case.execute(
        compounds=lineage,
        hierarchy=hierarchy,
        # All parameters use defaults from config
    )

    # Build peaks_dict for plotting
    # Must include ALL compounds in representative hierarchy (not just representatives)
    # because plotter needs descendant peaks to compute truncation boundaries
    peaks_dict = {}
    for compound in representative_hierarchy.compounds:
        peaks_dict[compound] = compound.detected_peaks

    total_peaks = sum(len(rep.detected_peaks) for rep in representatives)
    print(f"✓ Detected {total_peaks} total peaks")
    print(f"✓ Plotting {len(representatives)} traces (1 per equivalence class)")

    # Report correlation issues as warnings
    low_corr_classes = [
        eq for eq in equivalence_classes.values()
        if eq.consensus_status == ConsensusStatus.CONSENSUS_INVALID
    ]
    if low_corr_classes:
        print(f"⚠ Warning: {len(low_corr_classes)} classes have low correlation (< 0.8)")

    # STEP 6: Generate visualization
    print("\n" + "=" * 80)
    print("STEP 6: Generate Plot")
    print("=" * 80)
    plot_path = plots_dir / f"lineage_consensus_{reference.residue_sequence}.png"
    plotter = LineageOffsetPlotter()
    plotter.plot(
        representatives,  # Plot representatives (VirtualCompounds or real Compounds)
        peaks_dict,
        output_path=plot_path,
        reference=reference,  # Pass original reference for title (won't be highlighted since not in hierarchy)
        hierarchy_mode=HierarchyMode.MONOMER,
        hierarchy=representative_hierarchy,  # Use representative hierarchy for correct truncation boundaries
    )
    print(f"✓ Plot saved: {plot_path.name}")

    # STEP 7: Analyze results
    print("\n" + "=" * 80)
    print("STEP 7: Analyze Results")
    print("=" * 80)

    n_classes = len(equivalence_classes)
    n_variants = sum(len(eq.compounds) for eq in equivalence_classes.values())

    print(f"✓ Processed {n_classes} equivalence classes ({n_variants} variants)")
    print()

    # Count by status
    high_corr = sum(
        1 for eq in equivalence_classes.values()
        if eq.is_consensus_valid
    )
    low_corr = sum(
        1 for eq in equivalence_classes.values()
        if eq.consensus_status == ConsensusStatus.CONSENSUS_INVALID
    )
    single_variant = sum(
        1 for eq in equivalence_classes.values()
        if eq.consensus_status == ConsensusStatus.NOT_ATTEMPTED
    )

    print(f"  High correlation (≥0.8): {high_corr} classes")
    print(f"  Low correlation (<0.8): {low_corr} classes")
    print(f"  Single variant: {single_variant} classes")
    print()

    # Show details for each class
    print("\nEquivalence Class Details:")
    print("-" * 80)

    for residue_seq in sorted(equivalence_classes.keys()):
        eq_class = equivalence_classes[residue_seq]
        n_var = len(eq_class.compounds)

        # Status icon: ✓ for high correlation, ⚠ for low correlation, - for single variant
        if eq_class.consensus_status == ConsensusStatus.NOT_ATTEMPTED:
            status_icon = "-"
        elif eq_class.is_consensus_valid:
            status_icon = "✓"
        else:
            status_icon = "⚠"

        print(f"\n{status_icon} {residue_seq}")
        print(f"  Variants: {n_var}")

        if eq_class.correlation_min is not None:
            corr_str = f"{eq_class.correlation_min:.3f}"
            if eq_class.is_consensus_valid:
                print(f"  Correlation: {corr_str} (high)")
            else:
                print(f"  Correlation: {corr_str} (low - below 0.8 threshold)")

        if eq_class.consensus_peaks is not None:
            print(f"  Consensus peaks: {len(eq_class.consensus_peaks)}")

    # Export summary CSV (equivalence class level)
    summary_csv_path = args.output / f"consensus_summary_{reference.residue_sequence}.csv"
    with open(summary_csv_path, "w") as f:
        f.write("residue_sequence,n_variants,status,min_correlation,n_consensus_peaks\n")

        for residue_seq in sorted(equivalence_classes.keys()):
            eq_class = equivalence_classes[residue_seq]
            n_var = len(eq_class.compounds)
            status = eq_class.consensus_status.value
            corr = eq_class.correlation_min if eq_class.correlation_min else ""
            n_peaks = len(eq_class.consensus_peaks) if eq_class.is_consensus_valid else ""

            f.write(f"{residue_seq},{n_var},{status},{corr},{n_peaks}\n")

    print(f"\n✓ Summary CSV exported: {summary_csv_path.name}")

    # Export lineage CSV (representatives, same format as individual mode)
    lineage_csv_path = args.output / f"lineage_consensus_{reference.residue_sequence}.csv"
    with open(lineage_csv_path, "w") as f:
        f.write("sequence,monomer_level,n_peaks\n")

        # Sort representatives same way as visualization (maximal→minimal, grouped by canonical sequence)
        representatives_sorted = sorted(representatives, key=lambda c: (-c.monomer_level, c.residue_sequence))

        for cpd in representatives_sorted:
            peaks = peaks_dict.get(cpd, [])
            f.write(f"{cpd.positional_sequence},{cpd.monomer_level},{len(peaks)}\n")

    print(f"✓ Lineage CSV exported: {lineage_csv_path.name}")

    # Summary
    print("\n" + "=" * 80)
    print("✅ COMPLETE")
    print("=" * 80)
    print(f"\nReference: {reference.positional_sequence}")
    print(f"Lineage: {len(lineage)} compounds ({n_variants} total variants)")
    print(f"Equivalence classes: {n_classes}")
    print(f"Traces plotted: {len(representatives)} (1 per equivalence class)")
    print(f"Peaks: {total_peaks} detected")
    print(f"\nCorrelation Quality:")
    print(f"  High (≥0.8): {high_corr}/{n_classes}")
    print(f"  Low (<0.8): {low_corr}/{n_classes}")
    print(f"\nResults: {args.output}")
    print(f"  - {plot_path.name}")
    print(f"  - {lineage_csv_path.name}")
    print(f"  - {summary_csv_path.name}")


def main():
    """Main entry point for consensus mode demo."""
    parser = argparse.ArgumentParser(
        description="LC-Seq Consensus Mode Demo - Hybrid consensus strategy",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python examples/analyze_consensus.py --reference "Leu-Pro-Ala"

Consensus Mode (THEORY.md Section 4.2):
  - Groups positional variants by residue sequence
  - Detects peaks on consensus signal (expensive, once per class)
  - Integrates areas on individual variants (cheap, per variant)
  - Automatic fallback if correlation < 0.8
  - ~3-10× speedup for peak detection
        """,
    )

    parser.add_argument(
        "--reference",
        type=str,
        default="Leu-LA03-Pro-Leu-DLeuMe-DPro-Leu-Leu-DPro",
        help='Reference compound sequence',
    )
    parser.add_argument(
        "--data",
        type=Path,
        default=Path("test_data/processed_data.h5"),
        help="HDF5 data file",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("results/"),
        help="Output directory",
    )

    args = parser.parse_args()
    analyze_consensus(args)


if __name__ == "__main__":
    main()
