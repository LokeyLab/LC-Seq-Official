#!/usr/bin/env python3
"""LC-Seq Lineage Analysis - Config-driven analysis for reference compounds.

Analyzes a reference compound and its complete lineage using configuration from
configs/default.yaml to determine analysis mode (individual vs pooled), hierarchy
mode, and all parameters.

Analysis Modes (from config):
  individual - Process each variant separately (default)
  pooled     - Aggregate positional variants for ~3-10× speedup

Hierarchy Modes (from config):
  building_block - DAG with convergence at block granularity
  monomer        - DAG with convergence at monomer granularity

Usage:
    # Use default configuration (configs/default.yaml)
    python examples/analyze.py --reference "Phe-DNvl-DPhe"

    # Override configuration file
    python examples/analyze.py --reference "Phe-DNvl-DPhe" --config my_config.yaml

    # Override specific parameters
    python examples/analyze.py --reference "Phe-DNvl-DPhe" --variant-mode pooled
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
from lcseq.domain.models import HierarchyMode, AnalysisMode, PoolingStatus
from lcseq.application.use_cases import (
    ProcessChromatogramsUseCase,
    ProcessPooledChromatogramsUseCase,
)
from lcseq.presentation.visualization.plotters import LineageOffsetPlotter
from lcseq.infrastructure import HDF5CompoundLoader, JSONExporter
from lcseq.infrastructure.configuration.yaml_loader import ConfigurationLoader


def analyze_lineage(args, config):
    """Analyze reference compound and its lineage using configuration."""

    # Determine mode from config (can be overridden by args)
    variant_mode = args.variant_mode or config.analysis_mode
    hierarchy_mode = args.hierarchy_mode or config.hierarchy_mode

    mode_name = "Pooled" if variant_mode == AnalysisMode.POOLED else "Individual"

    print("=" * 80)
    print(f"LC-Seq Lineage Analysis ({mode_name} Mode)")
    print("=" * 80)
    print(f"\nReference Compound: {args.reference}")
    print(f"Data: {args.data}")
    print(f"Output: {args.output}")
    print(f"Analysis Mode: {variant_mode.value} (from {'CLI' if args.variant_mode else 'config'})")
    print(
        f"Hierarchy Mode: {hierarchy_mode.value} (from {'CLI' if args.hierarchy_mode else 'config'})"
    )

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
        print("\nAvailable sequences (first 20):")
        for i, cpd in enumerate(compounds[:20]):
            print(f"  {i+1}. {cpd.positional_block_sequence}")
        sys.exit(1)

    print(f"✓ Found: {reference.positional_block_sequence}")
    print(f"  Block Support Sequence: {reference.block_support_sequence}")

    # Display correct level based on hierarchy mode
    if hierarchy_mode == HierarchyMode.BUILDING_BLOCK:
        print(f"  Level: {reference.level} (Block Mode)")
    else:
        print(f"  Level: {reference.monomer_level} (Monomer Mode)")

    # STEP 3: Find lineage
    print("\n" + "=" * 80)
    print("STEP 3: Find Lineage (Principal Ideal ↓X)")
    print("=" * 80)
    print(f"Scanning {len(compounds):,} compounds to find lineage members...")
    lineage_finder = LineageFinderService()
    lineage = lineage_finder.find_principal_ideal(reference, compounds, hierarchy_mode)
    print(f"✓ Found: {len(lineage)} compounds (reference + {len(lineage)-1} descendants)")
    print(
        f"  Reduction: {len(compounds):,} → {len(lineage)} ({100*len(lineage)/len(compounds):.2f}%)"
    )

    # STEP 4: Build hierarchy
    print("\n" + "=" * 80)
    print("STEP 4: Build Minimal Hierarchy")
    print("=" * 80)
    print(f"Building hierarchy for lineage ({len(lineage)} compounds)...")
    print(f"  Mode: {hierarchy_mode.value}")
    builder = HierarchyBuilder()
    hierarchy = builder.build(lineage, hierarchy_mode)
    print(f"✓ Built: {hierarchy.size():,} compounds, {hierarchy.edge_count():,} edges")

    # Show structure
    use_monomer_level = hierarchy_mode == HierarchyMode.MONOMER
    by_level = lineage_finder.group_lineage_by_level(lineage, use_monomer_level=use_monomer_level)
    mode_label = "Monomer Mode" if use_monomer_level else "Block Mode"
    print(f"\n  Structure ({mode_label} levels):")
    for level in sorted(by_level.keys(), reverse=True):
        print(f"    Level {level}: {len(by_level[level])} compounds")

    # STEP 5: Process chromatograms (mode-dependent)
    print("\n" + "=" * 80)
    print("STEP 5: Process Chromatograms")
    print("=" * 80)

    if variant_mode == AnalysisMode.POOLED:
        print("Using hybrid pooled strategy (THEORY.md Section 4.2.3):")
        print("  Phase 1: Peak detection on pooled signal (expensive, once per class)")
        print("  Phase 2: Area integration on variants (cheap, per variant)")

        process_use_case = ProcessPooledChromatogramsUseCase()
        equivalence_classes, processed_compounds, used_hierarchy = process_use_case.execute(
            compounds=lineage,
            hierarchy=hierarchy,
        )

        # Build peaks_dict for plotting
        peaks_dict = {}
        for compound in used_hierarchy.compounds:
            peaks_dict[compound] = compound.detected_peaks

        total_peaks = sum(len(cpd.detected_peaks) for cpd in processed_compounds)
        n_classes = len(equivalence_classes)
        n_variants = sum(len(eq.members) for eq in equivalence_classes.values())

        print(f"✓ Processed: {n_classes} equivalence classes ({n_variants} variants)")
        print(f"  Peaks: {total_peaks} total")
        print(f"  Traces: {len(processed_compounds)} (1 per equivalence class)")

        # Report correlation warnings
        low_corr_classes = [
            eq
            for eq in equivalence_classes.values()
            if eq.pooling_status == PoolingStatus.POOLING_INVALID
        ]
        if low_corr_classes:
            print(f"  ⚠ Warning: {len(low_corr_classes)} classes have low correlation")
    else:
        # Individual mode
        print("Detecting peaks using Discrete Morse Theory + Poisson statistics...")

        process_use_case = ProcessChromatogramsUseCase()
        peaks_dict = process_use_case.execute(
            compounds=lineage,
            hierarchy=hierarchy,
        )

        processed_compounds = lineage
        used_hierarchy = hierarchy
        equivalence_classes = None  # Not used in individual mode

        total_peaks = sum(len(peaks) for peaks in peaks_dict.values())
        print(f"✓ Processed: {len(lineage)} compounds")
        print(f"  Peaks: {total_peaks} total, {total_peaks/len(lineage):.1f} avg/compound")

    # STEP 6: Generate plot
    print("\n" + "=" * 80)
    print("STEP 6: Generate Plot")
    print("=" * 80)

    # Build file suffix with both variant mode and hierarchy mode
    variant_suffix = "_pooled" if variant_mode == AnalysisMode.POOLED else ""
    hierarchy_suffix = "_block" if hierarchy_mode == HierarchyMode.BUILDING_BLOCK else "_monomer"
    mode_suffix = f"{variant_suffix}{hierarchy_suffix}"

    plot_path = plots_dir / f"lineage{mode_suffix}_{reference.block_support_sequence}.png"

    plotter = LineageOffsetPlotter()
    plotter.plot(
        processed_compounds,
        peaks_dict,
        output_path=plot_path,
        reference=reference if variant_mode == AnalysisMode.INDIVIDUAL else reference,
        hierarchy_mode=hierarchy_mode,
        hierarchy=used_hierarchy,
    )
    print(f"✓ Plot saved: {plot_path.name}")

    # STEP 7: Export results
    print("\n" + "=" * 80)
    print("STEP 7: Export Results")
    print("=" * 80)

    # Export comprehensive JSON with all analysis data
    json_exporter = JSONExporter()
    json_path = args.output / f"analysis{mode_suffix}_{reference.block_support_sequence}.json"

    json_exporter.export(
        reference=reference,
        lineage=lineage,
        hierarchy=used_hierarchy,
        peaks_dict=peaks_dict,
        output_path=json_path,
        variant_mode=variant_mode,
        hierarchy_mode=hierarchy_mode,
        equivalence_classes=equivalence_classes,
        data_file=args.data,
        total_library_size=len(compounds),
    )

    print(f"✓ Comprehensive JSON exported: {json_path.name}")
    print(f"  Contains: metadata, hierarchy, compounds, peaks, statistics")

    if variant_mode == AnalysisMode.POOLED and equivalence_classes:
        print(f"  + equivalence classes data")

        # Show pooling statistics
        print("\n  Pooling Statistics:")
        high_corr = sum(1 for eq in equivalence_classes.values() if eq.is_pooling_valid)
        low_corr = sum(
            1
            for eq in equivalence_classes.values()
            if eq.pooling_status == PoolingStatus.POOLING_INVALID
        )
        single = sum(
            1
            for eq in equivalence_classes.values()
            if eq.pooling_status == PoolingStatus.NOT_ATTEMPTED
        )
        print(f"    High correlation (≥0.8): {high_corr}/{n_classes}")
        print(f"    Low correlation (<0.8): {low_corr}/{n_classes}")
        print(f"    Single variant: {single}/{n_classes}")

    # Summary
    print("\n" + "=" * 80)
    print("✅ COMPLETE")
    print("=" * 80)
    print(f"\nReference: {reference.positional_block_sequence}")
    print(f"Lineage: {len(lineage)} compounds")

    if variant_mode == AnalysisMode.POOLED and equivalence_classes:
        print(f"Equivalence classes: {n_classes} ({n_variants} total variants)")

    print(f"Peaks: {total_peaks} detected")
    print(f"\nResults: {args.output}")
    print(f"  - {plot_path.name}")
    print(f"  - {json_path.name}")


def main():
    """Main entry point for lineage analysis."""
    parser = argparse.ArgumentParser(
        description="LC-Seq Lineage Analysis - Config-driven analysis",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Use default configuration (configs/default.yaml)
  python examples/analyze.py --reference "Phe-DNvl-DPhe"

  # Override with custom config file
  python examples/analyze.py --reference "Phe-DNvl-DPhe" --config my_config.yaml

  # Override analysis mode (config: individual → pooled)
  python examples/analyze.py --reference "Phe-DNvl-DPhe" --variant-mode pooled

  # Override hierarchy mode (config: monomer → building_block)
  python examples/analyze.py --reference "Phe-DNvl-DPhe" --hierarchy-mode building_block

Configuration:
  All parameters loaded from configs/default.yaml (Single Source of Truth).
  CLI arguments override config values when provided.

  Key config parameters:
    - analysis.variant_mode: individual or pooled
    - analysis.hierarchy_mode: building_block or monomer
    - detection.*: peak detection parameters
    - pooling.*: pooling parameters (for pooled mode)

Analysis Modes:
  individual - Process each variant separately (default)
  pooled     - Aggregate positional variants (~3-10× speedup)

Hierarchy Modes:
  monomer        - DAG with convergence at monomer granularity (default)
  building_block - DAG with convergence at block granularity
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
        help="HDF5 data file",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("results/"),
        help="Output directory",
    )
    parser.add_argument(
        "--config",
        type=Path,
        default=None,
        help="Configuration file (default: configs/default.yaml)",
    )
    parser.add_argument(
        "--variant-mode",
        choices=["individual", "pooled"],
        default=None,
        help="Override analysis mode from config",
    )
    parser.add_argument(
        "--hierarchy-mode",
        choices=["building_block", "monomer"],
        default=None,
        help="Override hierarchy mode from config",
    )

    args = parser.parse_args()

    # Load configuration
    if args.config:
        config = ConfigurationLoader.load_from_yaml(args.config)
        print(f"✓ Loaded configuration from: {args.config}")
    else:
        config = ConfigurationLoader.get_default_config()
        print("✓ Loaded configuration from: configs/default.yaml")

    # Convert CLI overrides to enums
    if args.variant_mode:
        args.variant_mode = (
            AnalysisMode.POOLED if args.variant_mode == "pooled" else AnalysisMode.INDIVIDUAL
        )

    if args.hierarchy_mode:
        args.hierarchy_mode = (
            HierarchyMode.BUILDING_BLOCK
            if args.hierarchy_mode == "building_block"
            else HierarchyMode.MONOMER
        )

    # Run analysis
    analyze_lineage(args, config)


if __name__ == "__main__":
    main()
