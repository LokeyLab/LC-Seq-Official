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
import matplotlib.pyplot as plt

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
from lcseq.presentation.visualization.plotters import (
    LineageOffsetPlotter,
    CompoundDiagnosticPlotter,
    LineageHeatmapPlotter,
    generate_diagnostic_plots,
)
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

    # STEP 4: Load cLPE reference data (optional)
    clpe_validator = None
    alogp_map = None
    scaffold_map = None

    if args.clpe_reference:
        print("\n" + "=" * 80)
        print("STEP 4: Load cLPE Reference Data")
        print("=" * 80)

        from lcseq.infrastructure.loaders.clpe_reference_loader import CLPEReferenceLoader
        from lcseq.domain.services.clpe_validator import CLPEValidator

        loader = CLPEReferenceLoader()
        ref_data = loader.load(args.clpe_reference)
        alogp_map = ref_data.alogp_map
        scaffold_map = ref_data.scaffold_map

        # Use CLI value if provided, otherwise use config value
        effective_clpe_threshold = (
            args.clpe_threshold if args.clpe_threshold is not None
            else config.validation_params['clpe_outlier_threshold']
        )
        effective_min_group_size = config.validation_params['clpe_min_group_size']

        # t0 will be set from L0 peak during classification
        clpe_validator = CLPEValidator(
            t0=1.0,  # Will be updated from L0
            outlier_threshold=effective_clpe_threshold,
            min_group_size=effective_min_group_size,
        )

        print(f"✓ Loaded: {len(alogp_map)} AlogP values, {len(set(scaffold_map.values()))} unique scaffolds")
        print(f"  Threshold: {effective_clpe_threshold} z-score")

    # STEP 5: Build hierarchy
    step_number = 5 if args.clpe_reference else 4
    print("\n" + "=" * 80)
    print(f"STEP {step_number}: Build Minimal Hierarchy")
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

    # STEP 6: Process chromatograms (mode-dependent)
    # Preprocessing is handled internally by use cases when preprocessing_params is passed
    step_number = 6 if args.clpe_reference else 5
    print("\n" + "=" * 80)
    print(f"STEP {step_number}: Process Chromatograms")
    print("=" * 80)

    # Show detection parameters being used
    print("\nDetection Parameters (from config):")
    print(f"  alpha: {config.peak_detection_params['alpha']}")
    print(f"  alpha_product: {config.peak_detection_params['alpha_product']}")
    print(f"  min_baseline_sds: {config.peak_detection_params['min_baseline_sds']}")
    print(f"  min_snr: {config.peak_detection_params['min_snr']}")
    print(f"  prominence_percentile: {config.peak_detection_params['prominence_percentile']}")
    print(f"  signal_variant: {config.peak_detection_params['signal_variant']}")

    if variant_mode == AnalysisMode.POOLED:
        print("Using hybrid pooled strategy (THEORY.md Section 4.2.3):")
        print("  Phase 1: Peak detection on pooled signal (expensive, once per class)")
        print("  Phase 2: Area integration on variants (cheap, per variant)")
        if args.clpe_reference:
            print("  cLPE validation: Enabled")

        process_use_case = ProcessPooledChromatogramsUseCase()
        equivalence_classes, processed_compounds, used_hierarchy, clpe_stats = process_use_case.execute(
            compounds=lineage,
            hierarchy=hierarchy,
            # Peak detection parameters (from config)
            alpha=config.peak_detection_params['alpha'],
            prominence_percentile=config.peak_detection_params['prominence_percentile'],
            min_snr=config.peak_detection_params['min_snr'],
            min_baseline_sds=config.peak_detection_params['min_baseline_sds'],
            signal_variant=config.peak_detection_params['signal_variant'],
            min_dispersion_r=config.peak_detection_params['min_dispersion_r'],
            sigma_clip_sigma=config.peak_detection_params['sigma_clip_sigma'],
            # Peak classification parameters (from config)
            alpha_product=config.peak_detection_params['alpha_product'],
            truncation_margin=config.classification_params['truncation_margin'],
            peak_matching_tolerance=config.classification_params['peak_matching_tolerance'],
            hungarian_min_threshold=config.classification_params['hungarian_min_threshold'],
            # Pooling parameters (from config)
            correlation_threshold=config.validation_params['correlation_threshold'],
            aggregation_method=config.validation_params['aggregation_method'],
            # Validation parameters (from config)
            include_rejected=config.peak_detection_params['include_rejected'],
            clpe_outlier_threshold=config.validation_params['clpe_outlier_threshold'],
            clpe_min_group_size=config.validation_params['clpe_min_group_size'],
            # Preprocessing parameters (from config)
            preprocessing_params=config.preprocessing_params,
            # Optional cLPE reference (user-provided)
            clpe_reference_csv=args.clpe_reference,
            clpe_t0=None,  # Will be derived from L0 peak RT
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
        if args.clpe_reference:
            print("  cLPE validation: Enabled")

        process_use_case = ProcessChromatogramsUseCase()
        peaks_dict, clpe_stats = process_use_case.execute(
            compounds=lineage,
            hierarchy=hierarchy,
            # Peak detection parameters (from config)
            alpha=config.peak_detection_params['alpha'],
            prominence_percentile=config.peak_detection_params['prominence_percentile'],
            min_snr=config.peak_detection_params['min_snr'],
            min_baseline_sds=config.peak_detection_params['min_baseline_sds'],
            signal_variant=config.peak_detection_params['signal_variant'],
            min_dispersion_r=config.peak_detection_params['min_dispersion_r'],
            sigma_clip_sigma=config.peak_detection_params['sigma_clip_sigma'],
            # Peak classification parameters (from config)
            alpha_product=config.peak_detection_params['alpha_product'],
            truncation_margin=config.classification_params['truncation_margin'],
            peak_matching_tolerance=config.classification_params['peak_matching_tolerance'],
            hungarian_min_threshold=config.classification_params['hungarian_min_threshold'],
            # Validation parameters
            include_rejected=config.peak_detection_params['include_rejected'],
            # Preprocessing parameters
            preprocessing_params=config.preprocessing_params,
            # Optional cLPE validator
            clpe_validator=clpe_validator,
            alogp_map=alogp_map,
            scaffold_map=scaffold_map,
        )

        processed_compounds = lineage
        used_hierarchy = hierarchy
        equivalence_classes = None  # Not used in individual mode

        total_peaks = sum(len(peaks) for peaks in peaks_dict.values())
        print(f"✓ Processed: {len(lineage)} compounds")
        print(f"  Peaks: {total_peaks} total, {total_peaks/len(lineage):.1f} avg/compound")

    # Report cLPE statistics if enabled
    if clpe_stats:
        print("\n  cLPE Statistics:")
        if clpe_stats.get("t0"):
            print(f"    t0 (dead time): {clpe_stats['t0']:.2f} min")

        levels_stats = clpe_stats.get("levels", {})
        if levels_stats:
            total_validated = sum(s.get("validated", 0) for s in levels_stats.values())
            total_outliers = sum(s.get("outliers", 0) for s in levels_stats.values())
            total_reselected = sum(s.get("reselected", 0) for s in levels_stats.values())

            print(f"    Validated: {total_validated} compounds")
            print(f"    Outliers: {total_outliers} ({100*total_outliers/total_validated:.1f}%)" if total_validated > 0 else "    Outliers: 0")
            print(f"    Re-selected: {total_reselected}")

            # Per-level details
            for level in sorted(levels_stats.keys()):
                level_stats = levels_stats[level]
                n_validated = level_stats.get("validated", 0)
                n_outliers = level_stats.get("outliers", 0)
                n_reselected = level_stats.get("reselected", 0)
                if n_validated > 0:
                    print(f"      L{level}: {n_validated} validated, {n_outliers} outliers, {n_reselected} reselected")

    # STEP 7: Generate plot
    step_number = 7 if args.clpe_reference else 6
    print("\n" + "=" * 80)
    print(f"STEP {step_number}: Generate Plot")
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

    # STEP 7b/8b: Generate Diagnostic Plots (optional)
    n_generated = 0
    if args.diagnostics:
        step_number_b = f"{step_number}b"
        print("\n" + "=" * 80)
        print(f"STEP {step_number_b}: Generate Diagnostic Plots")
        print("=" * 80)

        diagnostics_dir = plots_dir / f"diagnostics{mode_suffix}"
        print(f"Generating diagnostic plots for {len(processed_compounds)} compounds...")
        print(f"  Output directory: {diagnostics_dir}")

        n_generated = generate_diagnostic_plots(
            compounds=processed_compounds,
            hierarchy=used_hierarchy,
            output_dir=diagnostics_dir,
            filename_attr="block_support_sequence",
        )
        print(f"✓ Generated {n_generated} diagnostic plots")

    # STEP 7c/8c: Generate Purity Heatmap (optional)
    if args.heatmap:
        step_number_c = f"{step_number}c"
        print("\n" + "=" * 80)
        print(f"STEP {step_number_c}: Generate Purity Heatmap")
        print("=" * 80)

        heatmap_path = plots_dir / f"heatmap{mode_suffix}_{reference.block_support_sequence}.png"
        print(f"Building peak distribution matrix for {len(processed_compounds)} compounds...")

        heatmap_plotter = LineageHeatmapPlotter()
        heatmap_fig = heatmap_plotter.plot_from_compounds(
            compounds=processed_compounds,
            hierarchy=used_hierarchy,
            title=f"Peak Distribution: {reference.block_support_sequence}",
        )
        heatmap_fig.savefig(heatmap_path, dpi=150, bbox_inches='tight')
        plt.close(heatmap_fig)
        print(f"✓ Heatmap saved: {heatmap_path.name}")

    # STEP 8 or 9: Export results
    step_number_export = 8 if args.clpe_reference else 7
    print("\n" + "=" * 80)
    print(f"STEP {step_number_export}: Export Results")
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
    if args.diagnostics:
        print(f"  - diagnostics{mode_suffix}/ ({n_generated} plots)")
    if args.heatmap:
        print(f"  - heatmap{mode_suffix}_{reference.block_support_sequence}.png")


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

  # Enable cLPE validation with custom reference data
  python examples/analyze.py --reference "Phe-DNvl-DPhe" --clpe-reference test_data/raw_data.csv

  # Enable cLPE with custom outlier threshold
  python examples/analyze.py --reference "Phe-DNvl-DPhe" --clpe-reference test_data/raw_data.csv --clpe-threshold 3.0

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

cLPE Validation (optional):
  --clpe-reference - CSV with AlogP and scaffold data
  --clpe-threshold - Z-score threshold for outlier detection (default: from config)

  When enabled, validates product peaks using LogK ~ AlogP correlation
  and re-selects outlier peaks from UNKNOWN candidates. Requires reference
  CSV with columns: compound_id, AlogP, scaffold (or All Stereochem).
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
    parser.add_argument(
        "--diagnostics",
        action="store_true",
        help="Generate diagnostic plots for each compound",
    )
    parser.add_argument(
        "--heatmap",
        action="store_true",
        help="Generate purity heatmap showing peak distribution across lineage",
    )
    parser.add_argument(
        "--clpe-reference",
        type=Path,
        default=None,
        help="CSV with compound_id, AlogP, scaffold columns for cLPE validation"
    )
    parser.add_argument(
        "--clpe-threshold",
        type=float,
        default=None,
        help="Z-score threshold for cLPE outlier detection (default: from config)"
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
