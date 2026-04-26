#!/usr/bin/env python3
"""LC-Seq Full Library Analysis.

Processes entire library in single run with global hierarchy.
Each compound is classified exactly once via bottom-up processing.

Key Benefits:
- No redundant processing: Each compound classified exactly once
- Efficient hierarchy: O(k^2) instead of O(n^2) for edge detection
- Pooled peak detection: One detection per equivalence class
- Incremental output: Results written as processing completes
- Resumable: Checkpoint support for interrupted runs

Usage:
    # Basic usage - process all compounds, generate diagnostic plots
    python examples/analyze_library.py --data data.h5 --output results/

    # Skip diagnostic plots for faster processing
    python examples/analyze_library.py --data data.h5 --output results/ --no-diagnostics

    # Use 8 workers for faster diagnostic generation
    python examples/analyze_library.py --data data.h5 --output results/ --workers 8

    # Resume interrupted processing
    python examples/analyze_library.py --data data.h5 --output results/ --resume

    # Use monomer-level hierarchy
    python examples/analyze_library.py --data data.h5 --output results/ --hierarchy-mode monomer

Output Files:
    results/
    ├── library_analysis.jsonl     # One JSON record per compound (streaming)
    ├── summary.json               # Analysis summary and statistics
    └── diagnostics/               # Diagnostic plots (if enabled)
        ├── Leu-Pro-Val.png
        ├── Leu-Pro.png
        └── ...
"""

from pathlib import Path
import sys
import argparse

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from lcseq.application.use_cases import ProcessLibraryUseCase
from lcseq.domain.models import HierarchyMode
from lcseq.infrastructure.configuration.yaml_loader import ConfigurationLoader


def main():
    """Main entry point for library analysis."""
    parser = argparse.ArgumentParser(
        description="LC-Seq Full Library Analysis - Process entire library with global hierarchy",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # Basic usage with diagnostic plots
    python examples/analyze_library.py --data test_data/processed_data.h5 --output results/library/

    # Fast mode without diagnostics
    python examples/analyze_library.py --data data.h5 --output results/ --no-diagnostics

    # Resume interrupted processing
    python examples/analyze_library.py --data data.h5 --output results/ --resume

    # Use monomer-level hierarchy for finer granularity
    python examples/analyze_library.py --data data.h5 --output results/ --hierarchy-mode monomer

Output:
    The analysis produces:
    - library_analysis.jsonl: One JSON record per compound (streaming format)
    - summary.json: Overall statistics and any errors
    - diagnostics/: Per-compound diagnostic plots (if enabled)

Performance:
    For a typical 64k compound library:
    - Without diagnostics: ~10-15 minutes
    - With diagnostics: ~30-60 minutes (I/O bound)
    - Equivalence class reduction: 64k -> ~1.5k classes
        """,
    )

    parser.add_argument(
        "--data",
        type=Path,
        default=Path("test_data/processed_data.h5"),
        help="HDF5 data file containing compound chromatograms",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("results/library/"),
        help="Output directory for results and diagnostics",
    )
    parser.add_argument(
        "--config",
        type=Path,
        default=None,
        help="Configuration file (default: configs/default.yaml)",
    )
    parser.add_argument(
        "--hierarchy-mode",
        choices=["building_block", "monomer"],
        default=None,
        help="Hierarchy mode (default: from config)",
    )
    parser.add_argument(
        "--no-diagnostics",
        action="store_true",
        help="Skip generating diagnostic plots (faster)",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=4,
        help="Number of worker threads for parallel diagnostic generation (default: 4)",
    )
    parser.add_argument(
        "--resume",
        action="store_true",
        help="Resume from checkpoint if available",
    )
    parser.add_argument(
        "--clpe-reference",
        type=Path,
        default=None,
        help="Path to CSV with AlogP and scaffold data for cLPE peak validation",
    )
    parser.add_argument(
        "--clpe-t0",
        type=float,
        default=None,
        help="Dead time (t0) in minutes for cLPE LogK calculation (default: NULL peak RT)",
    )
    parser.add_argument(
        "--no-clpe-reselect",
        action="store_true",
        help="Don't re-select peaks for cLPE outliers (only report)",
    )

    args = parser.parse_args()

    # Load configuration
    if args.config:
        config = ConfigurationLoader.load_from_yaml(args.config)
        print(f"Loaded configuration from: {args.config}")
    else:
        config = ConfigurationLoader.get_default_config()
        print("Loaded configuration from: configs/default.yaml")

    # Determine hierarchy mode
    if args.hierarchy_mode:
        hierarchy_mode = (
            HierarchyMode.BUILDING_BLOCK
            if args.hierarchy_mode == "building_block"
            else HierarchyMode.MONOMER
        )
    else:
        hierarchy_mode = config.hierarchy_mode

    # Print header
    print("=" * 80)
    print("LC-Seq Full Library Analysis")
    print("=" * 80)
    print(f"\nData: {args.data}")
    print(f"Output: {args.output}")
    print(f"Hierarchy Mode: {hierarchy_mode.value}")
    print(f"Diagnostics: {'Disabled' if args.no_diagnostics else f'Enabled ({args.workers} workers)'}")
    print(f"Resume: {'Yes' if args.resume else 'No'}")
    if args.clpe_reference:
        print(f"cLPE Validation: {args.clpe_reference}")
        t0_str = f"{args.clpe_t0} min" if args.clpe_t0 else "NULL peak RT"
        print(f"  t0 = {t0_str}, reselect = {not args.no_clpe_reselect}")

    # Run analysis
    use_case = ProcessLibraryUseCase()
    result = use_case.execute(
        hdf5_path=args.data,
        output_dir=args.output,
        # Peak detection parameters (from config)
        alpha=config.peak_detection_params['alpha'],
        alpha_product=config.peak_detection_params['alpha_product'],
        prominence_percentile=config.peak_detection_params['prominence_percentile'],
        min_snr=config.peak_detection_params['min_snr'],
        min_baseline_sds=config.peak_detection_params['min_baseline_sds'],
        signal_variant=config.peak_detection_params['signal_variant'],
        min_dispersion_r=config.peak_detection_params['min_dispersion_r'],
        sigma_clip_sigma=config.peak_detection_params['sigma_clip_sigma'],
        # Peak classification parameters (from config)
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
        # Optional parameters
        generate_diagnostics=not args.no_diagnostics,
        resume=args.resume,
        hierarchy_mode=hierarchy_mode,
        n_workers=args.workers,
        # cLPE validation (optional)
        clpe_reference_csv=args.clpe_reference,
        clpe_t0=args.clpe_t0,
        clpe_reselect_peaks=not args.no_clpe_reselect,
    )

    # Print summary
    print("\n" + "=" * 80)
    print("ANALYSIS COMPLETE")
    print("=" * 80)
    print(f"\nTotal compounds: {result.total_compounds:,}")
    print(f"Equivalence classes: {result.equivalence_classes_count:,}")
    print(f"Hierarchy edges: {result.hierarchy_edges:,}")
    print(f"Total peaks: {result.total_peaks:,}")
    print(f"Processing time: {result.processing_time_seconds/60:.1f} minutes")
    print(f"Rate: {result.total_compounds/result.processing_time_seconds:.1f} compounds/second")

    print(f"\nPooling Statistics:")
    print(f"  High correlation: {result.high_correlation_classes}")
    print(f"  Low correlation: {result.low_correlation_classes}")
    print(f"  Single variant: {result.single_variant_classes}")

    print(f"\nOutput:")
    print(f"  Results: {result.output_dir / 'library_analysis.jsonl'}")
    print(f"  Summary: {result.output_dir / 'summary.json'}")
    if not args.no_diagnostics:
        print(f"  Diagnostics: {result.diagnostics_generated:,} plots")

    if result.errors:
        print(f"\n⚠ Errors: {len(result.errors)} (see summary.json for details)")

    return 0 if not result.errors else 1


if __name__ == "__main__":
    sys.exit(main())
