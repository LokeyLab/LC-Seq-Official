#!/usr/bin/env python3
"""LC-Seq Library Quality Filter.

Pre-filters a library based on signal quality metrics before analysis.
This is useful for removing low-quality compounds that might compromise
downstream analysis.

Quality Metrics:
    1. Replicate Correlation: Positional variants (same block support sequence)
       should have similar chromatograms. Low correlation indicates quality issues.

    2. Signal Intensity: Total integrated signal. Low counts indicate failed
       synthesis, poor ionization, or sample loss.

    3. Noise Level: Baseline noise relative to signal. High noise makes peak
       detection unreliable.

Key Points:
    - Single-variant classes are always included (can't assess correlation)
    - Filtering happens at equivalence class level (all variants pass or fail)
    - Outputs: filtered HDF5, inclusion CSV, and detailed QC report

Usage:
    # Basic filtering with defaults
    python examples/filter_library.py --data data.h5 --output filtered/

    # Strict correlation threshold
    python examples/filter_library.py --data data.h5 --output filtered/ --min-correlation 0.9

    # Filter bottom 10% by intensity
    python examples/filter_library.py --data data.h5 --output filtered/ --intensity-percentile 0.10

    # Generate only QC report (no filtered files)
    python examples/filter_library.py --data data.h5 --output qc/ --no-hdf5 --no-csv

Output Files:
    filtered/
    ├── filtered_library.h5      # HDF5 with only passing compounds
    ├── included_sequences.csv   # List of passing compound sequences
    └── qc_report.json           # Detailed quality metrics for all classes
"""

from pathlib import Path
import sys
import argparse

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from lcseq.application.use_cases import FilterLibraryUseCase


def main():
    """Main entry point for library filtering."""
    parser = argparse.ArgumentParser(
        description="LC-Seq Library Quality Filter - Pre-filter library based on signal quality",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # Basic filtering with defaults (0.8 correlation, bottom 5% intensity)
    python examples/filter_library.py --data test_data/processed_data.h5 --output filtered/

    # Strict correlation threshold
    python examples/filter_library.py --data data.h5 --output filtered/ --min-correlation 0.9

    # Filter bottom 10% by intensity
    python examples/filter_library.py --data data.h5 --output filtered/ --intensity-percentile 0.10

    # Use absolute intensity threshold instead of percentile
    python examples/filter_library.py --data data.h5 --output filtered/ --intensity-min 1000000

    # Enable noise filter
    python examples/filter_library.py --data data.h5 --output filtered/ --max-noise-ratio 0.5

    # Generate only QC report (inspect quality before filtering)
    python examples/filter_library.py --data data.h5 --output qc/ --no-hdf5 --no-csv

Output:
    The filter produces three output files:
    - filtered_library.h5: HDF5 with only passing compounds
    - included_sequences.csv: CSV list of passing compound sequences
    - qc_report.json: Detailed quality metrics for all equivalence classes

Quality Metrics:
    For each equivalence class (positional variants), the filter computes:
    - min_correlation: Minimum pairwise Pearson correlation between variants
    - mean_total_signal: Average total integrated signal across variants
    - mean_noise_ratio: Average (noise_std / median_signal) across variants
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
        default=Path("filtered/"),
        help="Output directory for filtered files",
    )

    # Correlation filter
    parser.add_argument(
        "--min-correlation",
        type=float,
        default=0.8,
        help="Minimum replicate correlation within equivalence class (default: 0.8)",
    )

    # Intensity filter
    parser.add_argument(
        "--intensity-percentile",
        type=float,
        default=0.05,
        help="Exclude compounds below this percentile (default: 0.05 = bottom 5%%)",
    )
    parser.add_argument(
        "--intensity-min",
        type=float,
        default=None,
        help="Absolute minimum total signal (overrides percentile)",
    )

    # Noise filter
    parser.add_argument(
        "--max-noise-ratio",
        type=float,
        default=None,
        help="Maximum noise_std / median_signal ratio (default: disabled)",
    )

    # Output options
    parser.add_argument(
        "--no-hdf5",
        action="store_true",
        help="Skip generating filtered HDF5 file",
    )
    parser.add_argument(
        "--no-csv",
        action="store_true",
        help="Skip generating inclusion CSV",
    )
    parser.add_argument(
        "--no-report",
        action="store_true",
        help="Skip generating QC report JSON",
    )

    args = parser.parse_args()

    # Print header
    print("=" * 80)
    print("LC-Seq Library Quality Filter")
    print("=" * 80)
    print(f"\nData: {args.data}")
    print(f"Output: {args.output}")
    print(f"\nFilter Settings:")
    print(f"  Min correlation: {args.min_correlation}")
    if args.intensity_min is not None:
        print(f"  Intensity threshold: {args.intensity_min:.0f} (absolute)")
    else:
        print(f"  Intensity threshold: {args.intensity_percentile*100:.0f}th percentile")
    if args.max_noise_ratio is not None:
        print(f"  Max noise ratio: {args.max_noise_ratio}")
    else:
        print(f"  Noise filter: disabled")

    # Run filtering
    use_case = FilterLibraryUseCase()
    result = use_case.execute(
        hdf5_path=args.data,
        output_dir=args.output,
        min_correlation=args.min_correlation,
        intensity_percentile=args.intensity_percentile if args.intensity_min is None else None,
        intensity_absolute=args.intensity_min,
        max_noise_ratio=args.max_noise_ratio,
        generate_hdf5=not args.no_hdf5,
        generate_csv=not args.no_csv,
        generate_report=not args.no_report,
    )

    # Print final summary
    print("\n" + "=" * 80)
    print("FILTERING SUMMARY")
    print("=" * 80)

    retention_rate = result.passed_compounds / result.total_compounds * 100
    print(f"\nRetention Rate: {retention_rate:.1f}%")
    print(f"  Passed: {result.passed_compounds:,} compounds")
    print(f"  Filtered: {result.filtered_compounds:,} compounds")

    print(f"\nEquivalence Class Breakdown:")
    print(f"  Total classes: {result.total_equivalence_classes:,}")
    print(f"  Passed: {result.passed_equivalence_classes:,}")
    print(f"  Single-variant (auto-passed): {result.single_variant_classes:,}")

    if result.filtered_equivalence_classes > 0:
        print(f"\nFilter Failure Breakdown:")
        print(f"  Failed correlation: {result.failed_correlation:,}")
        print(f"  Failed intensity: {result.failed_intensity:,}")
        print(f"  Failed noise: {result.failed_noise:,}")

    print(f"\nOutput Files:")
    if result.output_hdf5:
        print(f"  {result.output_hdf5}")
    if result.inclusion_csv:
        print(f"  {result.inclusion_csv}")
    if result.qc_report:
        print(f"  {result.qc_report}")

    print(f"\nProcessing time: {result.processing_time_seconds:.1f}s")

    return 0


if __name__ == "__main__":
    sys.exit(main())
