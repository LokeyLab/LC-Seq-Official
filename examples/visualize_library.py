#!/usr/bin/env python3
"""LC-Seq Library Visualization - Generate plots from JSONL analysis output.

Creates visualizations from library analysis JSONL output without requiring
raw chromatogram data.

Usage:
    python examples/visualize_library.py \
        --input results/library/library_analysis.jsonl \
        --output results/library/visualizations/ \
        --plots purity-histogram level-distribution peak-breakdown purity-by-level
"""

from pathlib import Path
import sys
import argparse
import json
from collections import defaultdict
from typing import Dict, List, Any

import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))
from lcseq.domain.services.purity_calculator import PurityCalculator


def load_jsonl(path: Path) -> List[Dict]:
    """Load JSONL records."""
    records = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line:
                records.append(json.loads(line))
    return records


def plot_purity_histogram(records: List[Dict], output_dir: Path) -> None:
    """Generate purity distribution histogram."""
    purities = [PurityCalculator.calculate_from_peaks(r.get("peaks", [])) for r in records]
    purities = [p for p in purities if p > 0]  # Exclude zero purity

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.hist(purities, bins=50, edgecolor='black', alpha=0.7)
    ax.set_xlabel("Purity", fontsize=12)
    ax.set_ylabel("Count", fontsize=12)
    ax.set_title(f"Purity Distribution (n={len(purities)})", fontsize=14)
    ax.axvline(np.median(purities), color='red', linestyle='--', label=f'Median: {np.median(purities):.2f}')
    ax.legend()

    plt.tight_layout()
    fig.savefig(output_dir / "purity_histogram.png", dpi=150)
    plt.close(fig)
    print(f"  Saved: purity_histogram.png")


def plot_level_distribution(records: List[Dict], output_dir: Path) -> None:
    """Generate compounds per level bar chart."""
    level_counts = defaultdict(int)
    for r in records:
        level = r.get("level", 0)
        level_counts[level] += 1

    levels = sorted(level_counts.keys())
    counts = [level_counts[l] for l in levels]

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.bar(levels, counts, edgecolor='black', alpha=0.7)
    ax.set_xlabel("Level", fontsize=12)
    ax.set_ylabel("Compound Count", fontsize=12)
    ax.set_title(f"Compounds by Level (total={sum(counts)})", fontsize=14)
    ax.set_xticks(levels)

    # Add count labels on bars
    for i, (level, count) in enumerate(zip(levels, counts)):
        ax.text(level, count + max(counts)*0.01, str(count), ha='center', va='bottom', fontsize=10)

    plt.tight_layout()
    fig.savefig(output_dir / "level_distribution.png", dpi=150)
    plt.close(fig)
    print(f"  Saved: level_distribution.png")


def plot_peak_type_breakdown(records: List[Dict], output_dir: Path) -> None:
    """Generate stacked bar chart of peak types by level."""
    # Count peak types per level
    level_peak_types = defaultdict(lambda: defaultdict(int))

    for r in records:
        level = r.get("level", 0)
        for peak in r.get("peaks", []):
            if peak.get("is_accepted", True):
                peak_type = peak.get("classification", "UNKNOWN")
                level_peak_types[level][peak_type] += 1

    if not level_peak_types:
        print("  Skipping peak_breakdown.png - no peak data")
        return

    levels = sorted(level_peak_types.keys())
    peak_types = ["PUTATIVE_PRODUCT", "TRUNCATION", "TRUNCATION_UNKNOWN", "NULL", "UNKNOWN"]
    colors = ["green", "orange", "darkorange", "gray", "red"]

    fig, ax = plt.subplots(figsize=(12, 6))

    bottom = np.zeros(len(levels))
    for peak_type, color in zip(peak_types, colors):
        counts = [level_peak_types[l].get(peak_type, 0) for l in levels]
        ax.bar(levels, counts, bottom=bottom, label=peak_type, color=color, edgecolor='black', alpha=0.8)
        bottom += np.array(counts)

    ax.set_xlabel("Level", fontsize=12)
    ax.set_ylabel("Peak Count", fontsize=12)
    ax.set_title("Peak Type Distribution by Level", fontsize=14)
    ax.set_xticks(levels)
    ax.legend(loc='upper right')

    plt.tight_layout()
    fig.savefig(output_dir / "peak_breakdown.png", dpi=150)
    plt.close(fig)
    print(f"  Saved: peak_breakdown.png")


def plot_purity_by_level(records: List[Dict], output_dir: Path) -> None:
    """Generate box plot of purity by level."""
    level_purities = defaultdict(list)

    for r in records:
        level = r.get("level", 0)
        purity = PurityCalculator.calculate_from_peaks(r.get("peaks", []))
        if purity > 0:
            level_purities[level].append(purity)

    levels = sorted(level_purities.keys())
    data = [level_purities[l] for l in levels]

    fig, ax = plt.subplots(figsize=(12, 6))
    bp = ax.boxplot(data, labels=levels, patch_artist=True)

    for patch in bp['boxes']:
        patch.set_facecolor('lightblue')

    ax.set_xlabel("Level", fontsize=12)
    ax.set_ylabel("Purity", fontsize=12)
    ax.set_title("Purity Distribution by Level", fontsize=14)
    ax.set_ylim(0, 1.05)

    # Add median labels
    for i, level in enumerate(levels):
        median = np.median(level_purities[level])
        ax.text(i+1, 1.02, f'{median:.2f}', ha='center', va='bottom', fontsize=9)

    plt.tight_layout()
    fig.savefig(output_dir / "purity_by_level.png", dpi=150)
    plt.close(fig)
    print(f"  Saved: purity_by_level.png")


def print_summary_statistics(records: List[Dict]) -> None:
    """Print summary statistics."""
    print("\n" + "="*60)
    print("SUMMARY STATISTICS")
    print("="*60)

    purities = [PurityCalculator.calculate_from_peaks(r.get("peaks", [])) for r in records]
    purities = [p for p in purities if p > 0]

    level_counts = defaultdict(int)
    for r in records:
        level_counts[r.get("level", 0)] += 1

    print(f"\nTotal compounds: {len(records)}")
    print(f"Compounds with purity > 0: {len(purities)}")
    print(f"Levels: {min(level_counts.keys())} to {max(level_counts.keys())}")

    print(f"\nPurity statistics:")
    print(f"  Mean:   {np.mean(purities):.3f}")
    print(f"  Median: {np.median(purities):.3f}")
    print(f"  Std:    {np.std(purities):.3f}")
    print(f"  Min:    {np.min(purities):.3f}")
    print(f"  Max:    {np.max(purities):.3f}")

    print(f"\nCompounds per level:")
    for level in sorted(level_counts.keys()):
        print(f"  Level {level}: {level_counts[level]:,}")


PLOT_FUNCTIONS = {
    "purity-histogram": plot_purity_histogram,
    "level-distribution": plot_level_distribution,
    "peak-breakdown": plot_peak_type_breakdown,
    "purity-by-level": plot_purity_by_level,
}


def main():
    parser = argparse.ArgumentParser(
        description="Generate visualizations from library analysis JSONL",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument("--input", type=Path, required=True, help="JSONL input file")
    parser.add_argument("--output", type=Path, default=Path("visualizations"), help="Output directory")
    parser.add_argument(
        "--plots",
        nargs="+",
        choices=list(PLOT_FUNCTIONS.keys()) + ["all"],
        default=["all"],
        help="Plots to generate"
    )
    parser.add_argument("--summary", action="store_true", help="Print summary statistics")

    args = parser.parse_args()

    if not args.input.exists():
        print(f"Error: Input file not found: {args.input}")
        sys.exit(1)

    args.output.mkdir(parents=True, exist_ok=True)

    print(f"Loading: {args.input}")
    records = load_jsonl(args.input)
    print(f"Loaded {len(records):,} records")

    # Determine which plots to generate
    if "all" in args.plots:
        plots = list(PLOT_FUNCTIONS.keys())
    else:
        plots = args.plots

    print(f"\nGenerating {len(plots)} plots...")
    for plot_name in plots:
        PLOT_FUNCTIONS[plot_name](records, args.output)

    if args.summary or "all" in args.plots:
        print_summary_statistics(records)

    print(f"\nDone! Output saved to: {args.output}")


if __name__ == "__main__":
    main()
