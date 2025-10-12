#!/usr/bin/env python3
"""Check if peaks exist in maximal compound signal at expected positions."""

import sys
import numpy as np
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent / "src"))

from lcseq.domain.services import CompoundSearchService, LineageFinderService
from lcseq.domain.models import HierarchyMode
from lcseq.infrastructure import HDF5CompoundLoader
from lcseq.domain.entities.peak import PeakType

# Load data
loader = HDF5CompoundLoader()
compounds = loader.load_all("test_data/processed_data.h5")

# Find compounds
search_service = CompoundSearchService()
maximal = search_service.find_by_sequence(compounds, "Leu-LA03-Pro-Leu-DLeuMe-DPro-Leu-Leu-DPro")
null = search_service.find_by_sequence(compounds, "AgxNull-AgxNull-AgxNull")
level1 = search_service.find_by_sequence(compounds, "AgxNull-AgxNull-DLeuMe")

print("Maximal compound signal analysis:")
print("=" * 80)

signal = maximal.chromatogram.get_signal("raw")
time_points = maximal.chromatogram.time_points

# Check signal at expected truncation positions
expected_positions = [645, 705, 735, 795]  # NULL and level-1 product peaks
print("\nSignal values at expected truncation positions:")
for expected_pos in expected_positions:
    # Find closest time point
    idx = np.argmin(np.abs(time_points - expected_pos))
    actual_pos = time_points[idx]
    value = signal[idx]
    
    # Check if it's a local maximum
    is_local_max = False
    if 0 < idx < len(signal) - 1:
        is_local_max = signal[idx] > signal[idx-1] and signal[idx] >= signal[idx+1]
    
    print(f"  Position {expected_pos}s -> actual {actual_pos:.0f}s: value={value:.1f}, local_max={is_local_max}")

# Compute background and Z-scores
background = np.percentile(signal, 10)
print(f"\nBackground (10th percentile): {background:.1f}")

print("\nZ-scores at expected positions:")
for expected_pos in expected_positions:
    idx = np.argmin(np.abs(time_points - expected_pos))
    value = signal[idx]
    z_score = (value - background) / np.sqrt(background + 1.0)
    print(f"  Position {time_points[idx]:.0f}s: value={value:.1f}, Z={z_score:.2f} {'✓ PASS' if z_score >= 3.0 else '✗ FAIL'}")

# Show signal statistics
print(f"\nSignal statistics:")
print(f"  Min: {np.min(signal):.1f}")
print(f"  Max: {np.max(signal):.1f}")
print(f"  Mean: {np.mean(signal):.1f}")
print(f"  Median: {np.median(signal):.1f}")
print(f"  10th percentile: {np.percentile(signal, 10):.1f}")

# Show all local maxima
print(f"\nAll local maxima in signal:")
local_maxima = []
for i in range(1, len(signal) - 1):
    if signal[i] > signal[i-1] and signal[i] >= signal[i+1]:
        local_maxima.append((time_points[i], signal[i]))

print(f"  Found {len(local_maxima)} local maxima")
print(f"\n  First 30 local maxima:")
for i, (pos, val) in enumerate(local_maxima[:30]):
    z_score = (val - background) / np.sqrt(background + 1.0)
    status = "✓" if z_score >= 3.0 else "✗"
    print(f"    {i+1}. {pos:.0f}s: height={val:.1f}, Z={z_score:.2f} {status}")
