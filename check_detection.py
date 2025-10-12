#!/usr/bin/env python3
"""Check what peaks are detected by different methods."""

import sys
import numpy as np
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent / "src"))

from lcseq.domain.services import CompoundSearchService, PeakDetector
from lcseq.infrastructure import HDF5CompoundLoader

# Load data
loader = HDF5CompoundLoader()
compounds = loader.load_all("test_data/processed_data.h5")

# Find maximal compound
search_service = CompoundSearchService()
maximal = search_service.find_by_sequence(compounds, "Leu-LA03-Pro-Leu-DLeuMe-DPro-Leu-Leu-DPro")

signal = maximal.chromatogram.get_signal("raw")
time_points = maximal.chromatogram.time_points

# Create detector
detector = PeakDetector()

# Get background
background = detector._estimate_background(signal)

# Find local maxima
local_maxima = detector._find_local_maxima(signal)
print(f"Local maxima: {len(local_maxima)}")
print("First 15 local maxima:")
for i, idx in enumerate(local_maxima[:15]):
    pos = time_points[idx]
    height = signal[idx]
    z_score = (height - background) / np.sqrt(background + 1.0)
    print(f"  {i+1}. {pos:.0f}s: height={height:.1f}, Z={z_score:.2f}")

# Find shoulder peaks
shoulder_peaks = detector._find_shoulder_peaks(signal, background, 3.0, local_maxima)
print(f"\nShoulder peaks: {len(shoulder_peaks)}")
if shoulder_peaks:
    print("Shoulder peaks:")
    for i, idx in enumerate(shoulder_peaks[:10]):
        pos = time_points[idx]
        height = signal[idx]
        z_score = (height - background) / np.sqrt(background + 1.0)
        print(f"  {i+1}. {pos:.0f}s: height={height:.1f}, Z={z_score:.2f}")

# Check specific positions
print("\nChecking specific positions (expected truncation product retention times):")
for expected_pos in [645, 705, 735, 795]:
    idx = np.argmin(np.abs(time_points - expected_pos))
    actual_pos = time_points[idx]
    height = signal[idx]
    z_score = (height - background) / np.sqrt(background + 1.0)
    
    is_local_max = idx in local_maxima
    is_shoulder = idx in shoulder_peaks
    
    print(f"  {expected_pos}s -> {actual_pos:.0f}s: height={height:.1f}, Z={z_score:.2f}, local_max={is_local_max}, shoulder={is_shoulder}")
