#!/usr/bin/env python3
"""Debug script to check peak classification for maximal compound."""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent / "src"))

from lcseq.domain.services import (
    HierarchyBuilder,
    CompoundSearchService,
    LineageFinderService,
)
from lcseq.domain.models import HierarchyMode
from lcseq.application.use_cases import ProcessChromatogramsUseCase
from lcseq.infrastructure import HDF5CompoundLoader
from lcseq.domain.entities.peak import PeakType

# Load data
print("Loading data...")
loader = HDF5CompoundLoader()
compounds = loader.load_all("test_data/processed_data.h5")
print(f"Loaded {len(compounds)} compounds")

# Find maximal compound
search_service = CompoundSearchService()
maximal = search_service.find_by_sequence(compounds, "Leu-LA03-Pro-Leu-DLeuMe-DPro-Leu-Leu-DPro")
print(f"\nMaximal compound: {maximal.positional_sequence}")
print(f"Level: {maximal.monomer_level}")

# Find lineage
lineage_finder = LineageFinderService()
lineage = lineage_finder.find_principal_ideal(maximal, compounds, HierarchyMode.MONOMER)
print(f"Lineage size: {len(lineage)} compounds")

# Build hierarchy
builder = HierarchyBuilder()
hierarchy = builder.build(lineage, HierarchyMode.MONOMER)
print(f"Hierarchy: {hierarchy.size()} compounds, {hierarchy.edge_count()} edges")

# Process chromatograms
print("\nProcessing chromatograms...")
process_use_case = ProcessChromatogramsUseCase()
peaks_dict = process_use_case.execute(
    compounds=lineage,
    hierarchy=hierarchy,
    z_threshold=3.0,
    prominence_percentile=0.2,
)

# Check maximal compound peaks
maximal_peaks = peaks_dict[maximal]
print(f"\nMaximal compound peaks: {len(maximal_peaks)}")
print("\nPeak classification breakdown:")
by_type = {}
for peak in maximal_peaks:
    peak_type = peak.peak_type.value
    by_type[peak_type] = by_type.get(peak_type, 0) + 1
    
for peak_type, count in sorted(by_type.items()):
    print(f"  {peak_type}: {count}")

# Show all peaks with details
print("\nAll peaks:")
for i, peak in enumerate(maximal_peaks):
    prom_str = f"{peak.prominence:.1f}" if peak.prominence is not None else "None"
    print(f"  {i+1}. Pos={peak.position:.1f}s, H={peak.height:.1f}, Type={peak.peak_type.value}, Prom={prom_str}")

# Check null peak
null_compound = search_service.find_by_sequence(lineage, "AgxNull-AgxNull-AgxNull")
if null_compound:
    null_peaks = peaks_dict[null_compound]
    print(f"\nNull compound (L0) peaks: {len(null_peaks)}")
    if null_peaks:
        # Find the global maximum (NULL peak)
        null_peak = max(null_peaks, key=lambda p: p.height)
        print(f"  NULL peak position: {null_peak.position:.1f}s, height={null_peak.height:.1f}")
        
        # Check if maximal has peak at this position
        print(f"\nChecking maximal compound for NULL peak at {null_peak.position:.1f}s:")
        for peak in maximal_peaks:
            if abs(peak.position - null_peak.position) < 10:  # Within 10s
                print(f"  Found peak at {peak.position:.1f}s (delta={abs(peak.position - null_peak.position):.1f}s), Type={peak.peak_type.value}")

# Check a few level-1 compounds
print("\n\nLevel-1 compounds (should have product peaks that become truncations in maximal):")
level1_compounds = [c for c in lineage if c.monomer_level == 1][:5]
for compound in level1_compounds:
    peaks = peaks_dict[compound]
    product_peaks = [p for p in peaks if p.peak_type == PeakType.PUTATIVE_PRODUCT]
    print(f"\n{compound.positional_sequence}:")
    print(f"  Total peaks: {len(peaks)}")
    print(f"  Product peaks: {len(product_peaks)}")
    if product_peaks:
        for pp in product_peaks:
            print(f"    Product at {pp.position:.1f}s (height={pp.height:.1f})")
            
            # Check if maximal has this as truncation
            found_in_maximal = False
            for mp in maximal_peaks:
                if abs(mp.position - pp.position) < 10:  # Within 10s
                    found_in_maximal = True
                    print(f"      -> Found in maximal at {mp.position:.1f}s as {mp.peak_type.value} (delta={abs(mp.position - pp.position):.1f}s)")
                    break
            if not found_in_maximal:
                print(f"      -> NOT FOUND in maximal!")
