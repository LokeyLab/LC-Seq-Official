#!/usr/bin/env python3
"""Debug script to investigate why consensus signals have fewer peaks."""

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from lcseq.domain.services import (
    HierarchyBuilder,
    CompoundSearchService,
    LineageFinderService,
    SignalAggregator,
    PeakDetector,
    EquivalenceClassBuilder,
)
from lcseq.domain.models import HierarchyMode
from lcseq.domain.entities import Compound
from lcseq.infrastructure import HDF5CompoundLoader
from lcseq.application.use_cases import ProcessChromatogramsUseCase
import numpy as np


def debug_consensus():
    """Debug consensus peak detection."""
    print("=" * 80)
    print("Debugging Consensus Peak Detection")
    print("=" * 80)

    # Load data
    loader = HDF5CompoundLoader()
    compounds = loader.load_all(Path("test_data/processed_data.h5"))
    print(f"\n✓ Loaded {len(compounds):,} compounds")

    # Find reference
    search_service = CompoundSearchService()
    reference = search_service.find_by_sequence(
        compounds, "Leu-LA03-Pro-Leu-DLeuMe-DPro-Leu-Leu-DPro"
    )
    print(f"✓ Found reference: {reference.positional_sequence}")

    # Find lineage
    lineage_finder = LineageFinderService()
    lineage = lineage_finder.find_principal_ideal(
        reference, compounds, HierarchyMode.MONOMER
    )
    print(f"✓ Found lineage: {len(lineage)} compounds")

    # Build equivalence classes
    class_builder = EquivalenceClassBuilder()
    equivalence_classes = class_builder.build(lineage)
    print(f"✓ Built {len(equivalence_classes)} equivalence classes")

    # Find a class with multiple variants to debug
    aggregator = SignalAggregator()
    detector = PeakDetector()

    # Build hierarchy for pipeline processing
    hierarchy_builder = HierarchyBuilder()
    hierarchy = hierarchy_builder.build(lineage, HierarchyMode.MONOMER)

    # Create pipeline use case
    process_use_case = ProcessChromatogramsUseCase()

    for eq_class in equivalence_classes:
        if len(eq_class.compounds) >= 3:
            variants = list(eq_class.compounds)

            print(f"\n" + "=" * 80)
            print(f"Debugging class: {eq_class.residue_sequence}")
            print(f"Variants: {len(variants)}")
            print("=" * 80)

            # Get individual peak counts
            print("\nIndividual variant peaks:")
            for variant in variants:
                peaks = detector.detect_peaks(variant.chromatogram)
                signal = variant.chromatogram.get_signal()
                print(f"  {variant.positional_sequence}: {len(peaks)} peaks")
                print(f"    Signal: min={signal.min():.1f}, max={signal.max():.1f}, mean={signal.mean():.1f}")

            # Aggregate and check consensus
            consensus_chrom, min_corr, is_valid, reason = aggregator.aggregate(
                variants, method="mean"
            )

            consensus_signal = consensus_chrom.get_signal()
            print(f"\nConsensus signal:")
            print(f"  min={consensus_signal.min():.1f}, max={consensus_signal.max():.1f}, mean={consensus_signal.mean():.1f}")
            print(f"  Correlation: {min_corr:.3f}")

            # Test 1: Direct peak detection
            consensus_peaks_direct = detector.detect_peaks(consensus_chrom)
            print(f"\n[Test 1] Direct peak detection: {len(consensus_peaks_direct)} peaks")
            if len(consensus_peaks_direct) > 0:
                for i, peak in enumerate(consensus_peaks_direct):
                    print(f"  Peak {i+1}: position = {peak.position:.1f}s, height = {peak.height:.1f}")

            # Test 2: Full pipeline processing
            template = variants[0]

            print(f"\nTemplate compound info:")
            print(f"  Template ID: {id(template)}")
            print(f"  Template sequence: {template.positional_sequence}")

            consensus_compound = Compound(
                building_blocks=template.building_blocks,
                chromatogram=consensus_chrom,
            )

            print(f"\nConsensus compound info:")
            print(f"  Consensus ID: {id(consensus_compound)}")
            print(f"  Consensus sequence: {consensus_compound.positional_sequence}")
            print(f"  Same object as template: {consensus_compound is template}")
            print(f"  Equal to template: {consensus_compound == template}")

            print(f"\n[Test 2] Running consensus compound through full pipeline...")
            print(f"  Consensus compound level: {consensus_compound.level}")
            print(f"  Consensus compound monomer level: {consensus_compound.monomer_level}")
            print(f"  Consensus compound in hierarchy: {consensus_compound in hierarchy.compounds}")

            # Check what compounds are in the hierarchy
            hierarchy_compounds = list(hierarchy.compounds)
            print(f"  Total compounds in hierarchy: {len(hierarchy_compounds)}")

            peaks_dict = process_use_case.execute(
                compounds=[consensus_compound],
                hierarchy=hierarchy,
            )

            print(f"  Returned peaks_dict keys: {len(peaks_dict)}")

            pipeline_peaks = consensus_compound.detected_peaks if consensus_compound.detected_peaks else []
            print(f"  Pipeline detected peaks: {len(pipeline_peaks)} peaks")

            # Check descendants' truncation positions
            descendants = hierarchy.get_descendants(consensus_compound)
            if descendants:
                print(f"\n  Consensus compound has {len(descendants)} descendants:")
                truncation_positions = []
                for desc in descendants:
                    if desc.selected_peak:
                        truncation_positions.append(desc.selected_peak.position)
                        print(f"    {desc.monomer_sequence}: selected peak at {desc.selected_peak.position:.1f}s")
                    else:
                        print(f"    {desc.monomer_sequence}: no selected peak")

                if truncation_positions:
                    max_trunc_pos = max(truncation_positions)
                    print(f"\n  Max truncation position: {max_trunc_pos:.1f}s")
                    from lcseq.config import DEFAULT_TRUNCATION_MARGIN
                    boundary = max_trunc_pos + DEFAULT_TRUNCATION_MARGIN
                    print(f"  Truncation boundary (max + {DEFAULT_TRUNCATION_MARGIN}s margin): {boundary:.1f}s")

                    if len(consensus_peaks_direct) > 0:
                        for i, peak in enumerate(consensus_peaks_direct):
                            if peak.position <= boundary:
                                print(f"  ⚠ Peak {i+1} at {peak.position:.1f}s is BEFORE boundary {boundary:.1f}s (filtered out)")
                            else:
                                print(f"  ✓ Peak {i+1} at {peak.position:.1f}s is AFTER boundary {boundary:.1f}s (should pass)")
            else:
                print(f"\n  Consensus compound has no descendants (L₀)")

            if len(consensus_peaks_direct) > 0 and len(pipeline_peaks) == 0:
                print("\n⚠ ROOT CAUSE: Pipeline filtering removed all peaks!")
                print("  This is likely due to truncation position constraints from descendants")
            elif len(consensus_peaks_direct) == 0:
                print("\n⚠ WARNING: No peaks detected in consensus!")
                print("  This suggests averaging is reducing peak prominence")

            # Only debug first multi-variant class
            break


if __name__ == "__main__":
    debug_consensus()
