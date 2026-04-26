"""
Use case for filtering LC-Seq library based on signal quality metrics.

This use case provides a preprocessing step that filters compounds based on
signal quality (intensity, noise) and replicate correlation within
equivalence classes.

References
----------
THEORY.md Section 4.2.8: Validity Requirements
"""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Any
from datetime import datetime
import json
import time

from lcseq.domain.entities import Compound
from lcseq.domain.models import EquivalenceClass
from lcseq.domain.services import (
    EquivalenceClassBuilder,
    QualityAssessor,
    SignalQualityMetrics,
    EquivalenceClassQuality,
)
from lcseq.infrastructure.loaders import HDF5CompoundLoader


@dataclass
class FilterResult:
    """Result of library filtering."""

    # Input statistics
    total_compounds: int
    total_equivalence_classes: int

    # Output statistics
    passed_compounds: int
    passed_equivalence_classes: int
    filtered_compounds: int
    filtered_equivalence_classes: int

    # Filter breakdown
    failed_correlation: int = 0
    failed_intensity: int = 0
    failed_noise: int = 0

    # Single-variant classes (always included, not filtered by correlation)
    single_variant_classes: int = 0
    single_variant_compounds: int = 0

    # Processing time
    processing_time_seconds: float = 0.0

    # Output paths
    output_hdf5: Optional[Path] = None
    inclusion_csv: Optional[Path] = None
    qc_report: Optional[Path] = None

    # Detailed QC data (for JSON report)
    class_quality: List[Dict[str, Any]] = field(default_factory=list)


class FilterLibraryUseCase:
    """
    Filter LC-Seq library based on signal quality metrics.

    This use case provides principled pre-filtering of low-quality data
    before downstream analysis. It computes quality metrics at both
    individual compound and equivalence class levels.

    Quality Filters
    ---------------
    1. Replicate Correlation: Min pairwise Pearson correlation within
       equivalence class (positional variants are pseudo-replicates)
    2. Signal Intensity: Total signal or percentile-based threshold
    3. Noise Level: Baseline noise / median signal ratio

    Output Options
    --------------
    - Filtered HDF5: New file with only passing compounds
    - Inclusion CSV: List of compound sequences to include
    - QC Report: Detailed JSON with filtering statistics

    Notes
    -----
    - Single-variant classes are always included (can't assess correlation)
    - Filtering happens at equivalence class level (all variants pass or fail)
    - Uses sigma-clipping for robust baseline/noise estimation

    References
    ----------
    THEORY.md Section 4.2.8: Validity Requirements
    """

    def __init__(
        self,
        loader: Optional[HDF5CompoundLoader] = None,
        class_builder: Optional[EquivalenceClassBuilder] = None,
        quality_assessor: Optional[QualityAssessor] = None,
    ):
        """
        Initialize with domain services (dependency injection).

        Parameters
        ----------
        loader : HDF5CompoundLoader, optional
            Compound loader (created if not provided)
        class_builder : EquivalenceClassBuilder, optional
            Equivalence class builder (created if not provided)
        quality_assessor : QualityAssessor, optional
            Quality assessor (created if not provided)
        """
        self.loader = loader or HDF5CompoundLoader()
        self.class_builder = class_builder or EquivalenceClassBuilder()
        self.quality_assessor = quality_assessor or QualityAssessor()

    def execute(
        self,
        hdf5_path: Path,
        output_dir: Path,
        # Replicate correlation filter
        min_correlation: float = 0.8,
        # Intensity filter
        intensity_percentile: Optional[float] = 0.05,  # Exclude bottom 5%
        intensity_absolute: Optional[float] = None,
        # Noise filter
        max_noise_ratio: Optional[float] = None,
        # Output options
        generate_hdf5: bool = True,
        generate_csv: bool = True,
        generate_report: bool = True,
    ) -> FilterResult:
        """
        Filter library and generate output files.

        Parameters
        ----------
        hdf5_path : Path
            Input HDF5 file with compound data
        output_dir : Path
            Output directory for filtered files
        min_correlation : float, optional
            Minimum replicate correlation threshold (default 0.8)
        intensity_percentile : float, optional
            Exclude compounds below this percentile (default 0.05 = bottom 5%)
        intensity_absolute : float, optional
            Absolute minimum total signal (overrides percentile)
        max_noise_ratio : float, optional
            Maximum noise_std / median_signal ratio (default None = disabled)
        generate_hdf5 : bool, optional
            Generate filtered HDF5 file (default True)
        generate_csv : bool, optional
            Generate inclusion CSV (default True)
        generate_report : bool, optional
            Generate QC report JSON (default True)

        Returns
        -------
        FilterResult
            Summary of filtering results and statistics
        """
        start_time = time.time()
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        # ============================================================
        # PHASE 1: Load & Structure
        # ============================================================
        print("\n" + "=" * 80)
        print("PHASE 1: Load & Structure")
        print("=" * 80)

        # Load all compounds
        compounds = self.loader.load_all(hdf5_path)
        total_compounds = len(compounds)

        # Build equivalence classes
        print("\nBuilding equivalence classes...")
        equivalence_classes = self.class_builder.build(compounds)
        n_classes = len(equivalence_classes)
        print(f"Built {n_classes:,} equivalence classes")

        # ============================================================
        # PHASE 2: Compute Quality Metrics
        # ============================================================
        print("\n" + "=" * 80)
        print("PHASE 2: Compute Quality Metrics")
        print("=" * 80)

        # Compute library-wide intensity percentile threshold
        intensity_threshold = None
        if intensity_percentile is not None:
            intensity_threshold = self.quality_assessor.compute_library_intensity_percentile(
                compounds, intensity_percentile
            )
            print(f"Intensity threshold ({intensity_percentile*100:.0f}th percentile): {intensity_threshold:.2f}")
        elif intensity_absolute is not None:
            intensity_threshold = intensity_absolute
            print(f"Intensity threshold (absolute): {intensity_threshold:.2f}")

        # Assess quality of each equivalence class
        class_quality_results: List[EquivalenceClassQuality] = []

        print("\nAssessing equivalence class quality...")
        for eq_class in equivalence_classes:
            quality = self.quality_assessor.assess_equivalence_class(
                variants=list(eq_class.members),
                block_support_sequence=eq_class.block_support_sequence,
                correlation_threshold=min_correlation,
                intensity_percentile_threshold=intensity_threshold,
                max_noise_ratio=max_noise_ratio,
            )
            class_quality_results.append(quality)

        # ============================================================
        # PHASE 3: Apply Filters
        # ============================================================
        print("\n" + "=" * 80)
        print("PHASE 3: Apply Filters")
        print("=" * 80)

        passed_classes: List[EquivalenceClass] = []
        failed_classes: List[EquivalenceClass] = []
        single_variant_classes: List[EquivalenceClass] = []

        failed_correlation = 0
        failed_intensity = 0
        failed_noise = 0

        for eq_class, quality in zip(equivalence_classes, class_quality_results):
            # Single-variant classes always pass (can't assess correlation)
            if quality.n_variants == 1:
                passed_classes.append(eq_class)
                single_variant_classes.append(eq_class)
                continue

            # Multi-variant class - apply all filters
            if quality.passes_all:
                passed_classes.append(eq_class)
            else:
                failed_classes.append(eq_class)

                # Count failure reasons (can have multiple)
                if not quality.passes_correlation:
                    failed_correlation += 1
                if not quality.passes_intensity:
                    failed_intensity += 1
                if not quality.passes_noise:
                    failed_noise += 1

        # Collect passing compounds
        passed_compounds: List[Compound] = []
        for eq_class in passed_classes:
            passed_compounds.extend(eq_class.members)

        # ============================================================
        # PHASE 4: Generate Outputs
        # ============================================================
        print("\n" + "=" * 80)
        print("PHASE 4: Generate Outputs")
        print("=" * 80)

        output_hdf5 = None
        inclusion_csv = None
        qc_report_path = None

        # Generate filtered HDF5
        if generate_hdf5 and passed_compounds:
            output_hdf5 = output_dir / "filtered_library.h5"
            self.loader.save_all(passed_compounds, output_hdf5)

        # Generate inclusion CSV
        if generate_csv:
            inclusion_csv = output_dir / "included_sequences.csv"
            self._write_inclusion_csv(passed_compounds, inclusion_csv)

        # Generate QC report
        class_quality_data = []
        if generate_report:
            qc_report_path = output_dir / "qc_report.json"
            class_quality_data = self._build_class_quality_data(
                equivalence_classes, class_quality_results
            )
            self._write_qc_report(
                total_compounds=total_compounds,
                total_classes=n_classes,
                passed_classes=passed_classes,
                failed_classes=failed_classes,
                single_variant_classes=single_variant_classes,
                class_quality_data=class_quality_data,
                min_correlation=min_correlation,
                intensity_threshold=intensity_threshold,
                max_noise_ratio=max_noise_ratio,
                output_path=qc_report_path,
            )

        # Build result
        elapsed = time.time() - start_time

        result = FilterResult(
            total_compounds=total_compounds,
            total_equivalence_classes=n_classes,
            passed_compounds=len(passed_compounds),
            passed_equivalence_classes=len(passed_classes),
            filtered_compounds=total_compounds - len(passed_compounds),
            filtered_equivalence_classes=len(failed_classes),
            failed_correlation=failed_correlation,
            failed_intensity=failed_intensity,
            failed_noise=failed_noise,
            single_variant_classes=len(single_variant_classes),
            single_variant_compounds=sum(len(eq.members) for eq in single_variant_classes),
            processing_time_seconds=elapsed,
            output_hdf5=output_hdf5,
            inclusion_csv=inclusion_csv,
            qc_report=qc_report_path,
            class_quality=class_quality_data,
        )

        # Print summary
        print(f"\n{'=' * 80}")
        print("FILTERING COMPLETE")
        print("=" * 80)
        print(f"\nInput:")
        print(f"  Total compounds: {total_compounds:,}")
        print(f"  Equivalence classes: {n_classes:,}")
        print(f"\nFilter Thresholds:")
        print(f"  Min correlation: {min_correlation}")
        if intensity_threshold:
            print(f"  Min intensity: {intensity_threshold:.2f}")
        if max_noise_ratio:
            print(f"  Max noise ratio: {max_noise_ratio}")
        print(f"\nResults:")
        print(f"  Passed compounds: {len(passed_compounds):,} ({100*len(passed_compounds)/total_compounds:.1f}%)")
        print(f"  Passed classes: {len(passed_classes):,} ({100*len(passed_classes)/n_classes:.1f}%)")
        print(f"  Single-variant (auto-passed): {len(single_variant_classes):,}")
        print(f"\nFiltered Out:")
        print(f"  Failed correlation: {failed_correlation:,} classes")
        print(f"  Failed intensity: {failed_intensity:,} classes")
        print(f"  Failed noise: {failed_noise:,} classes")
        print(f"\nProcessing time: {elapsed:.1f}s")

        if output_hdf5:
            print(f"\nOutput: {output_hdf5}")
        if inclusion_csv:
            print(f"Inclusion list: {inclusion_csv}")
        if qc_report_path:
            print(f"QC report: {qc_report_path}")

        return result

    def _write_inclusion_csv(
        self,
        compounds: List[Compound],
        output_path: Path,
    ) -> None:
        """Write inclusion CSV with compound sequences."""
        with open(output_path, "w") as f:
            f.write("sequence,block_support_sequence\n")
            for compound in compounds:
                f.write(f"{compound.positional_block_sequence},{compound.block_support_sequence}\n")

        print(f"Wrote inclusion list: {output_path}")

    def _build_class_quality_data(
        self,
        equivalence_classes: List[EquivalenceClass],
        quality_results: List[EquivalenceClassQuality],
    ) -> List[Dict[str, Any]]:
        """Build quality data for JSON export."""
        data = []
        for eq_class, quality in zip(equivalence_classes, quality_results):
            class_data = {
                "block_support_sequence": eq_class.block_support_sequence,
                "n_variants": quality.n_variants,
                "min_correlation": round(quality.min_correlation, 4),
                "mean_total_signal": round(quality.mean_total_signal, 2),
                "mean_noise_ratio": round(quality.mean_noise_ratio, 6),
                "passes_correlation": quality.passes_correlation,
                "passes_intensity": quality.passes_intensity,
                "passes_noise": quality.passes_noise,
                "passes_all": quality.passes_all,
                "variant_sequences": [
                    cpd.positional_block_sequence for cpd in eq_class.members
                ],
            }
            data.append(class_data)
        return data

    def _write_qc_report(
        self,
        total_compounds: int,
        total_classes: int,
        passed_classes: List[EquivalenceClass],
        failed_classes: List[EquivalenceClass],
        single_variant_classes: List[EquivalenceClass],
        class_quality_data: List[Dict[str, Any]],
        min_correlation: float,
        intensity_threshold: Optional[float],
        max_noise_ratio: Optional[float],
        output_path: Path,
    ) -> None:
        """Write detailed QC report JSON."""
        report = {
            "timestamp": datetime.now().isoformat(),
            "summary": {
                "total_compounds": total_compounds,
                "total_equivalence_classes": total_classes,
                "passed_compounds": sum(len(eq.members) for eq in passed_classes),
                "passed_equivalence_classes": len(passed_classes),
                "filtered_equivalence_classes": len(failed_classes),
                "single_variant_classes": len(single_variant_classes),
            },
            "thresholds": {
                "min_correlation": min_correlation,
                "intensity_threshold": intensity_threshold,
                "max_noise_ratio": max_noise_ratio,
            },
            "equivalence_classes": class_quality_data,
        }

        with open(output_path, "w") as f:
            json.dump(report, f, indent=2)

        print(f"Wrote QC report: {output_path}")
