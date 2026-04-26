"""
Use case for processing entire LC-Seq library with global hierarchy.

This use case processes all compounds in a single run, ensuring each compound
is classified exactly once via bottom-up processing.

Key Optimization: Uses equivalence classes to reduce hierarchy building from
O(n²) to O(k²) where k = number of unique block_support_sequence values.
For typical libraries: 64k compounds → ~1.5k classes = ~99.95% reduction.

References
----------
THEORY.md Section 4.2: Hierarchy Construction
THEORY.md Section 4.2.3: Hybrid Pooled Strategy
"""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Set, Any, Tuple
from datetime import datetime
from concurrent.futures import ThreadPoolExecutor, as_completed
import json
import time
import threading
from tqdm import tqdm

from lcseq.domain.entities import Compound, Peak
from lcseq.domain.models import (
    CompoundHierarchy,
    HierarchyMode,
    EquivalenceClass,
    PoolingStatus,
)
from lcseq.domain.services import (
    EquivalenceClassBuilder,
    HierarchyBuilder,
    PeakClassifier,
    PeakDetector,
    SignalAggregator,
)
from lcseq.application.use_cases import ProcessPooledChromatogramsUseCase
from lcseq.infrastructure.loaders import HDF5CompoundLoader


@dataclass
class LibraryAnalysisResult:
    """Result of full library analysis."""

    total_compounds: int
    equivalence_classes_count: int
    hierarchy_edges: int
    total_peaks: int
    processing_time_seconds: float
    output_dir: Path

    # Statistics by pooling status
    high_correlation_classes: int = 0
    low_correlation_classes: int = 0
    single_variant_classes: int = 0

    # Progress tracking
    compounds_processed: int = 0
    diagnostics_generated: int = 0

    # Error tracking
    errors: List[Dict[str, Any]] = field(default_factory=list)


@dataclass
class CheckpointData:
    """Checkpoint for resumable processing."""

    last_processed_index: int
    total_compounds: int
    equivalence_classes_processed: Set[str]
    statistics: Dict[str, int]
    errors: List[Dict[str, Any]]

    def to_dict(self) -> Dict[str, Any]:
        """Convert to JSON-serializable dictionary."""
        return {
            "last_processed_index": self.last_processed_index,
            "total_compounds": self.total_compounds,
            "equivalence_classes_processed": list(self.equivalence_classes_processed),
            "statistics": self.statistics,
            "errors": self.errors,
        }

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "CheckpointData":
        """Create from dictionary."""
        return cls(
            last_processed_index=data["last_processed_index"],
            total_compounds=data["total_compounds"],
            equivalence_classes_processed=set(data["equivalence_classes_processed"]),
            statistics=data["statistics"],
            errors=data["errors"],
        )


class ProgressReporter:
    """Report progress for long-running library analysis."""

    def __init__(self, total_classes: int, total_compounds: int):
        self.total_classes = total_classes
        self.total_compounds = total_compounds
        self.processed_classes = 0
        self.processed_compounds = 0
        self.start_time = time.time()
        self.last_report_time = self.start_time

    def update(self, classes_done: int, compounds_done: int) -> None:
        """Update progress and report if enough time has passed."""
        self.processed_classes = classes_done
        self.processed_compounds = compounds_done

        current_time = time.time()
        # Report every 5 seconds
        if current_time - self.last_report_time >= 5.0:
            self.report()
            self.last_report_time = current_time

    def report(self) -> None:
        """Print progress report."""
        elapsed = time.time() - self.start_time
        if elapsed > 0 and self.processed_compounds > 0:
            rate = self.processed_compounds / elapsed
            remaining = self.total_compounds - self.processed_compounds
            eta = remaining / rate if rate > 0 else 0

            print(
                f"\rProgress: {self.processed_classes}/{self.total_classes} classes "
                f"({self.processed_compounds:,}/{self.total_compounds:,} compounds) "
                f"[{rate:.1f} cpd/s, ETA: {eta/60:.1f}m]",
                end="", flush=True
            )

    def finish(self) -> None:
        """Print final report."""
        elapsed = time.time() - self.start_time
        print(
            f"\n✓ Completed: {self.processed_compounds:,} compounds "
            f"in {elapsed/60:.1f} minutes ({self.processed_compounds/elapsed:.1f} cpd/s)"
        )


class ProcessLibraryUseCase:
    """
    Process entire library with global hierarchy.

    Each compound is classified exactly once via bottom-up processing.
    Uses equivalence classes for efficient hierarchy building and pooled peak detection.

    Key Benefits
    ------------
    1. No redundant processing: Each compound classified exactly once
    2. Efficient hierarchy: O(k²) instead of O(n²) for edge detection
    3. Pooled peak detection: One detection per equivalence class
    4. Incremental output: Results written as processing completes
    5. Resumable: Checkpoint support for interrupted runs

    Notes
    -----
    Uses existing infrastructure:
    - EquivalenceClassBuilder: Group by block_support_sequence
    - ProcessPooledChromatogramsUseCase: Pooled detection pattern
    - PeakClassifier.classify_hierarchy(): Bottom-up classification
    - CompoundDiagnosticPlotter: Plot generation

    References
    ----------
    THEORY.md Section 4.2: Hierarchy Construction
    THEORY.md Section 4.2.3: Hybrid Pooled Strategy
    """

    def __init__(
        self,
        loader: Optional[HDF5CompoundLoader] = None,
        class_builder: Optional[EquivalenceClassBuilder] = None,
        hierarchy_builder: Optional[HierarchyBuilder] = None,
        pooled_processor: Optional[ProcessPooledChromatogramsUseCase] = None,
    ):
        """
        Initialize with domain services (dependency injection).

        Parameters
        ----------
        loader : HDF5CompoundLoader, optional
            Compound loader (created if not provided)
        class_builder : EquivalenceClassBuilder, optional
            Equivalence class builder (created if not provided)
        hierarchy_builder : HierarchyBuilder, optional
            Hierarchy builder (created if not provided)
        pooled_processor : ProcessPooledChromatogramsUseCase, optional
            Pooled processing use case (created if not provided)
        """
        self.loader = loader or HDF5CompoundLoader()
        self.class_builder = class_builder or EquivalenceClassBuilder()
        self.hierarchy_builder = hierarchy_builder or HierarchyBuilder()
        self.pooled_processor = pooled_processor or ProcessPooledChromatogramsUseCase()

    def execute(
        self,
        hdf5_path: Path,
        output_dir: Path,
        # Peak detection parameters (all required from config)
        alpha: float,
        alpha_product: float,
        prominence_percentile: float,
        min_snr: float,
        min_baseline_sds: float,
        signal_variant: str,
        min_dispersion_r: float,
        sigma_clip_sigma: float,
        # Peak classification parameters (all required from config)
        truncation_margin: float,
        peak_matching_tolerance: float,
        hungarian_min_threshold: float,
        # Pooling parameters (all required from config)
        correlation_threshold: float,
        aggregation_method: str,
        # Validation parameters (all required from config)
        include_rejected: bool,
        clpe_outlier_threshold: float,
        clpe_min_group_size: int,
        # Preprocessing parameters (all required from config)
        preprocessing_params: Dict[str, Any],
        # Optional parameters
        generate_diagnostics: bool = True,
        resume: bool = False,
        hierarchy_mode: HierarchyMode = HierarchyMode.BUILDING_BLOCK,
        n_workers: int = 4,
        # Optional cLPE reference (user-provided file path)
        clpe_reference_csv: Optional[Path] = None,
        # Optional dead time (None = derive from L0 peak RT)
        clpe_t0: Optional[float] = None,
        clpe_reselect_peaks: bool = True,
    ) -> LibraryAnalysisResult:
        """
        Process all compounds in library.

        Phases:
        1. Load & Structure: Load all, build equivalence classes, build global hierarchy
        2. Peak Detection: Pooled mode per equivalence class
        3. Classification: Bottom-up (L0 → Lmax), each compound once
        4. Output: Incremental diagnostic plots + JSONL results

        Parameters
        ----------
        hdf5_path : Path
            Path to HDF5 file with compound data
        output_dir : Path
            Output directory for results and diagnostics
        generate_diagnostics : bool, optional
            Whether to generate diagnostic plots (default True)
        resume : bool, optional
            Resume from checkpoint if available (default False)
        hierarchy_mode : HierarchyMode, optional
            Building block or monomer mode (default BUILDING_BLOCK)
        n_workers : int, optional
            Number of worker threads for parallel diagnostic generation (default 4)
        alpha : float, optional
            Significance level (false positive rate) for peak detection
        prominence_percentile : float, optional
            Prominence percentile threshold
        min_snr : float, optional
            Adaptive SNR threshold multiplier
        min_baseline_sds : float, optional
            Global baseline threshold in SDs
        signal_variant : str, optional
            Signal variant to use for detection
        truncation_margin : float, optional
            Margin beyond truncation positions (in seconds)
        correlation_threshold : float, optional
            Minimum correlation for pooling validity
        aggregation_method : str, optional
            Aggregation method ("mean" or "median")

        Returns
        -------
        LibraryAnalysisResult
            Summary of analysis results and statistics
        """
        start_time = time.time()
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        checkpoint_path = output_dir / ".checkpoint.json"
        checkpoint = None

        if resume and checkpoint_path.exists():
            with open(checkpoint_path, 'r') as f:
                checkpoint = CheckpointData.from_dict(json.load(f))
            print(f"✓ Resuming from checkpoint: {checkpoint.last_processed_index} compounds")

        # ============================================================
        # PHASE 1: Load & Structure
        # ============================================================
        print("\n" + "=" * 80)
        print("PHASE 1: Load & Structure")
        print("=" * 80)

        # Load all compounds
        compounds = self.loader.load_all(hdf5_path)
        total_compounds = len(compounds)
        print(f"✓ Loaded {total_compounds:,} compounds")

        # Build equivalence classes
        print("\nBuilding equivalence classes...")
        equivalence_classes = self.class_builder.build(compounds)
        n_classes = len(equivalence_classes)
        print(f"✓ Built {n_classes:,} equivalence classes")
        print(f"  Reduction: {total_compounds:,} → {n_classes:,} ({100*n_classes/total_compounds:.2f}%)")

        # Build global hierarchy
        print(f"\nBuilding global hierarchy (mode: {hierarchy_mode.value})...")
        hierarchy = self._build_global_hierarchy(
            compounds, equivalence_classes, hierarchy_mode
        )
        print(f"✓ Built hierarchy: {hierarchy.size():,} compounds, {hierarchy.edge_count():,} edges")

        # ============================================================
        # PHASE 2: Peak Detection (Pooled Mode)
        # ============================================================
        print("\n" + "=" * 80)
        print("PHASE 2: Peak Detection (Pooled Mode)")
        print("=" * 80)
        print("Using hybrid pooled strategy (THEORY.md Section 4.2.3):")
        print("  Phase 1: Peak detection on pooled signal (expensive, once per class)")
        print("  Phase 2: Area integration on variants (cheap, per variant)")

        # Pass cLPE params if provided - validation is now integrated into level-by-level processing
        eq_classes_dict, pooled_compounds, quotient_hierarchy, clpe_stats = self.pooled_processor.execute(
            compounds=compounds,
            hierarchy=hierarchy,
            # Peak detection parameters
            alpha=alpha,
            prominence_percentile=prominence_percentile,
            min_snr=min_snr,
            min_baseline_sds=min_baseline_sds,
            signal_variant=signal_variant,
            min_dispersion_r=min_dispersion_r,
            sigma_clip_sigma=sigma_clip_sigma,
            # Peak classification parameters
            alpha_product=alpha_product,
            truncation_margin=truncation_margin,
            peak_matching_tolerance=peak_matching_tolerance,
            hungarian_min_threshold=hungarian_min_threshold,
            # Pooling parameters
            correlation_threshold=correlation_threshold,
            aggregation_method=aggregation_method,
            # Validation parameters
            include_rejected=include_rejected,
            clpe_outlier_threshold=clpe_outlier_threshold,
            clpe_min_group_size=clpe_min_group_size,
            # Preprocessing parameters
            preprocessing_params=preprocessing_params,
            # cLPE validation (integrated into level-by-level processing)
            clpe_reference_csv=clpe_reference_csv,
            clpe_t0=clpe_t0,
        )

        total_peaks = sum(len(cpd.detected_peaks) for cpd in compounds)
        print(f"✓ Detected peaks: {total_peaks:,} total")

        # Count pooling statistics
        high_corr = sum(
            1 for eq in eq_classes_dict.values()
            if eq.pooling_status == PoolingStatus.POOLING_VALID
        )
        low_corr = sum(
            1 for eq in eq_classes_dict.values()
            if eq.pooling_status == PoolingStatus.POOLING_INVALID
        )
        single = sum(
            1 for eq in eq_classes_dict.values()
            if eq.pooling_status == PoolingStatus.NOT_ATTEMPTED
        )

        print(f"  Pooling statistics:")
        print(f"    High correlation (≥{correlation_threshold}): {high_corr}/{n_classes}")
        print(f"    Low correlation (<{correlation_threshold}): {low_corr}/{n_classes}")
        print(f"    Single variant: {single}/{n_classes}")

        # Print cLPE validation results if enabled
        if clpe_stats:
            print(f"\n  cLPE Validation (integrated, level-by-level):")
            print(f"    t0 from L0: {clpe_stats.get('t0', 'N/A'):.2f} min" if clpe_stats.get('t0') else "    t0: not determined")
            total_outliers = 0
            total_reselected = 0
            total_validated = 0
            for level, level_stats in clpe_stats.get("levels", {}).items():
                total_validated += level_stats.get("validated", 0)
                total_outliers += level_stats.get("outliers", 0)
                total_reselected += level_stats.get("reselected", 0)
            if total_validated > 0:
                outlier_rate = total_outliers / total_validated * 100
                print(f"    Validated: {total_validated:,}")
                print(f"    Outliers: {total_outliers} ({outlier_rate:.1f}%)")
                if clpe_reselect_peaks:
                    print(f"    Peaks re-selected: {total_reselected}")

        # ============================================================
        # PHASE 3: Output Generation (Incremental)
        # ============================================================
        print("\n" + "=" * 80)
        print("PHASE 3: Output Generation")
        print("=" * 80)

        # Create output directories
        results_file = output_dir / "library_analysis.jsonl"
        diagnostics_dir = output_dir / "diagnostics" if generate_diagnostics else None

        if diagnostics_dir:
            diagnostics_dir.mkdir(parents=True, exist_ok=True)

        # Process compounds and generate output
        processed_count = 0
        errors = []

        # Get compounds already processed if resuming
        processed_sequences = set()
        if checkpoint:
            processed_sequences = checkpoint.equivalence_classes_processed

        # Write JSONL results (sequential - fast)
        # OPTIMIZATION: Buffer writes for reduced I/O syscalls
        print("Writing JSONL results...")
        BUFFER_SIZE = 1000  # Write every 1000 records
        with open(results_file, 'a' if resume else 'w') as jsonl_file:
            write_buffer = []

            for eq_class in equivalence_classes:
                # Skip if already processed
                if eq_class.block_support_sequence in processed_sequences:
                    processed_count += len(eq_class.members)
                    continue

                try:
                    # Buffer compound results
                    for compound in eq_class.members:
                        record = self._compound_to_record(
                            compound,
                            hierarchy,
                            hierarchy_mode,
                        )
                        write_buffer.append(json.dumps(record) + "\n")
                        processed_count += 1

                        # Flush buffer when full
                        if len(write_buffer) >= BUFFER_SIZE:
                            jsonl_file.writelines(write_buffer)
                            write_buffer.clear()

                    processed_sequences.add(eq_class.block_support_sequence)

                except Exception as e:
                    errors.append({
                        "equivalence_class": eq_class.block_support_sequence,
                        "error": str(e),
                        "phase": "jsonl_writing"
                    })

            # Flush remaining buffer
            if write_buffer:
                jsonl_file.writelines(write_buffer)

        print(f"✓ Wrote {processed_count:,} compound records")

        # Generate diagnostic plots (parallel - slow, I/O bound)
        diagnostics_count = 0
        if diagnostics_dir:
            print(f"\nGenerating diagnostic plots ({n_workers} workers)...")
            diagnostics_count, diag_errors = self._generate_diagnostics_parallel(
                compounds=compounds,
                hierarchy=hierarchy,
                output_dir=diagnostics_dir,
                n_workers=n_workers,
            )
            errors.extend(diag_errors)

        # Clean up checkpoint on successful completion
        if checkpoint_path.exists():
            checkpoint_path.unlink()

        # Build result
        elapsed = time.time() - start_time
        result = LibraryAnalysisResult(
            total_compounds=total_compounds,
            equivalence_classes_count=n_classes,
            hierarchy_edges=hierarchy.edge_count(),
            total_peaks=total_peaks,
            processing_time_seconds=elapsed,
            output_dir=output_dir,
            high_correlation_classes=high_corr,
            low_correlation_classes=low_corr,
            single_variant_classes=single,
            compounds_processed=processed_count,
            diagnostics_generated=diagnostics_count,
            errors=errors,
        )

        # Write summary
        self._write_summary(result, output_dir / "summary.json")

        print(f"\n✓ Results written to: {results_file}")
        if diagnostics_dir:
            print(f"✓ Diagnostics: {diagnostics_count} plots in {diagnostics_dir}")
        if errors:
            print(f"⚠ Errors: {len(errors)} (see summary.json)")

        return result

    def _build_global_hierarchy(
        self,
        compounds: List[Compound],
        equivalence_classes: List[EquivalenceClass],
        mode: HierarchyMode,
    ) -> CompoundHierarchy:
        """
        Build global hierarchy efficiently using equivalence classes.

        Instead of O(n²) edge detection on all compounds, we:
        1. Build edges between representatives only (O(k²) where k << n)
        2. Expand to include all variants with same edges

        Parameters
        ----------
        compounds : List[Compound]
            All compounds
        equivalence_classes : List[EquivalenceClass]
            Equivalence classes grouping compounds by block_support_sequence
        mode : HierarchyMode
            Building block or monomer mode

        Returns
        -------
        CompoundHierarchy
            Complete hierarchy with all compounds
        """
        # Extract one representative per equivalence class
        representatives = []
        rep_to_class = {}

        for eq_class in equivalence_classes:
            rep = next(iter(eq_class.members))
            representatives.append(rep)
            rep_to_class[rep] = eq_class

        print(f"  Building edges on {len(representatives):,} representatives...")

        # Build hierarchy on representatives only
        rep_hierarchy = self.hierarchy_builder.build(representatives, mode)

        print(f"  Expanding to all {len(compounds):,} compounds...")

        # Create full hierarchy
        hierarchy = CompoundHierarchy(mode=mode)

        # Add all compounds
        for compound in compounds:
            hierarchy.add_compound(compound)

        # Copy edges from representative hierarchy to all variants
        # For each edge rep_a -> rep_b, add variant_a -> variant_b for all variants
        block_support_to_class = {
            eq.block_support_sequence: eq for eq in equivalence_classes
        }

        for rep in representatives:
            rep_class = rep_to_class[rep]

            # Get direct descendants of representative
            direct_descs = rep_hierarchy.get_direct_descendants(rep)

            for desc_rep in direct_descs:
                desc_class = rep_to_class[desc_rep]

                # Add edges from all variants of ancestor to all variants of descendant
                for ancestor_variant in rep_class.members:
                    for desc_variant in desc_class.members:
                        hierarchy.add_edge(ancestor_variant, desc_variant)

        return hierarchy

    def _compound_to_record(
        self,
        compound: Compound,
        hierarchy: CompoundHierarchy,
        hierarchy_mode: HierarchyMode,
    ) -> Dict[str, Any]:
        """Convert compound to JSONL record.

        Notes
        -----
        Retention times are converted from seconds to minutes for output.
        Building blocks are stored in cycle order (0 = C-terminus, synthesized first).
        """
        level_attr = "monomer_level" if hierarchy_mode == HierarchyMode.MONOMER else "level"

        peaks_data = []
        for peak in compound.detected_peaks:
            peaks_data.append({
                "retention_time_min": round(peak.position / 60.0, 4),
                "area": round(peak.area, 2),
                "height": round(peak.height, 2),
                "classification": peak.peak_type.value,
                "matched_compound": peak.matched_compound_sequence,
                "matched_position_min": (
                    round(peak.matched_peak_position / 60.0, 4)
                    if peak.matched_peak_position else None
                ),
            })

        # Export building blocks with cycle positions
        # Cycle 0 = C-terminus (synthesized first), higher cycles toward N-terminus
        building_blocks_data = []
        for bb in compound.building_blocks:
            building_blocks_data.append({
                "cycle": bb.cycle,
                "code": bb.code,
                "is_null": bb.is_null,
            })

        return {
            "sequence": compound.positional_block_sequence,
            "block_support_sequence": compound.block_support_sequence,
            "level": getattr(compound, level_attr),
            "building_blocks": building_blocks_data,
            "n_peaks": len(compound.detected_peaks),
            "peaks": peaks_data,
            "selected_peak_rt_min": (
                round(compound.selected_peak.position / 60.0, 4)
                if compound.selected_peak else None
            ),
        }

    def _generate_diagnostic(
        self,
        compound: Compound,
        hierarchy: CompoundHierarchy,
        output_dir: Path,
    ) -> None:
        """Generate diagnostic plot for a compound."""
        # Import here to avoid circular dependency
        from lcseq.presentation.visualization.plotters import CompoundDiagnosticPlotter

        plotter = CompoundDiagnosticPlotter()

        # Create filename from block support sequence (safe for filesystem)
        safe_name = compound.block_support_sequence.replace("/", "_").replace("\\", "_")
        output_path = output_dir / f"{safe_name}.png"

        plotter.plot(compound, hierarchy, output_path=output_path)

    def _generate_diagnostics_parallel(
        self,
        compounds: List[Compound],
        hierarchy: CompoundHierarchy,
        output_dir: Path,
        n_workers: int = 4,
    ) -> Tuple[int, List[Dict[str, Any]]]:
        """
        Generate diagnostic plots in parallel using thread pool.

        Uses ThreadPoolExecutor because:
        - Matplotlib releases GIL during rendering
        - No pickling overhead (hierarchy shared across threads)
        - I/O bound (writing PNGs to disk)

        Parameters
        ----------
        compounds : List[Compound]
            All compounds to generate plots for
        hierarchy : CompoundHierarchy
            Compound hierarchy (shared, read-only)
        output_dir : Path
            Output directory for plots
        n_workers : int
            Number of worker threads (default 4)

        Returns
        -------
        Tuple[int, List[Dict[str, Any]]]
            (success_count, list of error dicts)
        """
        # Import here to avoid circular dependency
        from lcseq.presentation.visualization.plotters import CompoundDiagnosticPlotter

        errors = []
        success_count = 0
        total = len(compounds)
        lock = threading.Lock()
        start_time = time.time()

        def generate_single(compound: Compound) -> Tuple[bool, Optional[Dict[str, Any]]]:
            """Generate a single diagnostic plot."""
            try:
                # Each thread gets its own plotter instance (matplotlib not thread-safe)
                plotter = CompoundDiagnosticPlotter()
                safe_name = compound.block_support_sequence.replace("/", "_").replace("\\", "_")
                output_path = output_dir / f"{safe_name}.png"
                plotter.plot(compound, hierarchy, output_path=output_path)
                return True, None
            except Exception as e:
                return False, {
                    "compound": compound.block_support_sequence,
                    "error": str(e),
                    "phase": "diagnostic"
                }

        with ThreadPoolExecutor(max_workers=n_workers) as executor:
            # Submit all tasks
            futures = {
                executor.submit(generate_single, cpd): cpd
                for cpd in compounds
            }

            # Process results as they complete
            for future in as_completed(futures):
                success, error = future.result()

                with lock:
                    if success:
                        success_count += 1
                    else:
                        errors.append(error)

                    # Progress reporting every 1000 plots
                    processed = success_count + len(errors)
                    if processed % 1000 == 0 or processed == total:
                        elapsed = time.time() - start_time
                        rate = processed / elapsed if elapsed > 0 else 0
                        remaining = total - processed
                        eta = remaining / rate if rate > 0 else 0
                        print(
                            f"\r  Plots: {processed:,}/{total:,} "
                            f"[{rate:.1f}/s, ETA: {eta/60:.1f}m]",
                            end="", flush=True
                        )

        print()  # Newline after progress
        print(f"✓ Generated {success_count:,} diagnostic plots")
        if errors:
            print(f"  ⚠ {len(errors)} plots failed")

        return success_count, errors

    def _save_checkpoint(
        self,
        checkpoint_path: Path,
        processed_count: int,
        total_compounds: int,
        processed_sequences: Set[str],
        diagnostics_count: int,
        errors: List[Dict[str, Any]],
    ) -> None:
        """Save checkpoint for resumable processing."""
        checkpoint = CheckpointData(
            last_processed_index=processed_count,
            total_compounds=total_compounds,
            equivalence_classes_processed=processed_sequences,
            statistics={
                "processed": processed_count,
                "diagnostics": diagnostics_count,
                "errors": len(errors),
            },
            errors=errors,
        )

        with open(checkpoint_path, 'w') as f:
            json.dump(checkpoint.to_dict(), f, indent=2)

    def _write_summary(
        self,
        result: LibraryAnalysisResult,
        output_path: Path,
    ) -> None:
        """Write analysis summary to JSON."""
        summary = {
            "timestamp": datetime.now().isoformat(),
            "total_compounds": result.total_compounds,
            "equivalence_classes": result.equivalence_classes_count,
            "hierarchy_edges": result.hierarchy_edges,
            "total_peaks": result.total_peaks,
            "processing_time_seconds": round(result.processing_time_seconds, 2),
            "processing_time_minutes": round(result.processing_time_seconds / 60, 2),
            "compounds_per_second": round(
                result.total_compounds / result.processing_time_seconds, 1
            ),
            "pooling_statistics": {
                "high_correlation": result.high_correlation_classes,
                "low_correlation": result.low_correlation_classes,
                "single_variant": result.single_variant_classes,
            },
            "output": {
                "compounds_processed": result.compounds_processed,
                "diagnostics_generated": result.diagnostics_generated,
                "errors": len(result.errors),
            },
            "errors": result.errors[:100] if result.errors else [],  # Limit error list
        }

        with open(output_path, 'w') as f:
            json.dump(summary, f, indent=2)

    # NOTE: _run_clpe_validation and _find_null_peak_rt methods have been removed.
    # cLPE validation is now integrated into level-by-level processing in PeakClassifier.
    # t0 is automatically extracted from the L0 NULL peak during classification.
