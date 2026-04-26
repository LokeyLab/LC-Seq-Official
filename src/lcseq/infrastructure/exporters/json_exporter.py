"""
JSON exporter for comprehensive lineage analysis results.

Exports all relevant analysis data in a structured JSON format including
metadata, hierarchy structure, compound details, peaks, and statistics.
"""

import json
from pathlib import Path
from typing import Dict, List, Optional, Any
from datetime import datetime
from ...domain.entities.compound import Compound
from ...domain.entities.peak import Peak
from ...domain.models import CompoundHierarchy, HierarchyMode, AnalysisMode, EquivalenceClass


class JSONExporter:
    """
    Export comprehensive lineage analysis results to JSON format.

    Creates structured JSON output with complete analysis metadata,
    hierarchy information, compound details, peak classifications,
    and equivalence class data (for pooled mode).

    Examples
    --------
    >>> exporter = JSONExporter()
    >>> exporter.export(
    ...     reference=reference,
    ...     lineage=lineage,
    ...     hierarchy=hierarchy,
    ...     peaks_dict=peaks_dict,
    ...     output_path=Path('results/analysis.json'),
    ...     variant_mode=AnalysisMode.INDIVIDUAL,
    ...     hierarchy_mode=HierarchyMode.MONOMER
    ... )
    """

    def export(
        self,
        reference: Compound,
        lineage: List[Compound],
        hierarchy: CompoundHierarchy,
        peaks_dict: Dict[Compound, List[Peak]],
        output_path: Path,
        variant_mode: AnalysisMode,
        hierarchy_mode: HierarchyMode,
        equivalence_classes: Optional[Dict[str, EquivalenceClass]] = None,
        data_file: Optional[Path] = None,
        total_library_size: Optional[int] = None,
    ) -> Path:
        """
        Export comprehensive analysis results to JSON.

        Parameters
        ----------
        reference : Compound
            Reference compound (root of lineage)
        lineage : List[Compound]
            All lineage members (reference + descendants)
        hierarchy : CompoundHierarchy
            Compound hierarchy structure
        peaks_dict : Dict[Compound, List[Peak]]
            Detected peaks for each compound
        output_path : Path
            Output JSON file path
        variant_mode : AnalysisMode
            Individual or pooled analysis mode
        hierarchy_mode : HierarchyMode
            Building block or monomer hierarchy mode
        equivalence_classes : Optional[Dict[str, EquivalenceClass]]
            Equivalence classes (pooled mode only)
        data_file : Optional[Path]
            Source data file path
        total_library_size : Optional[int]
            Total compounds in library before lineage filtering

        Returns
        -------
        Path
            Path to created JSON file

        Examples
        --------
        >>> json_path = exporter.export(
        ...     reference=reference,
        ...     lineage=lineage,
        ...     hierarchy=hierarchy,
        ...     peaks_dict=peaks_dict,
        ...     output_path=Path('results/analysis.json'),
        ...     variant_mode=AnalysisMode.POOLED,
        ...     hierarchy_mode=HierarchyMode.BUILDING_BLOCK,
        ...     equivalence_classes=equivalence_classes
        ... )
        """
        # Build comprehensive JSON structure
        data = {
            "metadata": self._build_metadata(
                reference,
                variant_mode,
                hierarchy_mode,
                len(lineage),
                data_file,
                total_library_size
            ),
            "hierarchy": self._build_hierarchy_info(
                hierarchy,
                hierarchy_mode
            ),
            "compounds": self._build_compounds_data(
                lineage,
                peaks_dict,
                hierarchy,
                hierarchy_mode
            ),
            "statistics": self._build_statistics(
                lineage,
                peaks_dict,
                equivalence_classes
            )
        }

        # Add equivalence classes if pooled mode
        if equivalence_classes:
            data["equivalence_classes"] = self._build_equivalence_classes_data(
                equivalence_classes,
                peaks_dict
            )

        # Write JSON file
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with open(output_path, 'w') as f:
            json.dump(data, f, indent=2)

        return output_path

    def _build_metadata(
        self,
        reference: Compound,
        variant_mode: AnalysisMode,
        hierarchy_mode: HierarchyMode,
        lineage_size: int,
        data_file: Optional[Path],
        total_library_size: Optional[int]
    ) -> Dict[str, Any]:
        """Build metadata section."""
        metadata = {
            "timestamp": datetime.now().isoformat(),
            "analysis_mode": variant_mode.value,
            "hierarchy_mode": hierarchy_mode.value,
            "reference_compound": {
                "positional_block_sequence": reference.positional_block_sequence,
                "block_support_sequence": reference.block_support_sequence,
                "block_level": reference.level,
                "monomer_level": reference.monomer_level
            },
            "lineage_size": lineage_size
        }

        if data_file:
            metadata["data_file"] = str(data_file)

        if total_library_size:
            metadata["total_library_size"] = total_library_size
            metadata["lineage_reduction_percent"] = round(
                100 * lineage_size / total_library_size, 2
            )

        return metadata

    def _build_hierarchy_info(
        self,
        hierarchy: CompoundHierarchy,
        hierarchy_mode: HierarchyMode
    ) -> Dict[str, Any]:
        """Build hierarchy structure information."""
        level_attr = "monomer_level" if hierarchy_mode == HierarchyMode.MONOMER else "level"

        # Count compounds at each level
        level_counts = {}
        for compound in hierarchy.compounds:
            level = getattr(compound, level_attr)
            level_counts[level] = level_counts.get(level, 0) + 1

        # Get maximal and minimal compounds
        maximal = hierarchy.get_maximal_compounds()
        minimal = hierarchy.get_minimal_compounds()

        return {
            "mode": hierarchy_mode.value,
            "compounds_count": hierarchy.size(),
            "edges_count": hierarchy.edge_count(),
            "levels": {str(k): v for k, v in sorted(level_counts.items(), reverse=True)},
            "maximal_compounds": [c.positional_block_sequence for c in maximal],
            "minimal_compounds": [c.positional_block_sequence for c in minimal]
        }

    def _build_compounds_data(
        self,
        lineage: List[Compound],
        peaks_dict: Dict[Compound, List[Peak]],
        hierarchy: CompoundHierarchy,
        hierarchy_mode: HierarchyMode
    ) -> List[Dict[str, Any]]:
        """Build comprehensive compound data with peaks and relationships."""
        level_attr = "monomer_level" if hierarchy_mode == HierarchyMode.MONOMER else "level"
        level_label = "monomer_level" if hierarchy_mode == HierarchyMode.MONOMER else "block_level"

        compounds_data = []

        for compound in lineage:
            peaks = peaks_dict.get(compound, [])

            # Count peaks by type (only accepted peaks for main counts)
            accepted_peaks = [p for p in peaks if p.is_accepted]
            rejected_peaks = [p for p in peaks if p.is_rejected]
            product_peaks = sum(1 for p in accepted_peaks if p.is_product_peak)
            truncation_peaks = sum(1 for p in accepted_peaks if p.peak_type.value == "TRUNCATION")
            truncation_unknown_peaks = sum(1 for p in accepted_peaks if p.peak_type.value == "TRUNCATION_UNKNOWN")
            null_peaks = sum(1 for p in accepted_peaks if p.peak_type.value == "NULL")
            unknown_peaks = sum(1 for p in accepted_peaks if p.peak_type.value == "UNKNOWN")

            # Get descendants and ancestors
            descendants = hierarchy.get_descendants(compound)
            ancestors = hierarchy.get_ancestors(compound)

            compound_data = {
                "sequence": compound.positional_block_sequence,
                "block_support_sequence": compound.block_support_sequence,
                level_label: getattr(compound, level_attr),
                "peaks": [self._peak_to_dict(p) for p in peaks],
                "statistics": {
                    "total_peaks": len(peaks),
                    "accepted_peaks": len(accepted_peaks),
                    "rejected_peaks": len(rejected_peaks),
                    "product_peaks": product_peaks,
                    "truncation_peaks": truncation_peaks,
                    "truncation_unknown_peaks": truncation_unknown_peaks,
                    "null_peaks": null_peaks,
                    "unknown_peaks": unknown_peaks
                },
                "descendants": [d.positional_block_sequence for d in descendants],
                "ancestors": [a.positional_block_sequence for a in ancestors]
            }

            # Add cLPE reference data if present (AlogP and scaffold only)
            if compound.alogp is not None:
                compound_data["clpe_data"] = {
                    "alogp": round(compound.alogp, 4),
                    "scaffold_group": compound.scaffold_group
                }

            compounds_data.append(compound_data)

        # Sort same way as visualization (highest level first)
        compounds_data.sort(
            key=lambda x: (-x[level_label], x["block_support_sequence"])
        )

        return compounds_data

    def _peak_to_dict(self, peak: Peak) -> Dict[str, Any]:
        """Convert Peak entity to dictionary."""
        result = {
            "retention_time": round(peak.position, 4),
            "area": round(peak.area, 2),
            "height": round(peak.height, 2),
            "classification": peak.peak_type.value,
            "validation_status": peak.validation_status.value,
            "is_accepted": peak.is_accepted,
            "rejection_reason": peak.rejection_reason.value if peak.is_rejected else None,
            "p_value": round(peak.p_value, 6) if peak.p_value is not None else None,
            "left_base": round(peak.left_base, 4),
            "right_base": round(peak.right_base, 4),
            "width": round(peak.width, 4),
            "prominence": round(peak.prominence, 2) if peak.prominence else None,
            # Match tracking fields - trace peak origin through hierarchy
            "matched_compound_sequence": peak.matched_compound_sequence,
            "matched_peak_position": (
                round(peak.matched_peak_position, 4)
                if peak.matched_peak_position is not None else None
            ),
            "matched_peak_type": (
                peak.matched_peak_type.value
                if peak.matched_peak_type is not None else None
            ),
        }

        # Add cLPE validation fields if present
        if peak.clpe_residual is not None:
            result["clpe_validation"] = {
                "residual": round(peak.clpe_residual, 4),
                "z_score": round(peak.clpe_z_score, 2) if peak.clpe_z_score is not None else None,
                "is_outlier": peak.clpe_is_outlier,
                "reselected": peak.clpe_reselected
            }

        return result

    def _build_equivalence_classes_data(
        self,
        equivalence_classes: Dict[str, EquivalenceClass],
        peaks_dict: Dict[Compound, List[Peak]]
    ) -> List[Dict[str, Any]]:
        """Build equivalence class data (pooled mode)."""
        eq_classes_data = []

        for block_support_seq in sorted(equivalence_classes.keys()):
            eq_class = equivalence_classes[block_support_seq]

            # Get correlation stats (if applicable)
            correlation_data = {}
            if eq_class.correlation_min is not None:
                correlation_data["correlation_min"] = round(eq_class.correlation_min, 4)

            # Get member sequences
            members = sorted([m.positional_block_sequence for m in eq_class.members])

            # Build pooled peaks data
            pooled_peaks_data = [
                self._peak_to_dict(p) for p in eq_class.pooled_peaks
            ]

            eq_class_data = {
                "block_support_sequence": block_support_seq,
                "n_variants": len(eq_class.members),
                "members": members,
                "pooling_status": eq_class.pooling_status.value,
                "pooled_peaks": pooled_peaks_data,
                **correlation_data
            }

            eq_classes_data.append(eq_class_data)

        return eq_classes_data

    def _build_statistics(
        self,
        lineage: List[Compound],
        peaks_dict: Dict[Compound, List[Peak]],
        equivalence_classes: Optional[Dict[str, EquivalenceClass]]
    ) -> Dict[str, Any]:
        """Build summary statistics."""
        # Basic peak statistics
        all_peaks = []
        for peaks in peaks_dict.values():
            all_peaks.extend(peaks)

        # Separate accepted vs rejected peaks
        accepted_peaks = [p for p in all_peaks if p.is_accepted]
        rejected_peaks = [p for p in all_peaks if p.is_rejected]

        total_peaks = len(all_peaks)
        compounds_with_peaks = sum(1 for peaks in peaks_dict.values() if len(peaks) > 0)
        compounds_with_accepted = sum(
            1 for peaks in peaks_dict.values()
            if any(p.is_accepted for p in peaks)
        )
        compounds_without_peaks = len(lineage) - compounds_with_peaks

        stats = {
            "total_peaks": total_peaks,
            "accepted_peaks": len(accepted_peaks),
            "rejected_peaks": len(rejected_peaks),
            "peaks_per_compound_avg": round(total_peaks / len(lineage), 2) if lineage else 0,
            "accepted_peaks_per_compound_avg": round(len(accepted_peaks) / len(lineage), 2) if lineage else 0,
            "compounds_with_peaks": compounds_with_peaks,
            "compounds_with_accepted_peaks": compounds_with_accepted,
            "compounds_without_peaks": compounds_without_peaks
        }

        # Add pooled mode statistics
        if equivalence_classes:
            from ...domain.models import PoolingStatus

            high_corr = sum(
                1 for eq in equivalence_classes.values()
                if eq.pooling_status == PoolingStatus.POOLING_VALID
            )
            low_corr = sum(
                1 for eq in equivalence_classes.values()
                if eq.pooling_status == PoolingStatus.POOLING_INVALID
            )
            single = sum(
                1 for eq in equivalence_classes.values()
                if eq.pooling_status == PoolingStatus.NOT_ATTEMPTED
            )
            total_variants = sum(len(eq.members) for eq in equivalence_classes.values())

            stats.update({
                "equivalence_classes_count": len(equivalence_classes),
                "total_variants": total_variants,
                "high_correlation_classes": high_corr,
                "low_correlation_classes": low_corr,
                "single_variant_classes": single
            })

        return stats
