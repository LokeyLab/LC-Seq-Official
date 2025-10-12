"""
Report generator for creating comprehensive analysis reports.

Generates human-readable reports summarizing analysis results.
"""

from pathlib import Path
from datetime import datetime
from ...application.dtos.analysis_response import AnalysisResponse


class ReportGenerator:
    """
    Generate analysis reports in various formats.

    Creates comprehensive reports summarizing analysis results,
    validation statistics, and quality metrics.

    Examples
    --------
    >>> generator = ReportGenerator()
    >>> response = AnalysisResponse(...)
    >>> generator.generate_markdown(response, Path('report.md'))
    """

    def generate_markdown(
        self,
        response: AnalysisResponse,
        output_path: Path
    ) -> Path:
        """
        Generate Markdown report.

        Parameters
        ----------
        response : AnalysisResponse
            Analysis results
        output_path : Path
            Output file path

        Returns
        -------
        Path
            Path to generated report

        Examples
        --------
        >>> report_path = generator.generate_markdown(response, Path('report.md'))
        """
        output_path.parent.mkdir(parents=True, exist_ok=True)

        report = []
        report.append(f"# LC-Seq Analysis Report\n")
        report.append(f"**Generated**: {response.timestamp.strftime('%Y-%m-%d %H:%M:%S')}\n")
        report.append(f"**Request ID**: {response.request_id}\n")
        report.append("\n---\n")

        # Validation Summary
        summary = response.validation_summary
        report.append("\n## Validation Summary\n")
        report.append(f"- **Total Compounds**: {summary.total_compounds}\n")
        report.append(f"- **Validation Rate**: {summary.validation_rate:.1%}\n")
        report.append(f"- **Median Purity**: {summary.median_purity:.3f}\n")
        report.append(f"- **Dataset Quality**: {summary.dataset_quality}\n")
        report.append("\n### Status Distribution\n")
        report.append(f"- **Validated**: {summary.validated_count} ({summary.validated_count/summary.total_compounds:.1%})\n")
        report.append(f"- **Likely Success**: {summary.likely_success_count} ({summary.likely_success_count/summary.total_compounds:.1%})\n")
        report.append(f"- **Uncertain**: {summary.uncertain_count} ({summary.uncertain_count/summary.total_compounds:.1%})\n")
        report.append(f"- **Likely Failure**: {summary.likely_failure_count} ({summary.likely_failure_count/summary.total_compounds:.1%})\n")
        report.append(f"- **Failed**: {summary.failed_count} ({summary.failed_count/summary.total_compounds:.1%})\n")

        # Dataset Statistics
        stats = response.dataset_stats
        report.append("\n## Dataset Statistics\n")
        report.append(f"- **Purity P25**: {stats.get('purity_p25', 0):.3f}\n")
        report.append(f"- **Purity P50** (median): {stats.get('purity_p50', 0):.3f}\n")
        report.append(f"- **Purity P75**: {stats.get('purity_p75', 0):.3f}\n")
        report.append(f"- **Purity P90**: {stats.get('purity_p90', 0):.3f}\n")
        report.append(f"- **Purity MAD**: {stats.get('purity_mad', 0):.3f}\n")
        report.append(f"- **Background**: {stats.get('background', 0):.2f}\n")

        # Processing Metadata
        metadata = response.processing_metadata
        report.append("\n## Processing Information\n")
        report.append(f"- **Runtime**: {metadata.get('runtime_seconds', 0):.2f} seconds\n")
        report.append(f"- **Hierarchy Mode**: {metadata.get('hierarchy_mode', 'unknown')}\n")
        report.append(f"- **Variant Mode**: {metadata.get('variant_mode', 'unknown')}\n")

        # Top Performers
        report.append("\n## Top Performers (by Purity)\n")
        top_compounds = sorted(
            response.compound_results,
            key=lambda x: x.purity,
            reverse=True
        )[:10]

        report.append("\n| Rank | Compound ID | Sequence | Purity | Status |\n")
        report.append("|------|-------------|----------|--------|--------|\n")
        for rank, compound in enumerate(top_compounds, 1):
            report.append(
                f"| {rank} | {compound.compound_id} | {compound.sequence} | "
                f"{compound.purity:.3f} | {compound.validation_status} |\n"
            )

        # Problematic Compounds
        failed = response.get_failed_compounds()
        if failed:
            report.append("\n## Failed Compounds\n")
            report.append(f"\n{len(failed)} compounds failed validation:\n\n")
            report.append("| Compound ID | Sequence | Purity | SNR |\n")
            report.append("|-------------|----------|--------|-----|\n")
            for compound in failed[:20]:  # Limit to first 20
                report.append(
                    f"| {compound.compound_id} | {compound.sequence} | "
                    f"{compound.purity:.3f} | {compound.snr:.2f} |\n"
                )

        # Errors and Warnings
        if response.errors:
            report.append("\n## Errors\n")
            for error in response.errors:
                report.append(f"- {error}\n")

        if response.warnings:
            report.append("\n## Warnings\n")
            for warning in response.warnings:
                report.append(f"- {warning}\n")

        # Write report
        with open(output_path, 'w') as f:
            f.writelines(report)

        return output_path

    def generate_text_summary(
        self,
        response: AnalysisResponse
    ) -> str:
        """
        Generate plain text summary.

        Parameters
        ----------
        response : AnalysisResponse
            Analysis results

        Returns
        -------
        str
            Text summary

        Examples
        --------
        >>> summary = generator.generate_text_summary(response)
        >>> print(summary)
        LC-Seq Analysis Summary
        Total Compounds: 1000
        Validation Rate: 75.2%
        ...
        """
        summary = response.validation_summary
        lines = []
        lines.append("=" * 60)
        lines.append("LC-Seq Analysis Summary")
        lines.append("=" * 60)
        lines.append(f"Total Compounds: {summary.total_compounds}")
        lines.append(f"Validation Rate: {summary.validation_rate:.1%}")
        lines.append(f"Median Purity: {summary.median_purity:.3f}")
        lines.append(f"Dataset Quality: {summary.dataset_quality}")
        lines.append("-" * 60)
        lines.append(f"Validated: {summary.validated_count}")
        lines.append(f"Likely Success: {summary.likely_success_count}")
        lines.append(f"Uncertain: {summary.uncertain_count}")
        lines.append(f"Likely Failure: {summary.likely_failure_count}")
        lines.append(f"Failed: {summary.failed_count}")
        lines.append("=" * 60)

        return "\n".join(lines)
