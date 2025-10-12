"""
LC-Seq Command Line Interface.

Provides user-friendly commands for analyzing DEL chromatography data.
"""

from pathlib import Path
from typing import Optional
import typer
from rich.console import Console
from rich.table import Table
from datetime import datetime
import json

from lcseq.application.pipelines.full_analysis_pipeline import FullAnalysisPipeline
from lcseq.application.dtos.analysis_request import AnalysisRequest
from lcseq.domain.models.analysis_configuration import (
    AnalysisConfiguration,
    AnalysisMode,
    HierarchyMode,
)

app = typer.Typer(
    name="lcseq",
    help="LC-Seq: DNA-encoded library chromatographic data analysis",
    add_completion=False,
)
console = Console()


@app.command()
def analyze(
    data_path: Path = typer.Argument(
        ...,
        exists=True,
        dir_okay=False,
        help="Path to chromatography data file (CSV, Excel, or HDF5)",
    ),
    library: Path = typer.Option(
        None,
        "--library",
        "-l",
        exists=True,
        dir_okay=False,
        help="Path to library design file (CSV or Excel)",
    ),
    output: Path = typer.Option(
        Path("results"),
        "--output",
        "-o",
        help="Output directory for results",
    ),
    config: Optional[Path] = typer.Option(
        None,
        "--config",
        "-c",
        exists=True,
        dir_okay=False,
        help="Path to configuration YAML file",
    ),
    variant_mode: str = typer.Option(
        "individual",
        "--variant-mode",
        "-v",
        help="Variant mode: 'individual' or 'consensus'",
    ),
    hierarchy_mode: str = typer.Option(
        "block",
        "--hierarchy-mode",
        "-h",
        help="Hierarchy mode: 'block' or 'monomer'",
    ),
    format: str = typer.Option(
        "csv",
        "--format",
        "-f",
        help="Output format: 'csv', 'excel', or 'json'",
    ),
    verbose: bool = typer.Option(
        False,
        "--verbose",
        "-V",
        help="Enable verbose output",
    ),
) -> None:
    """
    Analyze DEL chromatography data.

    This command performs complete analysis of DNA-encoded library
    chromatography data including:
    - Peak detection (Morse theory + persistent homology)
    - Peak classification (NULL/TRUNCATION/PUTATIVE_PRODUCT)
    - Synthesis validation (adaptive thresholds)

    Examples:
        # Basic analysis
        lcseq analyze data.csv --library design.csv

        # With custom output directory
        lcseq analyze data.csv --library design.csv --output my_results/

        # Using consensus mode
        lcseq analyze data.csv --library design.csv --variant-mode consensus

        # With configuration file
        lcseq analyze data.csv --library design.csv --config custom.yaml
    """
    console.print(f"\n[bold blue]LC-Seq Analysis[/bold blue]", justify="center")
    console.print("=" * 70 + "\n")

    # Validate inputs
    if not data_path.exists():
        console.print(f"[red]Error:[/red] Data file not found: {data_path}")
        raise typer.Exit(1)

    if library and not library.exists():
        console.print(f"[red]Error:[/red] Library file not found: {library}")
        raise typer.Exit(1)

    # Create output directory
    output.mkdir(parents=True, exist_ok=True)

    # Load configuration
    if config:
        console.print(f"[cyan]Loading configuration:[/cyan] {config}")
        # TODO: Implement YAML config loading
        analysis_config = _get_default_config(variant_mode, hierarchy_mode)
    else:
        analysis_config = _get_default_config(variant_mode, hierarchy_mode)

    # Display configuration
    if verbose:
        _display_config(analysis_config)

    # Create analysis request
    request = AnalysisRequest(
        data_path=str(data_path),
        library_path=str(library) if library else None,
        configuration=analysis_config,
        output_directory=str(output),
    )

    # Run analysis
    console.print(f"\n[cyan]Starting analysis...[/cyan]")
    start_time = datetime.now()

    try:
        pipeline = FullAnalysisPipeline()
        response = pipeline.execute(request)

        elapsed = (datetime.now() - start_time).total_seconds()

        # Display results summary
        console.print(f"\n[green]Analysis complete![/green] ({elapsed:.2f}s)\n")
        _display_summary(response)

        # Export results
        output_file = output / f"results_{response.request_id}.{format}"
        console.print(f"\n[cyan]Exporting results to:[/cyan] {output_file}")
        _export_results(response, output_file, format)

        # Display errors/warnings
        if response.errors:
            console.print("\n[yellow]Errors:[/yellow]")
            for error in response.errors:
                console.print(f"  • {error}")

        if response.warnings:
            console.print("\n[yellow]Warnings:[/yellow]")
            for warning in response.warnings:
                console.print(f"  • {warning}")

        console.print(f"\n[green]✓ Results saved to {output}/[/green]\n")

    except Exception as e:
        console.print(f"\n[red]Error during analysis:[/red] {e}")
        if verbose:
            console.print_exception()
        raise typer.Exit(1)


@app.command()
def validate(
    results_file: Path = typer.Argument(
        ...,
        exists=True,
        dir_okay=False,
        help="Path to results file to validate",
    ),
) -> None:
    """
    Validate analysis results file format.

    Check that a results file is properly formatted and contains
    all required fields.
    """
    console.print(f"\n[cyan]Validating results file:[/cyan] {results_file}")

    # TODO: Implement validation
    console.print("[yellow]Validation not yet implemented[/yellow]")


@app.command()
def info() -> None:
    """
    Display LC-Seq version and configuration information.
    """
    from lcseq import __version__

    console.print(f"\n[bold]LC-Seq Version:[/bold] {__version__}")
    console.print("[bold]DNA-Encoded Library Chromatographic Data Analysis[/bold]\n")

    console.print("[cyan]Capabilities:[/cyan]")
    console.print("  • Peak detection (Morse theory + persistent homology)")
    console.print("  • Peak classification (DAG constraints)")
    console.print("  • Synthesis validation (adaptive thresholds)")
    console.print("  • Consensus mode (automatic fallback)")
    console.print("  • Hierarchical compound organization\n")

    console.print("[cyan]Supported Formats:[/cyan]")
    console.print("  • Input: CSV, Excel, HDF5")
    console.print("  • Output: CSV, Excel, JSON\n")


def _get_default_config(
    variant_mode: str,
    hierarchy_mode: str,
) -> AnalysisConfiguration:
    """Create default analysis configuration."""
    # Map string to enum
    analysis_mode = AnalysisMode.INDIVIDUAL if variant_mode == "individual" else AnalysisMode.CONSENSUS
    hierarchy = HierarchyMode.BUILDING_BLOCK if hierarchy_mode == "block" else HierarchyMode.MONOMER

    return AnalysisConfiguration(
        analysis_mode=analysis_mode,
        hierarchy_mode=hierarchy,
    )


def _display_config(config: AnalysisConfiguration) -> None:
    """Display configuration in a table."""
    table = Table(title="Analysis Configuration")
    table.add_column("Parameter", style="cyan")
    table.add_column("Value", style="green")

    table.add_row("Analysis Mode", config.analysis_mode.value)
    table.add_row("Hierarchy Mode", config.hierarchy_mode.value)
    table.add_row("Min Persistence", f"{config.peak_detection_params.get('min_persistence', 0.05)}")
    table.add_row("Boundary Method", config.peak_detection_params.get("boundary_method", "valley_or_5pct"))
    table.add_row("Purity Threshold", config.validation_params.get("purity_threshold", "auto"))
    table.add_row("SNR Threshold", config.validation_params.get("snr_threshold", "auto"))

    console.print(table)


def _display_summary(response) -> None:
    """Display analysis results summary."""
    summary = response.validation_summary

    # Main results table
    table = Table(title="Analysis Results")
    table.add_column("Metric", style="cyan")
    table.add_column("Value", style="green", justify="right")

    table.add_row("Total Compounds", str(summary.total_compounds))
    table.add_row("Validated", str(summary.validated_count))
    table.add_row("Likely Success", str(summary.likely_success_count))
    table.add_row("Uncertain", str(summary.uncertain_count))
    table.add_row("Likely Failure", str(summary.likely_failure_count))
    table.add_row("Failed", str(summary.failed_count))
    table.add_row("Validation Rate", f"{summary.validation_rate * 100:.1f}%")
    table.add_row("Median Purity", f"{summary.median_purity:.3f}")
    table.add_row("Dataset Quality", summary.dataset_quality)

    console.print(table)


def _export_results(response, output_file: Path, format: str) -> None:
    """Export results to specified format."""
    if format == "json":
        # Export as JSON
        with open(output_file, "w") as f:
            json.dump(response.model_dump(), f, indent=2, default=str)

    elif format == "csv":
        # TODO: Implement CSV export using infrastructure layer
        console.print("[yellow]CSV export not yet implemented[/yellow]")

    elif format == "excel":
        # TODO: Implement Excel export using infrastructure layer
        console.print("[yellow]Excel export not yet implemented[/yellow]")

    else:
        console.print(f"[red]Unknown format:[/red] {format}")


if __name__ == "__main__":
    app()
