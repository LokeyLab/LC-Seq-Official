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
from lcseq.application.dtos.analysis_request import AnalysisRequest, CLPEParams
from lcseq.domain.models.analysis_configuration import (
    AnalysisConfiguration,
    AnalysisMode,
    HierarchyMode,
)
from lcseq.infrastructure.loaders.hdf5_compound_loader import HDF5CompoundLoader
from lcseq.infrastructure.configuration.yaml_loader import ConfigurationLoader
from lcseq.domain.services.baseline_estimator import MinimaSplineParams
from lcseq.presentation.visualization.plotters.baseline_debug_plotter import plot_minima_spline_baseline

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
        help="Variant mode: 'individual' or 'pooled'",
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
    clpe_reference: Optional[Path] = typer.Option(
        None,
        "--clpe-reference",
        exists=True,
        dir_okay=False,
        help="CSV with AlogP and scaffold data for cLPE validation",
    ),
    clpe_threshold: Optional[float] = typer.Option(
        None,
        "--clpe-threshold",
        help="Z-score threshold for cLPE outlier detection (default: from config)",
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

        # Using pooled mode
        lcseq analyze data.csv --library design.csv --variant-mode pooled

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
        analysis_config = ConfigurationLoader.load_from_yaml(config)
        # Override mode settings from CLI if different from config
        if variant_mode == "individual":
            analysis_config.analysis_mode = AnalysisMode.INDIVIDUAL
        elif variant_mode == "pooled":
            analysis_config.analysis_mode = AnalysisMode.POOLED
        if hierarchy_mode == "block":
            analysis_config.hierarchy_mode = HierarchyMode.BUILDING_BLOCK
        elif hierarchy_mode == "monomer":
            analysis_config.hierarchy_mode = HierarchyMode.MONOMER
    else:
        # Use default config from configs/default.yaml
        analysis_config = ConfigurationLoader.get_default_config()
        # Override mode settings from CLI
        if variant_mode == "individual":
            analysis_config.analysis_mode = AnalysisMode.INDIVIDUAL
        elif variant_mode == "pooled":
            analysis_config.analysis_mode = AnalysisMode.POOLED
        if hierarchy_mode == "block":
            analysis_config.hierarchy_mode = HierarchyMode.BUILDING_BLOCK
        elif hierarchy_mode == "monomer":
            analysis_config.hierarchy_mode = HierarchyMode.MONOMER

    # Display configuration
    if verbose:
        _display_config(analysis_config)

    # Create cLPE params if reference provided
    clpe_params = None
    if clpe_reference:
        # Use CLI value if provided, otherwise use config value
        effective_clpe_threshold = (
            clpe_threshold if clpe_threshold is not None
            else analysis_config.validation_params['clpe_outlier_threshold']
        )
        effective_min_group_size = analysis_config.validation_params['clpe_min_group_size']
        clpe_params = CLPEParams(
            enabled=True,
            reference_csv_path=clpe_reference,
            outlier_threshold=effective_clpe_threshold,
            min_group_size=effective_min_group_size,
            reselect_peaks=True,
        )
        console.print(f"[cyan]cLPE validation enabled:[/cyan] {clpe_reference}")
        console.print(f"  Outlier threshold: {effective_clpe_threshold} z-score")

    # Create analysis request with config params
    from lcseq.domain.services.signal_preprocessor import PreprocessingConfig
    request = AnalysisRequest(
        data_path=data_path,
        output_path=output,
        library_path=library,
        variant_mode=variant_mode,
        hierarchy_mode='monomer' if hierarchy_mode == 'monomer' else 'building_block',
        detection_params=analysis_config.peak_detection_params,
        validation_params=analysis_config.validation_params,
        preprocessing_params=PreprocessingConfig.from_dict(analysis_config.preprocessing_params),
        clpe_params=clpe_params,
        num_workers=analysis_config.performance_params.get('num_workers'),
    )

    # Load compounds from data file
    console.print(f"\n[cyan]Loading data from:[/cyan] {data_path}")
    loader = HDF5CompoundLoader()
    compounds = loader.load_all(data_path)
    console.print(f"[green]Loaded {len(compounds):,} compounds[/green]")

    # Run analysis
    console.print(f"\n[cyan]Starting analysis...[/cyan]")
    start_time = datetime.now()

    try:
        pipeline = FullAnalysisPipeline()
        response = pipeline.execute(compounds, request)

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


@app.command("filter")
def filter_library(
    data_path: Path = typer.Argument(
        ...,
        exists=True,
        dir_okay=False,
        help="Path to HDF5 data file containing compound chromatograms",
    ),
    output: Path = typer.Option(
        Path("filtered"),
        "--output",
        "-o",
        help="Output directory for filtered files",
    ),
    config: Optional[Path] = typer.Option(
        None,
        "--config",
        "-c",
        exists=True,
        dir_okay=False,
        help="Path to configuration YAML file (uses quality_filter section)",
    ),
    min_correlation: float = typer.Option(
        0.8,
        "--min-correlation",
        "-r",
        help="Minimum replicate correlation within equivalence class (0.0-1.0)",
    ),
    intensity_percentile: Optional[float] = typer.Option(
        0.05,
        "--intensity-percentile",
        "-p",
        help="Exclude compounds below this percentile (e.g., 0.05 = bottom 5%)",
    ),
    intensity_absolute: Optional[float] = typer.Option(
        None,
        "--intensity-min",
        "-m",
        help="Absolute minimum total signal (overrides percentile)",
    ),
    max_noise_ratio: Optional[float] = typer.Option(
        None,
        "--max-noise-ratio",
        "-n",
        help="Maximum noise_std / median_signal ratio",
    ),
    no_hdf5: bool = typer.Option(
        False,
        "--no-hdf5",
        help="Skip generating filtered HDF5 file",
    ),
    no_csv: bool = typer.Option(
        False,
        "--no-csv",
        help="Skip generating inclusion CSV",
    ),
    no_report: bool = typer.Option(
        False,
        "--no-report",
        help="Skip generating QC report JSON",
    ),
) -> None:
    """
    Filter library based on signal quality metrics.

    Pre-filters low-quality data before downstream analysis using:

    1. REPLICATE CORRELATION: Positional variants within each equivalence
       class should have similar chromatograms. Low correlation indicates
       quality issues (e.g., noise, failed synthesis, sample degradation).

    2. SIGNAL INTENSITY: Very low total signal indicates failed synthesis
       or sample loss. Filters by percentile or absolute threshold.

    3. NOISE LEVEL: High baseline noise relative to signal median makes
       peak detection unreliable.

    Single-variant equivalence classes are always included (can't assess
    replicate correlation for classes with only one member).

    Output files:
    - filtered_library.h5: HDF5 with only passing compounds
    - included_sequences.csv: List of passing compound sequences
    - qc_report.json: Detailed quality metrics for all classes

    Examples:
        # Basic filtering with defaults
        lcseq filter data.h5 --output filtered/

        # Strict correlation threshold
        lcseq filter data.h5 -r 0.9 --output filtered/

        # Filter bottom 10% by intensity
        lcseq filter data.h5 -p 0.10 --output filtered/

        # Use absolute intensity threshold
        lcseq filter data.h5 --intensity-min 1000000 --output filtered/

        # Generate only QC report (no filtered files)
        lcseq filter data.h5 --no-hdf5 --no-csv --output qc/
    """
    console.print(f"\n[bold blue]LC-Seq Quality Filter[/bold blue]", justify="center")
    console.print("=" * 70 + "\n")

    # Create output directory
    output.mkdir(parents=True, exist_ok=True)

    # Load configuration if provided
    if config:
        console.print(f"[cyan]Loading configuration:[/cyan] {config}")
        analysis_config = ConfigurationLoader.load_from_yaml(config)
        # Use config values as defaults, CLI args override
        qf_params = analysis_config.quality_filter_params
        if qf_params:
            # Only use config value if CLI arg wasn't explicitly changed from default
            if min_correlation == 0.8 and "min_correlation" in qf_params:
                min_correlation = qf_params["min_correlation"]
            if intensity_percentile == 0.05 and "intensity_percentile" in qf_params:
                intensity_percentile = qf_params["intensity_percentile"]
            if intensity_absolute is None and qf_params.get("intensity_absolute"):
                intensity_absolute = qf_params["intensity_absolute"]
            if max_noise_ratio is None and qf_params.get("max_noise_ratio"):
                max_noise_ratio = qf_params["max_noise_ratio"]
    else:
        console.print("[cyan]Using CLI parameters[/cyan]")

    # Display filter settings
    console.print(f"\n[cyan]Filter Settings:[/cyan]")
    console.print(f"  Min correlation: {min_correlation}")
    if intensity_absolute is not None:
        console.print(f"  Intensity threshold: {intensity_absolute:.0f} (absolute)")
    elif intensity_percentile is not None:
        console.print(f"  Intensity threshold: {intensity_percentile*100:.0f}th percentile")
    else:
        console.print(f"  Intensity threshold: disabled")
    if max_noise_ratio is not None:
        console.print(f"  Max noise ratio: {max_noise_ratio}")
    else:
        console.print(f"  Noise filter: disabled")

    # Run filtering
    from lcseq.application.use_cases import FilterLibraryUseCase

    use_case = FilterLibraryUseCase()

    try:
        result = use_case.execute(
            hdf5_path=data_path,
            output_dir=output,
            min_correlation=min_correlation,
            intensity_percentile=intensity_percentile if intensity_absolute is None else None,
            intensity_absolute=intensity_absolute,
            max_noise_ratio=max_noise_ratio,
            generate_hdf5=not no_hdf5,
            generate_csv=not no_csv,
            generate_report=not no_report,
        )

        # Display summary
        console.print(f"\n[green]✓ Filtering complete![/green]")
        console.print(f"  Passed: {result.passed_compounds:,} / {result.total_compounds:,} compounds")
        console.print(f"  Classes: {result.passed_equivalence_classes:,} / {result.total_equivalence_classes:,}")
        console.print(f"  Time: {result.processing_time_seconds:.1f}s\n")

        if result.output_hdf5:
            console.print(f"[cyan]Filtered data:[/cyan] {result.output_hdf5}")
        if result.inclusion_csv:
            console.print(f"[cyan]Inclusion list:[/cyan] {result.inclusion_csv}")
        if result.qc_report:
            console.print(f"[cyan]QC report:[/cyan] {result.qc_report}")

    except Exception as e:
        console.print(f"\n[red]Error during filtering:[/red] {e}")
        raise typer.Exit(1)


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
    console.print("  • Pooled mode (automatic fallback)")
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
    analysis_mode = AnalysisMode.INDIVIDUAL if variant_mode == "individual" else AnalysisMode.POOLED
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

    # Detection parameters
    table.add_row("─ Detection ─", "─" * 20)
    table.add_row("  Alpha (significance)", str(config.peak_detection_params.get("alpha", "N/A")))
    table.add_row("  Alpha Product", str(config.peak_detection_params.get("alpha_product", "N/A")))
    table.add_row("  Min Baseline SDs", str(config.peak_detection_params.get("min_baseline_sds", "N/A")))
    table.add_row("  Min SNR", str(config.peak_detection_params.get("min_snr", "N/A")))
    table.add_row("  Prominence Percentile", str(config.peak_detection_params.get("prominence_percentile", "N/A")))
    table.add_row("  Min Persistence", str(config.peak_detection_params.get("min_persistence", "N/A")))
    table.add_row("  Boundary Method", str(config.peak_detection_params.get("boundary_method", "N/A")))
    table.add_row("  Signal Variant", str(config.peak_detection_params.get("signal_variant", "N/A")))

    # Validation parameters
    table.add_row("─ Validation ─", "─" * 20)
    table.add_row("  Purity Threshold", str(config.validation_params.get("purity_threshold", "N/A")))
    table.add_row("  SNR Threshold", str(config.validation_params.get("snr_threshold", "N/A")))

    # Performance
    num_workers = config.performance_params.get("num_workers")
    workers_str = "sequential" if num_workers is None or num_workers == 1 else str(num_workers)
    table.add_row("─ Performance ─", "─" * 20)
    table.add_row("  Workers", workers_str)

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


@app.command("baseline-debug")
def baseline_debug(
    data_path: Path = typer.Argument(
        ...,
        exists=True,
        dir_okay=False,
        help="Path to HDF5 data file",
    ),
    compound: str = typer.Option(
        None,
        "--compound",
        "-c",
        help="Compound sequence to analyze (e.g., 'Leu-LA03-Pro'). If not specified, uses first compound.",
    ),
    index: int = typer.Option(
        None,
        "--index",
        "-i",
        help="Compound index to analyze (0-based). Alternative to --compound.",
    ),
    output: Path = typer.Option(
        Path("baseline_debug.png"),
        "--output",
        "-o",
        help="Output path for debug figure",
    ),
    order: int = typer.Option(
        3,
        "--order",
        help="Minima detection order (higher = fewer minima)",
    ),
    smoothing: float = typer.Option(
        None,
        "--smoothing",
        "-s",
        help="Spline smoothing factor (default: auto)",
    ),
    no_start: bool = typer.Option(
        False,
        "--no-start",
        help="Don't include start point as virtual minimum (default: included)",
    ),
    no_end: bool = typer.Option(
        False,
        "--no-end",
        help="Don't include end point as virtual minimum (default: included)",
    ),
) -> None:
    """
    Generate baseline correction debug figure for a compound.

    Fits a smoothed spline through local minima. Shows:
    - Original signal (green, with markers)
    - Baseline (red)
    - Corrected signal (black)

    Examples:
        # By compound sequence
        lcseq baseline-debug data.h5 --compound "Leu-LA03-Pro-Leu-DLeuMe"

        # By index
        lcseq baseline-debug data.h5 --index 100

        # With custom smoothing
        lcseq baseline-debug data.h5 -i 100 --smoothing 5.0
    """
    console.print(f"\n[bold blue]Baseline Debug[/bold blue]\n")

    # Load data
    console.print(f"[cyan]Loading data from:[/cyan] {data_path}")
    loader = HDF5CompoundLoader()
    compounds = loader.load_all(data_path)

    # Find compound
    if compound:
        target = next((c for c in compounds if c.positional_block_sequence == compound), None)
        if target is None:
            # Try partial match
            target = next((c for c in compounds if compound in c.positional_block_sequence), None)
        if target is None:
            console.print(f"[red]Error:[/red] Compound '{compound}' not found")
            raise typer.Exit(1)
    elif index is not None:
        if index < 0 or index >= len(compounds):
            console.print(f"[red]Error:[/red] Index {index} out of range (0-{len(compounds)-1})")
            raise typer.Exit(1)
        target = compounds[index]
    else:
        target = compounds[0]

    console.print(f"[green]Compound:[/green] {target.positional_block_sequence}")
    console.print(f"[green]Signal length:[/green] {len(target.chromatogram.counts)} points")

    # Fit baseline
    params = MinimaSplineParams(
        order=order,
        smoothing=smoothing,
        include_start=not no_start,
        include_end=not no_end,
    )
    fig, result = plot_minima_spline_baseline(
        signal=target.chromatogram.counts,
        time_axis=target.chromatogram.time_points,
        params=params,
        output_path=output,
        title=f"Baseline: {target.positional_block_sequence}",
    )

    console.print(f"\n[cyan]Results:[/cyan]")
    console.print(f"  Detected minima: {len(result.minima_indices)}")
    console.print(f"  Baseline range: [{result.baseline.min():.1f}, {result.baseline.max():.1f}]")
    console.print(f"  Corrected range: [{result.corrected.min():.1f}, {result.corrected.max():.1f}]")
    console.print(f"\n[green]✓ Saved to {output}[/green]\n")


if __name__ == "__main__":
    app()
