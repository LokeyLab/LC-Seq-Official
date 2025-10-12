# LC-Seq

**Production-ready DNA-Encoded Library chromatographic data analysis**

[![Python 3.11+](https://img.shields.io/badge/python-3.11+-blue.svg)](https://www.python.org/downloads/)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![Coverage](https://img.shields.io/badge/coverage-67%25-yellow.svg)]()

A complete Python package for analyzing chromatography data from DNA-encoded library (DEL) screens. LC-Seq provides robust peak detection, hierarchical truncation analysis, and adaptive synthesis validation based on rigorous mathematical foundations.

## Features

- 🔬 **Rigorous Peak Detection**: Morse theory + persistent homology for noise-robust peak identification
- 📊 **Adaptive Validation**: Dataset-relative thresholds with Bayesian confidence assessment
- 🌳 **Hierarchical Analysis**: DAG-based truncation analysis with automatic ancestor-descendant relationships
- 🎯 **Pooled Mode**: Automatic aggregation of positional variants with fallback handling
- 🏗️ **Clean Architecture**: Production-ready design with clear layer separation
- 🔒 **Type-Safe**: Full type hints and static checking with mypy strict mode
- ✅ **Well-Tested**: 590+ tests with 67% coverage (91% domain layer)
- 🚀 **User-Friendly CLI**: Simple commands with rich output formatting

## Installation

### From source (development)

```bash
# Clone the repository
git clone https://github.com/yourusername/LC-Seq-Official.git
cd LC-Seq-Official

# Create virtual environment
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install in editable mode with dev dependencies
pip install -e ".[dev]"

# Install pre-commit hooks
pre-commit install
```

### For users

```bash
pip install lcseq
```

## Quick Start

### Command Line Interface

```bash
# Install with CLI support
pip install -e ".[cli]"

# Run analysis with default settings
lcseq analyze data.csv --library design.csv --output results/

# Use custom configuration
lcseq analyze data.csv --library design.csv --config configs/default.yaml

# Enable pooled mode for positional variants
lcseq analyze data.csv --library design.csv --variant-mode pooled

# Show capabilities
lcseq info
```

### Programmatic Usage

```python
from lcseq import Compound, Chromatogram, Peak, AnalysisConfiguration
from lcseq.application.pipelines.full_analysis_pipeline import FullAnalysisPipeline
from lcseq.application.dtos.analysis_request import AnalysisRequest

# Create analysis request
request = AnalysisRequest(
    data_path="data.csv",
    library_path="design.csv",
    configuration=AnalysisConfiguration.get_default(),
)

# Run pipeline
pipeline = FullAnalysisPipeline()
response = pipeline.execute(request)

# Examine results
print(f"Validation rate: {response.validation_summary.validation_rate:.1%}")
print(f"Median purity: {response.validation_summary.median_purity:.3f}")
```

**See [QUICKSTART.md](docs/QUICKSTART.md) for complete usage guide.**

## Development

### Project Structure

```
LC-Seq-Official/
├── src/lcseq/              # Main package (Clean Architecture)
│   ├── domain/             # Pure domain logic (zero external dependencies)
│   │   ├── entities/       # Core entities (23 files)
│   │   ├── value_objects/  # Value objects (4 files)
│   │   ├── models/         # Graph models (5 files)
│   │   └── services/       # Domain services (13 files)
│   ├── application/        # Use cases and orchestration
│   │   ├── pipelines/      # Full analysis pipeline
│   │   └── dtos/           # Data transfer objects
│   ├── infrastructure/     # I/O adapters
│   │   ├── parsers/        # CSV/Excel parsing
│   │   ├── exporters/      # Result export
│   │   ├── repositories/   # Data persistence
│   │   └── configuration/  # YAML config loading
│   └── interfaces/         # External interfaces
│       └── cli/            # Command-line interface
├── tests/                  # Test suite (590+ tests, 67% coverage)
│   ├── test_domain/        # Domain layer tests (584 tests)
│   ├── test_integration/   # Integration tests
│   └── fixtures/           # Test data
├── configs/                # Configuration
│   └── default.yaml        # Single Source of Truth for all configuration
├── docs/                   # Documentation
│   ├── THEORY.md           # Mathematical foundations (3,772 lines)
│   ├── ARCHITECTURE.md     # Clean architecture guide
│   ├── IMPLEMENTATION_PLAN.md  # 8-week roadmap
│   ├── QUICKSTART.md       # User guide
│   └── PHASE_1_4_REPORT.md # Implementation report
└── examples/               # Example datasets and workflows
```

### Running Tests

```bash
# Run all tests
pytest

# Run with coverage
pytest --cov=lcseq --cov-report=html

# Run specific test file
pytest tests/test_domain/test_chromatogram.py

# Run with verbose output
pytest -v
```

### Code Quality

```bash
# Format code
black src/ tests/

# Lint code
ruff check src/ tests/

# Type check
mypy src/

# Run all checks (pre-commit)
pre-commit run --all-files
```

### Adding New Features

1. **Create a new branch**
   ```bash
   git checkout -b feature/my-new-feature
   ```

2. **Write tests first** (test-driven development)
   ```bash
   # Add tests in tests/test_*/
   pytest tests/test_algorithms/test_my_feature.py
   ```

3. **Implement the feature**
   - Add code to appropriate module in `src/lcseq/`
   - Follow type hints and docstring conventions
   - Keep functions focused and testable

4. **Verify code quality**
   ```bash
   pre-commit run --all-files
   pytest
   ```

5. **Create pull request**

## Architecture

The package follows **Clean Architecture** (Hexagonal Architecture) principles:

1. **Domain Layer**: Pure domain logic with zero external dependencies
   - Entities: `Compound`, `Chromatogram`, `Peak`, `BuildingBlock`
   - Algorithms: Morse theory peak detection, DAG-based classification
   - Services: Hierarchy building, validation, pooled mode analysis

2. **Application Layer**: Use cases and orchestration
   - Use cases coordinate domain services
   - Ports define I/O boundaries (file loading, storage, configuration)

3. **Infrastructure Layer**: I/O adapters
   - HDF5 persistence, YAML configuration, CSV/JSON export

4. **Interface Layer**: External adapters
   - CLI (Typer-based), visualization (matplotlib/plotly)

**Key Principle**: Dependencies always flow inward (Interface → Infrastructure → Application → Domain)

See [`docs/ARCHITECTURE.md`](docs/ARCHITECTURE.md) for complete architectural specifications.

## Implementation Status

**Current Status**: Production-ready v0.1.0 with complete core functionality

### Completed ✅

**Phase 1: Domain Layer** (23 components)
- ✅ Core entities: BuildingBlock, Compound, Chromatogram, Peak
- ✅ Value objects: RetentionTime, PeakBoundaries, BuildingBlockSequence, MonomerSequence
- ✅ Graph models: CompoundHierarchy, EquivalenceClass, PeakClassification
- ✅ Domain services: 13 services including HierarchyBuilder, PeakDetector, ValidationChecker
- ✅ 584 unit tests with 91% domain layer coverage

**Phase 2: Synthesis Validation Layer** (3 modules)
- ✅ BayesianValidator with prior probabilities
- ✅ Adaptive validation with dataset-relative thresholds
- ✅ Pooling validator for positional variants

**Phase 3: Application Layer** (4 modules)
- ✅ FullAnalysisPipeline orchestrating complete workflow
- ✅ DTOs for request/response boundaries
- ✅ Composition root with dependency injection

**Phase 4: Infrastructure Layer** (12 modules)
- ✅ CSV/Excel parsers and exporters
- ✅ Library repository for compound management
- ✅ Configuration infrastructure with YAML support
- ✅ Result export to multiple formats

**Phase 5: Interface Layer** (CLI)
- ✅ Typer-based CLI with rich formatting
- ✅ Comprehensive YAML configuration (single source of truth)
- ✅ Progress reporting and error handling
- ✅ QUICKSTART.md comprehensive user guide

### Statistics
- **Total Code**: ~9,938 lines (domain + application + infrastructure)
- **Total Tests**: 590+ tests with 67% overall coverage
- **Domain Coverage**: 91% (high confidence in core logic)
- **Configuration**: Single SSoT in configs/default.yaml
- **Documentation**: 3,772 lines of theory + architecture + guides

### Deferred to v0.2.0
- ⏳ Visualization layer (matplotlib/plotly adapters)
- ⏳ HDF5 data format support
- ⏳ Property-based testing with Hypothesis
- ⏳ Performance optimization for >10,000 compound libraries

## Contributing

Contributions are welcome! Please:

1. Fork the repository
2. Create a feature branch
3. Add tests for new functionality
4. Ensure all tests pass and code quality checks succeed
5. Submit a pull request

## Citation

If you use this software in your research, please cite:

```
[Citation to be added]
```

## License

MIT License - see LICENSE file for details

## Contact

For questions or issues, please open an issue on GitHub.
