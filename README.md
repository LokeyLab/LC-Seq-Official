# LC-Seq

**DNA-Encoded Library Chromatographic Data Analysis**

[![Python 3.11+](https://img.shields.io/badge/python-3.11+-blue.svg)](https://www.python.org/downloads/)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

A Python package for analyzing chromatography data from DNA-encoded library (DEL) screens. LC-Seq provides mathematically rigorous peak detection, hierarchical truncation analysis, and adaptive synthesis validation based on Discrete Morse Theory, DAG-based classification, and Bayesian validation.

## Features

- 🔬 **Rigorous Peak Detection**: Discrete Morse Theory for local maxima detection with Poisson statistics for significance testing
- 📊 **Adaptive Validation**: Dataset-relative thresholds with Bayesian confidence assessment (no magic numbers!)
- 🌳 **Hierarchical Analysis**: DAG-based truncation analysis with automatic ancestor-descendant relationships
- 🎯 **Pooled Mode**: Hybrid strategy for positional variant aggregation with automatic correlation-based fallback
- 🏗️ **Clean Architecture**: Clear separation of domain, application, infrastructure, and presentation layers
- 🔒 **Type-Safe**: Full type hints with mypy strict mode compatibility
- 📚 **Comprehensive Theory**: 2,270-line [THEORY.md](docs/THEORY.md) documenting mathematical foundations
- 🔧 **Working Examples**: Fully functional analysis pipeline in `examples/`

## Installation

### From source (development)

```bash
# Clone the repository
git clone https://github.com/LokeyLab/LC-Seq-Official.git
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

### Using the Example Analysis Script

The script supports **4 analysis modes** (2 variant modes × 2 hierarchy modes):

```bash
# 1. Individual + Monomer (default) - Full detail, chemical identity focus
python examples/analyze.py --reference "Phe-DNvl-DPhe" --hierarchy-mode monomer --variant-mode individual

# 2. Individual + Building-Block - Full detail, synthesis position focus
python examples/analyze.py --reference "Phe-DNvl-DPhe" --hierarchy-mode building_block

# 3. Pooled + Monomer - Fast processing, chemical identity focus
python examples/analyze.py --reference "Phe-DNvl-DPhe" --variant-mode pooled --hierarchy-mode monomer

# 4. Pooled + Building-Block - Fast processing, synthesis position focus
python examples/analyze.py --reference "Phe-DNvl-DPhe" --variant-mode pooled --hierarchy-mode building_block
```

**Data Requirements**: The example script works with the included test dataset (`test_data/processed_data.h5`). To use your own data, format it as HDF5 with the same structure (see [examples/README.md](examples/README.md) for details).

**Customizing Parameters**: All analysis parameters are controlled by `configs/default.yaml`. To customize:

1. Copy the default config: `cp configs/default.yaml configs/my_config.yaml`
2. Edit your copy to adjust parameters (peak detection thresholds, validation criteria, etc.)
3. Run with your config: `python examples/analyze.py --reference "Phe-DNvl-DPhe" --config configs/my_config.yaml`

See the **Configuration** section below for details on available parameters.

See [examples/README.md](examples/README.md) for complete usage guide and all mode combinations.

### Programmatic Usage

```python
from lcseq.domain.entities import Compound, Chromatogram, Peak
from lcseq.domain.services import HierarchyBuilder, PeakDetector
from lcseq.infrastructure.configuration import ConfigurationLoader
from lcseq.infrastructure.loaders import HDF5CompoundLoader

# Load configuration
config = ConfigurationLoader.get_default_config()

# Load data
loader = HDF5CompoundLoader()
compounds = loader.load_all("data.h5")

# Build hierarchy
builder = HierarchyBuilder()
hierarchy = builder.build(compounds, hierarchy_mode="monomer")

# Detect peaks
detector = PeakDetector(config)
for compound in compounds:
    peaks = detector.detect(compound.chromatogram)
```

See `examples/analyze.py` for complete working implementation.

## Development

### Project Structure

```
LC-Seq-Official/
├── src/lcseq/              # Main package (~15,800 lines)
│   ├── domain/             # Pure domain logic (zero external dependencies)
│   │   ├── entities/       # Core entities (6 files: Compound, Peak, Chromatogram, etc.)
│   │   ├── value_objects/  # Value objects (RetentionTime, PeakBoundaries, etc.)
│   │   ├── models/         # Graph models (CompoundHierarchy, EquivalenceClass, etc.)
│   │   └── services/       # Domain services (23 files: peak detection, validation, etc.)
│   ├── application/        # Use cases and orchestration
│   │   ├── pipelines/      # Full analysis pipeline
│   │   ├── use_cases/      # Individual/pooled chromatogram processing
│   │   └── dtos/           # Data transfer objects
│   ├── infrastructure/     # I/O adapters
│   │   ├── parsers/        # CSV/Excel parsing
│   │   ├── exporters/      # CSV/JSON/Excel export
│   │   ├── loaders/        # HDF5 compound loading
│   │   ├── repositories/   # Data persistence
│   │   └── configuration/  # YAML config loading
│   └── presentation/       # External interfaces
│       ├── cli/            # Command-line interface
│       └── visualization/  # Plotting and visualization
├── configs/                # Configuration files
│   └── default.yaml        # Default configuration (Single Source of Truth)
├── docs/                   # Documentation
│   └── THEORY.md           # Mathematical foundations (2,270 lines)
└── examples/               # Working examples
    ├── analyze.py          # Lineage analysis script (fully functional)
    └── README.md           # Usage guide
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

2. **Implement the feature**

   - Add code to appropriate module in `src/lcseq/`
   - Follow type hints and docstring conventions
   - Keep functions focused and testable
   - Update [THEORY.md](docs/THEORY.md) if adding new algorithms

3. **Verify code quality**

   ```bash
   pre-commit run --all-files
   mypy src/
   ```

4. **Test with examples**

   ```bash
   python examples/analyze.py --reference "Phe-DNvl-DPhe"
   ```

5. **Create pull request**

## Architecture

The package follows **Clean Architecture** (Hexagonal Architecture) principles:

1. **Domain Layer**: Pure domain logic with zero external dependencies

   - Entities: `Compound`, `Chromatogram`, `Peak`, `BuildingBlock`, `PooledCompound`
   - Value Objects: `RetentionTime`, `PeakBoundaries`, `BuildingBlockSequence`, `MonomerSequence`
   - Models: `CompoundHierarchy`, `EquivalenceClass`, `PeakClassification`
   - Services: 23 domain services including peak detection, hierarchy building, validation

2. **Application Layer**: Use cases and orchestration

   - `ProcessChromatogramsUseCase`: Individual variant processing
   - `ProcessPooledChromatogramsUseCase`: Pooled variant processing
   - `FullAnalysisPipeline`: Complete analysis workflow

3. **Infrastructure Layer**: I/O adapters

   - HDF5 compound loading
   - YAML configuration
   - CSV/JSON/Excel export
   - Library and result repositories

4. **Presentation Layer**: External interfaces
   - CLI (command-line interface)
   - Visualization (plotting and chromatogram display)

**Key Principle**: Dependencies always flow inward (Presentation → Infrastructure → Application → Domain)

See [THEORY.md](docs/THEORY.md) for mathematical foundations and `examples/analyze.py` for architectural implementation.

## Implementation Status

**Current Status**: Core algorithms implemented with working examples

### Completed ✅

**Domain Layer** (Pure business logic, ~15,800 lines)

- ✅ **Entities**: Compound, Chromatogram, Peak, BuildingBlock, PooledCompound (6 files)
- ✅ **Value Objects**: RetentionTime, PeakBoundaries, BuildingBlockSequence, MonomerSequence
- ✅ **Models**: CompoundHierarchy, EquivalenceClass, PeakClassification
- ✅ **Services**: 23 domain services implementing core algorithms:
  - Peak detection (Discrete Morse Theory + Poisson statistics)
  - Hierarchy building (DAG construction for building-block and monomer modes)
  - Peak classification (constraint propagation through hierarchy)
  - Validation (Bayesian confidence assessment)
  - Pooled mode (positional variant aggregation with correlation validation)

**Application Layer** (Use cases and orchestration)

- ✅ `ProcessChromatogramsUseCase`: Individual variant analysis
- ✅ `ProcessPooledChromatogramsUseCase`: Hybrid pooling strategy
- ✅ `FullAnalysisPipeline`: Complete workflow orchestration

**Infrastructure Layer** (I/O adapters)

- ✅ HDF5 compound loading
- ✅ YAML configuration system
- ✅ CSV/JSON/Excel exporters
- ✅ Library and result repositories

**Presentation Layer** (External interfaces)

- ✅ CLI framework
- ✅ Visualization plotters (chromatogram, lineage, hierarchy)

**Working Examples**

- ✅ `examples/analyze.py`: Fully functional lineage analysis script
- ✅ Supports both hierarchy modes (building-block, monomer)
- ✅ Supports both variant modes (individual, pooled)
- ✅ Config-driven with CLI overrides

**Documentation**

- ✅ **[THEORY.md](docs/THEORY.md)**: 2,270 lines documenting mathematical foundations
  - Discrete Morse Theory for peak detection
  - DAG structure and poset theory
  - Bayesian validation framework
  - Pooled mode mathematics
- ✅ **[examples/README.md](examples/README.md)**: Complete usage guide
- ✅ Configuration examples in `configs/default.yaml`

### In Progress / Planned

**Testing Infrastructure**

- ⏳ Comprehensive test suite (unit, integration, property-based)
- ⏳ Test coverage reporting
- ⏳ CI/CD pipeline

**Additional Documentation**

- ⏳ Architecture documentation (detailed design specifications)
- ⏳ User guide (step-by-step tutorials)
- ⏳ API reference documentation

**Features**

- ⏳ Performance optimization for large libraries (>10,000 compounds)
- ⏳ Additional export formats
- ⏳ Interactive visualization tools

## Contributing

Contributions are welcome! Please:

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/my-feature`)
3. Implement your feature following Clean Architecture principles
4. Add type hints and docstrings
5. Run code quality checks (`pre-commit run --all-files`)
6. Update [THEORY.md](docs/THEORY.md) if adding new algorithms or mathematical concepts
7. Test with `examples/analyze.py`
8. Submit a pull request with clear description

See the **Adding New Features** section above for detailed workflow.

## Configuration

All analysis parameters are controlled by `configs/default.yaml` (Single Source of Truth). The config file contains:

**Analysis Settings:**

- `variant_mode`: `individual` or `pooled` (default: `individual`)
- `hierarchy_mode`: `monomer` or `building_block` (default: `monomer`)

**Peak Detection:**

- `min_persistence`: Minimum peak persistence (prominence threshold)
- `z_threshold`: Z-score threshold for statistical significance (Poisson test)
- `prominence_percentile`: Percentile for adaptive prominence filtering
- `min_snr`: Minimum signal-to-noise ratio

**Validation:**

- Purity thresholds (percentile-based, dataset-relative)
- SNR thresholds (detection and quantitation limits)
- Retention time precision parameters

**Pooling (for pooled mode):**

- `correlation_threshold`: Minimum correlation for valid pooling
- `aggregation_method`: `mean` or `median` for signal pooling

### How to Customize

1. **Copy the default config:**

   ```bash
   cp configs/default.yaml configs/my_config.yaml
   ```

2. **Edit your config file:**

   ```yaml
   analysis:
     variant_mode: pooled # Enable pooled mode by default
     hierarchy_mode: building_block # Use building-block hierarchy

   detection:
     min_persistence: 0.03 # Lower threshold = more sensitive
     z_threshold: 2.5 # Lower threshold = detect weaker peaks

   pooling:
     correlation_threshold: 0.75 # More lenient pooling threshold
   ```

3. **Run with your config:**
   ```bash
   python examples/analyze.py --reference "Phe-DNvl-DPhe" --config configs/my_config.yaml
   ```

**Note**: CLI arguments override config values. For example, `--variant-mode pooled` overrides the config file setting.

See `configs/default.yaml` for complete parameter documentation with comments.

---

## Documentation

- **[THEORY.md](docs/THEORY.md)**: Comprehensive mathematical foundations (2,270 lines)
  - Discrete Morse Theory for peak detection
  - DAG-based classification theory
  - Bayesian validation framework
  - Pooled mode mathematics
  - Domain vocabulary and terminology
- **[examples/README.md](examples/README.md)**: Complete usage guide
  - All 4 mode combinations explained
  - Data format requirements
  - Troubleshooting guide
- **[configs/default.yaml](configs/default.yaml)**: Complete configuration reference with inline documentation

## Citation

If you use this software in your research, please cite:

```
[Citation to be added]
```

## License

MIT License - see LICENSE file for details

## Contact

For questions or issues, please open an issue on GitHub.
