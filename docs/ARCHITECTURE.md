# LC-Seq Architecture

**Date**: 2025-10-08
**Status**: Active
**Related Documents**:

- [THEORY.md](THEORY.md) - Mathematical foundations and domain vocabulary

---

## Overview

LC-Seq follows **Clean Architecture** (also known as **Hexagonal Architecture** or **Ports and Adapters**) principles to maintain clear separation of concerns, enable testability, and support long-term maintainability.

**Core Principle**: Dependencies always flow inward toward the domain.

```
┌─────────────────────────────────────────────┐
│        Presentation Layer (CLI, Viz)        │  ← External adapters
├─────────────────────────────────────────────┤
│       Infrastructure (Files, Storage)       │  ← I/O implementations
├─────────────────────────────────────────────┤
│     Application (Use Cases, Orchestration)  │  ← Business workflows
├─────────────────────────────────────────────┤
│   Domain (Entities, Services, Algorithms)   │  ← Pure business logic
└─────────────────────────────────────────────┘
        ↑         ↑         ↑         ↑
        All dependencies flow INWARD
```

---

## Architecture Layers

### Domain Layer (Core)

**Purpose**: Pure domain logic with zero external dependencies

**Contains**:

- **Entities**: `Compound`, `Chromatogram`, `Peak`, `BuildingBlock`
- **Value Objects**: `RetentionTime`, `PeakBoundaries`, `MonomerSequence`
- **Models**: `CompoundHierarchy` (DAG/poset), `EquivalenceClass`, `AnalysisConfiguration`
- **Domain Services**:
  - Algorithms: Peak detection, signal processing, classification
  - Graph operations: Hierarchy building, truncation analysis, compound ordering
  - Validation: Purity calculation, threshold computation, Bayesian validation

**Key Characteristics**:

- No dependencies on outer layers
- No I/O operations (no file access, no database calls)
- Pure functions and immutable data where possible
- 100% testable without mocks
- Technology-agnostic

**Why Algorithms Are Domain Services**:

Algorithms encode **chromatography domain knowledge**:

- What constitutes a valid peak (Morse theory, persistence)
- Signal processing parameters (smoothing, derivatives)
- Classification rules (retention time ordering, DAG constraints)

Using scipy/numpy is like using `math.sqrt()` - they're computational tools, not infrastructure. **Infrastructure is for crossing I/O boundaries** (files, databases, networks), not for domain logic.

See [THEORY.md Part 5-7](THEORY.md#part-5-peak-detection-mathematical-foundations) for mathematical foundations.

**Directory Structure**:

```
domain/
├── entities/
│   ├── compound.py           # Sequence, level, relationships
│   ├── chromatogram.py       # Time points, signal variants
│   ├── peak.py               # Position, boundaries, classification
│   └── building_block.py     # Cycle, code, monomer decomposition
│
├── value_objects/
│   ├── building_block_sequence.py
│   ├── monomer_sequence.py
│   ├── retention_time.py
│   └── peak_boundaries.py
│
├── models/
│   ├── compound_hierarchy.py      # DAG/poset (renamed from TruncationHierarchy)
│   ├── equivalence_class.py       # Positional variants (renamed from CanonicalGroup)
│   ├── peak_classification.py
│   ├── peak_type_classification.py
│   ├── analysis_configuration.py
│   └── analysis_result.py
│
└── services/
    ├── peak_detector.py                    # Morse theory peak detection
    ├── peak_classifier.py                  # Peak type classification
    ├── peak_integrator.py                  # Peak area calculation
    ├── hierarchy_builder.py                # DAG/poset construction
    ├── lineage_finder.py                   # Find all descendants
    ├── compound_search.py                  # Search and filtering
    ├── compound_ordering.py                # Hierarchical clustering
    ├── validation_checker.py               # DAG validation
    ├── path_finder.py                      # Graph path algorithms
    ├── compound_similarity.py              # Similarity metrics
    ├── sequence_similarity.py              # Levenshtein distance
    ├── level_analyzer.py                   # Level-based analysis
    ├── equivalence_class_builder.py        # Positional variant grouping
    └── validation/
        ├── adaptive_validator.py
        ├── bayesian_validator.py
        ├── consensus_validator.py
        └── validation_workflow.py
```

---

### Application Layer

**Purpose**: Orchestrate domain services to implement business workflows

**Contains**:

- **Use Cases**: High-level business operations
- **Ports**: Abstract interfaces for I/O boundaries ONLY
- **Services**: Coordinate domain services that need repository access
- **Factories**: Simple dependency injection

**Key Characteristics**:

- Depends on domain (inward)
- Defines ports (abstractions) for infrastructure
- No direct I/O operations (delegates to infrastructure via ports)
- Orchestrates domain services

**Ports Are Only for I/O Boundaries**:

✅ **DO create ports for**:

- File loading (`DataLoaderPort`)
- Storage access (`CompoundRepositoryPort`)
- Configuration loading (`ConfigurationPort`)
- Results export (`ExportPort`)

❌ **DO NOT create ports for**:

- Peak detection (domain service - call directly)
- Signal processing (domain service - call directly)
- Classification logic (domain service - call directly)
- Validation services (domain service - call directly)

**Why?** Ports are for crossing I/O boundaries (architecture boundaries with side effects). Domain services are pure logic within the same boundary.

**Directory Structure**:

```
application/
├── use_cases/
│   ├── analyze_compound_family.py
│   ├── build_hierarchy.py
│   ├── detect_peaks.py
│   ├── classify_peaks.py
│   ├── validate_synthesis.py
│   └── generate_report.py
│
├── ports/                          # I/O boundaries ONLY
│   ├── compound_repository_port.py
│   ├── data_loader_port.py
│   ├── configuration_port.py
│   └── export_port.py
│
├── services/
│   ├── family_service.py           # Needs repository (crosses I/O boundary)
│   ├── batch_processor.py
│   └── progress_tracker.py
│
└── factories/
    └── use_case_factory.py
```

---

### Infrastructure Layer

**Purpose**: Implement adapters for external systems (files, databases, networks)

**Contains**:

- Persistence implementations (HDF5, in-memory)
- Configuration loaders (YAML)
- Export implementations (CSV, JSON, HDF5)

**Key Characteristics**:

- Implements application ports
- Handles all I/O operations
- No domain logic
- Technology-specific implementations

**What Does NOT Belong in Infrastructure**:

❌ Peak detection algorithms → Domain service
❌ Signal processing → Domain service
❌ Compound ordering (hierarchical clustering) → Domain service
❌ Validation logic → Domain service

**Why?** These encode domain knowledge, not I/O operations.

**Directory Structure**:

```
infrastructure/
├── persistence/
│   ├── hdf5_compound_repository.py    # Implements CompoundRepositoryPort
│   ├── hdf5_loader.py
│   └── in_memory_repository.py
│
├── configuration/
│   ├── yaml_config_loader.py          # Implements ConfigurationPort
│   └── config_validator.py
│
└── export/
    ├── csv_exporter.py                # Implements ExportPort
    ├── json_exporter.py
    └── hdf5_exporter.py
```

---

### Presentation Layer

**Purpose**: Adapt external interfaces to application layer

**Contains**:

- CLI (Typer-based)
- Visualization adapters (view models, presenters, plotters)

**Key Characteristics**:

- Depends on application layer
- Transforms external requests → use case calls
- Transforms domain results → external formats
- No business logic

**Naming Note**:
This layer is called "Presentation" to avoid confusion with the term "interfaces"
in the abstract sense (protocols/ports). In Clean Architecture, "interfaces" refers
to abstract boundaries between layers (e.g., repository ports), NOT to concrete
adapters like CLI or visualization components.

**Directory Structure**:

```
presentation/
├── cli/
│   ├── main.py
│   ├── commands/
│   └── formatters/
│
└── visualization/
    ├── adapters/
    │   ├── compound_view_model.py     # Domain → Presentation
    │   └── peak_view_model.py
    │
    ├── presenters/
    │   ├── family_presenter.py
    │   └── chromatogram_presenter.py
    │
    └── plotters/
        ├── base_plotter.py            # Template method
        └── signal_plotter.py          # Concrete implementations
```

---

## Configuration Strategy

### YAML-Based Configuration

All configuration externalized to YAML files (no hardcoded parameters).

**Configuration Profiles**:

```yaml
# configs/standard.yaml
analysis:
  variant_mode: individual # or "consensus"
  hierarchy_mode: block # or "monomer"

detection:
  min_persistence: 0.05
  boundary_method: valley_or_5pct

validation:
  purity_threshold: auto # Use P₂₅ from dataset distribution
  snr_threshold: auto # Use P₅₀ from dataset distribution
```

**Key Principle**: No magic numbers. All thresholds derived from data (THEORY.md Section 6.2).

---

## Dependency Injection

### Simple Composition Root

**No elaborate DI framework** - simple factory functions:

```python
# application/factories/use_case_factory.py

class AnalysisCompositionRoot:
    def __init__(self, config: AnalysisConfiguration):
        self.config = config

    def build_analysis_workflow(self) -> AnalyzeCompoundFamilyUseCase:
        # Infrastructure (implements ports)
        repository = HDF5CompoundRepository(self.config.data_path)

        # Domain services (call directly - no ports!)
        hierarchy_builder = HierarchyBuilder()
        peak_detector = MorseDetector()
        classifier = PeakTypeClassifier()
        validator = ValidationClassifier()

        # Use case (orchestrates domain services)
        return AnalyzeCompoundFamilyUseCase(
            repository=repository,           # Port implementation
            hierarchy_builder=hierarchy_builder,  # Domain service
            peak_detector=peak_detector,          # Domain service
            classifier=classifier,                # Domain service
            validator=validator,                  # Domain service
            config=self.config
        )
```

---

## Domain Terminology

**Critical**: Use mathematical and ancestry-based terminology to avoid ambiguity.

See [THEORY.md Part 8: Domain Vocabulary](THEORY.md#part-8-domain-vocabulary-ubiquitous-language) for complete definitions.

**Key Concepts**:

| ❌ Avoid | ✅ Use Instead        | Definition                                      |
| -------- | --------------------- | ----------------------------------------------- |
| Parent   | Maximal compound      | Longest compound in dataset                     |
| Parent   | Reference compound    | The compound currently being analyzed           |
| Child    | Descendant            | Compound with fewer building blocks             |
| -        | Ancestor              | Compound with more building blocks              |
| -        | Lineage               | All ancestors + descendants + self              |
| -        | Principal Ideal (↓c)  | All descendants of compound c                   |
| -        | Principal Filter (↑c) | All ancestors of compound c                     |
| -        | DAG                   | Directed Acyclic Graph (mathematical structure) |
| -        | Poset                 | Partially Ordered Set (formal model)            |

**Sequence Conventions**:

- N→C order: `Leu-Ala-Val-Pro` (N-terminus to C-terminus)
- Right-align for visualization (C-terminus is synthesis anchor)
- Position 0 = C-terminus (rightmost, synthesized first)
- Position numbering increases toward N-terminus

See [THEORY.md Section 1.5.1](THEORY.md#151-peptide-sequence-convention) for details.

---

## Architectural Decisions

### Why Clean Architecture?

**Problems with legacy system** (see [backup/AUDIT.md](backup/AUDIT.md)):

- Domain → infrastructure dependencies (architectural violation)
- Circular dependencies
- 956-line visualization module with massive duplication
- Hardcoded configuration (320-line factory)
- Unclear layer boundaries

**Benefits of Clean Architecture**:

- ✅ Domain layer 100% testable without mocks
- ✅ Clear dependency flow (always inward)
- ✅ Easy to swap implementations (HDF5 → Parquet, CLI → Web API)
- ✅ No circular dependencies
- ✅ SOLID principles enforced by architecture

### Two-Track Implementation Strategy

**Track A**: Clean implementation in `src_refactor/`

- Build correct architecture from scratch
- Follow THEORY.md specifications exactly
- No compromises with legacy patterns

**Track B**: Legacy preservation in `src/`

- Keep existing code functional during transition
- Reference for validation (results must match)

**Why not refactor in place?**

- Audit identified 40-60% technical debt
- Circular dependencies difficult to untangle incrementally
- Cleaner to build correctly from start

**Migration Path**:

1. Implement complete system in `src_refactor/`
2. Validate results match legacy system
3. Tag `src/` as `v0.1.0-legacy`
4. Move `src_refactor/` → `src/`
5. Archive legacy to `archive/src_legacy/`

### Consensus Mode Design

**Key Insight**: Consensus mode is an **optional optimization**, not required.

**Implementation** (THEORY.md Section 4.2):

- Detect peaks on consensus signal (expensive, do once)
- Integrate areas on individual variants (cheap, do per-variant)
- Automatic correlation check: min(r) > 0.8
- Automatic fallback to individual mode if check fails

**Status Flags**:

- `CONSENSUS_VALID`: Correlation check passed
- `HETEROGENEOUS`: Correlation check failed, used individual mode
- `CONSENSUS_INVALID_BUT_SIMILAR`: Fallback successful

See [THEORY.md Section 4.2.8.1](THEORY.md#4281-operational-fallback-workflow) for detailed workflow.

---

## Testing Strategy

### Test Pyramid

**Unit Tests (60%)**:

- Domain entities and value objects
- Domain services (algorithms, validation)
- Pure functions with no mocks

**Integration Tests (30%)**:

- Use cases orchestrating domain services
- Repository implementations
- Configuration loading

**End-to-End Tests (10%)**:

- Full analysis pipeline
- CLI commands
- Results validation

### Property-Based Testing

Use Hypothesis to verify mathematical invariants:

```python
@given(library_design=st....)
def test_dag_is_acyclic(library_design):
    """DAG has no cycles - THEORY.md Section 1.4"""
    hierarchy = build_hierarchy(library_design)
    assert not has_cycles(hierarchy)

@given(compound=st....)
def test_descendant_levels_decrease(compound):
    """All descendants have lower level - THEORY.md Section 3.3"""
    for desc in compound.descendants:
        assert desc.level < compound.level
```

---

## Performance Optimizations

From [THEORY.md Part 7: Mathematical Optimizations](THEORY.md#part-7-mathematical-optimizations):

**Graph Algorithms**:

- Topological sort for bottom-up traversal: O(V + E)
- Transitive reduction: Store only direct descendants
- Memoize descendant sets: O(V) → O(1) for cached queries

**Data Structures**:

- Adjacency list: Fast neighbor lookup
- Level-order index: Fast same-level queries
- Union-Find for connected components: O(α(n)) ≈ O(1)

**Target Performance**:

- Process 1,000-compound library in < 5 minutes
- Memory efficient (streaming where possible)

---

## Code Quality Standards

**File Size Limits**:

- Maximum: 400 lines per file
- Average: ~150 lines per file
- Extract to separate files if exceeded

**Dependency Rules**:

- Domain: zero external dependencies
- Application: depends only on domain
- Infrastructure: implements application ports
- Interface: depends on application

**Documentation**:

- NumPy-style docstrings for all public APIs
- Reference THEORY.md sections in docstrings
- Type hints for all function signatures

**Static Analysis**:

- mypy strict mode (no type: ignore)
- ruff for linting
- black for formatting
- pre-commit hooks enforce all checks

---

## Success Metrics

**Architecture Metrics**:

- ✅ Zero circular dependencies
- ✅ Domain layer has zero external dependencies
- ✅ All dependencies flow inward
- ✅ Ports defined for I/O boundaries only

**Code Quality Metrics**:

- ✅ No file over 400 lines
- ✅ Test coverage >90%
- ✅ All SOLID principles followed
- ✅ Domain testable without mocks

**Performance Metrics**:

- ✅ Total codebase: ~8,600 lines (down from 13,500)
- ✅ 36% line reduction
- ✅ Results match legacy system
- ✅ Performance: < 5 min for 1,000-compound library

---

## References

- **[THEORY.md](THEORY.md)** - Complete theoretical foundations (3,772 lines)
- **Clean Architecture** - Robert C. Martin
- **Hexagonal Architecture** - Alistair Cockburn

---

**Last Updated**: 2025-10-08
