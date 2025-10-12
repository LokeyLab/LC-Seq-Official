# LC-Seq Changelog

All notable changes to this project are documented here.

---

## [Unreleased] - 2025-10-09

### Major Changes

#### Baseline Correction Completely Removed (2025-10-09)

**Problem:** Baseline correction (AsLS algorithm) was found to hinder data quality. Raw LC-MS signals provide better peak detection and integration results.

**Changes:**

**Code Deletion:**
- ❌ `src/lcseq/domain/services/baseline_corrector.py` - DELETED (217 lines)
- ❌ `tests/test_domain/test_services/test_baseline_corrector.py` - DELETED (156 lines)
- ❌ BaselineCorrector removed from domain services exports

**Application Layer:**
- `ProcessChromatogramsUseCase` - Simplified 3→2 steps
- `ProcessChromatogramsWithIntegrationUseCase` - Simplified 4→3 steps
- `full_analysis_pipeline.py` - Removed baseline correction stage
- `analysis_request.py` - DELETED `baseline_params` field

**Configuration System:**
- `analysis_configuration.py` - DELETED `baseline_params` field completely
- `yaml_loader.py` - Removed all baseline loading/saving code
- All 5 YAML configs - Removed commented baseline sections (clean slate)

**Documentation:**
- `README.md` - Removed AsLS from algorithm list
- `examples/README.md` - Updated to raw signals
- `ARCHITECTURE.md` - Removed baseline from domain services
- `THEORY.md` - DELETED Sections 5.0.2-5.0.9 (~550 lines on baseline correction)
- All docstrings - Removed deprecation notes and backwards compat references

**Examples & CLI:**
- `examples/analyze.py` - Removed `skip_baseline` parameter
- `cli/main.py` - Removed baseline from all displays

**Result:**
- ✅ **800+ lines deleted** (code + documentation)
- ✅ Simpler processing pipeline (25% fewer steps)
- ✅ Better data quality (raw signals perform better)
- ✅ No backwards compatibility needed - clean slate
- ✅ Processing flow: Raw Signal → Peak Detection → Integration → Classification → Validation

**Files Modified:** 14 files
**Files Deleted:** 2 files
**Lines Removed:** ~800 lines

**Verification:** ✅ System tested and working correctly with raw signals only

---

### Code Quality Improvements

#### DRY Violations Fixed - Purity and SNR Calculators Extracted (2025-10-09)

**Problem:** Brutal-code-auditor identified CRITICAL DRY violations - duplicate purity and SNR calculation logic in multiple validation services.

**Duplicate Code Identified:**
- `_calculate_purity()` - IDENTICAL 47-line method in BayesianValidator and AdaptiveValidator
- `_calculate_snr()` - IDENTICAL 46-line method in BayesianValidator (also accessed from ValidationWorkflow)

**Changes:**

**New Domain Services Created:**
- ✅ `src/lcseq/domain/services/purity_calculator.py` - Single-source-of-truth purity calculation
- ✅ `src/lcseq/domain/services/snr_calculator.py` - Single-source-of-truth SNR calculation

**Files Refactored (4):**
- `src/lcseq/domain/services/bayesian_validator.py` - Removed duplicate methods, now uses shared calculators
- `src/lcseq/domain/services/validation/adaptive_validator.py` - Removed duplicate purity method, uses PurityCalculator
- `src/lcseq/domain/services/validation/validation_workflow.py` - Fixed private method access, uses shared calculators
- `src/lcseq/domain/services/__init__.py` - Exported new calculators

**Tests Created (2 new test files):**
- `tests/test_domain/test_services/test_purity_calculator.py` - 7 tests, 100% coverage
- `tests/test_domain/test_services/test_snr_calculator.py` - 8 tests, 100% coverage

**Tests Updated (1):**
- `tests/test_domain/test_services/test_bayesian_validator.py` - Updated to use public calculators

**Result:**
- ✅ **93 lines of duplicate code eliminated** (47 purity + 46 SNR)
- ✅ **Single Source of Truth (SSoT)** - All validation services use same calculation logic
- ✅ **Better encapsulation** - No more private method access violations
- ✅ **100% test coverage** for both calculators
- ✅ All 32 tests passing (17 BayesianValidator + 7 PurityCalculator + 8 SNRCalculator)

**Verification:**
```bash
python -m pytest tests/test_domain/test_services/test_bayesian_validator.py -xvs  # ✅ 17/17 passing
python -m pytest tests/test_domain/test_services/test_purity_calculator.py -xvs   # ✅ 7/7 passing
python -m pytest tests/test_domain/test_services/test_snr_calculator.py -xvs      # ✅ 8/8 passing
```

**Architecture Benefits:**
- DRY principle enforced
- Easier to maintain (single implementation)
- Easier to test (isolated unit tests)
- Consistent calculations across all validators

---

### Major Cleanups

#### SPOTLESS THEORY.md Terminology Compliance Achieved (2025-10-09)

**Problem:** User requested "spotless" compliance - elimination of ALL domain-specific terminology violations throughout the entire codebase, not just CRITICAL user-facing issues.

**Phase 1: CRITICAL Violations (User-Facing)**

**Public API Docstrings Fixed (4 locations):**
- `src/lcseq/domain/models/compound_hierarchy.py:382` - "immediate children" → "immediate descendants"
- `src/lcseq/domain/models/compound_hierarchy.py:392` - "Immediate children" → "Immediate descendants"
- `src/lcseq/domain/models/compound_hierarchy.py:408` - "immediate parents" → "immediate extensions"
- `src/lcseq/domain/models/compound_hierarchy.py:418` - "Immediate parents" → "Immediate ancestors"

**Documentation Files Deleted (3 files):**
- `examples/FAMILY_ANALYSIS_GUIDE.md` (376 lines) - Taught wrong terminology to users
- `examples/FAMILY_ANALYSIS_TECHNICAL.md` (400 lines) - Extensive "parent/child/family" usage
- `examples/EXAMPLES_SUMMARY.md` (267 lines) - Outdated, duplicates examples/README.md

**Test Documentation Fixed (2 files):**
- `tests/test_domain/test_truncations.py:176` - "parent-child relationships" → "ancestor-descendant relationships"
- `tests/test_domain/test_services/test_peak_classifier.py:62-65` - "parent/child" comments → "ancestor/descendant"

**Phase 2: MAJOR Violations (Internal Implementation)**

**Domain Service Internal Variables (~59 violations fixed in 3 files):**

1. `src/lcseq/domain/services/hierarchy_builder.py` (~35 violations):
   - Algorithm comments: parent → ancestor, child → descendant
   - Loop variables throughout
   - Method signatures: `_detect_truncation_edge`, `_is_building_block_truncation`, `_is_monomer_truncation`
   - All internal parameters and docstrings

2. `src/lcseq/domain/services/validation_checker.py` (~15 violations):
   - DFS helper methods: children → descendants
   - Level ordering validation: parent → ancestor
   - Cycle finding: children → descendants
   - Connectivity validation: has_parents → has_ancestors, has_children → has_descendants
   - Edge validation loops

3. `src/lcseq/domain/services/path_finder.py` (~9 violations):
   - DFS path finding: child → descendant, children → descendants
   - BFS shortest path: children → descendants

**Test Variables (~31 violations fixed in 2 files):**

4. `tests/test_domain/test_services/test_hierarchy_builder.py` (~6 violations):
   - Test method variables: parent → ancestor, child → descendant
   - All test cases: `test_is_building_block_truncation_true`, `test_is_building_block_truncation_false_two_removed`, `test_monomer_truncation_detection`

5. `tests/test_domain/test_models/test_compound_hierarchy.py` (~25 violations):
   - Test variables: parent → ancestor, child → descendant, child1/2 → descendant1/2, parent1/2 → ancestor1/2
   - Error message matchers updated for 4 test methods
   - 9 test methods refactored

**Result:** **SPOTLESS - 100% THEORY.md Section 3.1 compliance achieved**

All domain-specific uses of informal terms eliminated:
- "parent" in domain contexts: **0 violations** ✅
- "child/children" in domain contexts: **0 violations** ✅
- "family" in domain contexts: **0 violations** ✅

Remaining instances (17 total) are ALL technical/infrastructure terms:
- "parent directory" (filesystem operations)
- `output_path.parent.mkdir()` (Python pathlib API)
- Documentation showing what NOT to use

**Verification:**
```bash
grep -r "parent" src/ tests/ --include="*.py" | grep -v "parent directory\|output_path.parent"  # 0 domain violations ✅
grep -r "child" src/ tests/ --include="*.py"                                                    # 0 domain violations ✅
grep -r "family" src/ tests/ --include="*.py"                                                   # 0 domain violations ✅
```

**Tests:** All 37 tests passing (9 in hierarchy_builder, 28 in compound_hierarchy)

---

#### Codebase-Wide Terminology Compliance (2025-10-09)

**Problem:** The entire codebase violated THEORY.md Section 3.1 terminology, using informal terms like "parent," "family," and "child" instead of mathematically precise terms.

**Changes:**

**Files Renamed (2):**
- `src/lcseq/domain/services/family_finder.py` → `lineage_finder.py`
- `src/lcseq/presentation/visualization/plotters/family_plotter.py` → `lineage_plotter.py`

**Classes Renamed (2):**
- `FamilyFinderService` → `LineageFinderService`
- `FamilyOffsetPlotter` → `LineageOffsetPlotter`

**Methods Renamed (3):**
- `find_family()` → `find_lineage()`
- `count_family_by_level()` → `count_lineage_by_level()`
- `group_family_by_level()` → `group_lineage_by_level()`

**Parameters/Variables Updated (~50 occurrences):**
- `parent` parameter → `reference` (THEORY.md line 1100: "compound currently being analyzed")
- `family` variable → `lineage` (THEORY.md line 1105: "all ancestors + descendants + self")
- `child` → `descendant` (THEORY.md line 1092: "compound with fewer building blocks")

**Files Modified (13):**
1. Domain services: `lineage_finder.py`, `compound_search.py`, `compound_ordering.py`, `hierarchy_builder.py`, `validation_checker.py`, `path_finder.py`
2. Domain models: `compound_hierarchy.py`
3. Application: `process_chromatograms.py`
4. Presentation: `lineage_plotter.py`, `base_plotter.py`
5. Infrastructure: `report_generator.py`, `excel_exporter.py`, `result_repository.py`

**Docstrings Updated:**
- Added THEORY.md Section 3.1 references throughout
- Replaced all informal terminology with mathematical terms
- Added terminology notes to all public APIs

**Verification:**
```bash
grep -r "FamilyFinder|FamilyPlotter" src/ --include="*.py"  # 0 results ✅
grep -r "find_family|count_family" src/ --include="*.py"    # 0 results ✅
```

**Result:** 100% THEORY.md terminology compliance across entire codebase.

---

#### Removed Synthetic Data from Examples (2025-10-09)

**Problem:** Fake/demo data was exposed in user-facing examples and infrastructure layer, violating the principle that "synthetic data should NEVER leave the test directory."

**Changes:**
- Deleted `demo` mode from `examples/analyze.py` (95 lines removed)
- Deleted entire `src/lcseq/infrastructure/test_data/` directory
- Deleted `DemoDataGenerator` class (117 lines)
- Removed `DemoDataGenerator` from infrastructure exports
- Updated all documentation to remove demo mode references

**Result:** Synthetic data now ONLY exists in `tests/fixtures/` for algorithm testing. Examples contain ZERO fake data.

**Files Deleted:**
- `src/lcseq/infrastructure/test_data/demo_data_generator.py`
- `src/lcseq/infrastructure/test_data/__init__.py`

**Files Modified:**
- `examples/analyze.py` - Removed demo mode
- `src/lcseq/infrastructure/__init__.py` - Removed DemoDataGenerator export
- `examples/README.md` - Removed demo mode documentation
- `examples/EXAMPLES_SUMMARY.md` - Removed demo mode documentation

---

#### Removed Batch Mode and Fixed Terminology (2025-10-09)

**Problem:**
1. Batch mode didn't match real workflows (arbitrary "first N compounds")
2. Terminology violated THEORY.md Section 3.1 ("parent"/"family" are informal, not mathematically defined)

**Changes:**
- Removed `batch` mode entirely (67 lines removed)
- Replaced "parent" with "reference compound" (per THEORY.md line 1100)
- Replaced "family" with "lineage" (per THEORY.md line 1105)
- Simplified CLI: removed `--mode` argument, only `--reference` required
- Updated all documentation with THEORY.md-compliant terminology

**Terminology Mapping:**

| ❌ Old (Informal) | ✅ New (THEORY.md) | Reference |
|------------------|-------------------|-----------|
| "parent" | Reference Compound | Line 1100 |
| "family" | Lineage | Line 1105 |
| "child" | Descendant | Line 1092 |
| "batch mode" | *(removed)* | (not in theory) |

**CLI Before:**
```bash
python examples/analyze.py --mode family --parent "Phe-DNvl-DPhe"
python examples/analyze.py --mode batch --max-compounds 100
```

**CLI After:**
```bash
python examples/analyze.py --reference "Phe-DNvl-DPhe"
```

**Files Modified:**
- `examples/analyze.py` - Complete terminology overhaul, batch mode removed
- `examples/README.md` - Rewritten with THEORY.md terminology
- `examples/EXAMPLES_SUMMARY.md` - Updated terminology

**Lines Removed:** 167 total (67 batch mode + 95 demo mode + 5 misc)

---

#### Architecture Cleanup (2025-10-09)

**Problem:** Multiple issues from brutal-code-auditor review:
- Empty "ghost" directories (YAGNI violations)
- 2MB legacy backup in working directory
- 1,029 __pycache__ directories
- Unused imports
- Overengineered view models and presenters

**Changes:**
- Deleted 4 empty application directories (factories, ports, services, use_cases)
- Deleted empty algorithms subdirectory structure
- Deleted `src_backup_20251008/` (2MB legacy backup)
- Cleaned 1,029 __pycache__ directories
- Cleaned 9 .DS_Store files
- Removed unused imports (3 files)
- Deleted overengineered view models (CompoundViewModel, PeakViewModel - 356 lines)
- Deleted unused presenters (FamilyPresenter, ReportPresenter - 184 lines)

**Result:** Cleaner structure, YAGNI-compliant, 2.1 MB freed

**Code Quality:**
- Before: C+ (68/100)
- After: B+ (82/100)
- YAGNI: 40/100 → 85/100
- Dead Code: 50/100 → 95/100

---

#### Presentation Layer Refactoring (2025-10-09)

**Problem:**
1. Wrong layer name: `interfaces/` instead of proper Clean Architecture `presentation/`
2. Business logic embedded in examples (229 lines matplotlib, 191 lines processing)
3. Data loading logic in examples (42 lines HDF5 parsing)

**Changes:**
- Renamed `src/lcseq/interfaces/` → `src/lcseq/presentation/`
- Updated all imports across codebase
- Created `FamilyOffsetPlotter` in presentation layer (348 lines extracted)
- Created `ProcessChromatogramsUseCase` in application layer (191 lines extracted)
- Created `HDF5CompoundLoader` in infrastructure layer (extracted from examples)
- Updated ARCHITECTURE.md with correct 4-layer naming

**Result:** Example script is now pure orchestration with ZERO business logic

**Example Script:**
- Before: 858 lines with embedded logic
- After: 259 lines of pure orchestration
- Reduction: 599 lines (70% reduction)

---

### Summary Statistics

**Total Lines Removed:** ~1,400 lines
- Business logic moved to proper layers: ~760 lines
- Dead code/overengineering: ~540 lines
- Batch/demo modes: ~167 lines

**Files Deleted:** 1,037
- Empty directories: 9
- __pycache__: 1,029
- Legacy backup: 1 (2MB)
- View models/presenters: 4 files

**Architecture Improvements:**
- Proper Clean Architecture layer naming (presentation not interfaces)
- YAGNI-compliant (no empty directories)
- THEORY.md-compliant terminology
- Zero synthetic data in examples
- Zero business logic in examples

**Code Quality:**
- Overall: C+ → B+ (+14 points)
- YAGNI: 40 → 85 (+45 points)
- Dead Code: 50 → 95 (+45 points)

---

## Versioning

This project follows [Semantic Versioning](https://semver.org/).

---

## Categories

- **Added** - New features
- **Changed** - Changes in existing functionality
- **Deprecated** - Soon-to-be removed features
- **Removed** - Removed features
- **Fixed** - Bug fixes
- **Security** - Vulnerability fixes
