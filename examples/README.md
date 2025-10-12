# LC-Seq Examples

This directory contains the example script demonstrating config-driven lineage analysis.

## Key Principle

**The example script contains ZERO business logic.** It only:
- Loads data (infrastructure)
- Calls domain services (orchestration)
- Handles visualization (presentation)

All algorithms (peak detection, hierarchy building, pooling) are in tested domain services (`src/lcseq/domain/services/`).

## analyze.py - Unified Config-Driven Analysis

**Single script that handles all analysis modes based on configuration.**

### Quick Start

```bash
# Use default configuration (individual mode, monomer hierarchy)
python examples/analyze.py --reference "Phe-DNvl-DPhe"

# Override to pooled mode
python examples/analyze.py --reference "Phe-DNvl-DPhe" --variant-mode pooled

# Use building-block hierarchy
python examples/analyze.py --reference "Phe-DNvl-DPhe" --hierarchy-mode building_block

# Use custom configuration file
python examples/analyze.py --reference "Phe-DNvl-DPhe" --config my_config.yaml
```

---

## Configuration-Driven Design

The script automatically loads **`configs/default.yaml`** (Single Source of Truth) and uses:
- **`analysis.variant_mode`** - `individual` (default) or `pooled`
- **`analysis.hierarchy_mode`** - `building_block` or `monomer` (default)
- All detection, classification, pooling, and visualization parameters

CLI arguments override config values when provided.

### To Change Defaults

**Edit `configs/default.yaml`** - changes propagate everywhere automatically!

```yaml
analysis:
  variant_mode: individual  # Change to pooled for ~3-10× speedup
  hierarchy_mode: monomer   # or building_block

detection:
  min_persistence: 0.05  # Adjust peak detection sensitivity
  z_threshold: 3.0
  prominence_percentile: 0.5

pooling:
  correlation_threshold: 0.8  # Minimum correlation for valid pooling
  aggregation_method: mean    # or median
```

---

## Analysis Modes

### Individual Mode (default)
- Processes each compound variant separately
- Full detail per variant
- Standard workflow

### Pooled Mode
- Groups positional variants into equivalence classes
- Detects peaks on pooled signal (mean/median)
- Integrates areas on individual variants
- **~3-10× speedup** for peak detection
- Automatic fallback if correlation < 0.8
- See THEORY.md Section 4.2

**When to use pooled mode:**
- Large datasets with many positional variants
- Variants expected to have similar signals (high correlation)
- Need faster processing without sacrificing purity measurements

---

## Hierarchy Modes

### Monomer Mode (default, per THEORY.md)
- DAG with convergence patterns
- Based on chemical identity after monomer decomposition
- More connections, richer structure
- Recommended for most analyses

### Building Block Mode
- Poset (partial order) structure
- Based on positional sequences
- Simpler, tree-like structure
- Useful for understanding sequence relationships

---

## Outputs

The script generates:

**Individual Mode:**
```
results/
├── lineage_<sequence>.csv      # Peak counts per compound
└── plots/
    └── lineage_<sequence>.png  # Offset chromatogram plot
```

**Pooled Mode:**
```
results/
├── lineage_pooled_<sequence>.csv       # Peak counts (1 per equivalence class)
├── pooled_summary_<sequence>.csv       # Pooling statistics per class
└── plots/
    └── lineage_pooled_<sequence>.png   # Offset plot (1 trace per class)
```

---

## Command-Line Options

```
--reference TEXT         Reference compound sequence (required)
                        Example: "Phe-DNvl-DPhe"

--variant-mode CHOICE   Override analysis mode from config
                        individual: Process each variant separately
                        pooled: Aggregate variants (~3-10× speedup)

--hierarchy-mode CHOICE Override hierarchy mode from config
                        monomer: DAG with convergence (default)
                        building_block: Poset structure

--config PATH           Configuration file (default: configs/default.yaml)

--data PATH             HDF5 data file (default: test_data/processed_data.h5)

--output PATH           Output directory (default: results/)
```

---

## Lineage Analysis Workflow

Per THEORY.md Section 3.1, the workflow analyzes a **reference compound** and its complete **lineage**:

1. **User specifies a reference compound** (e.g., "Phe-DNvl-DPhe")
2. **System finds the lineage** (reference + all descendants)
3. **Builds hierarchy** (monomer DAG or building-block poset)
4. **Processes all members** (peak detection ± pooling)
5. **Generates offset chromatogram plot** (stacked vertically)
6. **Exports results** (CSV + plots)

### Terminology (per THEORY.md Section 3.1)

**Use these terms:**
- ✅ **Reference Compound** - The compound currently being analyzed
- ✅ **Lineage** - All ancestors + descendants + self (Principal Ideal ↓X)
- ✅ **Descendant** - Compound with fewer building blocks
- ✅ **Ancestor** - Compound with more building blocks

**Do NOT use:**
- ❌ "Parent" (ambiguous in combinatorial libraries)
- ❌ "Child" (relative, not absolute)

---

## Architecture

### Clean Architecture Compliance

```
Example Script (Presentation Layer)
    ↓ loads config from
configs/default.yaml (Infrastructure)
    ↓ calls
Domain Services (Domain Layer)
    ↓ operates on
Domain Entities (Domain Layer)
```

**What's in the Example Script:**
- ✅ Configuration loading (infrastructure)
- ✅ Data loading (infrastructure: HDF5)
- ✅ Orchestration (calling domain services in sequence)
- ✅ Presentation (plotting, visualization)

**What's NOT in the Example Script:**
- ❌ Domain Logic (peak detection, pooling algorithms)
- ❌ Business Rules (validation criteria, classification logic)
- ❌ Hardcoded parameters (everything comes from config)

All algorithms are in `src/lcseq/domain/services/` with comprehensive test coverage.

---

## Finding Available Sequences

If you're not sure what sequences exist in your data:

```bash
# This will fail but show first 20 sequences
python examples/analyze.py --reference "INVALID"
```

---

## Examples

### Basic Usage

```bash
# Default: individual mode, monomer hierarchy
python examples/analyze.py --reference "Leu-LA03-Pro-Leu-DLeuMe-DPro-Leu-Leu-DPro"
```

### Pooled Mode

```bash
# Enable pooled mode for faster processing
python examples/analyze.py \
  --reference "Leu-LA03-Pro-Leu-DLeuMe-DPro-Leu-Leu-DPro" \
  --variant-mode pooled
```

### Building-Block Hierarchy

```bash
# Use building-block hierarchy instead of monomer
python examples/analyze.py \
  --reference "Phe-DNvl-DPhe" \
  --hierarchy-mode building_block
```

### Custom Configuration

```bash
# Use custom config file (overrides defaults)
python examples/analyze.py \
  --reference "Phe-DNvl-DPhe" \
  --config my_custom_config.yaml
```

---

## Troubleshooting

### "Reference compound not found"

**Problem:** The sequence you specified doesn't exist in the data.

**Solutions:**
1. Check sequence format (should be "BB1-BB2-BB3")
2. Check available sequences (run with "INVALID" to see first 20)
3. Verify HDF5 file has data

### Small lineage size

**Problem:** Lineage only has 1-2 compounds.

**Reasons:**
- Descendants not in dataset
- Sequence has no natural truncations

### Low correlation warnings (pooled mode)

**Problem:** Many equivalence classes have low correlation (< 0.8).

**Solutions:**
- This is expected when variants have different signals
- System automatically handles this (no action needed)
- Consider using individual mode if most classes have low correlation

---

## Questions?

See:
- `docs/QUICKSTART.md` - Full user guide
- `docs/THEORY.md` - Mathematical foundations
  - Section 3.1: Lineage terminology
  - Section 4.2: Pooled mode
- `docs/ARCHITECTURE.md` - System design
- `configs/default.yaml` - Single Source of Truth for all parameters
