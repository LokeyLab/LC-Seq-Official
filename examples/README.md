# LC-Seq Examples

This directory contains the example script demonstrating config-driven lineage analysis.

## Data Requirements

**IMPORTANT**: The example script is configured to work with the included test dataset (`test_data/processed_data.h5`). To use your own data, you must:

1. Format your data as HDF5 with the same structure as the test data
2. Include chromatogram signals (time points and scaled counts)
3. Include compound metadata (building block sequences, positional assignments)
4. Update the `--data` parameter to point to your file

See the HDF5 loader implementation in `src/lcseq/infrastructure/loaders/hdf5_compound_loader.py` for the expected format.

## Key Principle

**The example script contains ZERO business logic.** It only:

- Loads data (infrastructure)
- Calls domain services (orchestration)
- Handles visualization (presentation)

All algorithms (peak detection, hierarchy building, pooling) are in domain services (`src/lcseq/domain/services/`).

## analyze.py - Unified Config-Driven Analysis

**Single script that handles all analysis modes based on configuration.**

### Quick Start

```bash
# monomer hierarchy, individual mode
python examples/analyze.py --reference "Phe-DNvl-DPhe" --hierarchy-mode monomer --variant-mode individual

# building-block hierarchy, individual mode
python examples/analyze.py --reference "Phe-DNvl-DPhe" --hierarchy-mode building_block --variant-mode individual

# monomer hierarchy, pooled mode
python examples/analyze.py --reference "Phe-DNvl-DPhe" --variant-mode pooled --hierarchy-mode monomer

# building-block hierarchy, pooled mode
python examples/analyze.py --reference "Phe-DNvl-DPhe" --hierarchy-mode building_block --variant-mode individual


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

### How to Customize Parameters

**Option 1: Edit the default config directly** (affects all runs)

```bash
# Edit configs/default.yaml
vim configs/default.yaml
```

**Option 2: Create a custom config** (recommended for experiments)

```bash
# Copy default config
cp configs/default.yaml configs/my_experiment.yaml

# Edit your copy
vim configs/my_experiment.yaml

# Run with custom config
python examples/analyze.py --reference "Phe-DNvl-DPhe" --config configs/my_experiment.yaml
```

### Common Parameter Adjustments

**Make peak detection more sensitive:**

```yaml
detection:
  min_persistence: 0.03 # Lower = more sensitive (default: 0.05)
  z_threshold: 2.5 # Lower = detect weaker peaks (default: 3.0)
  prominence_percentile: 0.3 # Lower = retain more peaks (default: 0.5)
```

**Change pooling behavior:**

```yaml
pooling:
  correlation_threshold: 0.75 # More lenient (default: 0.8)
  aggregation_method: median # More robust to outliers (default: mean)
```

**Adjust validation stringency:**

```yaml
validation:
  purity_threshold: lenient # Options: auto, strict, lenient, or float (default: auto)
  snr_threshold: lenient # Options: auto, strict, lenient, or float (default: auto)
```

**Change default analysis modes:**

```yaml
analysis:
  variant_mode: pooled # Enable pooled mode by default (default: individual)
  hierarchy_mode: building_block # Use building-block hierarchy (default: monomer)
```

See `configs/default.yaml` for complete parameter documentation with inline comments explaining each setting.

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

- DAG with convergence at block granularity
- Based on block support sequences (position-independent)
- Positional variants with same blocks converge
- Useful for analyzing synthesis at block level

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

**Standard Terms:**

- **Reference Compound** - The compound currently being analyzed
- **Lineage** - All ancestors + descendants + self (Principal Ideal ↓X)
- **Descendant** - Compound with fewer building blocks
- **Ancestor** - Compound with more building blocks

**Note**: Terms like "parent" or "child" are ambiguous in combinatorial libraries.

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

All algorithms are in `src/lcseq/domain/services/` with full implementation.

---

## Finding Available Sequences

If you're not sure what sequences exist in your data:

```bash
# This will fail but show first 20 sequences
python examples/analyze.py --reference "INVALID"
```

---

## Examples - All 4 Mode Combinations

The script supports 4 analysis modes (2 variant modes × 2 hierarchy modes):

### 1. Individual + Monomer (Default)

```bash
# Process each variant separately with monomer-level hierarchy
python examples/analyze.py --reference "Phe-DNvl-DPhe"
# Equivalent to:
python examples/analyze.py --reference "Phe-DNvl-DPhe" --variant-mode individual --hierarchy-mode monomer
```

**Output**: `results/plots/lineage_monomer_Phe-DNvl-DPhe.png`

### 2. Individual + Building-Block

```bash
# Process each variant separately with building-block hierarchy
python examples/analyze.py --reference "Phe-DNvl-DPhe" --hierarchy-mode building_block
```

**Output**: `results/plots/lineage_block_Phe-DNvl-DPhe.png`

### 3. Pooled + Monomer

```bash
# Aggregate positional variants with monomer-level hierarchy
python examples/analyze.py --reference "Phe-DNvl-DPhe" --variant-mode pooled
```

**Output**: `results/plots/lineage_pooled_monomer_Phe-DNvl-DPhe.png`

### 4. Pooled + Building-Block

```bash
# Aggregate positional variants with building-block hierarchy
python examples/analyze.py --reference "Phe-DNvl-DPhe" --variant-mode pooled --hierarchy-mode building_block
```

**Output**: `results/plots/lineage_pooled_block_Phe-DNvl-DPhe.png`

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

- `docs/THEORY.md` - Mathematical foundations (2,270 lines)
  - Section 3.1: Lineage terminology
  - Section 4.2: Pooled mode analysis
  - Section 5: Peak detection (Discrete Morse Theory)
  - Section 6: Synthesis validation (Bayesian framework)
- `configs/default.yaml` - Single Source of Truth for all parameters
- Main `README.md` - Project overview and architecture
