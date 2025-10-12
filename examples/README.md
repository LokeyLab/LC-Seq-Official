# LC-Seq Examples

This directory contains a single example script demonstrating lineage analysis.

## Key Principle

**The example script contains ZERO business logic.** It only:
- Loads data (infrastructure)
- Calls domain services (orchestration)
- Handles visualization (presentation)

All algorithms (peak detection, hierarchy building) are in tested domain services (`src/lcseq/domain/services/`).

## Single Script

**`analyze.py`** - Lineage analysis for reference compounds

```bash
# Analyze reference compound and its descendants
python examples/analyze.py --reference "Phe-DNvl-DPhe"
```

---

## Lineage Analysis

**THE canonical workflow for LC-Seq data analysis.**

Per THEORY.md Section 3.1, the workflow analyzes a **reference compound** and its complete **lineage** (all descendants):

1. **User specifies a reference compound** (e.g., "Phe-DNvl-DPhe")
2. **System finds the lineage** (reference + all descendants)
3. **Processes all members** (peak detection on raw signals)
4. **Generates offset chromatogram plot** (stacked vertically)

### Terminology (per THEORY.md Section 3.1)

**Use these terms:**
- ✅ **Reference Compound** - The compound currently being analyzed
- ✅ **Lineage** - All ancestors + descendants + self
- ✅ **Descendant** - Compound with fewer building blocks
- ✅ **Ancestor** - Compound with more building blocks

**Do NOT use:**
- ❌ "Parent" (ambiguous in combinatorial libraries)
- ❌ "Child" (relative, not absolute)
- ❌ "Family" (not defined in THEORY.md)

### Usage Examples

```bash
# Basic usage (monomer hierarchy mode - default per THEORY.md)
python examples/analyze.py --reference "Phe-DNvl-DPhe"

# Use building-block hierarchy (poset structure)
python examples/analyze.py --reference "Phe-DNvl-DPhe" --hierarchy-mode building_block

# Custom data file
python examples/analyze.py --reference "Leu-Pro-Ala" --data my_data.h5
```

### Features

**Analysis:**
- ✅ Finds all descendants of reference compound
- ✅ Builds hierarchy (monomer DAG or building-block poset)
- ✅ Processes chromatograms with peak detection
- ✅ Generates offset chromatogram plots

**Architecture:**
- ❌ NO business logic in script
- ✅ Only orchestrates domain services (PeakDetector, HierarchyBuilder, LineageFinderService)
- ✅ Visualization is presentation concern (not domain)

### Output

```
examples/lineage_results/
├── lineage_<sequence>.csv      # Peak counts per compound
└── plots/
    └── lineage_<sequence>.png  # Offset chromatogram plot
```

### Command-Line Options

```
--reference TEXT         Reference compound sequence (required)
                        Example: "Phe-DNvl-DPhe"

--hierarchy-mode CHOICE  Hierarchy mode: monomer (default) or building_block
                        monomer: DAG with convergence (per THEORY.md)
                        building_block: Poset structure

--data PATH             HDF5 data file (default: test_data/processed_data.h5)

--output PATH           Output directory (default: examples/lineage_results)
```

---

## Finding Available Sequences

If you're not sure what sequences exist in your data:

```bash
# This will fail but show first 20 sequences
python examples/analyze.py --reference "INVALID"
```

---

## Architecture

### Clean Architecture Compliance

```
Example Script (Presentation Layer)
    ↓ calls
Domain Services (Domain Layer)
    ↓ operates on
Domain Entities (Domain Layer)
```

**What's in the Example Script:**
- ✅ Infrastructure: Loading HDF5 data
- ✅ Orchestration: Calling domain services in sequence
- ✅ Presentation: Plotting, visualization

**What's NOT in the Example Script:**
- ❌ Domain Logic: Peak detection algorithms
- ❌ Business Rules: Validation criteria, classification logic

All algorithms are in `src/lcseq/domain/services/` with 607 passing tests.

---

## Customization

### Change Detection Threshold

Edit `analyze.py` and search for `persistence_threshold`:

```python
persistence_threshold=0.05  # → 0.10 for stricter detection
```

---

## Typical Workflow

1. **Choose a reference sequence** you want to analyze
2. **Run lineage analysis:**
   ```bash
   python examples/analyze.py --reference "YourSequence"
   ```
3. **Check the plots** in `examples/lineage_results/plots/`
4. **Review the CSV** for peak counts
5. **Adjust if needed:**
   - Try different hierarchy mode: `--hierarchy-mode building_block`
   - Try different reference: `--reference "NewSequence"`

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

---

## Questions?

See:
- `docs/QUICKSTART.md` - Full user guide
- `docs/THEORY.md` - Mathematical foundations (Section 3.1 for terminology)
- `docs/ARCHITECTURE.md` - System design
