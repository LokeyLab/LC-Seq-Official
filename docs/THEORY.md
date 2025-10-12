# LC-Seq Theoretical Foundations

**Last Updated**: 2025-10-09
**Purpose**: Foundational concepts, mathematical structures, and domain vocabulary for the LC-Seq analysis system

---

## Table of Contents

**I. FOUNDATIONS**

1. [Mathematical Structure](#part-1-mathematical-structure)
2. [Domain Foundations](#part-2-domain-foundations)
3. [Hierarchical Relationships](#part-3-hierarchical-relationships)
4. [Graph Properties and Patterns](#part-4-graph-properties-and-patterns)

**II. ANALYSIS METHODS** 5. [Peak Detection Mathematical Foundations](#part-5-peak-detection-mathematical-foundations) 6. [Synthesis Validation Theory](#part-6-synthesis-validation-theory)

**III. COMPUTATIONAL IMPLEMENTATION** 7. [Mathematical Optimizations](#part-7-mathematical-optimizations)

**IV. REFERENCE** 8. [Domain Vocabulary](#part-8-domain-vocabulary-ubiquitous-language) 9. [Appendix: Quick Reference](#appendix-quick-reference)

---

# I. FOUNDATIONS

## PART 1: MATHEMATICAL STRUCTURE

### 1.1 Core Mathematical Model

The LC-Seq compound library forms a **Directed Acyclic Graph (DAG)** representing truncation relationships. More specifically, it is a **partially ordered set (poset)** with the following properties:

**Formal Definition:**

- **Set**: V = {all compounds in the library}
- **Relation**: ≼ = "is a truncation of" (directed)
- **Properties**:
  - Reflexive: A ≼ A (every compound truncates to itself)
  - Antisymmetric: If A ≼ B and B ≼ A, then A = B
  - Transitive: If A ≼ B and B ≼ C, then A ≼ C
  - Directed: Edges flow from longer → shorter sequences
  - Acyclic: No cycles possible (can't increase length via truncation)

### 1.2 Chemical Identity vs Positional Identity

**Critical Distinction:**

**Chemical Peptide (Molecular Identity)**

- What molecule actually exists chemically
- Independent of synthesis position
- Examples:
  - [Val, Null, Null] = chemically "Val" (single amino acid)
  - [Null, Val, Null] = chemically "Val" (same molecule!)
  - [Null, Null, Val] = chemically "Val" (same molecule!)

**Positional Sequence (Synthesis Path)**

- How/when the peptide was synthesized
- Position-dependent encoding
- Examples:
  - [Val, Null, Null] ≠ [Null, Val, Null] ≠ [Null, Null, Val] (different synthesis paths)

**Implication for Graph Structure:**

- Graph vertices represent **chemical peptides** (in monomer mode) or **positional sequences** (in building-block mode)
- Same chemical peptide can have multiple synthesis paths → convergence in monomer mode

### 1.3 Two Analysis Modes, Two Graph Structures

#### Building-Block Mode (Poset Structure)

**Vertices**: Positional sequences (unique by synthesis path)

**Example compounds**:

- [Val, Null, Null] = vertex A
- [Null, Val, Null] = vertex B (different from A)
- [Null, Null, Val] = vertex C (different from A and B)

**Structure**: Forest (multiple disconnected trees)

**Properties**:

- No convergence
- Each positional variant is separate vertex
- Edges: building block subset relationships

#### Monomer-Level Mode (DAG with Convergence)

**Vertices**: Chemical peptides (monomer sequences)

**Example compounds**:

- [Val, Null, Null] → vertex "Val"
- [Null, Val, Null] → vertex "Val" (SAME as above!)
- [Null, Null, Val] → vertex "Val" (SAME as above!)

**Structure**: DAG with massive convergence

**Properties**:

- Multiple paths to same vertex (diamonds)
- Positional variants collapse to chemical identity
- Edges: monomer subsequence relationships

**Example Diamond Pattern**:

Val-Phe-Leu and Leu-Phe-Val both truncate through intermediate compounds (Val-Phe, Phe-Leu, Leu-Phe, Phe-Val) which all eventually converge to the same single-monomer vertices (Val, Phe, Leu). Both maximal compounds share the same descendant "Val" vertex, creating a diamond convergence pattern.

### 1.4 Why DAG Structure Holds (No Cycles)

**Proof by Contradiction:**

Assume a cycle exists: A → B → C → A

1. A → B means A has fewer monomers than B (or equal length with subset)
2. B → C means B has fewer monomers than C
3. C → A means C has fewer monomers than A
4. By transitivity: A has fewer monomers than A (contradiction!)

**Therefore:** Cycles impossible. Convergence (diamonds) ≠ cycles.

- **Cycle**: A → B → C → A (forbidden - can't return to start)
- **Diamond**: A → C, B → C (allowed - multiple paths down, never up)

### 1.5 Truncation Mechanics and Decomposition

This section explains how truncation relationships are computed in block-level vs monomer-level modes, addressing the combinatorial complexity that arises from composite building blocks.

#### 1.5.1 Peptide Sequence Convention

**N→C Notation Standard**:

Peptide sequences are written from **N-terminus** (amino end) to **C-terminus** (carboxyl end):

```
Leu-Leu-Ala-Val-Pro-Leu-DLeuMe-DPro-Leu-Leu-DPro
N-term ←──────────────────────────────────────→ C-term
```

**DNA-Encoded Library Synthesis Direction**:

DEL synthesis typically proceeds **C→N** (reverse of written sequence):

- **Position 0**: C-terminus (rightmost), synthesized FIRST, attached to DNA tag
- **Position 8**: N-terminus (leftmost), synthesized LAST
- DNA tag serves as synthesis anchor at C-terminus

**Position Numbering**:

```
Building Block View (what was synthesized at each position):
Position:   8              7         6     5      4        3      2     1     0
Block:     Leu      Leu-Ala-Val     Pro   Leu   DLeuMe   DPro   Leu   Leu   DPro
           (BB)      (composite)    (BB)  (BB)   (BB)     (BB)   (BB)  (BB)   (BB)

Monomer Sequence (chemical reality):
Leu - Leu - Ala - Val - Pro - Leu - DLeuMe - DPro - Leu - Leu - DPro
N-term ←──────────────────────────────────────────────────────→ C-term

Direction: N-terminus ←──────────────────────→ C-terminus
Synthesis: Last ←────────────────────────────→ First
```

**Key Point**: Position 7 contains a composite building block "Leu-Ala-Val" which appears identical to a regular tripeptide sequence. The distinction is contextual - you must know which positions contain composite blocks from the library design.

**Chemical Significance**:

- C-terminus (position 0) is the **anchor point** (most stable, first synthesized)
- Truncation typically occurs **N→C** (failure to add N-terminal residues)
- Right-aligning sequences preserves this chemical reality

#### 1.5.2 Building-Block Level Truncation

**Truncation Rule**: Replace one building block with Null

**Example Library**:

```
Cycle 0: {Leu, Val, Ala}
Cycle 1: {Pro, Gly}
Cycle 2: {Phe, Trp}
```

**Maximal Compound**: `Leu-Pro-Phe` (3 building blocks)

**All Direct Descendants** (one position → Null):

```
Leu-Pro-Phe  (maximal, 3 blocks)
  ├─→ Null-Pro-Phe  (removed Leu)
  ├─→ Leu-Null-Phe  (removed Pro)
  └─→ Leu-Pro-Null  (removed Phe)
```

**Properties**:

- Each truncation step removes exactly **one building block**
- Unambiguous: Position information preserved
- Graph structure: **Forest** (no convergence)
- Total compounds: 4 × 3 × 3 = 36 (including L₀)

**Edge Generation Algorithm (Block Mode)**:

```
For each compound C = [B₀, B₁, B₂, ..., Bₙ]:
  For each position i where Bᵢ ≠ Null:
    Create descendant D by setting Dᵢ = Null, all other positions unchanged
    Add edge: C → D
```

#### 1.5.3 Monomer-Level Decomposition

**Decomposition Rule**: Composite building blocks decompose to individual monomers

**Example with Composite Block**:

**Building-block level compound**: `Leu-Leu-Ala-Val-Pro`

Where position assignment (from library design) is:
- Position 2: Leu (single monomer building block)
- Position 1: Leu-Ala-Val (composite building block - 3 monomers)
- Position 0: Pro (single monomer building block)

**Key Point**: The sequence "Leu-Leu-Ala-Val-Pro" looks identical whether:
- All 5 are individual monomer blocks, OR
- Position 1 contains composite block "Leu-Ala-Val"

You must know from the library design which positions contain composite blocks.

**Decomposition to monomer-level**:

```
Building-block view: Leu - Leu-Ala-Val - Pro
                     (pos 2) (pos 1)    (pos 0)
                     1 mono  3 mono     1 mono

Monomer-level view:  Leu-Leu-Ala-Val-Pro
                     (5 individual monomers)
```

In monomer mode, the composite block "Leu-Ala-Val" decomposes into three individual monomers, creating a 5-mer sequence.

#### 1.5.3.1 Position vs Sequence Indexing

**Key Distinction**: Building-block positions vs monomer sequence indices

**Building-Block Mode**:
- Positions are synthesis positions: 0, 1, 2, ... (from library design)
- Position 1 might contain composite block "Leu-Ala-Val" (3 monomers)
- These are distinct: Position ≠ Sequence Index

**Monomer Mode - Complete Restructuring**:
- Positions are discarded entirely
- Graph vertices identified by monomer sequence only
- Example: "Leu-Leu-Ala-Val-Pro" (5 monomers)
  - Building-block positions {2, 1, 0} are NOT preserved
  - Monomer mode knows only: "this is a 5-mer with sequence Leu-Leu-Ala-Val-Pro"
  - No concept of "position 1 contains Leu-Ala-Val"

**Critical Point**: When switching from building-block → monomer mode, positional information is LOST. The graph completely restructures based on chemical identity alone.

**Example Showing Loss of Position**:

Building-block mode (3 different vertices):
- Vertex A: Position assignment {2:Leu, 1:Leu-Ala-Val, 0:Pro}
- Vertex B: Position assignment {2:Pro, 1:Leu, 0:Leu-Ala-Val}
- Vertex C: Position assignment {2:Leu-Ala-Val, 1:Pro, 0:Leu}

All three have different positional sequences (different synthesis paths).

Monomer mode (1 vertex - convergence!):
- Single vertex: "Leu-Leu-Ala-Val-Pro"
- Position information gone
- All three building-block vertices converge to this one monomer vertex
- This is the "diamond convergence" pattern

**Why This Matters**:
- Alignment (Part 2.5) works differently in each mode
- Building-block mode: Use positional alignment (unambiguous)
- Monomer mode: Must use left/right-align with ambiguity flag

**Monomer-Level Truncation**:

Now truncation operates on **individual monomers**, not blocks:

```
Leu-Leu-Ala-Val-Pro (5 monomers)
  ├─→ Leu-Ala-Val-Pro (removed N-terminal Leu)
  ├─→ Leu-Leu-Val-Pro (removed Ala)
  ├─→ Leu-Leu-Ala-Pro (removed Val)
  └─→ Leu-Leu-Ala-Val (removed C-terminal Pro)
  ... (and many more truncations at various positions)
```

**Key Difference**: Monomer-level graph has **vastly more vertices and edges**

#### 1.5.4 Convergence in Monomer Mode

**Critical Insight**: Multiple positional sequences can represent the same chemical peptide.

**Example**:

Three different building-block sequences:

```
[Leu, Null, Null] at synthesis positions [2, 1, 0]
[Null, Leu, Null] at synthesis positions [2, 1, 0]
[Null, Null, Leu] at synthesis positions [2, 1, 0]
```

All three decompose to the same monomer sequence: **Leu** (single amino acid)

**In Monomer Mode Graph**:

- These three positional variants **converge** to the same vertex "Leu"
- Creates **diamond patterns** (multiple ancestors → one descendant)
- Graph structure: **DAG with convergence**, not forest

**Example Diamond Pattern**:

```
Building-block level:
  [Val-Phe, Null, Null]    [Null, Val-Phe, Null]
           ↘             ↙
            [Val, Null, Null]  [Null, Val, Null]
                   ↘         ↙
                    "Val" (monomer-level vertex)

Chemical interpretation: Same molecule (Val), different synthesis paths
```

#### 1.5.5 Combinatorial Explosion

**Example Library**:

```
3 positions × 3 building blocks per position
Each building block is a tripeptide (3 monomers)
```

**Building-Block Mode**:

```
Total vertices: 4³ = 64 (including nulls)
Max edges per vertex: 3 (one per position)
Graph structure: Forest (no convergence)
```

**Monomer-Level Mode**:

```
Total unique monomers: 3 × 3 × 3 = 27 different monomers
Max sequence length: 3 positions × 3 monomers/block = 9 monomers
Total possible monomer sequences: Exponential in sequence length
Graph structure: DAG with massive convergence
Vertices: Thousands to millions (depending on library)
```

**Why the Explosion**:

- Block mode: Position information prevents convergence
- Monomer mode: Chemical identity independent of synthesis path
- Same peptide synthesized via different paths → convergence → more complex graph

#### 1.5.6 Edge Generation Algorithm (Monomer Mode)

**Algorithm**:

```
Step 1: Decompose Building Blocks
For each compound C with building blocks [B₀, B₁, ..., Bₙ]:
  For each block Bᵢ:
    Decompose Bᵢ → [m₁, m₂, ..., mₖ] (individual monomers)
  Concatenate all monomer lists → monomer_sequence(C)

Step 2: Generate Monomer-Level Truncations
For each compound C:
  monomer_seq = monomer_sequence(C)
  For each monomer mⱼ in monomer_seq:
    Create truncation T by removing mⱼ
    T_seq = monomer_seq with mⱼ deleted
    Add edge: C → T (in monomer-level graph)

Step 3: Identify Convergent Vertices
Group compounds by monomer sequence (ignoring synthesis path):
  If multiple compounds have same monomer_seq:
    → They converge to same vertex in monomer graph
    → All share same descendants
```

**Result**: DAG with diamond convergence patterns everywhere

#### 1.5.7 Mode Selection Implications

**Use Building-Block Mode when**:

- Analyzing synthesis fidelity per position
- Troubleshooting position-specific issues
- Positional information is meaningful
- Simpler graph preferred (faster computation)

**Use Monomer-Level Mode when**:

- Chemical identity more important than synthesis path
- Comparing across different synthesis strategies
- Aggregating results by molecular composition
- Studying convergent synthesis pathways

**Computational Trade-offs**:

```
Block Mode:
  + Smaller graph (fewer vertices/edges)
  + Faster traversal
  + Unambiguous relationships
  - Position-dependent (same molecule appears multiple times)

Monomer Mode:
  + Chemical identity-based (one vertex per molecule)
  + Natural for chemistry-focused analysis
  - Larger graph (more vertices/edges)
  - More complex (convergence patterns)
```

#### 1.5.8 Example: Complete Decomposition

**Sequence**: `Leu-Leu-Ala-Val-Pro`

**Building-Block Level Interpretation**:
- Position 2: Leu (single monomer block)
- Position 1: Leu-Ala-Val (composite block containing 3 monomers)
- Position 0: Pro (single monomer block)

**Important**: The sequence "Leu-Leu-Ala-Val-Pro" is written the same way whether position 1 is a composite block or individual monomers. Only the library design metadata indicates which positions contain composite blocks.

**Step 1: Monomer-Level Decomposition**:

```
Building-block view (3 blocks total):
  Position 2: Leu (1 monomer)
  Position 1: Leu-Ala-Val (3 monomers)
  Position 0: Pro (1 monomer)

Monomer-level view (5 monomers total):
  Leu-Leu-Ala-Val-Pro
  (All monomers treated individually)
```

**Step 2: Building-Block Descendants** (replace entire block with Null):

```
At building-block level:
Leu-Leu-Ala-Val-Pro (where middle is composite block)
  ├─→ Null-Leu-Ala-Val-Pro  (removed position 2 Leu block)
  ├─→ Leu-Null-Pro          (removed position 1 Leu-Ala-Val composite block)
  └─→ Leu-Leu-Ala-Val-Null  (removed position 0 Pro block)
```

**Step 3: Monomer-Level Descendants** (remove any single monomer):

```
Leu-Leu-Ala-Val-Pro (5 monomers)
  ├─→ Leu-Ala-Val-Pro (removed position-0 Leu)
  ├─→ Leu-Ala-Val-Pro (removed position-1 Leu - SAME as above! Convergence!)
  ├─→ Leu-Leu-Val-Pro (removed Ala)
  ├─→ Leu-Leu-Ala-Pro (removed Val)
  └─→ Leu-Leu-Ala-Val (removed Pro)
```

**Notice**: Two different Leus removed → same descendant "Leu-Ala-Val-Pro"
**This is convergence**: Multiple truncation paths lead to same chemical product

**Full Monomer-Level Graph for This Compound** (partial):

```
                    Leu-Leu-Ala-Val-Pro (5-mer)
                   /  /    |    \     \
                  /  /     |     \     \
    Leu-Ala-Val-Pro Leu-Leu-Val-Pro ... (4-mers, many vertices)
              |   \  /   |
              |    \/    |  (convergence - diamonds form)
              |    /\    |
    Leu-Ala-Val  Leu-Val-Pro ... (3-mers)
         |    \  /   |
         |     \/    |
         |     /\    |
    Leu-Ala  Leu-Val ... (2-mers)
       |  \  /  |
       |   \/   |
       |   /\   |
      Leu Ala Val Pro (1-mers)
        \  |  |  /
         \ | | /
          \|_|/
           L₀ (null compound - all paths lead here)
```

---

## PART 2: DOMAIN FOUNDATIONS

### 2.1 Core Entities

**Compound**

- Entity defined by unique sequence identifier
- Contains: sequence, building blocks, chromatogram, detected peaks, selected peak
- May have hierarchical relationships (ancestors/descendants)
- Examples:
  - Leu-Val-Pro (maximal in dataset)
  - Leu-Null-Pro (descendant of above)

**Peak**

- Detected feature in chromatogram representing compound elution
- Defined by: position (retention time), boundaries, metrics
- May be classified as NULL, TRUNCATION, PUTATIVE_PRODUCT, or UNKNOWN
- Optionally has Gaussian fit parameters

**Chromatogram**

- Elution profile across fractions/time
- Contains: time_points, counts (raw signal)
- Supports multiple signal variants stored in dict:
  - "raw": Primary signal used for analysis
  - "derivative": First derivative (optional)
  - "derivative_2": Second derivative (optional)
  - Custom processing can add more

**BuildingBlock**

- Chemical building block used in library synthesis
- Defined by: cycle (synthesis round), code (identifier), is_null flag
- May be composite (e.g., "Leu-Ala-Val" as single building block = 3 monomers)
- **Key**: Composite blocks written identically to regular sequences (e.g., "Leu-Ala-Val")
- Context from library design determines if a sequence segment is one composite block vs multiple individual blocks
- Null detection: Any code containing "null" (case-insensitive) → null block
- Standard null code: "Null"

### 2.2 Sequence Representations

**Positional Sequence**

- Exact representation including null positions
- Example: Leu-Null-Pro
- Unique per compound (synthesis path)
- Attribute: Compound.positional_sequence or Compound.sequence

**Canonical Sequence** (or **Residue Sequence**)

- Non-null building blocks only
- Example: Leu-Pro (from Leu-Null-Pro)
- Multiple positional sequences → same canonical
- Used for grouping positional variants
- Attribute: Compound.residue_sequence

**Monomer Sequence**

- Fully decomposed individual monomers
- Example: "Leu-Leu-Ala-Val-Pro" where position 1 is composite block "Leu-Ala-Val"
  - Decomposes to: ["Leu", "Leu", "Ala", "Val", "Pro"] (5 individual monomers)
- Chemical peptide identity (what molecule exists)
- Composite blocks expand to their constituent monomers
- Method: BuildingBlock.decompose_to_monomers()

### 2.3 Temporal Representation and Alignment

#### 2.3.1 Absolute Time Representation

Chromatogram signals are defined over **absolute time** (seconds or minutes from injection). Peak retention times are scalar values representing the time at which a compound elutes.

**Core Principle**: All retention time comparisons operate on scalar time values, independent of signal boundaries or sampling rates.

**Properties**:

- Retention times are physically meaningful (not array indices)
- Peak position comparisons: direct scalar arithmetic
- Example: "Does peak at 45.2s match expected 45.0s?" → |45.2 - 45.0| = 0.2s

**Time Units**:

- Standard: seconds (s) or minutes (min)
- Must be consistent across dataset
- Attribute: Chromatogram.time (array of time points)

#### 2.3.2 Variable Signal Boundaries

Signals from different experimental runs may have different start and end times due to:

- Manual data collection timing variation
- Different run durations
- Preprocessing operations (trimming, artifact removal)

**Example**:

- Signal A: t ∈ [600, 3600] seconds (10:00 to 60:00 min)
- Signal B: t ∈ [610, 3630] seconds (10:10 to 60:30 min)
- Signal C: t ∈ [595, 3590] seconds (9:55 to 59:50 min)

**Impact on Analysis**:

**✅ No alignment required for**:

- **Peak detection**: Each signal analyzed independently in its own time domain
- **Peak classification**: Compares scalar retention times (position matching via absolute time)
- **Retention order validation**: Scalar comparisons (t_product > t_truncation)
- **Purity calculation**: Sum counts over peak regions per compound
- **SNR calculation**: Ratio of intensities per compound
- **Synthesis validation**: All metrics operate on peak positions or per-compound counts

**⚠️ Alignment required for**:

- **Similarity metrics**: Wasserstein distance, Jensen-Shannon divergence (signal-to-signal comparison)
- **Hierarchical clustering**: Requires pairwise similarity matrix
- **Similarity-based ordering**: For visualization (optional feature)

#### 2.3.3 Simple Alignment for Similarity Metrics

When similarity-based compound ordering is needed (hierarchical clustering for visualization), signals must be aligned to a common time grid.

**Algorithm**:

**Step 1: Define Global Time Window**

Compute global bounds across all signals in dataset:
t_min_global = min(signal.t_min for all signals)
t_max_global = max(signal.t_max for all signals)

This ensures all signals have real data coverage in the window (no padding needed).

**Step 2: Construct Common Grid**

Define uniform grid with finest sampling resolution:
dt = min(signal.sampling_interval for all signals)
t_grid = linspace(t_min_global, t_max_global, N)

where N = (t_max_global - t_min_global) / dt

**Step 3: Trim and Interpolate**

For each signal S:

1. **Trim**: Extract data where t ∈ [t_min_global, t_max_global]
2. **Interpolate**: Use linear interpolation to resample onto t_grid
3. **Store**: S.aligned_signal = interpolated values on t_grid

**Step 4: Compute Similarity**

For all signal pairs (A, B):
similarity(A, B) = metric(A.aligned_signal, B.aligned_signal)

where metric ∈ {Wasserstein distance, Jensen-Shannon divergence, etc.}

**Properties**:

- ✅ All signals have real data in [t_min_global, t_max_global] (by definition)
- ✅ No boundary padding needed (no assumptions about signal outside range)
- ✅ Minimal data loss (only extreme tails where one signal started earlier/ended later)
- ✅ Consistent treatment across entire dataset
- ✅ Preserves peak positions and relative intensities

**Typical Data Loss**:
For signals with 10-60 minute range and ±30 second variation:

- Global window: ≈ [10:30, 59:30] = 49 minutes
- Loss per signal: <1% of total duration
- All critical peaks retained (products typically elute in middle of gradient)

**When NOT to Use**:
If similarity-based ordering is not required, skip alignment entirely. Use alternative orderings:

- By hierarchy level (truncation rank)
- By sequence length (number of building blocks)
- By positional sequence (lexicographic)
- By validation status (VALIDATED first, etc.)

**Computational Complexity**:

- Finding global bounds: O(V) for V signals
- Interpolation: O(V × n) where n = average signal length
- Total alignment: O(V × n) per dataset

#### 2.3.4 Alignment in Building-Block vs Monomer Mode

**Critical Note**: The temporal alignment described above (for similarity metrics) is separate from **sequence alignment for visualization** (see examples in Part 1.5.3.1).

**Building-Block Mode** (no sequence alignment ambiguity):
- Position information preserved from synthesis
- When visualizing lineages: align by synthesis position
- Example: Position 7 in all variants aligns perfectly (even if composite block)
- **No ambiguity** - positions are explicit from library design
- Temporal alignment (this section) used only for similarity metrics

**Monomer Mode** (sequence alignment IS ambiguous):
- Position information discarded (graph restructured by chemical identity)
- When visualizing lineages: must use heuristic (right-align or left-align)
- Example: "Leu-Leu-Leu" could map to many position combinations
- **Always ambiguous** when repeated residues present
- Both temporal alignment (this section) AND sequence alignment needed

**Key Distinction**:
- **Temporal alignment** (Part 2.3): Aligning chromatogram time grids for similarity computation
  - Needed: When computing signal similarity metrics
  - Works on: Time points (seconds/minutes)
  - Result: Common time grid for correlation/distance calculation

- **Sequence alignment** (for visualization): Aligning sequences for display
  - Needed: When visualizing lineages or comparing sequences
  - Works on: Amino acid residues or building blocks
  - Result: Sequences with gaps showing truncation patterns

**Recommendation by Mode**:
- Building-block mode → No sequence alignment issues (use positional mapping)
- Monomer mode → Acknowledge alignment ambiguity (use right-align with flag)

Both modes may need temporal alignment if computing similarity metrics for clustering/ordering.

### 2.4 The Null Compound (L₀)

#### 2.4.1 Chemical Definition

**L₀** (the null compound) is defined as the compound with all building blocks set to null:
L₀ = [Null, Null, ..., Null]

**Physical Molecule**:

- DNA barcode tag (for sequencing identification)
- Linker chemistry (connects DNA to synthesis site)
- Any protecting groups or caps from synthesis
- **NO amino acids or building blocks attached**

**Chemical Interpretation**:

- L₀ represents the starting scaffold before peptide synthesis
- Every library member begins as L₀, then building blocks are sequentially added
- L₀ is the "bare tag" - the molecular entity present before any peptide construction

#### 2.4.2 Graph-Theoretic Properties

**Minimal Element**:

- L₀ is the unique minimal element in building-block mode DAG
- No compound has fewer building blocks than L₀ (zero is minimum)
- Formally: ∄C ∈ V such that C ≺ L₀

**Universal Descendant Property**:

- Every compound C can truncate to L₀ by removing all building blocks
- Formally: ∀C ∈ V, L₀ ≼ C (L₀ is descendant of all compounds)
- L₀ is reachable from every vertex via directed truncation edges

**Topological Role**:

- Bottom element in partial order (poset)
- Anchor point for bottom-up DAG traversal (topological sort starts here)
- First compound processed in constraint propagation algorithm

**Proof of Universal Descendant**:
For any compound C = [B₁, B₂, ..., Bₙ] where some Bᵢ ≠ Null:

- Truncate each Bᵢ → Null sequentially
- After n truncations: [Null, Null, ..., Null] = L₀
- Therefore C ≻ L₀ (C is ancestor of L₀)
- By definition of descendant: L₀ ≼ C
- QED

#### 2.4.3 Retention Time Ordering Principle

**Hydrophobicity Additivity**:

- DNA tag + linker has base hydrophobicity H₀
- Each building block adds hydrophobicity: H(C) = H₀ + Σ(H_building_block)
- Longer peptides → higher total hydrophobicity → longer retention

**Early Elution Hypothesis**:

- DNA backbone is polar/hydrophilic (phosphate groups, sugar-phosphate backbone)
- Linker chemistry typically hydrophilic or moderately hydrophobic
- Most peptide building blocks increase hydrophobicity
- **Expected**: t(L₀) ≤ t(compounds with peptides)

**Chemical Rationale**:

- L₀ = tag only → relatively hydrophilic
- Compounds with amino acids → additional hydrophobic character
- Result: L₀ should elute in early region of chromatogram

**Exceptions** (rare):

- Highly charged peptides (multiple Arg, Lys, Asp, Glu)
- Extremely hydrophilic building blocks
- Modified tags with high hydrophobicity
- These may violate t(L₀) < t(compound) ordering

#### 2.4.4 L₀ Peak Detection Algorithm

**Primary Method: Earliest Significant Peak in Early Region**

**Step 1: Define Early Elution Window**

Compute early region boundary:
t_early_max = t_min + 0.2 × (t_max - t_min)

This defines the first 20% of chromatogram as "early region" where L₀ is expected.

**Step 2: Detect Significant Peaks**

Find all statistically significant peaks in early window:
early_peaks = {peaks with Z > 3 where t ∈ [t_min, t_early_max]}

Poisson Z-score threshold ensures we only consider real peaks, not noise.

**Step 3: Select Earliest Peak**

Choose L₀ peak by minimum retention time:
L₀_peak = argmin(peak.time) for peak in early_peaks

**Rationale for "Earliest" Selection**:

- Earlier peak = more hydrophilic → more likely pure DNA tag
- Later peaks in early region could be very short peptides (e.g., single amino acid)
- Conservative choice: Ensures NULL reference is truly the minimal compound
- Chemically justified: DNA tag should elute before any peptide-containing compounds

**Multiple Peak Handling**:

- If multiple significant peaks detected in early region
- Select the one with **minimum retention time** (earliest elution)
- Later peaks likely represent short truncation products (L₁, L₂), not pure tag (L₀)
- Example: Peak at 10.2 min and peak at 12.5 min → Choose 10.2 min as L₀

**Fallback Method: Signal Start**

If no detectable peak in early region (no peaks pass Z > 3 threshold):
L₀_position = t_min (signal start time)
L₀_intensity = signal(t_min)

**Rationale**:

- Conservative: t_min is earliest possible position
- Safe NULL reference (any peak after start is valid)
- Enables analysis to proceed even with poor L₀ signal quality
- Worst case: NULL position = 0, all other peaks treated as non-NULL

#### 2.4.5 Quality Control and Validation

**Late L₀ Detection (Red Flag)**

If L₀_peak.time > t_min + 0.5 × (t_max - t_min):
Issue WARNING: "L₀ peak detected in late elution region (t = {L₀_peak.time})"

**Likely Causes**:

- Serious contamination (non-peptide impurity with high hydrophobicity)
- Incorrect dataset (not from expected library, wrong tag sequence)
- Severe signal distortion or artifact
- DNA degradation, modification, or tag loss
- Tag chemistry fundamentally different than expected

**Recommended Action**:

- **STOP automated analysis**
- Inspect raw chromatogram visually
- Verify dataset identity (correct library, correct experiment)
- Investigate signal processing issues
- Investigate synthesis or purification issues
- Consider excluding sample from analysis

**Rationale**: Late-eluting L₀ violates fundamental hydrophobicity principle and indicates data quality problems too severe for automated handling.

**Validation Checks**:

**Check 1: L₀ Exists**

- Verify L₀ = [Null, Null, ..., Null] is present in dataset
- If missing: ERROR - cannot proceed with NULL peak classification

**Check 2: Detectable Signal**

- Verify L₀ has signal above noise floor (SNR > 2)
- If below detection: WARNING - background estimation unreliable

**Check 3: Early Elution**

- Verify L₀_peak.time < 0.5 × (t_max - t_min)
- If violated: WARNING - review data quality

**Check 4: Retention Consistency**

- If multiple L₀ measurements across dataset, verify consistent retention
- Expected variation: < 5% of median retention time
- If high variation: WARNING - chromatography instability or tag heterogeneity

#### 2.4.6 Role as Universal Reference Point

**NULL Peak Definition**

**Definition**: A NULL peak is a peak at retention time t(L₀_peak)

**Universal Constraint**: Any compound C may have a peak at position t(L₀_peak)

- Reason: Any compound can undergo complete truncation (all BBs → Null)
- Result: Peak at L₀ position represents DNA tag only, no peptide
- Classification: These peaks are labeled as NULL type

**Used In**:

- Peak classification constraint propagation
- Distinguishing product from truncation peaks
- Identifying complete degradation to tag-only

**Background Estimation**

**Primary Method**: Use L₀ chromatogram off-peak regions
background = median(L₀_signal where not in peak region)

**Rationale**: L₀ has minimal complexity (tag only), so off-peak signal is pure noise

**Alternative Method** (if L₀ unavailable or unreliable):
background = median(bottom 10% of all scaled_count values across dataset)

**SNR Reference**: All other compounds compute SNR relative to this background:
SNR(C) = max(product_counts_C) / background

**Retention Anchor**

**Implicit Reference**: All retention times referenced to L₀ as reference point

- L₀ represents "zero peptide" retention
- Compounds measured by how much later they elute vs L₀
- Retention order: t(L₀) ≤ t(L₁) ≤ ... ≤ t(Lₙ) (ideally)

**Classification Anchor**

**Bottom-Up Traversal**: DAG constraint propagation starts at L₀

1. Process L₀ first (topological sort minimal element)
2. Extract L₀ peak position as NULL constraint
3. Propagate NULL constraint to all ancestors (all other compounds)
4. Continue processing compounds in topological order

**Foundation**: L₀ provides the base case for recursive peak classification

#### 2.4.7 Edge Cases and Failure Modes

**Case 1: L₀ Not in Dataset**

**Problem**: All-null compound [Null, Null, ..., Null] is missing

**Impact**:

- Cannot establish NULL peak reference
- Cannot perform NULL peak classification
- Background estimation compromised

**Resolution**:

- ERROR: Require L₀ compound in dataset
- User must add L₀ or skip NULL-based classification
- Fallback: Use alternative background method, omit NULL classification

**Case 2: L₀ Below Detection Limit**

**Problem**: L₀ exists but has no detectable peaks (all signal near noise floor)

**Causes**:

- DNA didn't amplify during PCR
- Tag loss during chromatography
- Signal outside collection time window
- Sequencing depth too low

**Impact**:

- No NULL peak position available
- Background estimation unreliable
- SNR calculations less accurate

**Resolution**:

- WARNING: Flag reduced confidence in validation
- Fallback L₀ position: use t_min (signal start)
- Use alternative background (low-count tail method)
- Reduce confidence in synthesis validation results

**Case 3: Multiple Peaks at L₀**

**Problem**: L₀ chromatogram has 2+ significant peaks

**Causes**:

- Tag degradation (partial fragments)
- Contaminant co-elution
- Tag sequence variants in library
- Oligomerization of tags

**Resolution**:

- Select earliest peak as NULL reference (most conservative)
- Alternative: Select peak with highest Z-score (most significant)
- Flag ambiguity in QC report
- If >3 peaks: WARNING - data quality issue

**Case 4: Mode-Dependent Minimality**

**Building-Block Mode**:

- L₀ is unique minimal element
- All compounds truncate to same L₀
- NULL peak unambiguous

**Monomer Mode**:

- Multiple minimal elements possible (individual amino acids: Leu, Val, Ala, etc.)
- Each single amino acid is minimal (cannot truncate further without → L₀)
- L₀ still exists as separate minimal element (tag only, no AAs)

**Resolution**:

- NULL peak definition still based on physical L₀ compound (tag only)
- Monomer-level minimal elements (single AAs) are distinct from L₀
- Classification unaffected (NULL = tag position, regardless of mode)

#### 2.4.8 Theoretical Constraints

**Universal NULL Constraint**

**Formal Statement**: ∀C ∈ V, C may have peak at position(L₀_peak)

**Justification**: Any compound can undergo complete truncation:
C = [B₁, B₂, ..., Bₙ] → [Null, Null, ..., Null] = L₀

**Used In**: Peak classification as global constraint

- When classifying peaks for compound C
- A peak at position ≈ t(L₀_peak) → candidate for NULL classification
- Applies to ALL compounds (universal)

**Retention Ordering Expectation**

**Expected Ordering**: t(C) ≥ t(L₀) for compounds with building blocks

**Physical Basis**: Hydrophobicity additivity (peptides increase retention)

**Validation**:

- Check: Does t(C) ≥ t(L₀) hold for >90% of compounds?
- If violated frequently: Data quality issue or unusual tag chemistry

**Violations**:

- Highly hydrophilic peptides (charged residues) may elute before L₀
- Flag as anomalies for review
- May indicate synthesis error, tag modification, or unusual chemistry

**Completeness Property**

**Definition**: L₀ represents maximal truncation (no building blocks remain)

**Chemical Interpretation**:

- L₀ peak = only DNA tag detected in sequencing
- No peptide product present
- Not a synthesis failure (L₀ is the starting material!)

**Implication for Validation**:

- Presence of L₀ peak ≠ synthesis failure
- L₀ is expected to have dominant peak (it's the "null" reference)
- Other compounds compared against L₀ as reference

---

## PART 3: HIERARCHICAL RELATIONSHIPS

### 3.1 Terminology: Ancestry and Lineage (Not Parent/Child)

**Problem with "Parent":**

- In combinatorial library, Val-Leu-Ala can be:
  - "Parent" of Val-Leu-Null (its truncation)
  - "Child" of Val-Leu-Ala-Pro (longer compound)
- No absolute "parent" - only relative relationships

**Correct Terminology:**

**Directional Relationships:**

- **Ancestor**: Compound B is ancestor of A if A truncates from B (B has more building blocks)
- **Descendant**: Compound A is descendant of B if A has subset of B's building blocks
- **Comparable**: Two compounds where one is ancestor/descendant of other
- **Incomparable**: Two compounds in different branches (neither truncates to other)

**Positions in Hierarchy:**

- **Maximal Compound** (not "parent"): No ancestors in dataset (longest available)
- **Minimal Compound**: No descendants (shortest, e.g., all-null)
- **Reference Compound** (not "parent"): The compound currently being analyzed
- **Query Compound**: Compound being studied in current analysis

**Analysis Scope:**

- **Lineage**: All ancestors + descendants + self
- **Descendant Set** (Principal Ideal): All compounds that truncate from X
- **Ancestor Set** (Principal Filter): All compounds that X truncates from
- **Connected Component**: All compounds reachable through truncation relationships

### 3.2 Mathematical Poset Terminology

From **Order Theory**:

- **Poset**: Partially ordered set (our hierarchy)
- **Chain**: Totally ordered sequence (A₁ ⊂ A₂ ⊂ A₃ ⊂ ... ⊂ Aₙ)
- **Antichain**: Set of mutually incomparable elements (all siblings)
- **Principal Ideal** ↓X: All descendants of compound X
- **Principal Filter** ↑X: All ancestors of compound X
- **Hasse Diagram**: Visualization of the poset (what we draw!)
- **Maximal Element**: Compound with no ancestors (in dataset)
- **Minimal Element**: Compound with no descendants

From **Graph Theory**:

- **Topological Sort**: Order compounds for bottom-up analysis
- **Transitive Closure**: All ancestor/descendant relationships
- **Reachability**: Can we get from A to B via truncation edges?
- **Longest Path**: Deepest truncation chain
- **Connected Component**: Isolated family/lineage

### 3.3 Hierarchy Properties

**Truncation Level** (or **Level** or **Rank**)

- Number of non-null building blocks
- Level 0 = all null (complete truncation)
- Level N = maximal in dataset (N building blocks)
- Attribute: Compound.level

**Monomer Level**

- Number of individual monomers present
- Example: Sequence "Leu-Leu-Ala-Val-Pro" where position 1 is composite = monomer level 5 (1+3+1 monomers)
- Must decompose building blocks to count total monomers
- Used in monomer-based hierarchy mode
- Attribute: Compound.monomer_level

**Hierarchy Modes:**

- **Block Mode**: Building blocks as atomic units

  - 3 building blocks → levels 0-3
  - Forest structure (separate trees per maximal compound)

- **Monomer Mode**: Individual monomers as atomic units
  - 3 trimeric blocks (9 monomers) → levels 0-9
  - DAG with convergence (diamonds everywhere)

---

## PART 4: GRAPH PROPERTIES AND PATTERNS

### 4.1 Convergence in Monomer-Level Analysis

**Key Insight:** Multiple positional variants represent the same chemical peptide

**Example**:

Chemical peptide "Val-Phe" can be synthesized via three different synthesis paths:

- [Val-Phe, Null, Null] (trimeric BB in cycle 0)
- [Val, Phe, Null] (monomeric BBs in cycles 0-1)
- [Null, Val-Phe, Null] (trimeric BB in cycle 1)

All three synthesis paths converge to the same chemical "Val-Phe" vertex in monomer mode.

**Convergence Patterns:**

1. **Same-length convergence**: Different synthesis paths, same peptide
2. **Multi-ancestor convergence**: Many ancestors → one descendant
3. **Diamond structures**: Multiple paths from different ancestors to shared descendant

**Why This Matters:**

- Consensus mode aggregates positional variants (same chemical peptide)
- Peak selection considers all synthesis paths for a peptide
- Similarity analysis groups by chemical identity, not synthesis

### 4.2 Equivalence Classes and Consensus Analysis

#### 4.2.1 EquivalenceClass Definition

**EquivalenceClass**: Collection of positional sequences that represent the same chemical peptide (same residue sequence).

**Formal Definition**: An equivalence class under the relation R: "has same residue sequence"

**Equivalence Relation Properties**:

- **Reflexive**: Compound C related to itself (C R C)
- **Symmetric**: If A R B, then B R A
- **Transitive**: If A R B and B R C, then A R C

**Example**:

```
Equivalence class for "Val":
  - [Val, Null, Null] (position 1)
  - [Null, Val, Null] (position 2)
  - [Null, Null, Val] (position 3)

All three positional variants represent the same chemical peptide: Val
```

**Chemical Interpretation**:

- Same residue sequence = same chemical molecule
- Different positional sequences = different synthesis paths to same molecule
- EquivalenceClass groups synthesis paths by chemical identity

**Entity**: EquivalenceClass (collection of positional variants)

#### 4.2.2 Consensus Mode: Optional Optimization

**Purpose**: Consensus mode is an **optional performance optimization** that can provide:

1. **Computational speedup**: Process N equivalence classes instead of 3N-10N positional variants (~3-10× faster)
2. **Noise reduction**: Average out positional encoding noise and experimental variability
3. **Molecular focus**: One result per chemical molecule, not per synthesis path

**Important**: Consensus mode is **NOT required** - individual mode (analyzing each positional variant separately) is fully valid and remains the default approach.

**When to Consider Consensus Mode**:

- ✅ Large datasets (many positional variants)
- ✅ Variants have similar chromatographic behavior
- ✅ Molecular-level analysis (focus on chemistry, not synthesis path)
- ✅ Exploratory analysis (overview of library)

**When to Use Individual Mode** (default):

- ✅ Position-specific diagnostics needed
- ✅ Synthesis troubleshooting (identify problematic positions)
- ✅ Small datasets (speedup not critical)
- ✅ Variants have different behavior (consensus invalid)
- ✅ Regulatory/validation requirements (need actual measurements)

**Performance Impact**:

Individual mode:

```
N_variants compounds to process (e.g., 3³ = 27 for 3-position library)
Each processed independently through full pipeline
```

Consensus mode:

```
N_classes equivalence classes (e.g., 7 unique residue sequences)
Peak detection on consensus (once per class)
Area integration on variants (per variant, but cheap)

Speedup ≈ N_variants / N_classes (typically 3-10×)
```

**Recommendation**: Start with individual mode (simpler, always valid); adopt consensus mode if speedup needed and variants validated as similar.

#### 4.2.3 Hybrid Consensus Strategy

**Key Insight**: Peak detection is expensive, area integration is cheap. Do expensive operations once (on consensus), cheap operations per variant (on individuals).

**Hybrid Approach**:

**Phase 1: Peak Detection on Consensus Signal** (expensive, do once)

```
For each equivalence class:
  1. Aggregate variant signals → consensus_signal
  2. Peak detection (Discrete Morse + Poisson + Prominence) on consensus
  3. Peak classification (DAG constraints) on consensus

Output: Representative peak positions and boundaries for the molecule
```

**Phase 2: Area Integration on Individual Variants** (cheap, do per-variant)

```
For each variant in equivalence class:
  1. Use consensus peak boundaries (from Phase 1)
  2. Integrate areas in variant's own signal (not consensus!)
  3. Compute purity for this variant
  4. Validate this variant independently

Output: Individual purities and validation categories
```

**Phase 3: Aggregate Statistics** (summary for reporting)

```
For equivalence class:
  mean_purity = mean(variant_purities)
  std_purity = std(variant_purities)
  min_purity = min(variant_purities)

Class-level validation:
  If all variants VALIDATED → Class: VALIDATED
  If any variant FAILED → Class: FAILED (flag problematic position)
  If mixed results → Class: UNCERTAIN (heterogeneous quality)
```

**Why This Works**:

- Speedup: Peak detection done once (most expensive step)
- Validity: Purity computed from real signals (not averaged)
- Information: Position-specific results available for diagnostics
- Summary: Molecular-level overview for reporting

**What This Avoids**:

- ❌ Averaging purities (loses position-specific info)
- ❌ Processing every variant through expensive peak detection
- ❌ Mixing consensus and individual purities in dataset (comparability issues)

#### 4.2.4 Consensus Signal Aggregation

**Prerequisite**: Signal alignment (Part 2.3) - all variants on common time grid

**Aggregation Algorithm**:

**Step 1: Align Variants to Common Grid**

```
For equivalence class with variants V = {v₁, v₂, ..., vₙ}:

Find global time window (Part 2.3):
  t_min = min(vᵢ.t_min for all i)
  t_max = max(vᵢ.t_max for all i)

Interpolate each variant to common grid:
  For each vᵢ:
    vᵢ_aligned = interpolate(vᵢ.signal, t_grid)
```

**Step 2: Compute Consensus Signal**

```
Aggregation method (recommend mean):
  consensus_signal(t) = mean(v₁_aligned(t), v₂_aligned(t), ..., vₙ_aligned(t))

Alternative (robust to outliers):
  consensus_signal(t) = median(v₁_aligned(t), v₂_aligned(t), ..., vₙ_aligned(t))

Recommendation: Mean (smoother, faster, standard)
```

**Step 3: Validity Check** (essential!)

```
Compute pairwise correlations between variant signals:
  For all pairs (vᵢ, vⱼ):
    corr(vᵢ, vⱼ) = correlation(vᵢ_aligned, vⱼ_aligned)

Check: min(all correlations) > threshold (e.g., 0.8)

If correlation too low:
  → Variants have different chromatographic behavior
  → Consensus aggregation invalid
  → Fall back to individual mode
  → Report: "Positional variants differ - individual analysis required"
```

**Validity Assumption**: Positional variants have similar chromatographic behavior (peaks at similar positions, similar shapes, similar intensities).

**Why Correlation Check Matters**:

- If variant 1 has peak at t=20, variant 2 at t=25 → averaging creates artificial broad "peak"
- If variant 1 clean, variant 2 very noisy → averaging doesn't represent either
- Low correlation = consensus invalid = must use individual mode

**Recommended Threshold**: correlation > 0.8 (strong similarity required)

#### 4.2.5 Peak Detection on Consensus Signal

**Apply full detection pipeline to consensus signal**:

**CRITICAL**: Consensus mode does NOT change processing parameters. The same uniform parameters (scale ranges) apply to both consensus signals and individual variant signals to ensure comparability.

**Step 1: Peak Detection**

```
peaks = detect_peaks(consensus_signal)
(Same parameters as individual mode - Part 5.1-5.2)
```

**Step 2: Peak Detection**

```
Discrete Morse theory (Part 5.1) + Poisson statistics (Part 5.2)
→ Detect significant peaks in consensus signal
Output: Peak positions, boundaries, prominence, Z-scores
```

**Step 3: Peak Classification**

```
DAG constraint propagation (Part 5.3-5.4)
→ Classify peaks: NULL, TRUNCATION, PUTATIVE_PRODUCT, UNKNOWN
Output: Peak labels for consensus
```

**Output**: Representative peak information for the equivalence class

- Peak positions (retention times)
- Peak boundaries [t_start, t_end] for each peak
- Peak classifications

**These boundaries will be used for ALL variants in the next step.**

**Rationale**:

- Peak positions should be similar across variants (if correlation valid)
- Consensus signal has better SNR (noise averaged out)
- Detecting once is much faster than detecting N times

#### 4.2.6 Area Integration on Individual Variants

**Use consensus peak boundaries, but integrate on individual variant signals** (NOT consensus).

**Algorithm**:

```
For each variant vᵢ in equivalence class:

  Input:
    - vᵢ_signal_raw (variant's raw signal)
    - peak_boundaries from consensus (Phase 1)

  Step 1: Use raw variant signal directly
    vᵢ_signal = vᵢ_signal_raw

  Step 2: Integrate areas using consensus boundaries
    For each peak with boundaries [t_start, t_end]:
      Area_peak_vᵢ = Σ vᵢ_signal(t) for t ∈ [t_start, t_end]

  Step 3: Compute purity for variant vᵢ
    Purity(vᵢ) = Area(product_vᵢ) / [Area(product_vᵢ) + Area(truncations_vᵢ) + Area(unknowns_vᵢ) + Area(null_vᵢ)]

  Step 4: Validate variant vᵢ (Part 6)
    Compare Purity(vᵢ) to dataset percentiles
    → Validation category: VALIDATED, LIKELY_SUCCESS, UNCERTAIN, LIKELY_FAILURE, FAILED

  Step 5: Store individual results
    Store: Individual purity and validation for this variant
```

**Critical**: Each variant's purity computed from its OWN signal, not from consensus.

**Why This Approach**:

- ✅ Purity is real measurement (from actual variant signal)
- ✅ Position-specific information preserved
- ✅ Can identify problematic synthesis paths
- ✅ Valid for cross-compound comparison (all measured consistently)
- ✅ Cheap operation (simple integration, no expensive peak detection)

**What This Avoids**:

- ❌ Averaging purities (loses position-specific info)
- ❌ Purity from consensus signal (not a real measurement)
- ❌ Mixing consensus and individual purities (comparability issues)

#### 4.2.7 Aggregate Statistics and Reporting

**Class-Level Summary**:

From individual variant purities, compute aggregate statistics:

```
For equivalence class with variants {v₁, v₂, ..., vₙ}:

Purity statistics:
  mean_purity = mean(Purity(v₁), Purity(v₂), ..., Purity(vₙ))
  std_purity = std(Purity(v₁), Purity(v₂), ..., Purity(vₙ))
  min_purity = min(Purity(v₁), Purity(v₂), ..., Purity(vₙ))
  max_purity = max(Purity(v₁), Purity(v₂), ..., Purity(vₙ))

Validation summary:
  If all variants: VALIDATED → Class: VALIDATED
  If all variants: FAILED → Class: FAILED
  If any variant: FAILED → Class: FAILED (flag specific variant)
  If mixed results → Class: UNCERTAIN or HETEROGENEOUS
```

**Reporting**:

**Molecular-level report** (high-level overview):

```
EquivalenceClass "Val":
  Mean purity: 0.83 ± 0.15 (std)
  Range: [0.65, 0.95]
  Validation: HETEROGENEOUS (mixed results)
```

**Position-specific report** (diagnostic detail):

```
EquivalenceClass "Val":
  [Val, Null, Null] (position 1): Purity 0.90, VALIDATED
  [Null, Val, Null] (position 2): Purity 0.65, UNCERTAIN ← FLAG
  [Null, Null, Val] (position 3): Purity 0.95, VALIDATED

Interpretation: Val synthesis good at positions 1&3, problematic at position 2
```

**Use Case**:

- Molecular report for high-level screening
- Position report for synthesis troubleshooting
- Both available from hybrid approach

#### 4.2.8 Validity Requirements

**Consensus mode is valid ONLY IF positional variants have similar chromatographic behavior.**

**Validity Checks**:

**Check 1: Signal Correlation**

```
For each equivalence class:
  Compute pairwise correlations between aligned variant signals

  correlation_threshold = 0.8

  If min(correlations) < correlation_threshold:
    → Variants too different
    → Consensus invalid
    → Fall back to individual mode
    → Flag: "Position-dependent chromatography detected"
```

**Check 2: Peak Position Consistency**

```
If using individual mode first (to validate):
  Detect peaks in each variant independently
  Check: Are peak positions similar across variants?

  For each peak type (product, truncations):
    Δt = max(variant_positions) - min(variant_positions)

    If Δt > retention_precision × 3:
      → Peak positions vary significantly
      → Consensus may merge/blur peaks
      → Consider individual mode
```

**Check 3: Purity Variance**

```
After computing individual purities:
  std_purity = std(variant_purities)

  If std_purity > 0.2:
    → High variance between variants
    → Positions have different synthesis quality
    → Report heterogeneous class
    → Flag for investigation
```

**Assumption**: Chromatographic behavior independent of synthesis position

**Physical Basis**: Same molecule, same retention (hydrophobicity determined by sequence, not synthesis path)

**When Assumption Fails**:

- Positional effects on folding/conformation
- Incomplete deprotection (position-dependent)
- Tag-peptide interactions vary by position
- In these cases: Individual mode mandatory

#### 4.2.8.1 Operational Fallback Workflow

**What Happens When Validation Fails**:

**During Processing** (automated):

```
For each equivalence class:

  Step 1: Attempt consensus mode
    1a. Aggregate signals → consensus
    1b. Run correlation check (Check 1)

    If min(correlation) < 0.8:
      → Log warning: "EquivalenceClass [sequence] failed correlation check"
      → Automatically fall back to individual mode
      → Flag class as "CONSENSUS_INVALID" in metadata
      → Proceed with Step 2

  Step 2: Individual mode processing
    For each variant in class:
      - Peak detection (on raw signals)
      - Classification
      - Purity calculation
      - Validation

  Step 3: Aggregate individual results
    - Compute mean_purity, std_purity
    - Run Check 2 (peak position variance)
    - Run Check 3 (purity variance)
    - Assign class-level status
```

**Class-Level Status**:

```
If consensus mode succeeded:
  status = "CONSENSUS_VALID"

If fell back to individual mode:
  If Check 2 or Check 3 failed:
    status = "HETEROGENEOUS" (variants have different behavior)
  Else:
    status = "CONSENSUS_INVALID_BUT_SIMILAR" (low correlation but stable)
```

**Output Flags**:

```python
class_metadata = {
    "consensus_attempted": True,
    "consensus_valid": False,  # Fell back
    "correlation_min": 0.67,  # Below threshold
    "fallback_reason": "Low signal correlation between positional variants",
    "mode_used": "individual",
    "heterogeneous": False,  # Checks 2&3 passed
    "recommendation": "Review positional variants - may have position-specific effects"
}
```

**User Reporting**:

In summary output:
```
EquivalenceClass: Val (3 variants)
  Mode: Individual (consensus invalid - correlation 0.67 < 0.8)
  Mean purity: 0.85 ± 0.08
  Recommendation: Review variants - possible position-dependent effects
```

**Key Principle**: Fallback is automatic and transparent. User sees which classes used consensus vs individual mode and why.

#### 4.2.9 Individual vs Consensus Mode Comparison

**Individual Mode** (default, always valid):

**Approach**: Analyze each positional variant independently

```
For each of N_variants:
  1. Peak detection (Discrete Morse + Poisson + Prominence)
  2. Peak classification (DAG constraints)
  4. Area integration
  5. Purity calculation
  6. Validation

Output: N_variants individual results
```

**Advantages**:
✅ Always valid (no assumptions about variant similarity)
✅ Full position-specific information
✅ No validity checks required
✅ Simpler workflow
✅ Regulatory-compliant (actual measurements)

**Disadvantages**:
❌ Slower (process each variant through expensive pipeline)
❌ More results to review (N_variants instead of N_classes)
❌ Noise not averaged (lower SNR per variant)

**When to use**: Default choice, synthesis diagnostics, regulatory requirements

---

**Consensus Mode** (optional optimization, hybrid approach):

**Approach**: Detect peaks on consensus, quantify on individuals

```
For each of N_classes equivalence classes:
  1. Aggregate variants → consensus signal (cheap)
  2. Peak detection on consensus (expensive, done once)
  3. Peak classification on consensus
  5. For each variant:
     - Area integration (cheap, per-variant)
     - Purity calculation (per-variant)
     - Validation (per-variant)

Output: N_variants individual results + N_classes summaries
```

**Advantages**:
✅ Faster (~3-10× speedup for peak detection)
✅ Noise reduction (consensus has higher SNR)
✅ Molecular-level summaries (mean_purity, class validation)
✅ Still get individual purities (validity maintained)
✅ Position-specific diagnostics available

**Disadvantages**:
❌ Requires variant similarity (correlation check)
❌ More complex workflow
❌ Additional validation steps required
❌ Consensus peak boundaries may not be optimal for all variants

**When to use**: Large datasets, similar variants, molecular-level analysis, speedup needed

---

**Recommendation**:

**Start with Individual Mode**: Simpler, always valid, complete information

**Consider Consensus Mode** if:

1. Dataset is large (many positional variants)
2. Speedup is needed (peak detection bottleneck)
3. Variants validated as similar (correlation > 0.8)
4. Molecular focus appropriate (chemical identity over synthesis path)

**Best Practice**: Run correlation check on sample of classes before committing to consensus mode for full dataset

### 4.3 Truncation Hierarchy Structure

**CompoundHierarchy** (formerly TruncationHierarchy):

- Maps levels to compounds: {level: [compounds]}
- Maps compounds to direct descendants: {compound: [descendants]}
- Decomposition type: "block" or "monomer"

**Direct Descendants:**

- A compound C is a **direct descendant** of P if:
  1. C is a descendant of P (C ⊂ P)
  2. No other compound B exists where C ⊂ B ⊂ P

**Properties:**

- Direct descendants = edges in DAG
- All descendants = transitive closure
- Hierarchy can have gaps (level 9 → level 7, skipping level 8)

---

# II. ANALYSIS METHODS

## PART 5: PEAK DETECTION MATHEMATICAL FOUNDATIONS

### 5.0 Signal Properties and Preprocessing

#### 5.0.1 LC-Seq Signal Characteristics

LC-Seq chromatograms have unique properties derived from the fractionation and DNA sequencing workflow:

**Data Type: Discrete Fraction Counts (Pre-Scaled)**

- Library injected onto LC column and separated into **96 discrete fractions** (~30-37 second bins)
- Each fraction PCR-amplified and sequenced via NGS
- Raw signal = **DNA barcode sequencing counts per fraction**
- **Working signal = scaled counts** (normalized for sequencing depth, UMI deduplication, amplification bias)
- Can be fractional after scaling/normalization
- Underlying distribution: Poisson (raw counts), approximately Poisson-like after scaling

**Spatial Resolution**

- **Pre-binned data**: ~96 fractions over ~3000 seconds = 30s/fraction
- No sub-fraction resolution (discrete time points)
- Adjacent fractions independent samples (no interpolation needed)
- Molecular diffusion already averaged within fraction collection

**Count Statistics**

- Typical experiment: 3 × 10¹³ molecules injected, 5 × 10⁸ total sequencing reads
- Average: 5.3 × 10⁶ counts per fraction (good representation)
- Peak maximum: ~250 counts (10× above baseline)
- Baseline: ~25 counts (background)
- Noise: σ ≈ √c (Poisson property)

**Signal-to-Noise Ratio (SNR)**

- Variable across compounds (high for abundant, low for rare)
- Abundant compounds: Clear peaks, SNR >10
- Rare compounds: Noisy signals, SNR <3
- Detection limit: SNR >3 (statistical significance)
- Quantitation limit: SNR >10 (reproducible measurement)

---

### 5.1 Discrete Morse Theory Framework

This section establishes the rigorous mathematical framework for peak detection in discrete count data.

**Peak detection for LC-Seq is fundamentally about finding local maxima in discrete sequences.**

#### Mathematical Setup

An LC-Seq chromatogram is a **discrete sequence** c = {c₁, c₂, ..., c₉₆} where:

- cᵢ = scaled sequencing counts in fraction i
- Time points tᵢ correspond to fraction midpoints
- No smoothing applied (data already pre-binned)

#### Discrete Morse Theory

For a discrete sequence, a **local maximum** at index i is defined as:

```
c[i] > c[i-1]  AND  c[i] ≥ c[i+1]
```

**Properties:**
- Well-defined for discrete data (no derivatives needed)
- Complete: Finds ALL local maxima
- Mathematically rigorous: Extension of continuous Morse theory to discrete spaces
- Handles ties consistently (≥ on right allows plateau detection)

**Morse Index (Discrete):**
- Index 1: Local maximum (peak)
- Index 0: Local minimum (valley)

#### Why Discrete Morse Theory?

Discrete Morse theory is the **natural mathematical language** for LC-Seq:

- **Appropriate**: Discrete data requires discrete mathematics
- **Complete**: Finds all critical points
- **No smoothing**: Preserves chromatographic resolution
- **Rigorous**: Based on discrete topology (Forman 1998)

#### Algorithm

```
For each index i in [2, n-1]:
    If c[i] > c[i-1] AND c[i] ≥ c[i+1]:
        peak[i] = True
```

**Advantages:**
- **Direct**: No numerical derivatives
- **Fast**: O(n) single pass
- **Exact**: No approximation or smoothing
- **Complete**: Guaranteed to find all maxima

### 5.2 Statistical Significance Testing for Peak Detection

**Key Question**: Which peaks are "real signal" vs "statistical noise"?

#### Poisson Count Statistics

LC-Seq data follows Poisson-like statistics (after scaling):

**Distribution**: c[i] ~ Poisson-like with variance proportional to mean

**Noise model**: σ[i] ≈ √(c[i] + ε) where ε prevents division by zero at low counts

**Background estimation**: μ_bg = percentile(all counts, 10) - captures low-count baseline

#### Statistical Hypothesis Testing

For each detected local maximum at position i:

**Null Hypothesis H₀**: Peak is random Poisson fluctuation of background
**Alternative H₁**: Peak is real signal above background

**Test statistic**:
```
Z = (c[i] - μ_bg) / √(μ_bg + ε)
```

**Decision rule**:
- Z > 3: Reject H₀ (detection threshold, p < 0.001)
- Z > 10: High confidence (quantitation threshold)

#### Prominence (Local Significance)

Statistical testing alone is insufficient - must also measure peak height relative to **local baseline**:

**Definition**: Prominence measures how much a peak rises above surrounding valleys

```
prominence[i] = c[i] - max(valley_left[i], valley_right[i])

where:
- valley_left = min(c[j]) for j ∈ (prev_peak, i)
- valley_right = min(c[k]) for k ∈ (i, next_peak)
```

**Interpretation**:
- High prominence = peak well-separated from neighbors (distinct chemical entity)
- Low prominence = shoulder or unresolved peak (may be real but overlapping)

#### Why Prominence Instead of Persistence?

**Persistent Homology Assumes** (incorrect for LC-Seq):
- Multiple close peaks = noise/over-segmentation
- Smoothing reveals "true" structure
- Broad stable features = significant

**LC-Seq Reality**:
- Multiple close peaks = multiple chemical entities (product + truncations)
- Each peak = distinct synthesis outcome
- Prominence captures chromatographic significance without destroying resolution

**Prominence Properties**:
- Scale-invariant (no smoothing parameter)
- Respects valley separation (natural chemistry)
- Fast to compute (O(n) single pass)
- Standard in analytical chemistry

#### Adaptive Threshold (No Magic Numbers!)

Filter peaks using data-derived prominence threshold:

**Option 1: Percentile-based** - Compute percentile threshold (e.g., 20th percentile retains top 80% most prominent peaks)

**Option 2: Gap-based** - Sort prominence values in descending order, find largest gap between consecutive values, use the value after the gap as threshold. This identifies natural separation between signal and noise clusters.

#### Peak Boundary Determination

**After detecting significant peaks**, we need to define boundaries [t_start, t_end] for area integration.

**Algorithm**:

**Step 1: Find Peak Maximum**

```
t_peak = position of critical point (local maximum from Morse theory)
height_peak = signal(t_peak)
```

**Step 2: Find Left Boundary (t_start)**

```
Starting from t_peak, scan leftward (decreasing time):
  For each point t_i moving left from t_peak:
    If signal(t_i) < threshold_fraction × height_peak:
      t_start = t_i
      break
    If reach valley (local minimum where f''(t) > 0):
      t_start = t_valley
      break
    If reach signal start:
      t_start = t_min
      break

threshold_fraction = 0.05 (5% of peak height)
```

**Step 3: Find Right Boundary (t_end)**

```
Starting from t_peak, scan rightward (increasing time):
  For each point t_i moving right from t_peak:
    If signal(t_i) < threshold_fraction × height_peak:
      t_end = t_i
      break
    If reach valley (local minimum):
      t_end = t_valley
      break
    If reach signal end:
      t_end = t_max
      break
```

**Step 4: Store Peak with Boundaries**

```
Peak:
  position = t_peak
  height = height_peak
  left_boundary = t_start
  right_boundary = t_end
  prominence = height_peak - max(valley_left, valley_right)
  z_score = (height_peak - μ_bg) / √(μ_bg + ε)
```

**Rationale**:

- **Valley detection**: Natural separation between adjacent peaks (best case)
- **Threshold-based**: When peaks overlap or no clear valley (height-based cutoff)
- **Signal edges**: Boundary case handling (peaks at start/end of signal)
- **5% threshold**: Conservative (captures full peak including tails)

**Use in Integration**:
Peak area = Σ corrected_signal(t) for t ∈ [t_start, t_end]

This provides the peak boundaries needed for purity calculation (Part 5.0.7) and consensus mode area integration (Part 4.2.6).

### 5.3 Peak Type Classification

#### Updated Peak Types

Peak classification is based on topological position in the lineage DAG:

**IMPORTANT**: These labels describe POSITIONAL classification, NOT chemical identity or synthesis validation!

**Peak Types:**

- **NULL**: Peak at L₀ (minimal element) position
- **TRUNCATION**: Matches ancestor product position or null position
- **PUTATIVE_PRODUCT**: Positionally consistent with expected product (NOT chemically validated!)
- **UNKNOWN**: Peak that doesn't match any expected position (see [Part 5.8](#58-classification-limitations-and-unknown-peaks))

#### Classification Logic

Sequential classification based on position:

**Priority 0: NULL Peak**

- Position: Matches L₀ (full-null compound) global maximum
- Significance: Universal truncation reference point
- Cardinality: Exactly 1 at L₀; 0 or 1 in all other compounds
- Constraint: All compounds may have peak at this position

**Priority 1: TRUNCATION Peaks**

- Position: Matches ancestor product positions OR null position
- Sources: (1) Specific ancestor products, (2) Universal null position
- Cardinality: 0 or more per compound
- Constraint: position < putative_product position

**Priority 2: PUTATIVE_PRODUCT Peak**

- Position: First significant peak after truncations
- Significance: **Positionally consistent** with product elution
- Cardinality: 0 or 1 per compound
- **Critical**: NOT validated as pure product! (see [Part 6](#part-6-synthesis-validation-theory))
- Constraint: position > max(truncation positions)

**Priority 3: UNKNOWN Peaks**

- Position: Peaks not matching any expected position
- Includes: Late-eluting peaks after product, unmatched early peaks, ambiguous assignments
- Cardinality: 0 or more per compound
- Cannot identify chemical identity without orthogonal data (MS, NMR)

#### Example Classification

For detected peaks at positions [10, 15, 25, 35, 50] with NULL position 10 and ancestor product positions [15, 25]:

- Peak at position 10 is classified as NULL (matches L₀)
- Peaks at positions 15 and 25 are classified as TRUNCATION (match ancestors)
- Peak at position 35 is classified as PUTATIVE_PRODUCT (first significant peak after truncations)
- Peak at position 50 is classified as UNKNOWN (no expected position match - could be oligomer, contaminant, etc.)

### 5.4 Global Classification via Constraint Propagation

**Key Insight**: Peak classification must respect the **entire lineage DAG**, not just individual compounds.

#### Why Global Classification?

**Problem with local classification:**

- Each compound classified independently
- Ignores relationships between compounds
- Inconsistent across lineage

**Solution with global classification:**

- Process DAG in topological order
- Propagate constraints through edges
- Ensure consistency across lineage

#### Algorithm: Bottom-Up Propagation

**Step 1**: Find and process L₀ (minimal element) FIRST. Detect NULL peak (global maximum at L₀).

**Step 2**: Extract NULL peak position as universal constraint that applies to ALL compounds.

**Step 3**: Perform topological sort of remaining compounds (bottom-up: minimal → maximal elements).

**Step 4**: Process each compound in topological order:

- Retrieve detected peaks for this compound
- Get constraints from already-processed descendants (descendant product positions)
- Build expected truncation set: descendant products + NULL position
- Classify peaks using constraints and ordering rules
- Store peak labels for use as constraints by ancestors

#### Constraint Types

**1. NULL Constraint** (universal): All compounds may have peak at position(L₀_peak)

**2. Lineage Constraint** (ancestor-descendant): If descendant d has product at position p_d, then ancestor a may have truncation peak at position p_d

**3. Ordering Constraint** (within compound): position(truncation) < position(putative_product) < position(unknown)

**4. Downstream Constraint** (hierarchical): If a is ancestor of d, then position(product_a) > position(product_d)

### 5.5 Optimal Assignment and Matching

#### The Matching Problem

**Given:**

- Detected peaks D = {d₁, d₂, d₃, ...}
- Expected positions E = {e₁, e₂, e₃, ...} (from descendants)

**Find:** Optimal assignment D → E minimizing total cost

This is the **Linear Assignment Problem** (LAP).

#### Cost Function (Scale-Invariant!)

The cost matrix is computed by calculating positional distance between each detected peak and each expected position, normalized by the total signal length to achieve scale invariance. Cost(i,j) = |peak_i.position - expected_j| / signal_length.

#### Hungarian Algorithm

The optimal assignment is found using the Hungarian algorithm (linear sum assignment) which solves the assignment problem in O(n³) time. An adaptive acceptance threshold is computed as the median peak spacing divided by signal length (avoiding magic numbers). Only assignments with cost below this threshold are accepted, ensuring assigned peaks are genuinely close to expected positions.

**Key advantages:**

- **Optimal**: Minimizes total assignment cost (not greedy)
- **Adaptive threshold**: Based on peak spacing, not fixed seconds
- **Scale-invariant**: Normalized by signal characteristics

### 5.6 Classification Scope and Limitations

**CRITICAL DISTINCTION: Peak Classification ≠ Synthesis Validation**

#### What Peak Classification CAN Determine

✅ **Positional Consistency**

- Peak location relative to expected elution order
- Satisfies constraints from DAG structure
- Matches expected retention time pattern

✅ **Statistical Significance**

- Peak Z-score for Poisson counts (real vs noise)
- High prominence (not baseline fluctuation)
- Mathematically well-defined feature

✅ **Relational Constraints**

- Peak ordering within compound (truncation < product < oligomer)
- Peak relationships across lineage (descendant → ancestor)
- Hierarchical consistency in DAG

✅ **Hypothesis Ranking**

- Which peak is **most likely** putative product (by position)
- Confidence scores based on constraint satisfaction
- Alternative interpretations

#### What Peak Classification CANNOT Determine

❌ **Chemical Identity**

- Is this the intended molecule? (Unknown!)
- Could be correct product, truncation mixture, modified product, contaminant
- Requires mass spectrometry, NMR, or other orthogonal methods

❌ **Purity**

- Is this a single compound or mixture? (Unknown!)
- Co-eluting compounds appear as single peak
- Requires chromatographic resolution analysis or spectroscopy

❌ **Synthesis Success**

- Did the synthesis reaction work? (Unknown!)
- Peak presence ≠ successful synthesis
- Requires chemical validation (MS, NMR, standards)

❌ **Quantitation**

- How much product was made? (Unknown!)
- Peak area ≠ absolute quantity without calibration
- Requires quantitative standards and validated methods

#### Putative Product: What It Means

**A peak labeled "PUTATIVE_PRODUCT" indicates:**

✓ Positionally consistent with expected product elution
✓ Appears after all known truncation positions
✓ Satisfies hierarchical constraints from DAG
✓ Statistically significant (Z > 3, high prominence)
✓ **Hypothesis** that this _might_ be the product

✗ **NOT** confirmed as pure product
✗ **NOT** validated synthesis success
✗ **NOT** chemical identity confirmed
✗ May be mixture, modified product, or contaminant

**Synthesis validation requires additional evidence:**

- Mass spectrometry (molecular weight confirmation)
- NMR spectroscopy (structure confirmation)
- Chromatographic purity analysis (single compound)
- Comparison to authentic standard (identity confirmation)

#### Appropriate Use Cases

**This algorithm IS suitable for:**

- **Library screening**: Identify candidate compounds for follow-up
- **Comparative analysis**: Which compounds elute where in family
- **Quality control**: Detect anomalies in elution patterns
- **Hypothesis generation**: What might be product (requires validation)
- **Process monitoring**: Track consistency across batches

**This algorithm is NOT suitable for:**

- **Definitive synthesis validation**: Requires MS/NMR/standards
- **Regulatory compliance**: Requires validated analytical methods
- **Quantitative yield determination**: Requires calibration curves
- **Chemical identity proof**: Requires spectroscopic confirmation
- **Publication claims**: "Synthesis successful" requires orthogonal proof

### 5.7 Mathematical Algorithm Summary

#### Three-Stage Pipeline

**Stage 1: Local Detection** (Discrete Morse Theory + Poisson Statistics)

- Input: Discrete fraction counts {c₁, c₂, ..., c₉₆}
- Process:
  1. Find local maxima (discrete Morse theory) - O(n)
  2. Find peak boundaries (valley detection) - O(n)
  3. Compute prominence (height above local baseline) - O(n)
  4. Filter by Poisson significance (Z-score > 3) - O(n)
  5. Filter by prominence (percentile or gap-based) - O(n log n)
- Output: All significant peaks (position, prominence, height, boundaries)

**Stage 2: Global Classification** (Constraint Propagation on DAG)

- Input: Detected peaks + Lineage DAG
- Process:
  1. Process L₀ → extract NULL position
  2. Topological sort of DAG (bottom-up)
  3. For each compound (in order):
     - Get descendant constraints
     - Match peaks to expected positions (Hungarian)
     - Classify: NULL, TRUNCATION, PUTATIVE_PRODUCT, UNKNOWN
  4. Propagate product positions to ancestors
- Output: Peak labels for entire lineage

**Stage 3: Validation** (Constraint Satisfaction Check)

- Input: Peak labels + DAG constraints
- Process:
  1. Verify ordering constraints (truncation < product, unknown after product)
  2. Verify cardinality (≤1 putative product per compound)
  3. Verify downstream constraint (ancestor product > descendant product)
  4. Verify lineage constraint (descendant product → ancestor truncation)
- Output: Validated classifications + confidence scores + caveats

#### Computational Complexity

- **Local maxima detection**: O(n) for signal of length n
- **Prominence computation**: O(n) single pass
- **Statistical filtering**: O(n) per peak
- **Topological sort**: O(V + E) for DAG with V vertices, E edges
- **Hungarian algorithm**: O(p³) for p peaks per compound
- **Total pipeline**: O(V × p³) dominated by peak matching (no multi-scale overhead)

#### Parameters (Uniform Across Dataset!)

**CRITICAL**: All parameters must be **uniform across the entire dataset** for comparability.

**No Smoothing Parameters Needed:**

Since LC-Seq data is pre-binned into fractions, no smoothing is applied. Data used as-is.

**Statistical Significance Parameters:**

```
Poisson significance threshold:
  Z_threshold = 3.0  (3σ detection, p < 0.001)

Background estimation:
  μ_bg = percentile(all_counts, 10)  (10th percentile)
  ε = 1.0  (regularization for low counts)
```

**Rationale**:
- ✅ **Standard**: 3σ is universal detection threshold in analytical chemistry
- ✅ **Data-driven**: Background from actual count distribution
- ✅ **Conservative**: 10th percentile captures baseline without signal contamination

**Prominence Threshold:**

Adaptive, data-derived approach (no magic numbers):

```
Option 1: Percentile-based
  threshold = percentile(all_prominences, 20)
  (Retains top 80% most prominent peaks)

Option 2: Gap-based
  threshold = value_after_largest_gap(sorted_prominences)
  (Identifies natural separation between signal and noise)
```

**Rationale**:
- ✅ **Adaptive**: Threshold derived from data distribution
- ✅ **Relative**: Works across different signal intensities
- ✅ **No magic numbers**: Data determines cutoff

**Position Tolerance (for classification):**

Adaptive based on peak spacing:

```
tolerance = median(peak_spacings) / signal_length

where peak_spacings computed from sorted expected positions
```

**Rationale**:
- ✅ **Scale-invariant**: Normalized by signal length
- ✅ **Data-driven**: Based on actual peak separations
- ✅ **Chemical**: Reflects typical retention time variability

**Key Principle**: Uniformity > per-compound optimization.

### 5.8 Classification Limitations and Unknown Peaks

#### 5.8.1 What Can Be Classified

Peak classification based on chromatographic retention time and DAG constraints can identify:

**NULL Peak** (High Confidence):

- Peak at retention time t(L₀_peak)
- Definitive identification based on L₀ reference
- Chemical meaning: DNA tag only, complete truncation

**TRUNCATION Peaks** (High Confidence):

- Peaks matching ancestor product positions
- Peaks matching L₀ (null) position
- Confident identification via positional constraints from DAG

**PUTATIVE_PRODUCT Peak** (Positional Hypothesis):

- Peak positionally consistent with expected product elution
- First significant peak after all known truncation positions
- **NOT chemical confirmation** - positional hypothesis only
- See Part 5.6 for scope and limitations

**Basis for Classification**:

- Absolute retention time (scalar values)
- DAG constraint propagation (ancestor/descendant relationships)
- Statistical significance (Poisson Z-score, prominence)
- Ordering constraints (truncation < product)

#### 5.8.2 What Cannot Be Classified (UNKNOWN)

**Peaks labeled UNKNOWN** include:

**Late-Eluting Peaks**:

- Peaks appearing after putative product position
- Retention time > t(putative_product)
- Could be multiple distinct species

**Unmatched Peaks**:

- Peaks that don't match any expected position
- Outside tolerance windows for truncations
- No ancestor constraint satisfaction

**Ambiguous Assignments**:

- Multiple candidate peaks for same classification
- Position matches multiple hypotheses
- Insufficient confidence to assign label

**Fundamental Limitation**: Without orthogonal analytical data (mass spectrometry, NMR, etc.), chemical identity cannot be determined from retention time alone.

#### 5.8.3 Possible Identities of Unknown Peaks

**Late-eluting peaks could be**:

**Oligomers** (n-mers):

- Dimers (2 copies of compound): Expected retention ≈ 2 × t(monomer) - t(L₀)
- Trimers (3 copies): Expected retention ≈ 3 × t(monomer) - 2 × t(L₀)
- Higher oligomers: n-mers following hydrophobicity additivity
- Formation: Aggregation, cross-linking, non-covalent association

**Contaminants**:

- Synthesis reagents (unreacted building blocks, coupling agents)
- Degradation products (hydrolysis, oxidation)
- Column bleed (stationary phase breakdown)
- Carry-over from previous samples

**Modified Products**:

- Incomplete reactions (partial deprotection, incomplete coupling)
- Side reactions (epimerization, racemization, unwanted cyclization)
- Post-synthesis modifications (oxidation, hydrolysis during storage)

**Artifacts**:

- Signal distortions
- Ghost peaks from instrumentation
- Air bubbles or solvent effects

**Cannot distinguish without orthogonal data**: Mass spectrometry (molecular weight), NMR (structure), comparison to authentic standards.

#### 5.8.4 Handling Unknown Peaks in Analysis

**Classification Strategy**:

- Label all unidentifiable peaks as UNKNOWN
- Honest acknowledgment of analytical limitations
- No speculation about chemical identity without evidence

**Purity Calculation**:

All non-product peaks reduce purity, regardless of identity:

Purity(C) = Σ(counts_putative_product) / [Σ(counts_putative_product) + Σ(counts_truncation) + Σ(counts_unknown) + Σ(counts_null)]

**Key Point**: Unknown peaks count as impurities. Whether a late peak is an oligomer, contaminant, or artifact doesn't matter for purity assessment - it's NOT the intended product.

**Synthesis Validation Impact**:

- High UNKNOWN fraction → reduces purity → affects validation category
- VALIDATED requires purity > P₇₅ (includes all impurities)
- Unknown peaks weighted equally with truncations in purity calculation

**Reporting**:

- Report presence of UNKNOWN peaks
- Report fraction of total signal from unknowns
- Flag compounds with high UNKNOWN content (>20% of signal)
- Recommend follow-up analysis if unknowns dominate

#### 5.8.5 Oligomer Hypothesis Generation (Not Classification)

While oligomers cannot be definitively identified from retention time, we can generate hypotheses for follow-up investigation:

**Pattern-Based Oligomer Hypothesis**:

**Evidence suggesting oligomerization**:

1. **Position matching**: Late peak at retention ≈ 2 × t(product) - t(L₀)
2. **Ladder pattern**: Multiple peaks at regular intervals (t, 2t, 3t)
3. **Intensity decay**: Peak heights decrease geometrically (M >> M₂ >> M₃)
4. **Family consistency**: Same pattern across related compounds

**Hypothesis Strength**:

- Single late peak at ~2× position: WEAK hypothesis (could be coincidence)
- Ladder pattern (2×, 3×, 4×): MODERATE hypothesis (systematic behavior)
- Ladder + intensity decay + family-wide: STRONG hypothesis (likely oligomers)

**Use Case**:

- Flag compounds for follow-up MS analysis
- Identify systematic aggregation issues in library
- Guide synthesis/purification optimization
- **NOT for automated classification** (remains UNKNOWN)

**Recommended Workflow**:

1. Classify late peaks as UNKNOWN (conservative)
2. Generate oligomer hypothesis scores (informative)
3. Report hypothesis to user
4. If hypothesis strong → recommend MS confirmation
5. After MS analysis → update labels with confident identification

#### 5.8.6 Quality Metrics Based on Unknown Peaks

**Per-Compound Metrics**:

Unknown fraction:
unknown_fraction(C) = Σ(counts_unknown) / Σ(counts_all_peaks)

Categories:

- Low unknown: <5% (clean synthesis)
- Moderate unknown: 5-20% (acceptable)
- High unknown: >20% (quality concern)

**Dataset-Wide Metrics**:

Median unknown fraction:
median_unknown = median(unknown_fraction_i for all compounds)

Compounds with high unknowns:
high_unknown_rate = count(unknown_fraction > 0.2) / total_compounds

**Quality Flags**:

- If median_unknown > 0.15: WARNING "Dataset has high unknown peak content - review synthesis/purification"
- If high_unknown_rate > 0.3: WARNING "30% of compounds have high unknown peaks - systematic issue suspected"

**Interpretation**:

- Moderate unknowns: Normal (some aggregation/impurities expected)
- High unknowns dataset-wide: Synthesis problem, purification issue, or chromatography problem

#### 5.8.7 When Definitive Identification IS Possible

**Mass Spectrometry (MS)**:

- Confirms molecular weight
- Distinguishes oligomers (2×, 3× mass) from contaminants
- Identifies modified products (mass shifts)
- **After MS**: Update UNKNOWN → OLIGOMER_DIMER (confident label)

**NMR Spectroscopy**:

- Confirms chemical structure
- Identifies side products, impurities
- Verifies intended product structure

**Authentic Standards**:

- Synthesize reference compounds
- Match retention times with high precision
- Confirm peak identity by co-elution

**Orthogonal Chromatography**:

- Different separation mode (ion exchange, size exclusion)
- Complementary retention mechanism
- Resolves co-eluting species

**With Orthogonal Data Available**:

- Reclassify UNKNOWN peaks with confident labels
- Update purity calculations (same formula, more informative labels)
- Synthesis validation unchanged (impurity is impurity regardless of identity)
- More actionable for synthesis optimization (know what to fix)

#### 5.8.8 Classification Decision Summary

**Decision Tree**:

Question: Does peak match L₀ position?

- YES → NULL

Question: Does peak match ancestor product or null position?

- YES → TRUNCATION

Question: Is peak first significant peak after truncations?

- YES → PUTATIVE_PRODUCT

Otherwise:

- → UNKNOWN

**Conservative Principle**: When in doubt, label UNKNOWN. Better to acknowledge uncertainty than to speculate without evidence.

**Purity Impact**: All classifications except PUTATIVE_PRODUCT reduce purity (UNKNOWN, TRUNCATION, NULL all count as impurities).

---

## PART 6: SYNTHESIS VALIDATION THEORY

This section establishes a rigorous, data-driven framework for validating synthesis success using DNA-encoded library sequencing data, DAG structure constraints, and chromatographic physics.

### 6.1 DNA-Encoded Library Context

**Signal Type:** DNA barcode scaled counts vs elution time (NOT direct chemical measurement)

**Critical Assumptions:**

- Scaled counts normalize for sequencing depth and amplification bias
- Relative abundances are comparable across compounds
- DNA tag tracks compound through chromatography
- Tag-compound linkage remains intact during analysis

**Key Limitation:** DNA tag abundance ≠ chemical purity

- DNA can survive truncation
- PCR amplification introduces bias
- Sequencing is statistical sampling
- Tag could dissociate from compound

**Primary Validation Signal:** Chromatographic physics (retention time order) is MORE reliable than scaled counts alone.

### 6.2 Adaptive Validation Principle

**Core Principle:** All validation metrics must be **dataset-relative** and **scale-invariant** to accommodate:

- Different sequencing depths
- Different scaling methods
- Different library quality levels
- Different experimental conditions

**NO fixed thresholds.** All parameters derived from data distribution.

**Approach:**

1. Learn dataset distribution (bootstrap phase)
2. Define adaptive thresholds (percentiles, MAD)
3. Classify relative to dataset characteristics
4. Adjust stringency based on data quality

### 6.3 Purity Definition and Calculation

**Purity for Compound C:**

Purity(C) = Σ(scaled_counts_product) / [Σ(scaled_counts_product) + Σ(scaled_counts_truncations) + Σ(scaled_counts_unknowns) + Σ(scaled_counts_null)]

Where summation is over all elution fractions.

**Interpretation:**

- Purity = 1.0 → Only product peak present
- Purity = 0.5 → Product and impurities equally abundant
- Purity → 0 → Dominated by truncations/unknowns

**Statistical Uncertainty:**

Standard error of purity estimate:
SE(purity) = √[purity × (1-purity) / total_scaled_counts]

95% Confidence Interval:
CI = purity ± 1.96 × SE(purity)

**Minimum count threshold:**
total_scaled_counts > 100 for CI width < 0.2

### 6.4 Distribution-Based Thresholds

For dataset D with compounds C₁, C₂, ..., Cₙ:

**Step 1: Characterize Dataset**

1. Compute purity(Cᵢ) for all compounds
2. Extract percentiles: P₁₀, P₂₅, P₅₀, P₇₅, P₉₀, P₉₅
3. Compute Median Absolute Deviation: MAD = median(|purity - P₅₀|)
4. Estimate background from L₀ or low-count tail

**Step 2: Define Adaptive Categories**

- **Exceptional purity**: purity > P₉₀ (top 10%)
- **High purity**: P₇₅ < purity ≤ P₉₀ (75th-90th percentile)
- **Moderate purity**: P₅₀ < purity ≤ P₇₅ (50th-75th percentile)
- **Low purity**: P₂₅ < purity ≤ P₅₀ (25th-50th percentile)
- **Very low purity**: purity ≤ P₂₅ (bottom 25%)

**Step 3: Adjust for Dataset Quality**

If MAD(purity) < 0.1: # High-quality library
→ Use strict thresholds (P₇₅ for validation)

If MAD(purity) > 0.2: # Variable library
→ Use lenient thresholds (P₅₀ for validation)

This ensures fair evaluation regardless of library-wide synthesis quality.

### 6.5 Signal-to-Noise Ratio (Universal Metric)

**Background Estimation:**

**Option 1:** From L₀ (full-null compound)
background = median(scaled_counts) across L₀ signal

**Option 2:** From low-count tail
background = median(bottom 10% of all scaled_count values)

**Signal-to-Noise Ratio:**

SNR(C) = max(scaled_counts_product_C) / background

**Interpretation:**

- SNR > 10: High confidence detection (clear signal)
- 3 ≤ SNR ≤ 10: Moderate confidence (detectable)
- SNR < 3: Near noise floor (unreliable)

**Why SNR is universal:**

- Scale-invariant (ratio eliminates absolute scale)
- Works regardless of count magnitudes
- Independent of sequencing depth
- Robust to experimental variation

### 6.6 Retention Time Constraints (Chromatographic Physics)

**Physical Law:** Hydrophobicity additivity

For compound C (n building blocks) and truncation T (n-1 building blocks):
retention_time(C) > retention_time(T)

This is a **hard constraint** - violation implies synthesis failed OR peak assignment incorrect.

**Adaptive Retention Precision:**

Learn minimum resolvable time difference from data:
Δt*min = min(retention_time*{i+1} - retention_time_i) across all peaks
retention_precision = Δt_min / 2

**Retention Order Validation:**

For confident ordering:
t_product - t_truncation > 2 × retention_precision

If difference < precision → ambiguous (cannot distinguish)

**Why retention time > counts:**

- DNA tags can lie (dissociation, amplification bias)
- Physics doesn't lie (longer peptides elute later)
- Retention violations are impossible if chemistry is correct

### 6.7 Bayesian Validation Framework

**Synthesis Success Probability:**

P(synthesis_succeeded | purity, retention_order, descendants) ∝
P(purity | synthesis_succeeded) ×
P(retention_order | synthesis_succeeded) ×
P(descendants | synthesis_succeeded) ×
P(synthesis_succeeded)

**Likelihood Functions:**

**Purity likelihood:**

- P(purity=p | succeeded) ~ Beta(α=19, β=1) [mode ≈ 0.95]
- P(purity=p | failed) ~ Beta(α=2, β=8) [mode ≈ 0.20]

**Retention order likelihood:**

- P(order_correct | succeeded) = 0.95
- P(order_correct | failed) = 0.05

**Descendant evidence:**

- P(descendants_validated | succeeded) = 0.90^n for n descendants
- P(descendants_validated | failed) = 0.10^n

**Prior probability:**

- P(synthesis_succeeded) = dataset_success_rate (learned from data)

### 6.8 DAG Constraint Propagation

**Processing Order:** Topological sort (L₀ first, then bottom-up)

**Algorithm:**

**Step 1:** Process L₀ (minimal element)

- Detect NULL peak (global maximum)
- Estimate background noise
- Set universal NULL constraint

**Step 2:** For each compound in topological order:

1. Retrieve detected peaks
2. Calculate purity from scaled counts
3. Calculate SNR relative to background
4. Check retention_time(C) > retention_time(all descendants)
5. Get descendant validation status
6. Update P(C succeeded | all_evidence) using Bayes rule
7. Propagate C's validation status to ancestors

**Propagation Rules:**

If ALL descendants validated → strongly increases P(C succeeded)
If ANY descendant failed → decreases P(C succeeded)
If retention order violated → P(C succeeded) ≈ 0 (override all other evidence)

### 6.9 Robust Statistical Methods

**Median Absolute Deviation (MAD):**

MAD = median(|purity(Cᵢ) - median_purity|)

**Why MAD over standard deviation:**

- Robust to outliers (unlike σ)
- Works with skewed distributions
- No normality assumption
- More appropriate for small datasets

**Outlier Detection:**

Exceptionally clean: purity > median + 2×MAD
Exceptionally poor: purity < median - 2×MAD

### 6.10 Validation Classification

**Decision Framework:**

**VALIDATED** (synthesis confirmed)

- ✅ Retention time order correct (Δt > 2×precision)
- ✅ Purity > P₇₅ (dataset 75th percentile)
- ✅ SNR > 5
- ✅ All descendants validated
- ✅ Confidence interval excludes low purity
- **Confidence: Very High (>95%)**

**LIKELY_SUCCESS** (high confidence)

- ✅ Retention time order correct
- ✅ Purity > P₅₀ (dataset median)
- ✅ SNR > 3
- ✅ Majority (>50%) descendants validated
- **Confidence: High (80-95%)**

**UNCERTAIN** (ambiguous)

- ⚠️ P₂₅ < Purity < P₇₅ (middle range) OR
- ⚠️ SNR ≈ 3 (near detection limit) OR
- ⚠️ Retention difference < 3×precision (ambiguous order) OR
- ⚠️ Mixed descendant results OR
- ⚠️ Wide confidence interval on purity
- **Confidence: Moderate (50-80%)**

**LIKELY_FAILURE** (low confidence)

- ❌ Purity < P₂₅ (bottom quartile) OR
- ❌ SNR < 3 (too weak) OR
- ❌ Retention order suspicious (marginal) OR
- ❌ Multiple descendants failed
- **Confidence: Low (20-50%)**

**FAILED** (synthesis confirmed failed)

- ❌ Retention time order violated (Δt < 0) OR
- ❌ No putative product peak detected (SNR < 2) OR
- ❌ All descendants failed
- **Confidence: Very High failure (>95%)**

### 6.11 Uncertainty Quantification

**Sources of Uncertainty:**

1. **Sampling variance** (Poisson counting statistics)
2. **Retention time precision** (fraction width)
3. **Background estimation** (noise level)
4. **Amplification bias** (unknown for most compounds)
5. **Tag-compound dissociation** (unmeasured)

**Confidence Score Calculation:**

confidence_score = w₁×P(purity|success) + w₂×P(retention|success) + w₃×P(descendants|success)

Where weights sum to 1 and are adjusted based on data quality.

**Reporting Standard:**

Always report:

- Validation category (VALIDATED, LIKELY_SUCCESS, etc.)
- Confidence score (0-1)
- Purity with 95% CI
- SNR value
- Number of validated descendants / total descendants
- Any caveats (low signal, retention ambiguous, etc.)

### 6.12 Validation Decision Tree

**Level 1: Retention Time Check** (HARD CONSTRAINT)

Question: retention_time(product) > retention_time(truncations)?

- NO → FAILED (physics violation)
- YES → Proceed to Level 2

**Level 2: Signal Strength** (DETECTION)

Question: SNR > 3?

- NO → LIKELY_FAILURE (below detection)
- YES → Proceed to Level 3

**Level 3: Purity Assessment** (RELATIVE QUALITY)

Question: Purity > P₇₅?

- YES → Proceed to Level 4 (potential VALIDATED)

Question: Purity > P₅₀?

- YES → Proceed to Level 4 (potential LIKELY_SUCCESS)

Question: Purity > P₂₅?

- YES → UNCERTAIN

Otherwise:

- NO → LIKELY_FAILURE

**Level 4: Descendant Consistency** (DAG PROPAGATION)

Question: All descendants validated?

- YES + High purity → VALIDATED
- YES + Moderate purity → LIKELY_SUCCESS

Question: Majority descendants validated?

- YES → LIKELY_SUCCESS

Otherwise:

- NO → Downgrade classification by one level

### 6.13 Comparison to Peak Classification

**Peak Classification** ([Part 5](#part-5-peak-detection-mathematical-foundations)):

- Assigns labels: NULL, TRUNCATION, PUTATIVE_PRODUCT, UNKNOWN
- Based on topological position in DAG
- Does NOT determine synthesis success
- Reports: "This peak is positionally consistent with product"

**Synthesis Validation** (Part 6):

- Evaluates synthesis success probability
- Based on purity, SNR, retention physics, DAG consistency
- Determines: VALIDATED, LIKELY_SUCCESS, UNCERTAIN, LIKELY_FAILURE, FAILED
- Reports: "Synthesis likely succeeded with 85% confidence"

**Relationship:**

1. First run Peak Classification → identify PUTATIVE_PRODUCT
2. Then run Synthesis Validation → assess if synthesis actually worked
3. PUTATIVE_PRODUCT + High purity + Retention correct → VALIDATED
4. PUTATIVE_PRODUCT + Low purity → LIKELY_FAILURE

### 6.14 Mathematical Summary

**Core Equations:**

Purity: Purity(C) = counts(product) / [counts(product) + counts(impurities)]

SNR: SNR(C) = max(counts_product) / background

Bayesian posterior: P(success|data) ∝ P(data|success) × P(success) / Z

Retention constraint: t_product - t_truncation > 2 × Δt_min/2

**Adaptive Thresholds:**

- High purity threshold: P₇₅(dataset)
- Moderate purity threshold: P₅₀(dataset)
- SNR threshold: 3.0 (universal)
- Retention precision: Δt_min / 2 (dataset-specific)

**Computational Complexity:**

- Purity calculation: O(peaks) per compound
- SNR calculation: O(1) per compound
- DAG propagation: O(V + E) topological sort
- Total: O(V × peaks + E) for full validation

---

# III. COMPUTATIONAL IMPLEMENTATION

## PART 7: MATHEMATICAL OPTIMIZATIONS

### 7.1 Algorithms from Graph Theory

**Topological Sort** (O(V + E))

Instead of manual level-by-level iteration (O(V²)), use Kahn's algorithm or DFS-based topological sort to process compounds in guaranteed correct bottom-up order, handling each compound exactly once and automatically managing gaps in levels.

**Benefits:**

- Guaranteed correct bottom-up order
- Processes each compound exactly once
- Handles gaps in levels automatically

**Transitive Reduction** (Remove redundant edges)

Store only DIRECT descendants. Do not store redundant edges where A → C if A → B → C already exists. This saves memory and provides faster traversal.

**Transitive Closure** (Precompute all relationships)

Perform one-time O(V³) computation of transitive closure matrix to enable O(1) ancestor/descendant relationship checks via simple matrix lookup.

### 7.2 Dynamic Programming and Memoization

**Memoize Descendant Sets**

Use caching to compute descendant sets once and reuse for all subsequent queries. For minimal compounds, return empty set immediately. For others, recursively compute direct descendants and their descendants, storing results in cache.

**Benefits:**

- Avoid recomputation (O(V) → O(1) for cached queries)
- Critical for hierarchical peak selection (many ancestor queries)

### 7.3 Graph Indexing Structures

**Adjacency List** (Current Implementation)

Maps each compound to its list of descendants, enabling fast neighbor lookup in O(degree) time.

**Level-Order Index** (Current Implementation)

Maps each level to list of compounds at that level, enabling fast same-level queries in O(1) time.

**Additional Optimizations:**

**Union-Find for Connected Components** (O(α(n)) amortized)

Use Union-Find data structure to identify all disconnected families in the forest structure. Build by iterating through compounds and unioning each with its descendants. Query whether two compounds belong to the same family in near-constant time O(α(n)) ≈ O(1).

**Interval Tree for Position Queries**

Build interval tree indexed by peak boundaries (left_base, right_base) to enable fast "find peaks near position X" queries in O(log n + k) time, where k is the number of results.

### 7.4 Domain Service Operations

**Proposed CompoundHierarchy API** (replaces TruncationHierarchy):

A DAG representing truncation relationships as a partially ordered set (poset).

**Core Poset Operations:**

- Check if compound A is ancestor/descendant of compound B
- Check if two compounds are comparable (one is ancestor/descendant of other)

**Principal Ideal/Filter Operations:**

- Get all descendants of compound C (↓c - principal ideal)
- Get all ancestors of compound C (↑c - principal filter)
- Get complete lineage of compound C (↓c ∪ c ∪ ↑c)

**Graph Properties:**

- Get maximal elements (compounds with no ancestors in dataset)
- Get minimal elements (compounds with no descendants)
- Compute longest chain length

**Efficient Traversal:**

- Topological sort for correct processing order
- Iterate bottom-up through hierarchy
- Iterate top-down through hierarchy

**Analysis:**

- Get level (truncation rank) of a compound
- Get all compounds at a specific level
- Get connected components (disconnected families/lineages)

---

# IV. REFERENCE

## PART 8: DOMAIN VOCABULARY (UBIQUITOUS LANGUAGE)

This section defines precise terminology used throughout the codebase. Consistent vocabulary is critical for maintainability.

### 8.1 Core Entities (Summary)

See [Part 2.1](#21-core-entities) for detailed definitions.

- **Compound**: Synthesis product with sequence, building blocks, chromatogram
- **Peak**: Detected chromatogram feature with position, boundaries, metrics
- **Chromatogram**: Elution profile with time_points, counts, signal variants
- **BuildingBlock**: Chemical block with cycle, code, null flag
- **EquivalenceClass**: Collection of positional variants with same residue sequence (see [Part 4.2](#42-equivalence-classes-and-consensus))

### 8.2 Hierarchical Terminology

**AVOID** ❌:

- "Parent compound" (ambiguous in combinatorial library)
- "Child compound" (relative, not absolute)
- "Top compound" (suggests unique root)

**USE** ✅:

- **Reference compound**: The compound being analyzed
- **Maximal compound**: Longest compound in dataset (no ancestors)
- **Minimal compound**: Shortest compound (no descendants)
- **Ancestor**: Compound with more building blocks
- **Descendant**: Compound with fewer building blocks
- **Lineage**: All related compounds (ancestors + self + descendants)

### 8.3 Peak Classification

**PeakType Enum:**

- NULL: Peak at L₀ (minimal element) position - universal truncation reference
- TRUNCATION: Matches ancestor product position or null position
- PUTATIVE_PRODUCT: Positionally consistent with expected product (NOT chemically validated!)
- UNKNOWN: Peak not matching any expected position (late-eluting, unmatched, ambiguous) - see Part 5.8

**Classification Approach:**

- **Local Detection**: Detect peaks per-compound using Discrete Morse Theory + Poisson Statistics
- **Global Classification**: Classify ALL peaks using DAG constraint propagation
- **Processing Order**: Topological sort (L₀ first, then bottom-up through lineage)

**Mathematical Framework** (see [Part 5](#part-5-peak-detection-mathematical-foundations)):

- **Discrete Morse Theory**: Local maxima in discrete sequences (rigorous peak definition)
- **Poisson Statistics**: Statistical significance testing for count data (Z-score)
- **Prominence**: Chromatographic significance (height above surrounding valleys)
- **Optimal Assignment**: Hungarian algorithm for peak-to-position matching
- **Constraint Propagation**: Descendant constraints flow up through DAG edges

**CRITICAL Distinction - Classification vs Validation:**

- **Peak Classification**: Positional consistency, topological significance (what we DO)
- **Synthesis Validation**: Chemical identity, purity via orthogonal methods (see [Part 6](#part-6-synthesis-validation-theory))
- PUTATIVE_PRODUCT means: "Positionally consistent with expected product"
- PUTATIVE_PRODUCT does NOT mean: "Synthesis was successful" (could be truncation mixtures!)
- See [Part 5.6](#56-classification-scope-and-limitations) for complete scope discussion

### 8.4 Signal Processing

**Signal Variants** (stored in Chromatogram.signals dict):

- "raw": Raw signal (primary signal used for analysis)
- "derivative": First derivative (optional, for advanced analysis)
- "derivative_2": Second derivative (optional, for advanced analysis)
- Custom variants: Add via chromatogram.with_signal(name, array)

**Access:** Retrieve signals via chromatogram.get_signal(name) and add custom signals via chromatogram.with_signal(name, array).

### 8.5 Analysis Modes

**Hierarchy Mode:**

- BLOCK: Building blocks as atomic units (default)
- MONOMER: Individual monomers as atomic units
- NONE: No hierarchical analysis

**Variant Mode:**

- INDIVIDUAL: Analyze each positional variant separately
- CONSENSUS: Average variants by canonical sequence

**Detection Method:**

- MORSE_THEORY: Use discrete Morse theory local maxima (recommended)
- POISSON_PROMINENCE: Use Poisson statistics + prominence filtering (recommended)
- SECOND_DERIVATIVE: Use 2nd derivative zero-crossings (legacy)
- SCIPY: Use scipy.signal.find_peaks (legacy)

**Decision Strategy:**

- HIERARCHICAL_HYPOTHESIS: Hypothesis-based selection
- MAX_SCORE: Maximum score decision
- GAUSSIAN_MEAN: Use Gaussian fit mean

### 8.6 Similarity and Ordering

**Chromatographic Similarity**

- Wasserstein distance (Earth Mover's Distance)
- Jensen-Shannon divergence
- Range: [0, 1] where 1 = identical
- Service: CompoundSimilarityAnalyzer (DOMAIN service)

**Similarity-Based Ordering**

- Hierarchical clustering on distance matrix
- Algorithm: Average linkage with optimal leaf ordering
- **THIS IS DOMAIN LOGIC** (not presentation!)
- Service: CompoundOrderingService (domain/services/)
- Output: Ordered list for visualization

### 8.7 Peak Detection Terminology

**Mathematical Foundations** (see [Part 5](#part-5-peak-detection-mathematical-foundations)):

- **Critical point**: Position where local maximum occurs (discrete: c[i] > c[i-1] and c[i] ≥ c[i+1])
- **Valley**: Local minimum between peaks (defines peak boundaries)
- **Prominence**: Height of peak above surrounding valleys (chromatographic significance)
- **Poisson Z-score**: (c[i] - μ_bg) / √(μ_bg + ε) - statistical significance for count data
- **Background (μ_bg)**: Low-count baseline (10th percentile of all counts)
- **Statistical significance**: Peaks with Z > 3 are real signal (p < 0.001)

**Classification Terminology:**

- **Local detection**: Per-compound peak finding (Morse theory)
- **Global classification**: Lineage-wide constraint propagation
- **L₀ (minimal element)**: Full-null compound, descendant of ALL compounds
- **NULL peak**: Global maximum in L₀ chromatogram (universal reference)
- **Putative product**: Positionally consistent peak (NOT synthesis validation!)

### 8.8 Synthesis Validation Terminology

**Validation Metrics** (see [Part 6](#part-6-synthesis-validation-theory)):

- **Purity**: Fraction of product signal vs total signal
- **SNR (Signal-to-Noise Ratio)**: Product peak height / background
- **Retention order**: Chromatographic physics constraint (product elutes after truncations)
- **Dataset percentile**: Relative ranking within library (P₂₅, P₅₀, P₇₅, P₉₀)
- **MAD (Median Absolute Deviation)**: Robust spread measure

**Validation Categories:**

- **VALIDATED**: High confidence synthesis succeeded (purity > P₇₅, SNR > 5, retention correct)
- **LIKELY_SUCCESS**: Moderate-high confidence (purity > P₅₀, SNR > 3)
- **UNCERTAIN**: Ambiguous result (mixed signals, low counts)
- **LIKELY_FAILURE**: Low confidence (purity < P₂₅ or SNR < 3)
- **FAILED**: High confidence synthesis failed (retention violated, no peak)

**Key Distinction:**

- **Peak Classification** → Positional label (PUTATIVE_PRODUCT)
- **Synthesis Validation** → Success probability (VALIDATED, LIKELY_SUCCESS, etc.)

### 8.9 Anti-Patterns to Avoid

**DO NOT say:**

- ❌ "Successful synthesis" when you mean "found putative product peak"
- ❌ "Product peak" when you mean "putative product" (implies validation)
- ❌ "The peak is the product" (we don't know chemical identity!)
- ❌ "Parent" when context is ambiguous
- ❌ "Sorting" when you mean "hierarchical clustering ordering"
- ❌ "Visualization sorting" (sorting is domain logic)
- ❌ "Infrastructure algorithm" (algorithms are domain services)
- ❌ "Magic numbers" or "hardcoded thresholds" (use adaptive parameters!)
- ❌ "Synthesis succeeded" based only on PUTATIVE_PRODUCT classification

**DO say:**

- ✅ "Putative product peak identified" (honest about limitations)
- ✅ "Positionally consistent with expected product"
- ✅ "Synthesis validated with 90% confidence" (after validation analysis)
- ✅ "Maximal compound" or "reference compound"
- ✅ "Detected via discrete Morse theory local maxima"
- ✅ "Statistically significant peak" (Z > 3, high prominence)
- ✅ "Similarity-based ordering" (domain service)
- ✅ "Compound ordering service" (domain logic)
- ✅ "Peak detection algorithm" (domain service)
- ✅ "Adaptive threshold derived from data" (no magic numbers)
- ✅ "Purity 85% (75th percentile of dataset)" (dataset-relative)

---

## APPENDIX: Quick Reference

### Mathematical Model

- **Structure**: Directed Acyclic Graph (DAG), Partially Ordered Set (Poset)
- **Vertices**: Chemical peptides (monomer mode) or Positional sequences (block mode)
- **Edges**: "is a truncation of" (directed, length-decreasing)
- **Properties**: Reflexive, Antisymmetric, Transitive, Acyclic

### Graph Patterns

- **Building-Block Mode**: Forest (no convergence)
- **Monomer-Level Mode**: DAG with convergence (diamonds)
- **Convergence**: Multiple positional variants → same chemical peptide

### Key Algorithms

- **Topological Sort**: O(V + E) bottom-up traversal
- **Transitive Closure**: O(V³) precompute, O(1) queries
- **Memoization**: Cache descendant sets, O(1) reuse
- **Union-Find**: O(α(n)) connected components
- **Discrete Morse Theory**: O(n) local maxima detection
- **Poisson Significance Testing**: O(n) statistical filtering
- **Prominence Calculation**: O(n) chromatographic significance
- **Hungarian Assignment**: O(p³) optimal peak matching

### Peak Detection Pipeline

1. **Local Detection**: Discrete Morse theory + Poisson statistics + Prominence → significant peaks
2. **Global Classification**: DAG constraint propagation → NULL, TRUNCATION, PUTATIVE_PRODUCT, UNKNOWN
3. **Synthesis Validation**: Purity, SNR, retention order → VALIDATED, LIKELY_SUCCESS, UNCERTAIN, LIKELY_FAILURE, FAILED

### Terminology

- **Maximal compound** (not "parent")
- **Reference compound** (current analysis target)
- **Ancestor/Descendant** (directed relationships)
- **Lineage** (all related compounds)
- **Principal Ideal** (all descendants)
- **Principal Filter** (all ancestors)
- **EquivalenceClass** (positional variants with same residue sequence)
- **L₀** (minimal element, full-null)
- **PUTATIVE_PRODUCT** (positionally consistent, NOT validated)
- **VALIDATED** (synthesis succeeded with high confidence)

---

**END OF THEORETICAL FOUNDATIONS**
