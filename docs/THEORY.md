# LC-Seq Theoretical Foundations

**Last Updated**: 2025-10-12
**Purpose**: Foundational concepts, mathematical structures, and domain vocabulary for the LC-Seq analysis system

---

## Table of Contents

**I. FOUNDATIONS**

1. [Mathematical Structure](#part-1-mathematical-structure)
2. [Domain Foundations](#part-2-domain-foundations)
3. [Hierarchical Relationships](#part-3-hierarchical-relationships)
4. [Graph Properties and Patterns](#part-4-graph-properties-and-patterns)

**II. ANALYSIS METHODS**

5. [Peak Detection Mathematical Foundations](#part-5-peak-detection-mathematical-foundations)
6. [Synthesis Validation Theory](#part-6-synthesis-validation-theory)

**III. COMPUTATIONAL IMPLEMENTATION**

7. [Mathematical Optimizations](#part-7-mathematical-optimizations)

**IV. REFERENCE**

8. [Domain Vocabulary](#part-8-domain-vocabulary-ubiquitous-language)
9. [Appendix: Quick Reference](#appendix-quick-reference)

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

**Positional Block Sequence (Synthesis Path)**

- How/when the peptide was synthesized
- Position-dependent encoding
- Examples:
  - [Val, Null, Null] ≠ [Null, Val, Null] ≠ [Null, Null, Val] (different synthesis paths)

**Implication for Graph Structure:**

- Graph vertices represent **chemical peptides** (in monomer mode) or **positional sequences** (in building-block mode)
- Same chemical peptide can have multiple synthesis paths → convergence in monomer mode

### 1.3 Two Analysis Modes, Two Graph Structures

#### Building-Block Mode (DAG with Convergence at Block Granularity)

**Vertices**: Equivalence classes (block support sequences)

**Example compounds**:

- [Val, Null, Null] = vertex "Val"
- [Null, Val, Null] = vertex "Val" (SAME as above!)
- [Null, Null, Val] = vertex "Val" (SAME as above!)

**Structure**: DAG with convergence (multiple synthesis paths → same block composition)

**Properties**:

- Convergence at block granularity
- Positional variants with same blocks converge to same vertex
- Edges: block support subsequence relationships

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
- Unambiguous: Position information preserved from synthesis
- Graph structure: **DAG with convergence at block granularity**
- Total compounds: 4 × 3 × 3 = 36 (including L₀)

**Edge Generation Algorithm (Block Mode)**:

```
For each compound C with block support sequence S:
  For each compound D in hierarchy:
    If D's block support sequence is a proper subsequence of S:
      If D is at nearest level below C:
        Add edge: C → D

Note: Positional variants with same block support sequence will have
identical descendants because edges are based on subsequence relationships.
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
Graph structure: DAG with convergence at block granularity
  - Positional variants with same block support converge to same equivalence class
```

**Monomer-Level Mode**:

```
Total unique monomers: 3 × 3 × 3 = 27 different monomers
Max sequence length: 3 positions × 3 monomers/block = 9 monomers
Total possible monomer sequences: Exponential in sequence length
Graph structure: DAG with massive convergence
Vertices: Thousands to millions (depending on library)
```

**Why the Complexity Difference**:

- Block mode: Convergence at block granularity (coarser), fewer total sequences
- Monomer mode: Convergence at monomer granularity (finer), many more sequences possible
- Same peptide synthesized via different paths → convergence in both modes
- Monomer mode has exponentially more possible sequences due to finer granularity

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

The LC-Seq system uses three distinct sequence representations at different granularities:

**Positional Block Sequence**

- Full synthesis path including null positions at block granularity
- Example: Leu-Null-Pro
- Unique per compound (encodes synthesis path)
- Attribute: Compound.positional_block_sequence
- Mathematical: σ: Cycles → BuildingBlocks ∪ {Null}

**Block Support Sequence**

- Restriction to support (non-null positions) at block granularity
- Example: Leu-Pro (from Leu-Null-Pro)
- Multiple positional block sequences → same block support sequence (positional variants)
- Used as equivalence class identifier
- Attribute: Compound.block_support_sequence
- Mathematical: π(σ) where π projects to support domain

**Monomer Support Sequence**

- Fully decomposed to individual monomers, support only (no nulls)
- Example: "Leu-Leu-Ala-Val-Pro" where position 1 is composite block "Leu-Ala-Val"
  - Decomposes to: ["Leu", "Leu", "Ala", "Val", "Pro"] (5 individual monomers)
- Chemical peptide identity at monomer granularity
- Composite blocks expand to their constituent monomers
- Attribute: Compound.monomer_support_sequence
- Method: BuildingBlock.decompose_to_monomers()

**Note on Terminology:**

- "Support" refers to the mathematical concept: the domain where the function is non-zero/non-null
- Nulls exist only at block level (synthesis positions); at monomer level there are no nulls by definition

### 2.3 Temporal Representation and Alignment

#### 2.3.1 Absolute Time Representation

Chromatogram signals are defined over **absolute time** (seconds or minutes from injection). All retention time comparisons operate on scalar time values, independent of signal boundaries or sampling rates.

**Properties**:

- Retention times are physically meaningful scalar values (not array indices)
- Peak comparisons use direct scalar arithmetic: |t_observed - t_expected|
- Time units must be consistent across dataset (standard: seconds or minutes)

#### 2.3.2 When Alignment Is Needed

Signals from different runs may have variable start/end times due to manual timing variation or preprocessing.

**No alignment required for core analysis**:

- Peak detection (independent per signal)
- Peak classification (scalar retention time matching)
- Validation metrics (per-compound peak properties)

**Alignment required only for**:

- Similarity metrics (Wasserstein distance, Jensen-Shannon divergence)
- Hierarchical clustering for visualization ordering

#### 2.3.3 Alignment Algorithm (When Needed)

For similarity-based ordering, align signals to common time grid:

1. **Global bounds**: t_min = min(all signals), t_max = max(all signals)
2. **Common grid**: Use finest sampling resolution across signals
3. **Interpolate**: Resample each signal onto common grid via linear interpolation
4. **Compute**: Calculate similarity metrics on aligned signals

**Properties**: Minimal data loss, no padding assumptions, preserves peak positions.

**Alternative**: Skip alignment if using non-similarity orderings (by level, sequence length, validation status).

#### 2.3.4 Mode-Specific Considerations

**Temporal alignment** (this section) aligns time grids for similarity metrics. This is distinct from **sequence alignment** for visualization.

- **Building-block mode**: Position information preserved; use positional mapping for visualization
- **Monomer mode**: Position information lost; use heuristic alignment (right/left-align) for visualization with ambiguity flag

Both modes may need temporal alignment for similarity-based clustering. See Part 1.5.3.1 for sequence alignment details.

### 2.4 The Null Compound (L₀)

#### 2.4.1 Chemical and Graph-Theoretic Definition

**Definition**: L₀ = [Null, Null, ..., Null] is the compound with all building blocks set to null.

**Physical Molecule**: DNA barcode tag + linker chemistry with **no peptide building blocks attached**. L₀ represents the starting scaffold before synthesis.

**Graph Properties**:

- **Minimal element**: Unique bottom element in building-block mode DAG (∄C ∈ V such that C ≺ L₀)
- **Universal descendant**: Every compound C can truncate to L₀ (∀C ∈ V, L₀ ≼ C)
- **Topological role**: Anchor for bottom-up DAG traversal; constraint propagation starts here

#### 2.4.2 Retention Time Principle

**Hydrophobicity Ordering**: DNA tag is hydrophilic; peptides add hydrophobicity. H(C) = H₀ + Σ(H_building_block) → longer peptides elute later.

**Expected**: t(L₀) ≤ t(compounds with peptides). L₀ should elute in early chromatogram region.

**Exceptions** (rare): Highly charged peptides may violate ordering.

#### 2.4.3 L₀ Peak Detection Algorithm

**Primary Method** - Earliest significant peak in early region:

1. **Early window**: Define early chromatogram region (configurable fraction of total time range)
2. **Significant peaks**: Find statistically significant peaks in early window
3. **Select earliest**: L₀_peak = argmin(peak.time) for peaks in early window

**Rationale**: Earliest peak = most hydrophilic = most likely pure DNA tag. Later peaks likely short truncation products.

**Fallback**: If no peak detected, L₀_position = t_min (signal start).

#### 2.4.4 Quality Control

**Warning Thresholds**:

- **Late elution**: If L₀_peak elutes unusually late in chromatogram → WARNING (likely data quality issue - STOP automated analysis)
- **Missing L₀**: If L₀ not in dataset → ERROR (cannot proceed with NULL classification)
- **Weak signal**: If SNR below detection threshold → WARNING (unreliable background estimation)
- **Retention inconsistency**: High variation across dataset → WARNING (chromatography instability)

**Validation Checks**:

1. L₀ exists in dataset
2. Detectable signal (SNR above threshold)
3. Early elution (within expected hydrophilic region)
4. Consistent retention across dataset

#### 2.4.5 Role as Universal Reference

**NULL Peak Classification**: Any compound C may have peak at t(L₀_peak) due to universal truncation constraint (∀C ∈ V, C → L₀ possible).

**Background Estimation**: Use L₀ off-peak regions as noise reference:

- Primary: background = median(L₀_signal where not in peak region)
- Alternative: background = median(bottom 10% of all scaled_count values)

**SNR Reference**: All compounds compute SNR relative to L₀ background.

**Constraint Propagation**: Bottom-up DAG traversal starts at L₀ (topological minimal element); NULL constraint propagates to all ancestors.

#### 2.4.6 Edge Cases

**Case 1 - Missing L₀**: ERROR - require L₀ in dataset or skip NULL classification; use alternative background method.

**Case 2 - Weak/No Signal**: WARNING - fallback position (t_min), alternative background method, reduced validation confidence.

**Case 3 - Multiple Peaks**: Select earliest peak (most conservative); flag if many peaks detected (data quality issue).

**Case 4 - Monomer Mode**: Multiple minimal elements (single AAs) distinct from L₀; NULL peak still based on physical L₀ compound.

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
  - DAG with convergence at block granularity
  - Positional variants with same blocks converge to same equivalence class

- **Monomer Mode**: Individual monomers as atomic units
  - 3 trimeric blocks (9 monomers) → levels 0-9
  - DAG with convergence at monomer granularity
  - Positional variants with same monomers converge to same equivalence class

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

- Pooled mode aggregates positional variants (same chemical peptide)
- Peak selection considers all synthesis paths for a peptide
- Similarity analysis groups by chemical identity, not synthesis

### 4.2 Equivalence Classes and Pooled Mode Analysis

#### 4.2.1 EquivalenceClass Definition

**EquivalenceClass**: Collection of positional variants that represent the same chemical peptide at block granularity.

**Formal Definition**: An equivalence class under the relation R: "has same block support sequence"

**Equivalence Relation Properties**:

- **Reflexive**: Compound C related to itself (C R C)
- **Symmetric**: If A R B, then B R A
- **Transitive**: If A R B and B R C, then A R C

**Example**:

```
Equivalence class for block support sequence "Val":
  - [Val, Null, Null] (positional block sequence - synthesis path 1)
  - [Null, Val, Null] (positional block sequence - synthesis path 2)
  - [Null, Null, Val] (positional block sequence - synthesis path 3)

All three positional variants (members) represent the same chemical peptide: Val
```

**Chemical Interpretation**:

- Same block support sequence = same chemical molecule at block granularity
- Different positional block sequences = different synthesis paths to same molecule
- EquivalenceClass groups synthesis paths by chemical identity
- Members of an equivalence class are positional variants

**Entity**: EquivalenceClass (collection of members)
**Identifier**: block_support_sequence (the equivalence class identifier)

#### 4.2.2 Pooled Mode Overview

**Purpose**: Optional performance optimization providing significant speedup by processing N equivalence classes instead of multiple positional variants per class.

**Benefits**: Computational speedup, noise reduction, molecular-level focus.

**Important**: Individual mode is fully valid; pooled mode is optional.

**Use Pooled Mode When**: Large datasets, similar variant behavior, molecular-level analysis.

**Use Individual Mode When**: Position-specific diagnostics, synthesis troubleshooting, small datasets, dissimilar variants, regulatory requirements.

**Recommendation**: Start with individual mode; adopt pooled mode if speedup needed and variants validated as similar.

#### 4.2.3 Hybrid Pooled Strategy

**Key Insight**: Peak detection is expensive; area integration is cheap. Detect peaks once on pooled signal, integrate areas per variant.

**Three-Phase Workflow**:

1. **Peak Detection** (expensive, once per class): Aggregate signals → detect peaks on pooled signal
2. **Area Integration** (cheap, per variant): Use pooled peak boundaries on individual signals → compute purity per variant
3. **Aggregate Statistics** (reporting): mean_purity, std_purity, class-level validation status

**Result**: Speedup from shared peak detection + individual purities from real signals + position-specific diagnostics available.

#### 4.2.4 Pooled Signal Aggregation

**Algorithm**:

1. **Align**: Interpolate variants to common time grid (see Part 2.3)
2. **Aggregate**: pooled_signal(t) = mean(v₁(t), v₂(t), ..., vₙ(t)) [or median if outliers present]
3. **Validate**: Compute pairwise correlations; check min(corr) exceeds similarity threshold

**Validity Check** (essential): If minimum correlation below threshold → variants differ → fall back to individual mode.

**Threshold**: High correlation required for valid pooling (strong similarity assumption).

#### 4.2.5 Peak Detection on Pooled Signal

Apply full detection pipeline (Discrete Morse + Poisson + Prominence, see Part 5.1-5.2) to pooled signal using same parameters as individual mode.

**Output**: Peak positions, boundaries [t_start, t_end], and classifications (NULL, TRUNCATION, PUTATIVE_PRODUCT, UNKNOWN) for equivalence class.

**These boundaries used for ALL variants in next step.**

**Rationale**: Pooled signal has better SNR (noise averaged); detecting once faster than N detections.

#### 4.2.6 Area Integration on Individual Variants

**Critical**: Use pooled peak boundaries, but integrate on each variant's OWN signal (not pooled signal).

**Algorithm**:

```
For each variant vᵢ:
  1. Use vᵢ's raw signal
  2. Integrate areas using pooled boundaries [t_start, t_end]
  3. Compute Purity(vᵢ) from vᵢ's integrated areas
  4. Validate vᵢ individually (see Part 6)
  5. Store individual purity and validation status
```

**Result**: Real measurements from actual signals + position-specific information preserved + cheap operation.

#### 4.2.7 Aggregate Statistics and Reporting

Compute class-level statistics: mean_purity, std_purity, min/max_purity.

**Validation Summary**: If all variants VALIDATED → Class VALIDATED; if any FAILED → Class FAILED; if mixed → HETEROGENEOUS.

**Two Reporting Levels**:

- **Molecular**: Mean purity ± std, validation status (high-level screening)
- **Position-specific**: Per-variant purity and status (synthesis troubleshooting)

#### 4.2.8 Validity Requirements and Fallback

**Pooled mode valid ONLY IF variants have similar chromatographic behavior.**

**Validity Checks**:

1. **Signal Correlation**: min(pairwise correlations) exceeds threshold; if fails → fall back to individual mode
2. **Peak Position Consistency** (optional): Δt within acceptable retention precision window
3. **Purity Variance** (diagnostic): High std_purity → flag as heterogeneous

**Assumption**: Same molecule → same retention (position-independent chromatography).

**Fails When**: Positional folding effects, incomplete deprotection, tag-peptide interactions.

**Automatic Fallback**: If correlation below threshold, automatically use individual mode and flag class with reason.

#### 4.2.9 Individual vs Pooled Mode Comparison

**Individual Mode**: Full pipeline per variant. **Pros**: Always valid, simpler, full position info. **Cons**: Slower, more results, lower SNR.

**Pooled Mode** (optional): Detect peaks once on pooled signal, integrate per variant. **Pros**: Significant speedup, noise reduction, molecular summaries. **Cons**: Requires high similarity, more complex.

**Recommendation**: Start with individual mode. Consider pooled if large dataset, speedup needed, variants highly similar.

#### 4.2.10 PooledCompound: Delegation Proxy

**Pattern**: Immutable proxy delegating most attributes to real compound, overriding chromatogram with pooled signal.

**Key Methods**: `__getattr__` delegates properties, `__eq__/__hash__` use real compound for hierarchy compatibility.

**Usage**: Create PooledCompound(real_compound, pooled_chromatogram) → process through pipeline → transfer detected_peaks/selected_peak to all variants.

**Properties**: Temporary (not stored), immutable (preserves original data), transparent (duck typing).

**Implementation**: `src/lcseq/domain/entities/pooled_compound.py`

#### 4.2.11 Quotient Hierarchy: Edge Projection

**Problem**: Need hierarchy where nodes are equivalence classes, not individual compounds.

**Solution**: Project edges from original hierarchy onto equivalence classes using direct descendants only.

**Algorithm**:

1. Create mappings: block_support → pooled_compound, compound → block_support
2. For each edge (ancestor, descendant) in original hierarchy, project to (ancestor_class, descendant_class)
3. Deduplicate edges (same edge may arise from multiple variants)
4. Skip self-loops (variants in same class)

**Result**: Quotient poset where nodes are equivalence classes, edges are truncation relationships.

**Why**: Preserves all relationships from all variants; mathematically correct quotient structure.

**Implementation**: `ProcessPooledChromatogramsUseCase._build_quotient_hierarchy()`

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

- Library injected onto LC column and separated into discrete fractions
- Each fraction PCR-amplified and sequenced via NGS
- Raw signal = **DNA barcode sequencing counts per fraction**
- **Working signal = scaled counts** (normalized for sequencing depth, UMI deduplication, amplification bias)
- Can be fractional after scaling/normalization
- Underlying distribution: Poisson (raw counts), approximately Poisson-like after scaling

**Spatial Resolution**

- **Pre-binned data**: Discrete fractions over total elution time
- No sub-fraction resolution (discrete time points)
- Adjacent fractions independent samples (no interpolation needed)
- Molecular diffusion already averaged within fraction collection

**Count Statistics**

- Large molecular input typically processed through sequencing workflow
- Scaled counts after normalization for sequencing depth and amplification bias
- Peak maximum typically multiple times above baseline
- Baseline represents background level
- Noise: σ ≈ √c (Poisson property)

**Signal-to-Noise Ratio (SNR)**

- Variable across compounds (high for abundant, low for rare)
- Abundant compounds: Clear peaks, high SNR
- Rare compounds: Noisy signals, low SNR
- Detection limit: Minimum SNR threshold for statistical significance
- Quantitation limit: Higher SNR threshold for reproducible measurement

---

### 5.1 Discrete Morse Theory Framework

This section establishes the rigorous mathematical framework for peak detection in discrete count data.

**Peak detection for LC-Seq is fundamentally about finding local maxima in discrete sequences.**

#### Mathematical Setup

An LC-Seq chromatogram is a **discrete sequence** c = {c₁, c₂, ..., cₙ} where:

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

**Algorithm** (O(n) single pass, no smoothing, exact):

```
For each index i in [2, n-1]:
    If c[i] > c[i-1] AND c[i] ≥ c[i+1]:
        peak[i] = True
```

Discrete Morse theory is appropriate for discrete count data - finds all local maxima without derivatives, smoothing, or approximation.

### 5.2 Statistical Significance Testing for Peak Detection

**Key Question**: Which peaks are "real signal" vs "statistical noise"?

#### Poisson Count Statistics

LC-Seq data follows Poisson-like statistics (after scaling):

**Distribution**: c[i] ~ Poisson-like with variance proportional to mean

**Noise model**: σ[i] ≈ √(c[i] + ε) where ε prevents division by zero at low counts

**Background estimation**: μ_bg = low percentile(all counts) - captures low-count baseline

#### Statistical Hypothesis Testing

For each detected local maximum at position i:

**Null Hypothesis H₀**: Peak is random Poisson fluctuation of background
**Alternative H₁**: Peak is real signal above background

**Test statistic**:

```
Z = (c[i] - μ_bg) / √(μ_bg + ε)
```

**Decision rule**:

- Z exceeds detection threshold: Reject H₀ (statistical significance)
- Z exceeds quantitation threshold: High confidence for quantitation

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

**Prominence vs Persistence**: Persistent homology assumes multiple close peaks = noise, but in LC-Seq they represent distinct synthesis outcomes (product + truncations). Prominence respects valley separation, computes in O(n), and is standard in analytical chemistry.

#### Adaptive Threshold (No Fixed Values!)

Filter peaks using data-derived prominence threshold:

**Option 1: Percentile-based** - Compute percentile threshold from distribution (retains most prominent peaks)

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

threshold_fraction = configurable baseline fraction (small value to capture peak base)
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

Boundary detection uses valley separation when available, otherwise height thresholds (small fraction of peak height). Peak area = Σ corrected_signal(t) for t ∈ [t_start, t_end].

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

#### 5.3.1 Truncation Boundary: Retention Time Margin

To prevent late-eluting peaks from being incorrectly classified as truncations, we define a **truncation boundary** beyond which peaks cannot be TRUNCATION:

```
truncation_boundary = max(L₀_position, max(descendant_products)) + margin

where margin accommodates retention time variability
```

**Algorithm**: (1) Compute boundary, (2) Validate TRUNCATION assignments (peak.position ≤ boundary), (3) Identify PRODUCT candidates (unassigned peaks with position > boundary), (4) Select first candidate as PUTATIVE_PRODUCT.

**Margin Configuration**: Adjust based on LC stability - tight margin for stable chromatography, loose margin for high variability. Trade-off: too small misclassifies late truncations as products; too large misclassifies early products as truncations.

### 5.4 Global Classification via Constraint Propagation

**Key Insight**: Peak classification must respect the **entire lineage DAG**, not just individual compounds.

#### Algorithm: Bottom-Up Propagation

Global classification processes the DAG in topological order, propagating constraints through edges to ensure consistency across the entire lineage.

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

#### What Classification Can and Cannot Determine

**CAN determine**: Positional consistency (peak location vs DAG constraints, ordering), statistical significance (Z-score, prominence), relational constraints (lineage consistency), hypothesis ranking (most likely product candidate).

**CANNOT determine**: Chemical identity (requires MS/NMR), purity (co-elution possible), synthesis success (peak presence ≠ successful reaction), quantitation (requires calibration).

**PUTATIVE_PRODUCT means**: Positionally consistent peak appearing after truncations, statistically significant, satisfies DAG constraints. This is a hypothesis, NOT chemical confirmation. May be mixture, modified product, or contaminant.

**Synthesis validation requires:**

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

### 5.7 Classification Limitations and Unknown Peaks

#### 5.7.1 What Can Be Classified

**NULL** (high confidence): Peak at t(L₀_peak) - DNA tag only, complete truncation.

**TRUNCATION** (high confidence): Peaks matching ancestor product or L₀ positions - DAG constraint propagation.

**PUTATIVE_PRODUCT** (positional hypothesis): First significant peak after truncations - NOT chemical confirmation, see Part 5.6.

**Basis**: Retention time, DAG constraints, statistical significance, ordering constraints.

#### 5.7.2 What Cannot Be Classified (UNKNOWN)

**UNKNOWN peaks** include:

- **Late-eluting**: Retention time > t(putative_product)
- **Unmatched**: Outside tolerance windows, no ancestor constraint satisfaction
- **Ambiguous**: Multiple hypotheses, insufficient confidence

**Fundamental Limitation**: Without orthogonal data (MS, NMR), chemical identity cannot be determined from retention time alone.

#### 5.7.3 Possible Identities of Unknown Peaks

**Late-eluting peaks could be**:

- **Oligomers** (n-mers): Higher-order aggregates at integer multiples of monomer retention time
- **Contaminants**: Synthesis reagents, degradation products, column bleed, carry-over
- **Modified products**: Incomplete reactions, side reactions, post-synthesis modifications
- **Artifacts**: Signal distortions, instrumentation issues

**Cannot distinguish without MS, NMR, or authentic standards.**

#### 5.7.4 Handling Unknown Peaks

**Strategy**: Label all unidentifiable peaks as UNKNOWN (honest acknowledgment of limitations).

**Purity**: All non-product peaks reduce purity:
Purity(C) = counts(putative_product) / counts(all_peaks)

**Key Point**: Unknown peaks count as impurities regardless of identity - NOT intended product.

**Reporting**: Flag compounds with high unknown fraction; recommend follow-up analysis if unknowns dominate.

#### 5.7.5 Oligomer Hypothesis (Not Classification)

**Evidence for oligomers**: Position at integer multiples of monomer retention time, ladder pattern, intensity decay, family consistency.

**Hypothesis Strength**: Single peak without pattern (weak evidence); ladder pattern (moderate); ladder with decay and family consistency (strong).

**Use**: Flag for MS follow-up, identify aggregation issues, guide optimization. **NOT automated classification** - remains UNKNOWN until MS confirmation.

#### 5.7.6 Quality Metrics

**Per-compound**: unknown_fraction = counts(unknown) / counts(all_peaks). Categories: low (clean), moderate (acceptable), high (concern).

**Dataset**: Flag if median_unknown or high_unknown_rate exceeds threshold (systematic issue).

#### 5.7.7 When Definitive Identification Possible

**MS**: Confirms molecular weight, distinguishes oligomers from contaminants. **NMR**: Confirms structure. **Authentic standards**: Match retention times. **After orthogonal data**: Reclassify UNKNOWN with confident labels.

#### 5.7.8 Decision Summary

**Decision Tree**: L₀ position → NULL; ancestor position → TRUNCATION; first peak after truncations → PUTATIVE_PRODUCT; otherwise → UNKNOWN.

**Conservative Principle**: When in doubt, label UNKNOWN.

**Purity Impact**: All non-product peaks (UNKNOWN, TRUNCATION, NULL) count as impurities.

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
- Purity ≈ 0.5 → Product and impurities comparable
- Purity → 0 → Dominated by truncations/unknowns

**Statistical Uncertainty:**

Standard error of purity estimate:
SE(purity) = √[purity × (1-purity) / total_scaled_counts]

95% Confidence Interval:
CI = purity ± 1.96 × SE(purity)

**Minimum count threshold:**
Sufficient total_scaled_counts required for narrow confidence interval

### 6.4 Distribution-Based Thresholds

For dataset D with compounds C₁, C₂, ..., Cₙ:

**Step 1: Characterize Dataset**

1. Compute purity(Cᵢ) for all compounds
2. Extract percentiles across distribution
3. Compute Median Absolute Deviation: MAD = median(|purity - median|)
4. Estimate background from L₀ or low-count tail

**Step 2: Define Adaptive Categories**

- **Exceptional purity**: purity in top percentile tier
- **High purity**: purity in upper quartile range
- **Moderate purity**: purity near median
- **Low purity**: purity in lower quartile range
- **Very low purity**: purity in bottom percentile tier

**Step 3: Adjust for Dataset Quality**

If MAD(purity) is small: # High-quality library
→ Use strict thresholds (upper quartile for validation)

If MAD(purity) is large: # Variable library
→ Use lenient thresholds (median for validation)

This ensures fair evaluation regardless of library-wide synthesis quality.

### 6.5 Signal-to-Noise Ratio (Universal Metric)

**Background Estimation:**

**Option 1:** From L₀ (full-null compound)
background = median(scaled_counts) across L₀ signal

**Option 2:** From low-count tail
background = median(low percentile tail of all scaled_count values)

**Signal-to-Noise Ratio:**

SNR(C) = max(scaled_counts_product_C) / background

**Interpretation:**

- SNR above quantitation threshold: High confidence detection (clear signal)
- SNR between detection and quantitation thresholds: Moderate confidence (detectable)
- SNR below detection threshold: Near noise floor (unreliable)

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
Δt_min = min(retention_time_{i+1} - retention_time_i) across all peaks
retention_precision = Δt_min / 2

**Retention Order Validation:**

For confident ordering:
t_product - t_truncation > multiple × retention_precision

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

- P(purity=p | succeeded) ~ Beta distribution with high mode (successful synthesis yields high purity)
- P(purity=p | failed) ~ Beta distribution with low mode (failed synthesis yields low purity)

**Retention order likelihood:**

- P(order_correct | succeeded) = high probability
- P(order_correct | failed) = low probability

**Descendant evidence:**

- P(descendants_validated | succeeded) = high_prob^n for n descendants
- P(descendants_validated | failed) = low_prob^n

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

Exceptionally clean: purity > median + multiple×MAD
Exceptionally poor: purity < median - multiple×MAD

### 6.10 Validation Classification

**Decision Framework:**

**VALIDATED** (synthesis confirmed)

- ✅ Retention time order correct (Δt exceeds minimum multiple of precision)
- ✅ Purity in upper percentile tier
- ✅ SNR above quantitation threshold
- ✅ All descendants validated
- ✅ Confidence interval excludes low purity
- **Confidence: Very High**

**LIKELY_SUCCESS** (high confidence)

- ✅ Retention time order correct
- ✅ Purity above median
- ✅ SNR above detection threshold
- ✅ Majority descendants validated
- **Confidence: High**

**UNCERTAIN** (ambiguous)

- ⚠️ Purity in middle range OR
- ⚠️ SNR near detection limit OR
- ⚠️ Retention difference ambiguous OR
- ⚠️ Mixed descendant results OR
- ⚠️ Wide confidence interval on purity
- **Confidence: Moderate**

**LIKELY_FAILURE** (low confidence)

- ❌ Purity in lower percentile tier OR
- ❌ SNR below detection threshold OR
- ❌ Retention order suspicious (marginal) OR
- ❌ Multiple descendants failed
- **Confidence: Low**

**FAILED** (synthesis confirmed failed)

- ❌ Retention time order violated (Δt < 0) OR
- ❌ No putative product peak detected OR
- ❌ All descendants failed
- **Confidence: Very High failure**

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

Question: SNR above detection threshold?

- NO → LIKELY_FAILURE (below detection)
- YES → Proceed to Level 3

**Level 3: Purity Assessment** (RELATIVE QUALITY)

Question: Purity in upper tier?

- YES → Proceed to Level 4 (potential VALIDATED)

Question: Purity above median?

- YES → Proceed to Level 4 (potential LIKELY_SUCCESS)

Question: Purity in middle range?

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

- High purity threshold: Upper quartile (dataset-derived)
- Moderate purity threshold: Median (dataset-derived)
- SNR threshold: Detection limit (dataset-derived)
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

Use Union-Find data structure to identify connected components in the DAG. Build by iterating through compounds and unioning each with its descendants. Query whether two compounds belong to the same connected component in near-constant time O(α(n)) ≈ O(1). Note: Within a lineage, all compounds form a single connected component.

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
- **EquivalenceClass**: Collection of positional variants with same block support sequence (see [Part 4.2](#42-equivalence-classes-and-pooled-mode-analysis))
- **PooledCompound**: Immutable proxy for processing pooled chromatograms (see [Part 4.2.10](#4210-pooledcompound-pooled-signal-processing-proxy))

#### 8.1.1 Sequence Representations

LC-Seq uses a systematic naming convention for sequences based on two dimensions: **granularity** (blocks vs monomers) and **support** (with or without nulls).

**The Three Sequence Types:**

1. **positional_block_sequence** (Compound property)

   - Full synthesis path at block granularity INCLUDING null positions
   - Example: `"Val-Null-Leu"` (3 cycles: Val at position 1, skipped position 2, Leu at position 3)
   - Use case: Identify exact synthesis path, position-specific diagnostics
   - Access: `compound.positional_block_sequence`

2. **block_support_sequence** (Compound property, EquivalenceClass identifier)

   - Non-null building blocks only (SUPPORT = non-zero domain)
   - Example: `"Val-Leu"` (same chemistry as above, ignoring positional encoding)
   - Use case: Chemical identity at block granularity, equivalence class identifier
   - Access: `compound.block_support_sequence`
   - **Key property**: Equivalence classes group by this sequence

3. **monomer_support_sequence** (Compound property)
   - Fully decomposed to individual monomers, no nulls
   - Example: `"Val"` (single block) → `"Val"` (single monomer); `"ValLeu"` (dipeptide block) → `"Val-Leu"` (two monomers)
   - Use case: Chemical identity at monomer granularity, finest resolution
   - Access: `compound.monomer_support_sequence`

**Mathematical Terminology: "Support"**

- **Support** of a function = domain where function is non-zero/non-null
- In LC-Seq: building blocks that are NOT null (skipped cycles)
- **block_support_sequence** = projection onto support (removes nulls)
- **positional_block_sequence** = full domain (includes nulls)

**Why This Naming?**

- **Consistent**: Two-dimensional taxonomy (granularity × support)
- **Unambiguous**: No conflict with "canonical amino acids" terminology
- **Mathematical**: "Support" is standard mathematical concept
- **Precise**: Clearly specifies what each sequence represents

### 8.2 Hierarchical Terminology

**Standard Terms:**

- **Reference compound**: The compound being analyzed
- **Maximal compound**: Longest compound in dataset (no ancestors)
- **Minimal compound**: Shortest compound (no descendants)
- **Ancestor**: Compound with more building blocks
- **Descendant**: Compound with fewer building blocks
- **Lineage**: All related compounds (ancestors + self + descendants)

**Note**: Terms like "parent", "child", or "top" are ambiguous in combinatorial libraries and should not be used.

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

- BLOCK: Building blocks as atomic units
- MONOMER: Individual monomers as atomic units
- NONE: No hierarchical analysis

**Variant Mode:**

- INDIVIDUAL: Analyze each positional variant separately
- POOLED: Hybrid approach - peak detection on aggregated signal, purity on individual variants (optional optimization)

**Detection Method:**

- MORSE_THEORY: Use discrete Morse theory local maxima (recommended)
- POISSON_PROMINENCE: Use Poisson statistics + prominence filtering (recommended)
- SECOND_DERIVATIVE: Use 2nd derivative zero-crossings (legacy)
- SCIPY: Use scipy.signal.find_peaks (legacy)

**Decision Strategy:**

- HIERARCHICAL_HYPOTHESIS: Hypothesis-based selection
- MAX_SCORE: Maximum score decision
- GAUSSIAN_MEAN: Use Gaussian fit mean

#### 8.5.1 Pooled Mode Terminology

**Key Terms:**

- **Equivalence class**: Collection of positional variants with same block support sequence (chemical identity at block level)
- **Members**: Positional variants within an equivalence class (formerly "compounds" in EquivalenceClass)
- **Pooled signal**: Aggregated signal (mean or median) across all members of equivalence class
- **Pooled compound**: Processing entity - either real Compound (single variant) or PooledCompound (aggregate)
- **Quotient hierarchy**: Quotient structure where nodes are equivalence classes (via edge projection from original hierarchy)
- **Correlation threshold**: Minimum pairwise correlation for pooling validity (high similarity required)

**IMPORTANT Distinctions:**

❌ **NOT a "representative"**: Pooled compound is NOT selected from members
✅ **Synthetic aggregate**: Pooled signal = mean/median of ALL members
❌ **NOT stored**: PooledCompound is temporary processing artifact, discarded after pipeline execution
✅ **Results transferred**: Detected peaks copied from pooled compound to all members

**Pooled Mode Workflow:**

1. Group compounds → equivalence classes (by block support sequence)
2. Aggregate signals → pooled_chromatogram (mean/median)
3. Validate correlation (automatic fallback if < threshold)
4. Create PooledCompound (proxy with pooled chromatogram)
5. Process through pipeline (peak detection, classification)
6. Transfer results to all members (detected_peaks, selected_peak)
7. Integrate areas on individual members (NOT pooled signal!)
8. Compute individual purities and validation categories

**Why "Pooled" not "Representative":**

- **Representative** implies selection (choose one member to represent class)
- **Pooled** implies aggregation (synthesize new signal from all members)
- LC-Seq uses aggregation (mean/median), not selection
- PooledCompound is synthetic entity, not a real member

**References:**

- THEORY.md Section 4.2: Equivalence Classes and Pooled Mode Analysis
- THEORY.md Section 4.2.10: PooledCompound pattern
- THEORY.md Section 4.2.11: Quotient hierarchy construction
- `src/lcseq/application/use_cases/process_pooled_chromatograms.py`

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
- **Background (μ_bg)**: Low-count baseline (low percentile of all counts)
- **Statistical significance**: Peaks exceeding Z-score threshold are real signal

**Classification Terminology:**

- **Local detection**: Per-compound peak finding (Morse theory)
- **Global classification**: Lineage-wide constraint propagation
- **L₀ (minimal element)**: Full-null compound, descendant of ALL compounds
- **NULL peak**: Global maximum in L₀ chromatogram (universal reference)
- **Putative product**: Positionally consistent peak (NOT synthesis validation!)
- **Truncation boundary**: max(truncation_positions) + retention_margin - hard cutoff for TRUNCATION classification
- **Truncation margin**: Time margin accounting for retention time variability

**Truncation Boundary:**

- **Purpose**: Prevent Hungarian algorithm from incorrectly assigning late-eluting peaks as truncations
- **Definition**: `truncation_boundary = max(L₀_position, max(descendant_products)) + retention_margin`
- **Usage**: Peaks must elute BEFORE boundary to be classified as TRUNCATION
- **Product constraint**: Peaks must elute AFTER boundary to be considered for PUTATIVE_PRODUCT
- **References**: THEORY.md Section 5.3.1

### 8.8 Synthesis Validation Terminology

**Validation Metrics** (see [Part 6](#part-6-synthesis-validation-theory)):

- **Purity**: Fraction of product signal vs total signal
- **SNR (Signal-to-Noise Ratio)**: Product peak height / background
- **Retention order**: Chromatographic physics constraint (product elutes after truncations)
- **Dataset percentile**: Relative ranking within library (quartiles, median, upper/lower tiers)
- **MAD (Median Absolute Deviation)**: Robust spread measure

**Validation Categories:**

- **VALIDATED**: High confidence synthesis succeeded (purity in upper tier, high SNR, retention correct)
- **LIKELY_SUCCESS**: Moderate-high confidence (purity above median, sufficient SNR)
- **UNCERTAIN**: Ambiguous result (mixed signals, low counts)
- **LIKELY_FAILURE**: Low confidence (purity in lower tier or insufficient SNR)
- **FAILED**: High confidence synthesis failed (retention violated, no peak)

**Key Distinction:**

- **Peak Classification** → Positional label (PUTATIVE_PRODUCT)
- **Synthesis Validation** → Success probability (VALIDATED, LIKELY_SUCCESS, etc.)

---

## APPENDIX: Quick Reference

### Mathematical Model

- **Structure**: Directed Acyclic Graph (DAG), Partially Ordered Set (Poset)
- **Vertices**: Chemical peptides (monomer mode) or Positional sequences (block mode)
- **Edges**: "is a truncation of" (directed, length-decreasing)
- **Properties**: Reflexive, Antisymmetric, Transitive, Acyclic

### Graph Patterns

- **Building-Block Mode**: DAG with convergence at block granularity
  - Positional variants with same block support converge
- **Monomer-Level Mode**: DAG with convergence at monomer granularity
  - Positional variants with same monomer sequence converge
- **Convergence**: Multiple synthesis paths → same atomic unit composition

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
- **EquivalenceClass** (positional variants with same block support sequence)
- **L₀** (minimal element, full-null)
- **PUTATIVE_PRODUCT** (positionally consistent, NOT validated)
- **VALIDATED** (synthesis succeeded with high confidence)

---

**END OF THEORETICAL FOUNDATIONS**
