#!/usr/bin/env python3
"""LC-Seq Coupling Heatmap Generator - Log-Linear Regression Analysis.

Estimates the contribution of each block-block interface to compound purity
using a principled log-linear regression model.

The model:
    log(purity) = β₀ + Σ β_interface × I(interface present) + ε

Each β coefficient estimates the log-contribution of that interface to purity,
controlling for all other interfaces present. This handles the multiplicative
nature of coupling (purity ≈ p₁ × p₂ × ... × pₖ) and avoids the bias of
simple averaging (which conflates "bad interface" with "many steps").

The heatmap shows:
- Rows: prev_block (on chain) - sorted by LAST monomer (N-terminus, receives coupling)
- Columns: next_block (being coupled) - sorted by FIRST monomer (C-terminus, does coupling)
- Cell color: β coefficient (negative = reduces purity, positive = increases purity)

Uses cross-validated ridge regression (RidgeCV) to automatically handle
potential collinearity in the design matrix - no manual tuning required.

Usage:
    python examples/coupling_heatmap.py \\
        --input results/library/library_analysis.jsonl \\
        --output results/library/coupling_heatmap.png \\
        --sort-by-interface
"""

from pathlib import Path
import sys
import argparse
import json
from dataclasses import dataclass
from typing import Dict, List, Tuple
from collections import defaultdict

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from sklearn.linear_model import RidgeCV, ElasticNetCV
from sklearn.kernel_ridge import KernelRidge
from sklearn.model_selection import GridSearchCV

# Import PurityCalculator from lcseq.domain.services
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))
from lcseq.domain.services.purity_calculator import PurityCalculator


# =============================================================================
# Monomer Sorting Order
# =============================================================================

# Amino acid base order (within each group)
AMINO_ACID_ORDER = ["Ala", "Nvl", "Val", "Leu", "Phe", "Pro"]

# Group priority (lower = earlier in sort)
GROUP_ORDER = {
    "L_plain": 0,      # L non-NMe (e.g., Leu, Pro)
    "D_plain": 1,      # D non-NMe (e.g., DLeu, DPro)
    "L_NMe": 2,        # L NMe (e.g., LeuMe)
    "D_NMe": 3,        # D NMe (e.g., DLeuMe)
    "L_bHomo": 4,      # L βHomo (e.g., βHomoleu)
    "D_bHomo": 5,      # D βHomo (e.g., DβHomoleu)
    "peptoid": 6,      # Peptoids (LA03, LA18)
    "unknown": 7,      # Fallback for unknown monomers
}


def classify_monomer(monomer: str) -> Tuple[int, int]:
    """
    Classify monomer and return (group_priority, aa_priority).

    Returns tuple for sorting: lower values sort first.

    Group order:
    1. L non-NMe amino acids (Leu, Pro, etc.)
    2. D non-NMe amino acids (DLeu, DPro, etc.)
    3. L NMe amino acids (LeuMe, etc.)
    4. D NMe amino acids (DLeuMe, etc.)
    5. L βHomo amino acids (βHomoleu)
    6. D βHomo amino acids (DβHomoleu)
    7. Peptoids (LA03, LA18)

    Within each group, sort by amino acid: Ala → Nvl → Val → Leu → Phe → Pro
    """
    # Peptoids (LA03, LA18)
    if monomer.startswith("LA"):
        # LA03 before LA18
        return (GROUP_ORDER["peptoid"], int(monomer[2:]) if monomer[2:].isdigit() else 99)

    # βHomo variants (note: data uses lowercase, e.g., βHomoleu)
    if monomer.startswith("DβHomo") or monomer.startswith("Dβhomo"):
        base = monomer[6:].capitalize()  # strip DβHomo, capitalize
        aa_idx = AMINO_ACID_ORDER.index(base) if base in AMINO_ACID_ORDER else 99
        return (GROUP_ORDER["D_bHomo"], aa_idx)
    if monomer.startswith("βHomo") or monomer.startswith("βhomo"):
        base = monomer[5:].capitalize()  # strip βHomo, capitalize
        aa_idx = AMINO_ACID_ORDER.index(base) if base in AMINO_ACID_ORDER else 99
        return (GROUP_ORDER["L_bHomo"], aa_idx)

    # Check D prefix (D followed by uppercase letter)
    is_D = monomer.startswith("D") and len(monomer) > 1 and monomer[1].isupper()
    base = monomer[1:] if is_D else monomer

    # Check NMe suffix
    is_NMe = base.endswith("Me")
    aa_base = base[:-2] if is_NMe else base

    # Get amino acid priority
    aa_idx = AMINO_ACID_ORDER.index(aa_base) if aa_base in AMINO_ACID_ORDER else 99

    # Determine group
    if is_D and is_NMe:
        group = GROUP_ORDER["D_NMe"]
    elif is_D:
        group = GROUP_ORDER["D_plain"]
    elif is_NMe:
        group = GROUP_ORDER["L_NMe"]
    else:
        group = GROUP_ORDER["L_plain"]

    return (group, aa_idx)


# =============================================================================
# Data Classes
# =============================================================================

@dataclass
class InterfaceEffect:
    """Regression coefficient for a single block-block interface."""
    coefficient: float      # β value (log-scale effect)
    count: int              # Number of compounds with this interface

    @property
    def multiplicative_effect(self) -> float:
        """exp(β) - multiplicative effect on purity."""
        return np.exp(self.coefficient)


# =============================================================================
# Purity Calculation
# =============================================================================
# Note: Purity calculation is now provided by PurityCalculator.calculate_from_peaks()
# from lcseq.domain.services.purity_calculator


# =============================================================================
# Boundary Extraction
# =============================================================================

def extract_boundaries(record: Dict) -> List[Tuple[str, str]]:
    """
    Extract all block boundaries from a compound.

    Returns list of (prev_block, next_block) tuples where:
    - prev_block is already on chain
    - next_block is being coupled

    Building blocks are in cycle order (0 = C-term, synthesized first).
    """
    if "building_blocks" not in record:
        return []

    blocks_data = record["building_blocks"]

    # Get non-null blocks in cycle order
    non_null_blocks = []
    for bb in sorted(blocks_data, key=lambda x: x["cycle"]):
        if not bb.get("is_null", False):
            non_null_blocks.append(bb["code"])

    # Extract boundaries between consecutive blocks
    boundaries = []
    for i in range(len(non_null_blocks) - 1):
        prev_block = non_null_blocks[i]      # Already on chain
        next_block = non_null_blocks[i + 1]  # Being coupled
        boundaries.append((prev_block, next_block))

    return boundaries


def extract_boundaries_with_step(record: Dict) -> List[Tuple[int, str, str]]:
    """
    Extract all block boundaries from a compound with step information.

    Returns list of (step, prev_block, next_block) tuples where:
    - step is the synthesis step (0-indexed: step 0 = first block to support)
    - prev_block is already on chain (or "SUPPORT" for step 0)
    - next_block is being coupled

    Building blocks are in cycle order (0 = C-term, synthesized first).

    Steps:
    - Step 0: SUPPORT -> first_block (initial coupling to resin)
    - Step 1: first_block -> second_block
    - Step 2: second_block -> third_block
    - etc.
    """
    if "building_blocks" not in record:
        return []

    blocks_data = record["building_blocks"]

    # Get non-null blocks in cycle order
    non_null_blocks = []
    for bb in sorted(blocks_data, key=lambda x: x["cycle"]):
        if not bb.get("is_null", False):
            non_null_blocks.append(bb["code"])

    if not non_null_blocks:
        return []

    # Extract boundaries with step number
    boundaries = []

    # Step 0: N (resin/support) -> first block (initial coupling)
    boundaries.append((0, "N", non_null_blocks[0]))

    # Steps 1+: Block-to-block couplings
    for i in range(len(non_null_blocks) - 1):
        step = i + 1
        prev_block = non_null_blocks[i]      # Already on chain
        next_block = non_null_blocks[i + 1]  # Being coupled
        boundaries.append((step, prev_block, next_block))

    return boundaries


# =============================================================================
# Compositional Encoding
# =============================================================================

def build_vocabulary(records: List[Dict]) -> Tuple[List[str], int]:
    """
    Collect all unique monomers and determine max block length.

    Returns
    -------
    vocabulary : List[str]
        Sorted list of all unique monomers
    max_block_length : int
        Maximum number of monomers in any block
    """
    monomers = set()
    max_length = 0

    for record in records:
        for bb in record.get("building_blocks", []):
            if not bb.get("is_null", False):
                parts = bb["code"].split("-")
                monomers.update(parts)
                max_length = max(max_length, len(parts))

    return sorted(monomers), max_length


def encode_block(
    block_code: str,
    vocab_to_idx: Dict[str, int],
    vocab_size: int,
    max_length: int
) -> np.ndarray:
    """
    Positional one-hot encoding of a block.

    Parameters
    ----------
    block_code : str
        Block code like "Leu-DLeu-Pro"
    vocab_to_idx : Dict[str, int]
        Mapping from monomer name to index
    vocab_size : int
        Size of vocabulary
    max_length : int
        Maximum block length (for padding)

    Returns
    -------
    np.ndarray
        Vector of length max_length × vocab_size
    """
    encoding = np.zeros(max_length * vocab_size)

    monomers = block_code.split("-")
    for pos, monomer in enumerate(monomers):
        if pos >= max_length:
            break
        if monomer in vocab_to_idx:
            idx = pos * vocab_size + vocab_to_idx[monomer]
            encoding[idx] = 1.0

    return encoding


def encode_block_integer(
    block_code: str,
    vocab_to_idx: Dict[str, int],
    max_length: int
) -> List[int]:
    """
    Integer positional encoding of a block.

    Parameters
    ----------
    block_code : str
        Block code like "Leu-DLeu-Pro"
    vocab_to_idx : Dict[str, int]
        Mapping from monomer name to index (1-based, 0 = padding)
    max_length : int
        Maximum block length (for padding)

    Returns
    -------
    List[int]
        List of length max_length with monomer indices (0 = padding)
    """
    encoding = [0] * max_length
    monomers = block_code.split("-")
    for pos, monomer in enumerate(monomers):
        if pos >= max_length:
            break
        encoding[pos] = vocab_to_idx.get(monomer, 0)
    return encoding


def build_integer_design_matrix(
    records: List[Dict],
) -> Tuple[np.ndarray, np.ndarray, List[str], Dict]:
    """
    Build design matrix with integer positional encoding for tree models.

    Creates one row per interface (not per compound):
    [prev_pos0, prev_pos1, prev_pos2, next_pos0, next_pos1, next_pos2]

    Each interface represents "coupling block [i,j,k] to block [l,m,n]".

    Returns
    -------
    y : np.ndarray
        log(purity) for each interface (duplicated for multi-interface compounds)
    X : np.ndarray
        Integer-encoded features [n_interfaces × n_features]
    feature_names : List[str]
        Names for each feature
    metadata : Dict
        vocabulary, vocab_to_idx, max_block_length, n_compounds
    """
    # Build vocabulary and determine max block length
    vocabulary, max_block_length = build_vocabulary(records)
    vocab_to_idx = {m: i + 1 for i, m in enumerate(vocabulary)}  # 0 = padding

    n_features = 2 * max_block_length  # prev + next positions

    # Build feature names
    feature_names = []
    for role in ["prev", "next"]:
        for pos in range(max_block_length):
            feature_names.append(f"{role}_pos{pos}")

    # Build design matrix - one row per interface
    y_list = []
    X_list = []
    n_compounds = 0

    for record in records:
        purity = PurityCalculator.calculate_from_peaks(record.get("peaks", []))
        if purity <= 0:
            continue

        log_purity = np.log(purity)
        boundaries = extract_boundaries(record)

        if not boundaries:
            continue

        n_compounds += 1

        # One row per interface
        for prev_block, next_block in boundaries:
            prev_enc = encode_block_integer(prev_block, vocab_to_idx, max_block_length)
            next_enc = encode_block_integer(next_block, vocab_to_idx, max_block_length)
            row = prev_enc + next_enc
            X_list.append(row)
            y_list.append(log_purity)

    metadata = {
        "vocabulary": vocabulary,
        "vocab_to_idx": vocab_to_idx,
        "max_block_length": max_block_length,
        "n_features": n_features,
        "n_compounds": n_compounds,
    }

    return np.array(y_list), np.array(X_list), feature_names, metadata


def build_compositional_design_matrix(
    records: List[Dict],
) -> Tuple[np.ndarray, np.ndarray, List[str], Dict]:
    """
    Build design matrix using compositional encoding of blocks.

    Instead of block identity features, encodes each block as positional
    one-hot vectors of its constituent monomers.

    Returns
    -------
    y : np.ndarray
        log(purity) for each compound
    X : np.ndarray
        Compositional features [n_compounds × n_features]
    feature_names : List[str]
        Interpretable names for each feature
    metadata : Dict
        vocabulary, max_block_length, etc.
    """
    # Build vocabulary and determine max block length
    vocabulary, max_block_length = build_vocabulary(records)
    vocab_size = len(vocabulary)
    vocab_to_idx = {m: i for i, m in enumerate(vocabulary)}

    n_features = 2 * max_block_length * vocab_size  # prev + next

    # Build feature names
    feature_names = []
    for role in ["prev", "next"]:
        for pos in range(max_block_length):
            for monomer in vocabulary:
                feature_names.append(f"{role}_pos{pos}_{monomer}")

    # Build design matrix
    y_list = []
    X_list = []

    for record in records:
        purity = PurityCalculator.calculate_from_peaks(record.get("peaks", []))
        if purity <= 0:
            continue

        y_list.append(np.log(purity))

        row = np.zeros(n_features)
        half = max_block_length * vocab_size

        for prev_block, next_block in extract_boundaries(record):
            prev_enc = encode_block(prev_block, vocab_to_idx, vocab_size, max_block_length)
            next_enc = encode_block(next_block, vocab_to_idx, vocab_size, max_block_length)
            row[:half] += prev_enc
            row[half:] += next_enc

        X_list.append(row)

    metadata = {
        "vocabulary": vocabulary,
        "vocab_size": vocab_size,
        "max_block_length": max_block_length,
        "n_features": n_features,
    }

    return np.array(y_list), np.array(X_list), feature_names, metadata


# =============================================================================
# Design Matrix Construction (Block Identity)
# =============================================================================

def build_design_matrix_per_step(
    records: List[Dict],
    step: int,
) -> Tuple[np.ndarray, np.ndarray, List[Tuple[str, str]], Dict[Tuple[str, str], int]]:
    """
    Build design matrix for a specific synthesis step.

    Only includes compounds where the specified step is the TERMINAL step
    (i.e., step k analysis uses only level-(k+1) compounds). This ensures
    that purity measurements reflect only steps 0..k without contamination
    from later synthesis steps.

    Parameters
    ----------
    records : List[Dict]
        JSONL records with peaks and building_blocks
    step : int
        Synthesis step to filter (0-indexed: step 0 = N→first block)

    Returns
    -------
    y : np.ndarray
        log(purity) for each compound (excluding purity=0)
    X : np.ndarray
        Binary matrix [n_compounds × n_interfaces]
    interface_list : List[Tuple[str, str]]
        List of (prev_block, next_block) interface tuples for this step
    interface_counts : Dict
        Count of compounds containing each interface
    """
    # First pass: collect all unique interfaces for this step and count them
    # Only from compounds where this step is the TERMINAL step
    interface_counts: Dict[Tuple[str, str], int] = defaultdict(int)
    for record in records:
        boundaries = extract_boundaries_with_step(record)
        if not boundaries:
            continue
        max_step = max(s for s, _, _ in boundaries)
        if max_step != step:
            continue  # Only terminal compounds
        for s, prev_block, next_block in boundaries:
            if s == step:
                interface_counts[(prev_block, next_block)] += 1

    interface_list = sorted(interface_counts.keys())
    interface_to_idx = {iface: i for i, iface in enumerate(interface_list)}

    # Second pass: build matrix
    y_list = []
    X_list = []

    for record in records:
        purity = PurityCalculator.calculate_from_peaks(record.get("peaks", []))
        if purity <= 0:
            continue  # Can't take log of zero

        # Check if this compound has the specified step as its FINAL step
        # This ensures purity reflects only steps 0..k, not later steps
        boundaries = extract_boundaries_with_step(record)
        if not boundaries:
            continue

        max_step = max(s for s, _, _ in boundaries)
        if max_step != step:
            continue  # Only include compounds where this is the terminal step

        step_boundaries = [(p, n) for s, p, n in boundaries if s == step]
        if not step_boundaries:
            continue  # Compound doesn't have this step

        y_list.append(np.log(purity))

        row = np.zeros(len(interface_list))
        for boundary in step_boundaries:
            if boundary in interface_to_idx:
                row[interface_to_idx[boundary]] = 1.0
        X_list.append(row)

    return (
        np.array(y_list),
        np.array(X_list),
        interface_list,
        dict(interface_counts)
    )


def build_design_matrix(
    records: List[Dict],
) -> Tuple[np.ndarray, np.ndarray, List[Tuple[str, str]], Dict[Tuple[str, str], int]]:
    """
    Build design matrix for log-linear regression.

    Parameters
    ----------
    records : List[Dict]
        JSONL records with peaks and building_blocks

    Returns
    -------
    y : np.ndarray
        log(purity) for each compound (excluding purity=0)
    X : np.ndarray
        Binary matrix [n_compounds × n_interfaces]
    interface_list : List[Tuple[str, str]]
        List of (prev_block, next_block) interface tuples
    interface_counts : Dict
        Count of compounds containing each interface
    """
    # First pass: collect all unique interfaces and count them
    interface_counts: Dict[Tuple[str, str], int] = defaultdict(int)
    for record in records:
        for boundary in extract_boundaries(record):
            interface_counts[boundary] += 1

    interface_list = sorted(interface_counts.keys())
    interface_to_idx = {iface: i for i, iface in enumerate(interface_list)}

    # Second pass: build matrix
    y_list = []
    X_list = []

    for record in records:
        purity = PurityCalculator.calculate_from_peaks(record.get("peaks", []))
        if purity <= 0:
            continue  # Can't take log of zero

        y_list.append(np.log(purity))

        row = np.zeros(len(interface_list))
        for boundary in extract_boundaries(record):
            row[interface_to_idx[boundary]] = 1.0
        X_list.append(row)

    return (
        np.array(y_list),
        np.array(X_list),
        interface_list,
        dict(interface_counts)
    )


# =============================================================================
# Regression Fitting
# =============================================================================

def fit_interface_effects(
    y: np.ndarray,
    X: np.ndarray,
    regularization: str = "ridge"
) -> Tuple[np.ndarray, dict, float]:
    """
    Fit log-linear model using cross-validated regularized regression.

    Model: log(purity) = β₀ + Σ βᵢ × Iᵢ + ε

    Parameters
    ----------
    y : np.ndarray
        log(purity) values
    X : np.ndarray
        Binary design matrix
    regularization : str
        "ridge" for L2 only, "elastic" for L1+L2, "kernel" for Kernel Ridge

    Returns
    -------
    coefficients : np.ndarray
        β values for each interface (None for kernel)
    reg_params : dict
        Regularization parameters selected by CV
    r_squared : float
        R² score of the fitted model
    """
    # Search across wide range of λ values
    alphas = np.logspace(-6, 6, 50)

    if regularization == "kernel":
        # Kernel Ridge with polynomial kernel (implicit interactions)
        param_grid = {
            "alpha": np.logspace(-3, 3, 10),
            "kernel": ["poly"],
            "degree": [2, 3],
            "coef0": [1],
        }
        model = GridSearchCV(
            KernelRidge(),
            param_grid,
            cv=5,
            scoring='r2',
            n_jobs=-1
        )
        model.fit(X, y)
        best_model = model.best_estimator_
        reg_params = {
            "alpha": model.best_params_["alpha"],
            "kernel": model.best_params_["kernel"],
            "degree": model.best_params_["degree"],
            "cv_r2": model.best_score_,
            "n_nonzero": None  # Not applicable for kernel methods
        }
        # Compute training R²
        y_pred = best_model.predict(X)
        ss_res = np.sum((y - y_pred) ** 2)
        ss_tot = np.sum((y - np.mean(y)) ** 2)
        r_squared = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0.0
        # Kernel Ridge doesn't have explicit coefficients
        return None, reg_params, r_squared

    elif regularization == "elastic":
        # Elastic Net: L1 + L2 regularization
        # l1_ratio: 0 = pure L2 (Ridge), 1 = pure L1 (Lasso)
        l1_ratios = [0.1, 0.5, 0.7, 0.9, 0.95, 0.99]
        model = ElasticNetCV(
            l1_ratio=l1_ratios,
            alphas=alphas,
            cv=5,
            fit_intercept=True,
            max_iter=10000,
            random_state=42
        )
        model.fit(X, y)
        reg_params = {
            "alpha": model.alpha_,
            "l1_ratio": model.l1_ratio_,
            "n_nonzero": np.sum(model.coef_ != 0)
        }
    else:
        # Ridge: L2 only
        model = RidgeCV(alphas=alphas, cv=5, fit_intercept=True, scoring='r2')
        model.fit(X, y)
        reg_params = {
            "alpha": model.alpha_,
            "l1_ratio": 0.0,
            "n_nonzero": len(model.coef_)
        }

    # Compute R² on training data
    y_pred = model.predict(X)
    ss_res = np.sum((y - y_pred) ** 2)
    ss_tot = np.sum((y - np.mean(y)) ** 2)
    r_squared = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0.0

    return model.coef_, reg_params, r_squared


# =============================================================================
# Monomer Property Parsing (for sorting)
# =============================================================================

def extract_base_amino_acid(monomer: str) -> str:
    """Extract base amino acid from modified monomer name."""
    base = monomer

    if base.startswith("D") and len(base) > 1 and base[1].isupper():
        base = base[1:]
    if base.endswith("Me"):
        base = base[:-2]
    if "βHomo" in base:
        base = base.replace("βHomo", "")

    return base or monomer


def parse_monomer_properties(monomer: str) -> dict:
    """Extract chemical properties from monomer name."""
    is_peptoid = (
        monomer.startswith("LA") and
        len(monomer) > 2 and
        monomer[2:].isdigit()
    )

    return {
        "is_d": monomer.startswith("D") and len(monomer) > 1 and monomer[1].isupper(),
        "is_nme": monomer.endswith("Me"),
        "is_peptoid": is_peptoid,
        "is_beta": "βHomo" in monomer,
        "base": extract_base_amino_acid(monomer),
    }


def block_sort_key_by_first(block: str) -> tuple:
    """Sort key based on FIRST monomer (C-terminus, rightmost in string).

    For block "Leu-DLeu-Pro": FIRST = Pro (C-terminus, coupled first in synthesis)

    Uses chemical classification: L → D → L-NMe → D-NMe → L-βHomo → D-βHomo → peptoids
    Within each group, sorts by amino acid: Ala → Nvl → Val → Leu → Phe → Pro
    """
    if block == "N":
        return (-1, -1, "")  # N first

    monomers = block.split("-") if "-" in block else [block]
    group, aa = classify_monomer(monomers[-1])  # FIRST = rightmost = C-terminus
    return (group, aa, block)


def block_sort_key_by_last(block: str) -> tuple:
    """Sort key based on LAST monomer (N-terminus, leftmost in string).

    For block "Leu-DLeu-Pro": LAST = Leu (N-terminus, coupled last in synthesis)

    Uses chemical classification: L → D → L-NMe → D-NMe → L-βHomo → D-βHomo → peptoids
    Within each group, sorts by amino acid: Ala → Nvl → Val → Leu → Phe → Pro
    """
    if block == "N":
        return (-1, -1, "")  # N first

    monomers = block.split("-") if "-" in block else [block]
    group, aa = classify_monomer(monomers[0])  # LAST = leftmost = N-terminus
    return (group, aa, block)


# =============================================================================
# Main Processing
# =============================================================================

def load_and_compute(
    jsonl_path: Path,
    sort_by_interface: bool = False,
    regularization: str = "ridge",
) -> Tuple[Dict[Tuple[str, str], InterfaceEffect], List[str], List[str], dict]:
    """
    Load JSONL and compute interface effects via log-linear regression.

    Parameters
    ----------
    jsonl_path : Path
        Path to JSONL file
    sort_by_interface : bool
        Whether to sort by interface monomer properties
    regularization : str
        "ridge" for L2 only, "elastic" for L1+L2 (Elastic Net)

    Returns
    -------
    interface_effects : Dict
        Mapping (prev_block, next_block) -> InterfaceEffect
    row_blocks : List[str]
        Ordered list of prev blocks (for rows)
    col_blocks : List[str]
        Ordered list of next blocks (for columns)
    metadata : dict
        Additional info (reg_params, r_squared, n_compounds, n_interfaces)
    """
    # Load records
    records = []
    skipped = 0

    with open(jsonl_path) as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if not line:
                continue
            try:
                records.append(json.loads(line))
            except json.JSONDecodeError as e:
                print(f"Warning: Skipping line {line_num}: {e}")
                skipped += 1

    if skipped > 0:
        print(f"Skipped {skipped} records due to JSON errors")

    print(f"Loaded {len(records)} records")

    # Build design matrix
    y, X, interface_list, interface_counts = build_design_matrix(records)

    print(f"Built design matrix: {X.shape[0]} compounds × {X.shape[1]} interfaces")
    print(f"Compounds with purity > 0: {len(y)}")

    if len(y) == 0:
        raise ValueError("No compounds with purity > 0")

    # Fit regression
    coefficients, reg_params, r_squared = fit_interface_effects(y, X, regularization)

    if regularization == "kernel":
        print(f"Fitted model: kernel={reg_params['kernel']}, degree={reg_params['degree']}, α={reg_params['alpha']:.2e}")
        print(f"  CV R²: {reg_params['cv_r2']:.3f}, Training R²: {r_squared:.3f}")
    elif regularization == "elastic":
        print(f"Fitted model: α = {reg_params['alpha']:.2e}, l1_ratio = {reg_params['l1_ratio']:.2f}, R² = {r_squared:.3f}")
        print(f"  Non-zero coefficients: {reg_params['n_nonzero']} / {len(coefficients)}")
    else:
        print(f"Fitted model: λ = {reg_params['alpha']:.2e}, R² = {r_squared:.3f}")

    # Build interface effects dict (None for kernel methods)
    interface_effects = {}
    if coefficients is not None:
        for i, iface in enumerate(interface_list):
            interface_effects[iface] = InterfaceEffect(
                coefficient=coefficients[i],
                count=interface_counts[iface]
            )
    else:
        # For kernel methods, no explicit coefficients
        for iface in interface_list:
            interface_effects[iface] = InterfaceEffect(
                coefficient=0.0,  # Placeholder
                count=interface_counts[iface]
            )

    # Get unique prev and next blocks
    all_prev_blocks = set(iface[0] for iface in interface_list)
    all_next_blocks = set(iface[1] for iface in interface_list)

    # Sort blocks by interface-facing monomer
    # prev_block (on chain): interface = N-terminus = LAST monomer (leftmost)
    # next_block (being coupled): interface = C-terminus = FIRST monomer (rightmost)
    if sort_by_interface:
        row_blocks = sorted(all_prev_blocks, key=block_sort_key_by_last)   # N-terminus faces interface
        col_blocks = sorted(all_next_blocks, key=block_sort_key_by_first)  # C-terminus faces interface
    else:
        row_blocks = sorted(all_prev_blocks)
        col_blocks = sorted(all_next_blocks)

    metadata = {
        "reg_params": reg_params,
        "r_squared": r_squared,
        "n_compounds": len(y),
        "n_interfaces": len(interface_list),
        "regularization": regularization,
    }

    return interface_effects, row_blocks, col_blocks, metadata


def load_and_compute_per_step(
    jsonl_path: Path,
    step: int,
    sort_by_interface: bool = False,
    regularization: str = "ridge",
) -> Tuple[Dict[Tuple[str, str], InterfaceEffect], List[str], List[str], dict]:
    """
    Load JSONL and compute interface effects for a specific synthesis step.

    Parameters
    ----------
    jsonl_path : Path
        Path to JSONL file
    step : int
        Synthesis step (1-indexed)
    sort_by_interface : bool
        Whether to sort by interface monomer properties
    regularization : str
        "ridge" for L2 only, "elastic" for L1+L2

    Returns
    -------
    interface_effects : Dict
        Mapping (prev_block, next_block) -> InterfaceEffect
    row_blocks : List[str]
        Ordered list of prev blocks (for rows)
    col_blocks : List[str]
        Ordered list of next blocks (for columns)
    metadata : dict
        Additional info (reg_params, r_squared, n_compounds, n_interfaces)
    """
    # Load records
    records = []
    skipped = 0

    with open(jsonl_path) as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if not line:
                continue
            try:
                records.append(json.loads(line))
            except json.JSONDecodeError as e:
                print(f"Warning: Skipping line {line_num}: {e}")
                skipped += 1

    if skipped > 0:
        print(f"Skipped {skipped} records due to JSON errors")

    print(f"Loaded {len(records)} records")

    # Build design matrix for specific step
    y, X, interface_list, interface_counts = build_design_matrix_per_step(records, step)

    print(f"Step {step}: {X.shape[0]} compounds × {X.shape[1]} interfaces")
    print(f"Compounds with purity > 0 at step {step}: {len(y)}")

    if len(y) == 0:
        raise ValueError(f"No compounds with purity > 0 at step {step}")

    if X.shape[1] == 0:
        raise ValueError(f"No interfaces found at step {step}")

    # Fit regression
    coefficients, reg_params, r_squared = fit_interface_effects(y, X, regularization)

    if regularization == "elastic":
        print(f"Fitted model: α = {reg_params['alpha']:.2e}, l1_ratio = {reg_params['l1_ratio']:.2f}, R² = {r_squared:.3f}")
        print(f"  Non-zero coefficients: {reg_params['n_nonzero']} / {len(coefficients)}")
    else:
        print(f"Fitted model: λ = {reg_params['alpha']:.2e}, R² = {r_squared:.3f}")

    # Build interface effects dict
    interface_effects = {}
    if coefficients is not None:
        for i, iface in enumerate(interface_list):
            interface_effects[iface] = InterfaceEffect(
                coefficient=coefficients[i],
                count=interface_counts[iface]
            )

    # Get unique prev and next blocks
    all_prev_blocks = set(iface[0] for iface in interface_list)
    all_next_blocks = set(iface[1] for iface in interface_list)

    # Sort blocks by interface-facing monomer
    if sort_by_interface:
        row_blocks = sorted(all_prev_blocks, key=block_sort_key_by_last)
        col_blocks = sorted(all_next_blocks, key=block_sort_key_by_first)
    else:
        row_blocks = sorted(all_prev_blocks)
        col_blocks = sorted(all_next_blocks)

    metadata = {
        "reg_params": reg_params,
        "r_squared": r_squared,
        "n_compounds": len(y),
        "n_interfaces": len(interface_list),
        "regularization": regularization,
        "step": step,
    }

    return interface_effects, row_blocks, col_blocks, metadata


def get_max_steps(jsonl_path: Path) -> int:
    """
    Determine the maximum number of synthesis steps in the dataset.

    Returns
    -------
    int
        Maximum step number (1-indexed)
    """
    max_step = 0
    with open(jsonl_path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            try:
                record = json.loads(line)
                boundaries = extract_boundaries_with_step(record)
                for step, _, _ in boundaries:
                    max_step = max(max_step, step)
            except json.JSONDecodeError:
                continue
    return max_step


def load_and_compute_compositional(
    jsonl_path: Path,
) -> Tuple[np.ndarray, List[str], dict]:
    """
    Load JSONL and compute compositional feature effects via log-linear regression.

    Parameters
    ----------
    jsonl_path : Path
        Path to JSONL file

    Returns
    -------
    coefficients : np.ndarray
        Regression coefficients for each feature
    feature_names : List[str]
        Names for each feature (e.g., "prev_pos0_Leu")
    metadata : dict
        Additional info (lambda_chosen, r_squared, vocabulary, etc.)
    """
    # Load records
    records = []
    skipped = 0

    with open(jsonl_path) as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if not line:
                continue
            try:
                records.append(json.loads(line))
            except json.JSONDecodeError as e:
                print(f"Warning: Skipping line {line_num}: {e}")
                skipped += 1

    if skipped > 0:
        print(f"Skipped {skipped} records due to JSON errors")

    print(f"Loaded {len(records)} records")

    # Build compositional design matrix
    y, X, feature_names, comp_meta = build_compositional_design_matrix(records)

    print(f"Built compositional design matrix: {X.shape[0]} compounds × {X.shape[1]} features")
    print(f"  Vocabulary: {comp_meta['vocab_size']} monomers")
    print(f"  Max block length: {comp_meta['max_block_length']}")
    print(f"Compounds with purity > 0: {len(y)}")

    if len(y) == 0:
        raise ValueError("No compounds with purity > 0")

    # Fit regression (always use ridge for compositional)
    coefficients, reg_params, r_squared = fit_interface_effects(y, X, "ridge")

    print(f"Fitted model: λ = {reg_params['alpha']:.2e}, R² = {r_squared:.3f}")

    metadata = {
        "lambda_chosen": reg_params['alpha'],
        "r_squared": r_squared,
        "n_compounds": len(y),
        "n_features": len(feature_names),
        "vocabulary": comp_meta["vocabulary"],
        "vocab_size": comp_meta["vocab_size"],
        "max_block_length": comp_meta["max_block_length"],
    }

    return coefficients, feature_names, metadata


def load_and_compute_tree(
    jsonl_path: Path,
) -> Tuple[object, List[str], dict]:
    """
    Load JSONL and compute effects using gradient boosting with block identity encoding.

    Uses HistGradientBoostingRegressor on the same design matrix as Ridge
    (1 row per compound, binary interface indicators).

    Parameters
    ----------
    jsonl_path : Path
        Path to JSONL file

    Returns
    -------
    model : HistGradientBoostingRegressor
        Fitted model
    interface_list : List[Tuple[str, str]]
        List of interface tuples
    metadata : dict
        CV scores, interface counts, etc.
    """
    from sklearn.ensemble import HistGradientBoostingRegressor
    from sklearn.model_selection import cross_val_score

    # Load records
    records = []
    skipped = 0

    with open(jsonl_path) as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if not line:
                continue
            try:
                records.append(json.loads(line))
            except json.JSONDecodeError as e:
                print(f"Warning: Skipping line {line_num}: {e}")
                skipped += 1

    if skipped > 0:
        print(f"Skipped {skipped} records due to JSON errors")

    print(f"Loaded {len(records)} records")

    # Build block identity design matrix (same as Ridge)
    y, X, interface_list, interface_counts = build_design_matrix(records)

    print(f"Built design matrix: {X.shape[0]} compounds × {X.shape[1]} interfaces")

    if len(y) == 0:
        raise ValueError("No compounds with purity > 0")

    # Create model - more iterations for sparse features
    model = HistGradientBoostingRegressor(
        max_iter=500,
        max_depth=6,
        learning_rate=0.05,
        min_samples_leaf=20,
        random_state=42
    )

    # Cross-validation for R²
    print("Running 5-fold cross-validation...")
    cv_scores = cross_val_score(model, X, y, cv=5, scoring='r2')
    cv_r2_mean = cv_scores.mean()
    cv_r2_std = cv_scores.std()

    print(f"CV R²: {cv_r2_mean:.3f} ± {cv_r2_std:.3f}")

    # Fit final model on all data
    model.fit(X, y)

    # Compute training R²
    y_pred = model.predict(X)
    ss_res = np.sum((y - y_pred) ** 2)
    ss_tot = np.sum((y - np.mean(y)) ** 2)
    train_r2 = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0.0

    print(f"Training R²: {train_r2:.3f}")

    metadata = {
        "cv_r2_mean": cv_r2_mean,
        "cv_r2_std": cv_r2_std,
        "train_r2": train_r2,
        "n_compounds": len(y),
        "n_interfaces": len(interface_list),
        "interface_counts": interface_counts,
    }

    return model, interface_list, metadata


# =============================================================================
# Plotting
# =============================================================================

def plot_coefficient_heatmap(
    interface_effects: Dict[Tuple[str, str], InterfaceEffect],
    row_blocks: List[str],
    col_blocks: List[str],
    title: str = None,
    figsize: Tuple[int, int] = None,
    show_counts: bool = True,
    min_count: int = 1,
    show_exp: bool = False,
    warning_text: str = None,
) -> Figure:
    """
    Plot block × block heatmap of regression coefficients.

    Parameters
    ----------
    interface_effects : Dict
        Mapping (prev_block, next_block) -> InterfaceEffect
    row_blocks : List[str]
        Ordered list of prev blocks
    col_blocks : List[str]
        Ordered list of next blocks
    title : str, optional
        Plot title
    figsize : Tuple[int, int], optional
        Figure size
    show_counts : bool
        Show observation counts in cells
    min_count : int
        Minimum observations to show a cell
    show_exp : bool
        If True, show exp(β) instead of β

    Returns
    -------
    Figure
    """
    n_rows = len(row_blocks)
    n_cols = len(col_blocks)

    if figsize is None:
        figsize = (max(10, n_cols * 0.5), max(8, n_rows * 0.4))

    # Build matrix
    matrix = np.full((n_rows, n_cols), np.nan)
    counts = np.zeros((n_rows, n_cols), dtype=int)

    row_idx = {block: i for i, block in enumerate(row_blocks)}
    col_idx = {block: j for j, block in enumerate(col_blocks)}

    for (prev_block, next_block), effect in interface_effects.items():
        if prev_block in row_idx and next_block in col_idx:
            i = row_idx[prev_block]
            j = col_idx[next_block]
            if effect.count >= min_count:
                if show_exp:
                    matrix[i, j] = effect.multiplicative_effect
                else:
                    matrix[i, j] = effect.coefficient
                counts[i, j] = effect.count

    # Create figure
    fig, ax = plt.subplots(figsize=figsize)

    # Determine colormap and scale
    if show_exp:
        # Multiplicative scale: center at 1.0
        vmin, vmax = 0.5, 1.5
        cmap = 'RdYlGn'
        cbar_label = 'Multiplicative Effect exp(β)'
    else:
        # Log scale: center at 0
        max_abs = np.nanmax(np.abs(matrix)) if not np.all(np.isnan(matrix)) else 1.0
        vmin, vmax = -max_abs, max_abs
        cmap = 'RdYlGn'
        cbar_label = 'Coefficient β (log-scale)'

    im = ax.imshow(matrix, cmap=cmap, vmin=vmin, vmax=vmax, aspect='auto')

    # Colorbar
    cbar = plt.colorbar(im, ax=ax, shrink=0.8)
    cbar.set_label(cbar_label, fontsize=10)

    # Labels
    ax.set_xticks(range(n_cols))
    ax.set_xticklabels(col_blocks, rotation=45, ha='right', fontsize=8)
    ax.set_yticks(range(n_rows))
    ax.set_yticklabels(row_blocks, fontsize=8)

    ax.set_xlabel('Next Block (being coupled) - sorted by FIRST monomer (C-term)', fontsize=10)
    ax.set_ylabel('Prev Block (on chain) - sorted by LAST monomer (N-term)', fontsize=10)

    if title:
        ax.set_title(title, fontsize=12, fontweight='bold')
    else:
        ax.set_title('Interface Effect on Purity', fontsize=12, fontweight='bold')

    # Add warning text if provided
    if warning_text:
        fig.text(
            0.5, 0.02, warning_text,
            ha='center', va='bottom',
            fontsize=9, color='red', fontweight='bold',
            bbox=dict(boxstyle='round', facecolor='lightyellow', edgecolor='red', alpha=0.9)
        )

    # Annotations
    if show_counts:
        for i in range(n_rows):
            for j in range(n_cols):
                if not np.isnan(matrix[i, j]):
                    val = matrix[i, j]
                    count = counts[i, j]

                    # Text color based on value
                    if show_exp:
                        text_color = 'white' if val < 0.7 or val > 1.3 else 'black'
                    else:
                        text_color = 'white' if abs(val) > max_abs * 0.5 else 'black'

                    ax.text(
                        j, i, f'n={count}',
                        ha='center', va='center',
                        fontsize=5, color=text_color
                    )

    plt.tight_layout()
    return fig


# =============================================================================
# Main
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Generate coupling heatmap using log-linear regression",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Model:
    log(purity) = β₀ + Σ β_interface × I(interface present) + ε

Each cell shows the β coefficient for that interface:
    β < 0: Interface reduces purity (bad coupling)
    β > 0: Interface increases purity (good coupling)
    β = 0: Interface has no effect

Use --show-exp to display exp(β) instead (multiplicative effect on purity).

Examples:
    python examples/coupling_heatmap.py \\
        --input results/library/library_analysis.jsonl \\
        --output coupling_heatmap.png \\
        --sort-by-interface
        """,
    )

    parser.add_argument("--input", type=Path, required=True, help="Path to JSONL file")
    parser.add_argument("--output", type=Path, default=Path("coupling_heatmap.png"), help="Output path")
    parser.add_argument("--title", type=str, default=None, help="Custom title")
    parser.add_argument("--no-counts", action="store_true", help="Hide observation counts")
    parser.add_argument("--min-count", type=int, default=1, help="Minimum observations per cell")
    parser.add_argument("--figsize", type=float, nargs=2, default=None, metavar=("W", "H"))
    parser.add_argument("--sort-by-interface", action="store_true", help="Sort by interface monomer properties")
    parser.add_argument("--show-exp", action="store_true", help="Show exp(β) instead of β")
    parser.add_argument(
        "--compositional",
        action="store_true",
        help="Use compositional encoding (positional one-hot) instead of block identity"
    )
    parser.add_argument(
        "--model",
        choices=["ridge", "tree"],
        default="ridge",
        help="Model type: ridge (linear) or tree (gradient boosting)"
    )
    parser.add_argument(
        "--regularization",
        choices=["ridge", "elastic", "kernel"],
        default="ridge",
        help="Regularization: ridge (L2), elastic (L1+L2), or kernel (polynomial kernel)"
    )
    parser.add_argument(
        "--per-step",
        action="store_true",
        help="Generate separate heatmaps for each synthesis step (output becomes directory)"
    )

    args = parser.parse_args()

    if not args.input.exists():
        print(f"Error: Input file not found: {args.input}")
        sys.exit(1)

    # Load and compute
    print(f"Loading: {args.input}")

    # Per-step mode: generate separate heatmaps for each synthesis step
    if args.per_step:
        if args.model == "tree" or args.compositional or args.regularization == "kernel":
            print("Error: --per-step is only supported with block identity mode (ridge or elastic)")
            sys.exit(1)

        # Determine max steps
        max_steps = get_max_steps(args.input)
        print(f"Found {max_steps} synthesis steps in dataset")

        if max_steps == 0:
            print("Error: No synthesis steps found in dataset")
            sys.exit(1)

        # Create output directory
        output_dir = args.output.parent / args.output.stem
        output_dir.mkdir(parents=True, exist_ok=True)

        # Generate heatmap for each step (1+ = block->block couplings)
        # Skip step 0 (N→first_block) as it's not a true block-block interface
        # and typically has too few samples for reliable regression
        for step in range(1, max_steps + 1):
            print(f"\n{'='*60}")
            print(f"Processing Step {step} / {max_steps}")
            print(f"{'='*60}")

            try:
                interface_effects, row_blocks, col_blocks, metadata = load_and_compute_per_step(
                    args.input,
                    step=step,
                    sort_by_interface=args.sort_by_interface,
                    regularization=args.regularization,
                )
            except ValueError as e:
                print(f"Warning: {e}")
                continue

            # Check for insufficient data
            n_compounds = metadata['n_compounds']
            n_interfaces = metadata['n_interfaces']
            samples_per_feature = n_compounds / n_interfaces if n_interfaces > 0 else 0
            MIN_SAMPLES_PER_FEATURE = 10

            is_underpowered = samples_per_feature < MIN_SAMPLES_PER_FEATURE

            if is_underpowered:
                print(f"  WARNING: Insufficient data ({samples_per_feature:.1f} samples/interface, need ≥{MIN_SAMPLES_PER_FEATURE})")

            # Build title
            model_name = "Elastic Net" if args.regularization == "elastic" else "Ridge"
            if args.title:
                plot_title = f"{args.title} - Step {step}"
            else:
                plot_title = f"Step {step} Interface Effect ({model_name}, R²={metadata['r_squared']:.3f})"

            # Add warning subtitle if underpowered
            warning_text = None
            if is_underpowered:
                warning_text = (
                    f"⚠️ UNDERPOWERED: {samples_per_feature:.1f} samples/interface "
                    f"(need ≥{MIN_SAMPLES_PER_FEATURE} for reliable estimates)"
                )

            # Generate plot
            fig = plot_coefficient_heatmap(
                interface_effects=interface_effects,
                row_blocks=row_blocks,
                col_blocks=col_blocks,
                title=plot_title,
                figsize=tuple(args.figsize) if args.figsize else None,
                show_counts=not args.no_counts,
                min_count=args.min_count,
                show_exp=args.show_exp,
                warning_text=warning_text,
            )

            output_file = output_dir / f"step_{step}_coupling_heatmap.png"
            fig.savefig(output_file, dpi=300, bbox_inches='tight')
            plt.close(fig)
            print(f"Saved: {output_file}")

            # Summary for this step
            if interface_effects:
                sorted_effects = sorted(
                    interface_effects.items(),
                    key=lambda x: x[1].coefficient,
                    reverse=True
                )
                print(f"\n  Best interfaces for step {step}:")
                for iface, effect in sorted_effects[:2]:
                    print(f"    {iface[0]} -> {iface[1]}: β={effect.coefficient:.3f}, n={effect.count}")
                print(f"  Worst interfaces for step {step}:")
                for iface, effect in sorted_effects[-2:]:
                    print(f"    {iface[0]} -> {iface[1]}: β={effect.coefficient:.3f}, n={effect.count}")

        print(f"\n{'='*60}")
        print(f"Per-step heatmaps saved to: {output_dir}/")
        print(f"{'='*60}")
        return

    if args.model == "tree":
        # Gradient boosting with block identity encoding
        model, interface_list, metadata = load_and_compute_tree(args.input)

        print(f"\nRegression results (gradient boosting, block identity):")
        print(f"  CV R²: {metadata['cv_r2_mean']:.3f} ± {metadata['cv_r2_std']:.3f}")
        print(f"  Training R²: {metadata['train_r2']:.3f}")
        print(f"  Compounds: {metadata['n_compounds']}")
        print(f"  Interfaces: {metadata['n_interfaces']}")

        # Note about feature importance with many sparse features
        print(f"\nNote: With {metadata['n_interfaces']} sparse binary features,")
        print(f"individual feature importances are less meaningful.")
        print(f"The CV R² indicates model quality.")

        return

    if args.compositional:
        # Compositional encoding mode
        coefficients, feature_names, metadata = load_and_compute_compositional(args.input)

        print(f"\nRegression results (compositional encoding):")
        print(f"  λ (regularization): {metadata['lambda_chosen']:.2e}")
        print(f"  R²: {metadata['r_squared']:.3f}")
        print(f"  Compounds: {metadata['n_compounds']}")
        print(f"  Features: {metadata['n_features']}")

        # Print top/bottom features
        sorted_idx = np.argsort(coefficients)

        print(f"\nTop features (highest β):")
        for i in sorted_idx[-10:]:
            if coefficients[i] != 0:
                print(f"  {feature_names[i]}: β={coefficients[i]:.3f}, exp(β)={np.exp(coefficients[i]):.3f}")

        print(f"\nBottom features (lowest β):")
        for i in sorted_idx[:10]:
            if coefficients[i] != 0:
                print(f"  {feature_names[i]}: β={coefficients[i]:.3f}, exp(β)={np.exp(coefficients[i]):.3f}")

        print(f"\nNote: Compositional mode does not generate heatmap (features are not block pairs)")
        return

    # Block identity mode
    interface_effects, row_blocks, col_blocks, metadata = load_and_compute(
        args.input,
        sort_by_interface=args.sort_by_interface,
        regularization=args.regularization,
    )

    reg_params = metadata['reg_params']
    print(f"\nRegression results (block identity, {args.regularization}):")
    if args.regularization == "kernel":
        print(f"  Kernel: {reg_params['kernel']}, degree: {reg_params['degree']}")
        print(f"  α: {reg_params['alpha']:.2e}")
        print(f"  CV R²: {reg_params['cv_r2']:.3f}")
        print(f"  Training R²: {metadata['r_squared']:.3f}")
    elif args.regularization == "elastic":
        print(f"  α: {reg_params['alpha']:.2e}, l1_ratio: {reg_params['l1_ratio']:.2f}")
        print(f"  Non-zero: {reg_params['n_nonzero']} / {metadata['n_interfaces']} interfaces")
        print(f"  R²: {metadata['r_squared']:.3f}")
    else:
        print(f"  λ (regularization): {reg_params['alpha']:.2e}")
        print(f"  R²: {metadata['r_squared']:.3f}")
    print(f"  Compounds: {metadata['n_compounds']}")
    print(f"  Interfaces: {metadata['n_interfaces']}")

    # Kernel mode doesn't have explicit coefficients, so skip heatmap
    if args.regularization == "kernel":
        print(f"\nNote: Kernel Ridge doesn't produce explicit interface coefficients.")
        print(f"The CV R² of {reg_params['cv_r2']:.3f} indicates model quality.")
        print(f"Use --regularization ridge or elastic to generate heatmap.")
        return

    # Generate plot
    print(f"\nGenerating heatmap...")

    # Build title with model info
    if args.title:
        plot_title = args.title
    else:
        model_name = "Elastic Net" if args.regularization == "elastic" else "Ridge"
        plot_title = f"Interface Effect on Purity ({model_name}, R²={metadata['r_squared']:.3f})"

    fig = plot_coefficient_heatmap(
        interface_effects=interface_effects,
        row_blocks=row_blocks,
        col_blocks=col_blocks,
        title=plot_title,
        figsize=tuple(args.figsize) if args.figsize else None,
        show_counts=not args.no_counts,
        min_count=args.min_count,
        show_exp=args.show_exp,
    )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved: {args.output}")

    # Summary
    if interface_effects:
        sorted_effects = sorted(
            interface_effects.items(),
            key=lambda x: x[1].coefficient,
            reverse=True
        )
        best = sorted_effects[0]
        worst = sorted_effects[-1]

        print(f"\nTop interfaces (highest β):")
        for iface, effect in sorted_effects[:3]:
            print(f"  {iface[0]} -> {iface[1]}: β={effect.coefficient:.3f}, exp(β)={effect.multiplicative_effect:.3f}, n={effect.count}")

        print(f"\nBottom interfaces (lowest β):")
        for iface, effect in sorted_effects[-3:]:
            print(f"  {iface[0]} -> {iface[1]}: β={effect.coefficient:.3f}, exp(β)={effect.multiplicative_effect:.3f}, n={effect.count}")


if __name__ == "__main__":
    main()
