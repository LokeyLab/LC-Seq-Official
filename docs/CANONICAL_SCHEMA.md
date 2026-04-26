# LC-Seq Canonical Data Schema

This document defines the internal data schema used by LC-Seq. All input data formats must be transformed to this canonical schema before processing.

## Schema Overview

```json
{
  "metadata": {
    "source_file": "string",
    "source_format": "string",
    "libraries": ["DEL-0044", "DEL-0045"],
    "building_block_order": "N_to_C | C_to_N",
    "time_unit": "seconds | minutes",
    "version": "1.0"
  },
  "compounds": [
    {
      "identity": { ... },
      "building_blocks": [ ... ],
      "properties": { ... },
      "retention_times": { ... },
      "chromatogram": { ... }
    }
  ]
}
```

---

## Compound Schema

### 1. Identity

Compound identification and structural information.

```json
{
  "identity": {
    "name": "DβHomoleu-AgxNull-AgxNull",
    "id": "LDEL00002",
    "library_id": "DEL-0045",
    "smiles": "C([C@H](CC(=O)O)NC(OCC1...)=O)C(C)C",
    "product_id": "DEL-0045;BB17416;BB17646;BB21019",
    "stereochemistry": "DB,x,x"
  }
}
```

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| `name` | string | Yes | Human-readable compound name |
| `id` | string | Yes | Unique compound identifier |
| `library_id` | string | Yes | Library this compound belongs to |
| `smiles` | string | No | Full compound SMILES |
| `product_id` | string | No | Extended product ID (library + BBs) |
| `stereochemistry` | string | No | Stereochemistry string |

---

### 2. Building Blocks

Ordered list of building blocks (in synthesis order, typically N→C for peptides).

```json
{
  "building_blocks": [
    {
      "position": 1,
      "name": "DβHomoleu",
      "code": "BB17416",
      "smiles": "C([C@H](CC(=O)O)NC...)C(C)C",
      "stereochemistry": "DB",
      "properties": {
        "n_substitution": 0,
        "aa_count": 0
      }
    },
    {
      "position": 2,
      "name": "AgxNull",
      "code": "BB17646",
      "smiles": "AGXnull-0001",
      "stereochemistry": "x",
      "properties": {
        "n_substitution": 0,
        "aa_count": 0
      }
    }
  ]
}
```

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| `position` | int | Yes | Position in synthesis order (1-indexed) |
| `name` | string | Yes | Building block name |
| `code` | string | No | Building block catalog code |
| `smiles` | string | No | Building block SMILES |
| `stereochemistry` | string | No | Stereochemistry designation |
| `properties` | object | No | Additional BB-specific properties |

---

### 3. Chemical/Physical Properties

Computed or measured molecular properties.

```json
{
  "properties": {
    "exact_mass": {
      "cooh": 456.23,
      "dma": 469.25
    },
    "logp": 2.34,
    "h_donors": 3,
    "h_acceptors": 5,
    "aromatic_rings": 2,
    "tpsa": 89.5,
    "peptide_size": 3,
    "n_substitution_total": 1,
    "custom": {
      "any_other_property": "value"
    }
  }
}
```

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| `exact_mass` | object | No | Exact mass under different conditions |
| `logp` | float | No | Calculated LogP (AlogP, ClogP, etc.) |
| `h_donors` | int | No | Number of H-bond donors |
| `h_acceptors` | int | No | Number of H-bond acceptors |
| `aromatic_rings` | int | No | Number of aromatic rings |
| `tpsa` | float | No | Topological polar surface area |
| `peptide_size` | int | No | Number of amino acid equivalents |
| `n_substitution_total` | int | No | Total N-substitutions |
| `custom` | object | No | Any additional custom properties |

---

### 4. Retention Times

Reference retention time measurements.

```json
{
  "retention_times": {
    "linear": {
      "rt_minutes": 12.5,
      "logk": 0.85
    },
    "cyclized": {
      "rt_minutes": 13.1,
      "logk": 0.72
    },
    "delta_rt": 0.6,
    "weighted": {
      "combined": 1057.2,
      "by_library": {
        "DEL-0044": 962.5,
        "DEL-0045": 1057.2
      }
    }
  }
}
```

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| `linear` | object | No | Linear form retention data |
| `cyclized` | object | No | Cyclized form retention data |
| `delta_rt` | float | No | RT difference (cyclized - linear) |
| `weighted` | object | Yes* | Weighted RT from LC-Seq analysis |

*Required if chromatogram data is provided

---

### 5. Chromatogram Data

Time series data from LC-MS analysis.

```json
{
  "chromatogram": {
    "total_fractions": 96,
    "max_count": 1330,
    "time_unit": "seconds",
    "channels": ["DEL-0044", "DEL-0045"],
    "data": [
      {"time": 615.0, "signals": [295, 420]},
      {"time": 645.0, "signals": [370, 512]},
      {"time": 675.0, "signals": [522, 781]},
      {"time": 705.0, "signals": [722, 1133]}
    ]
  }
}
```

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| `total_fractions` | int | No | Total number of fractions collected |
| `max_count` | int | No | Maximum signal count |
| `time_unit` | string | Yes | "seconds" or "minutes" |
| `channels` | list | Yes | Names of signal channels |
| `data` | list | Yes | Time series: time + signals per channel |

---

## Data Mapping Configuration

To transform external data formats to the canonical schema, use a mapping configuration:

```yaml
# Example: mapping raw_data.csv to canonical schema
source:
  format: csv
  encoding: utf-8

mappings:
  identity:
    name: "Common_Name (N-->C)"
    id: "Identifier"
    library_id: "lid"
    smiles: "enumerated_smiles"
    product_id: "intended_product_id"
    stereochemistry: "All Stereochem (N-->C)"

  building_blocks:
    order: "N_to_C"  # BB1 is N-terminal
    blocks:
      - position: 1
        name: "BB1 Name"
        smiles: "bb1_smiles"
        stereochemistry: "BB1 Stereochem"
        properties:
          n_substitution: "BB1 N-Sub"
          aa_count: "BB1 AA Count"
      - position: 2
        name: "BB2 Name"
        smiles: "bb2_smiles"
        stereochemistry: "BB2 Stereochem"
        properties:
          n_substitution: "BB2 N-Sub"
          aa_count: "BB2 AA Count"
      - position: 3
        name: "BB3 Name"
        smiles: "bb3_smiles"
        stereochemistry: "BB3 Stereochem"
        properties:
          n_substitution: "BB3 N-Sub"
          aa_count: "BB3 AA Count"

  properties:
    exact_mass:
      cooh: "Exact Mass (COOH)"
      dma: "Exact Mass (DMA)"
    logp: "AlogP"
    h_donors: "NumHDonors"
    h_acceptors: "NumHAcceptors"
    aromatic_rings: "NumAromaticRings"
    tpsa: "TPSA"
    peptide_size: "Peptide Size"
    n_substitution_total: "SUM Nsub"

  retention_times:
    linear:
      rt_minutes: "Linear RT (min)"
      logk: "LogK10-40 (Linear)"
    cyclized:
      rt_minutes: "Cyclized RT (min)"
      logk: "LogK10-40 (Cyclized)"
    delta_rt: "Delta RT (min) (Cyclized-Linear)"
    weighted:
      combined: "wtd_rt"
      by_library:
        DEL-0044: "wtd_rt DEL-0044"
        DEL-0045: "wtd_rt DEL-0045"

  chromatogram:
    total_fractions: "total_fractions"
    max_count: "max_count"
    time_series:
      column: "all_datapoints"
      format: "time:sig1;sig2, time:sig1;sig2, ..."
      time_unit: seconds
      channels: ["DEL-0044", "DEL-0045"]
```

---

## Time Series Format Variants

The canonical schema supports multiple input time series formats:

### Format 1: Semicolon-Separated Channels (current test data)
```
time:signal1;signal2, time:signal1;signal2, ...
```
Example: `735.00:852;564, 885.00:352;194`

### Format 2: JSON Array
```json
[[time1, [sig1, sig2]], [time2, [sig1, sig2]], ...]
```
Example: `[[735.0, [852, 564]], [885.0, [352, 194]]]`

### Format 3: Separate Columns per Channel
Each channel in its own column with format `time:signal, time:signal, ...`

### Format 4: Long Format (one row per timepoint)
| compound_id | time | channel | signal |
|-------------|------|---------|--------|
| LDEL00002 | 735.0 | DEL-0044 | 852 |
| LDEL00002 | 735.0 | DEL-0045 | 564 |

---

## Validation Rules

1. **Required fields**: Every compound must have `identity.name`, `identity.id`, and `identity.library_id`
2. **Building block order**: Must be consistent across all compounds
3. **Time series alignment**: All channels must have the same time points
4. **Numeric fields**: Must be valid numbers (null/empty allowed for optional fields)
5. **Library consistency**: `library_id` must match one of the channels

---

## Example: Complete Compound Record

```json
{
  "identity": {
    "name": "DβHomoleu-AgxNull-AgxNull",
    "id": "LDEL00002",
    "library_id": "DEL-0045",
    "smiles": "C([C@H](CC(=O)O)NC(OCC1c(c(c(c2)c1ccc2)ccc3)c3)=O)C(C)C",
    "product_id": "DEL-0045;BB17416;BB17646;BB21019",
    "stereochemistry": "DB,x,x"
  },
  "building_blocks": [
    {
      "position": 1,
      "name": "DβHomoleu",
      "code": "BB17416",
      "smiles": "C([C@H](CC(=O)O)NC...)C(C)C",
      "stereochemistry": "DB",
      "properties": {"n_substitution": 0, "aa_count": 0}
    },
    {
      "position": 2,
      "name": "AgxNull",
      "code": "BB17646",
      "smiles": "AGXnull-0001",
      "stereochemistry": "x",
      "properties": {"n_substitution": 0, "aa_count": 0}
    },
    {
      "position": 3,
      "name": "AgxNull",
      "code": "BB21019",
      "smiles": "AGXnull-0002",
      "stereochemistry": "x",
      "properties": {"n_substitution": 0, "aa_count": 1}
    }
  ],
  "properties": {
    "exact_mass": {"cooh": -13.03, "dma": 14.01},
    "logp": null,
    "h_donors": null,
    "h_acceptors": null,
    "aromatic_rings": null,
    "tpsa": null,
    "peptide_size": 3,
    "n_substitution_total": 0
  },
  "retention_times": {
    "linear": {"rt_minutes": null, "logk": null},
    "cyclized": {"rt_minutes": 13.12, "logk": 0.72},
    "delta_rt": 0.41,
    "weighted": {
      "combined": 1057.21,
      "by_library": {
        "DEL-0044": 962.51,
        "DEL-0045": 1057.21
      }
    }
  },
  "chromatogram": {
    "total_fractions": 96,
    "max_count": 1330,
    "time_unit": "seconds",
    "channels": ["DEL-0044", "DEL-0045"],
    "data": [
      {"time": 615.0, "signals": [295, 420]},
      {"time": 645.0, "signals": [370, 512]},
      {"time": 675.0, "signals": [522, 781]}
    ]
  }
}
```
