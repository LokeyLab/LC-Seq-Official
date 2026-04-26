"""
Data transformer for converting external data formats to canonical schema.

Transforms CSV/Excel files with arbitrary column structures into the
LC-Seq canonical data format using YAML mapping configurations.
"""

import re
import json
import yaml
import numpy as np
import pandas as pd
from pathlib import Path
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Tuple, Union


@dataclass
class TimeSeriesConfig:
    """Configuration for parsing time series data."""
    column: str
    format: str  # "semicolon", "json", "separate_columns"
    time_unit: str = "seconds"
    channels: List[str] = field(default_factory=list)
    delimiter: str = ", "  # delimiter between timepoints


@dataclass
class BuildingBlockMapping:
    """Mapping for a single building block position."""
    position: int
    name: str
    smiles: Optional[str] = None
    code: Optional[str] = None
    stereochemistry: Optional[str] = None
    properties: Dict[str, str] = field(default_factory=dict)


@dataclass
class DataMapping:
    """Complete mapping configuration for transforming external data."""
    # Source configuration
    source_format: str = "csv"
    encoding: str = "utf-8"

    # Identity mappings (column name -> canonical field)
    identity: Dict[str, str] = field(default_factory=dict)

    # Building block configuration
    building_block_order: str = "N_to_C"
    building_blocks: List[BuildingBlockMapping] = field(default_factory=list)

    # Property mappings
    properties: Dict[str, Any] = field(default_factory=dict)

    # Retention time mappings
    retention_times: Dict[str, Any] = field(default_factory=dict)

    # Chromatogram configuration
    chromatogram: Dict[str, Any] = field(default_factory=dict)

    @classmethod
    def from_yaml(cls, path: Path) -> 'DataMapping':
        """Load mapping configuration from YAML file."""
        with open(path, 'r') as f:
            config = yaml.safe_load(f)

        mapping = cls()

        # Source config
        source = config.get('source', {})
        mapping.source_format = source.get('format', 'csv')
        mapping.encoding = source.get('encoding', 'utf-8')

        # Mappings
        mappings = config.get('mappings', {})

        # Identity
        mapping.identity = mappings.get('identity', {})

        # Building blocks
        bb_config = mappings.get('building_blocks', {})
        mapping.building_block_order = bb_config.get('order', 'N_to_C')

        for bb in bb_config.get('blocks', []):
            mapping.building_blocks.append(BuildingBlockMapping(
                position=bb['position'],
                name=bb['name'],
                smiles=bb.get('smiles'),
                code=bb.get('code'),
                stereochemistry=bb.get('stereochemistry'),
                properties=bb.get('properties', {})
            ))

        # Properties
        mapping.properties = mappings.get('properties', {})

        # Retention times
        mapping.retention_times = mappings.get('retention_times', {})

        # Chromatogram
        mapping.chromatogram = mappings.get('chromatogram', {})

        return mapping


class TimeSeriesParser:
    """Parse various time series formats."""

    @staticmethod
    def parse_semicolon_format(
        value: str,
        channels: List[str],
        delimiter: str = ", "
    ) -> List[Dict[str, Any]]:
        """
        Parse semicolon-separated format: time:sig1;sig2, time:sig1;sig2, ...

        Parameters
        ----------
        value : str
            Time series string
        channels : List[str]
            Channel names (e.g., ["DEL-0044", "DEL-0045"])
        delimiter : str
            Delimiter between timepoints

        Returns
        -------
        List[Dict]
            List of {"time": float, "signals": [int, ...]}
        """
        if not value or pd.isna(value):
            return []

        data = []
        # Split by delimiter, handling potential whitespace
        parts = [p.strip() for p in str(value).split(delimiter) if p.strip()]

        for part in parts:
            if ':' not in part:
                continue

            time_str, signals_str = part.split(':', 1)

            try:
                time_val = float(time_str)

                # Parse signals (semicolon-separated)
                signal_parts = signals_str.split(';')
                signals = [int(float(s.strip())) for s in signal_parts if s.strip()]

                # Pad with zeros if fewer signals than channels
                while len(signals) < len(channels):
                    signals.append(0)

                data.append({
                    "time": time_val,
                    "signals": signals[:len(channels)]  # Truncate if more
                })
            except (ValueError, IndexError):
                continue

        # Sort by time
        data.sort(key=lambda x: x["time"])
        return data

    @staticmethod
    def parse_json_format(value: str, channels: List[str]) -> List[Dict[str, Any]]:
        """
        Parse JSON array format: [[time, [sig1, sig2]], ...]
        """
        if not value or pd.isna(value):
            return []

        try:
            parsed = json.loads(value)
            data = []
            for item in parsed:
                if isinstance(item, (list, tuple)) and len(item) >= 2:
                    time_val = float(item[0])
                    signals = item[1] if isinstance(item[1], list) else [item[1]]
                    signals = [int(s) for s in signals]
                    data.append({"time": time_val, "signals": signals})
            return sorted(data, key=lambda x: x["time"])
        except (json.JSONDecodeError, ValueError):
            return []


class DataTransformer:
    """
    Transform external data formats to LC-Seq canonical schema.

    Examples
    --------
    >>> mapping = DataMapping.from_yaml(Path('config/mapping.yaml'))
    >>> transformer = DataTransformer(mapping)
    >>> records = transformer.transform(Path('data/raw_data.csv'))
    >>> len(records)
    1000
    """

    def __init__(self, mapping: DataMapping):
        """
        Initialize transformer with mapping configuration.

        Parameters
        ----------
        mapping : DataMapping
            Configuration for transforming data
        """
        self.mapping = mapping
        self.ts_parser = TimeSeriesParser()

    def transform(self, file_path: Path) -> List[Dict[str, Any]]:
        """
        Transform external data file to canonical format.

        Parameters
        ----------
        file_path : Path
            Path to input data file (CSV or Excel)

        Returns
        -------
        List[Dict]
            List of compound records in canonical format
        """
        # Read data file
        if self.mapping.source_format == 'csv':
            df = pd.read_csv(file_path, encoding=self.mapping.encoding)
        elif self.mapping.source_format in ('xlsx', 'excel'):
            df = pd.read_excel(file_path)
        else:
            raise ValueError(f"Unsupported format: {self.mapping.source_format}")

        records = []
        for _, row in df.iterrows():
            record = self._transform_row(row)
            if record:
                records.append(record)

        return records

    def _transform_row(self, row: pd.Series) -> Optional[Dict[str, Any]]:
        """Transform a single row to canonical format."""
        record = {}

        # Identity
        record['identity'] = self._extract_identity(row)
        if not record['identity'].get('id'):
            return None  # Skip rows without valid ID

        # Building blocks
        record['building_blocks'] = self._extract_building_blocks(row)

        # Properties
        record['properties'] = self._extract_properties(row)

        # Retention times
        record['retention_times'] = self._extract_retention_times(row)

        # Chromatogram
        record['chromatogram'] = self._extract_chromatogram(row)

        return record

    def _extract_identity(self, row: pd.Series) -> Dict[str, Any]:
        """Extract identity fields from row."""
        identity = {}

        id_map = self.mapping.identity
        for canonical_field, source_col in id_map.items():
            if source_col in row.index:
                val = row[source_col]
                identity[canonical_field] = None if pd.isna(val) else str(val)

        return identity

    def _extract_building_blocks(self, row: pd.Series) -> List[Dict[str, Any]]:
        """Extract building block information from row."""
        blocks = []

        for bb_map in self.mapping.building_blocks:
            block = {
                'position': bb_map.position,
                'properties': {}
            }

            # Name (required)
            if bb_map.name in row.index:
                val = row[bb_map.name]
                block['name'] = None if pd.isna(val) else str(val).strip()
            else:
                block['name'] = None

            # Optional fields
            if bb_map.smiles and bb_map.smiles in row.index:
                val = row[bb_map.smiles]
                block['smiles'] = None if pd.isna(val) else str(val)

            if bb_map.code and bb_map.code in row.index:
                val = row[bb_map.code]
                block['code'] = None if pd.isna(val) else str(val)

            if bb_map.stereochemistry and bb_map.stereochemistry in row.index:
                val = row[bb_map.stereochemistry]
                block['stereochemistry'] = None if pd.isna(val) else str(val)

            # Properties
            for prop_name, prop_col in bb_map.properties.items():
                if prop_col in row.index:
                    val = row[prop_col]
                    block['properties'][prop_name] = None if pd.isna(val) else val

            blocks.append(block)

        return blocks

    def _extract_properties(self, row: pd.Series) -> Dict[str, Any]:
        """Extract chemical/physical properties from row."""
        props = {'custom': {}}

        prop_map = self.mapping.properties

        # Handle nested mappings (like exact_mass: {cooh: "col", dma: "col"})
        for prop_name, source in prop_map.items():
            if prop_name == 'custom':
                continue

            if isinstance(source, dict):
                # Nested mapping
                nested = {}
                for sub_name, sub_col in source.items():
                    if sub_col in row.index:
                        val = row[sub_col]
                        nested[sub_name] = None if pd.isna(val) else float(val)
                props[prop_name] = nested
            else:
                # Direct mapping
                if source in row.index:
                    val = row[source]
                    if pd.isna(val):
                        props[prop_name] = None
                    elif prop_name in ('h_donors', 'h_acceptors', 'aromatic_rings',
                                      'peptide_size', 'n_substitution_total'):
                        props[prop_name] = int(val) if not pd.isna(val) else None
                    else:
                        props[prop_name] = float(val) if not pd.isna(val) else None

        return props

    def _extract_retention_times(self, row: pd.Series) -> Dict[str, Any]:
        """Extract retention time data from row."""
        rt = {}

        rt_map = self.mapping.retention_times

        def extract_nested(mapping: Dict, target: Dict):
            for key, value in mapping.items():
                if isinstance(value, dict):
                    target[key] = {}
                    extract_nested(value, target[key])
                else:
                    if value in row.index:
                        val = row[value]
                        target[key] = None if pd.isna(val) else float(val)

        extract_nested(rt_map, rt)
        return rt

    def _extract_chromatogram(self, row: pd.Series) -> Dict[str, Any]:
        """Extract chromatogram/time series data from row."""
        chrom = {
            'total_fractions': None,
            'max_count': None,
            'time_unit': 'seconds',
            'channels': [],
            'data': []
        }

        chrom_config = self.mapping.chromatogram

        # Metadata
        if 'total_fractions' in chrom_config:
            col = chrom_config['total_fractions']
            if col in row.index:
                val = row[col]
                chrom['total_fractions'] = int(val) if not pd.isna(val) else None

        if 'max_count' in chrom_config:
            col = chrom_config['max_count']
            if col in row.index:
                val = row[col]
                chrom['max_count'] = int(val) if not pd.isna(val) else None

        # Time series
        ts_config = chrom_config.get('time_series', {})
        if ts_config:
            col = ts_config.get('column')
            fmt = ts_config.get('format', 'semicolon')
            channels = ts_config.get('channels', [])
            time_unit = ts_config.get('time_unit', 'seconds')
            delimiter = ts_config.get('delimiter', ', ')

            chrom['time_unit'] = time_unit
            chrom['channels'] = channels

            if col and col in row.index:
                value = row[col]

                if 'semicolon' in fmt or ':' in fmt:
                    chrom['data'] = self.ts_parser.parse_semicolon_format(
                        value, channels, delimiter
                    )
                elif fmt == 'json':
                    chrom['data'] = self.ts_parser.parse_json_format(value, channels)

        return chrom

    def transform_to_jsonl(
        self,
        input_path: Path,
        output_path: Path
    ) -> int:
        """
        Transform data file and write to JSONL format.

        Parameters
        ----------
        input_path : Path
            Input data file
        output_path : Path
            Output JSONL file path

        Returns
        -------
        int
            Number of records written
        """
        records = self.transform(input_path)

        with open(output_path, 'w') as f:
            for record in records:
                f.write(json.dumps(record) + '\n')

        return len(records)


def create_default_mapping_for_test_data() -> DataMapping:
    """
    Create a DataMapping configured for the test_data/raw_data.csv format.

    This serves as a reference implementation for the test data format.
    """
    mapping = DataMapping(
        source_format='csv',
        encoding='utf-8'
    )

    # Identity mappings
    mapping.identity = {
        'name': 'Common_Name (N-->C)',
        'id': 'Identifier',
        'library_id': 'lid',
        'smiles': 'enumerated_smiles',
        'product_id': 'intended_product_id',
        'stereochemistry': 'All Stereochem (N-->C)'
    }

    # Building blocks (N->C order: BB1 is first/N-terminal)
    mapping.building_block_order = 'N_to_C'
    mapping.building_blocks = [
        BuildingBlockMapping(
            position=1,
            name='BB1 Name',
            smiles='bb1_smiles',
            stereochemistry='BB1 Stereochem',
            properties={'n_substitution': 'BB1 N-Sub', 'aa_count': 'BB1 AA Count'}
        ),
        BuildingBlockMapping(
            position=2,
            name='BB2 Name',
            smiles='bb2_smiles',
            stereochemistry='BB2 Stereochem',
            properties={'n_substitution': 'BB2 N-Sub', 'aa_count': 'BB2 AA Count'}
        ),
        BuildingBlockMapping(
            position=3,
            name='BB3 Name',
            smiles='bb3_smiles',
            stereochemistry='BB3 Stereochem',
            properties={'n_substitution': 'BB3 N-Sub', 'aa_count': 'BB3 AA Count'}
        )
    ]

    # Properties
    mapping.properties = {
        'exact_mass': {
            'cooh': 'Exact Mass (COOH)',
            'dma': 'Exact Mass (DMA)'
        },
        'logp': 'AlogP',
        'h_donors': 'NumHDonors',
        'h_acceptors': 'NumHAcceptors',
        'aromatic_rings': 'NumAromaticRings',
        'tpsa': 'TPSA',
        'peptide_size': 'Peptide Size',
        'n_substitution_total': 'SUM Nsub'
    }

    # Retention times
    mapping.retention_times = {
        'linear': {
            'rt_minutes': 'Linear RT (min)',
            'logk': 'LogK10-40 (Linear)'
        },
        'cyclized': {
            'rt_minutes': 'Cyclized RT (min)',
            'logk': 'LogK10-40 (Cyclized)'
        },
        'delta_rt': 'Delta RT (min) (Cyclized-Linear)',
        'weighted': {
            'combined': 'wtd_rt',
            'by_library': {
                'DEL-0044': 'wtd_rt DEL-0044',
                'DEL-0045': 'wtd_rt DEL-0045'
            }
        }
    }

    # Chromatogram
    mapping.chromatogram = {
        'total_fractions': 'total_fractions',
        'max_count': 'max_count',
        'time_series': {
            'column': 'all_datapoints',
            'format': 'time:sig1;sig2, time:sig1;sig2, ...',
            'time_unit': 'seconds',
            'channels': ['DEL-0044', 'DEL-0045'],
            'delimiter': ', '
        }
    }

    return mapping


def generate_mapping_template(output_path: Path) -> None:
    """
    Generate a template YAML mapping configuration.

    Parameters
    ----------
    output_path : Path
        Path to write the template
    """
    template = """# LC-Seq Data Mapping Configuration
# Maps external data columns to canonical schema fields

source:
  format: csv  # csv, xlsx, excel
  encoding: utf-8

mappings:
  # Compound identity fields
  identity:
    name: "Common_Name"          # Human-readable name
    id: "Identifier"             # Unique identifier
    library_id: "Library_ID"     # Library this compound belongs to
    smiles: "SMILES"             # Full compound SMILES (optional)
    product_id: "Product_ID"     # Extended product ID (optional)
    stereochemistry: "Stereo"    # Stereochemistry string (optional)

  # Building block definitions
  building_blocks:
    order: "N_to_C"  # or "C_to_N"
    blocks:
      - position: 1
        name: "BB1_Name"
        smiles: "BB1_SMILES"      # optional
        code: "BB1_Code"          # optional
        stereochemistry: "BB1_Stereo"  # optional
        properties:               # optional, any additional BB properties
          n_substitution: "BB1_NSub"
          aa_count: "BB1_AACount"

      - position: 2
        name: "BB2_Name"
        smiles: "BB2_SMILES"

      - position: 3
        name: "BB3_Name"
        smiles: "BB3_SMILES"

  # Chemical/physical properties
  properties:
    exact_mass:
      cooh: "Mass_COOH"
      dma: "Mass_DMA"
    logp: "AlogP"
    h_donors: "HBD"
    h_acceptors: "HBA"
    aromatic_rings: "ArRings"
    tpsa: "TPSA"
    peptide_size: "Size"
    n_substitution_total: "TotalNSub"
    # Add any custom properties:
    # custom_property: "Column_Name"

  # Retention time data
  retention_times:
    linear:
      rt_minutes: "Linear_RT"
      logk: "Linear_LogK"
    cyclized:
      rt_minutes: "Cyclized_RT"
      logk: "Cyclized_LogK"
    delta_rt: "Delta_RT"
    weighted:
      combined: "Weighted_RT"
      by_library:
        LIB1: "Weighted_RT_LIB1"
        LIB2: "Weighted_RT_LIB2"

  # Chromatogram/time series data
  chromatogram:
    total_fractions: "Total_Fractions"
    max_count: "Max_Count"
    time_series:
      column: "All_Datapoints"
      # Supported formats:
      # - "time:sig1;sig2, time:sig1;sig2, ..." (semicolon format)
      # - "json" for [[time, [sig1, sig2]], ...] format
      format: "time:sig1;sig2, time:sig1;sig2, ..."
      time_unit: seconds  # or minutes
      channels: ["LIB1", "LIB2"]
      delimiter: ", "
"""

    with open(output_path, 'w') as f:
        f.write(template)
