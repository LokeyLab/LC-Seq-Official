"""
Library design parser.

Parses library design files containing compound sequences and building blocks.
"""

import pandas as pd
from pathlib import Path
from typing import List, Dict, Any
from ...domain.entities.building_block import BuildingBlock
from ...domain.entities.compound import Compound


class LibraryParser:
    """
    Parse library design files.

    Expected format (CSV/Excel):
    - Column 'compound_id': Unique compound identifier
    - Columns 'pos_0', 'pos_1', ...: Building blocks at each position
    - Optional columns for metadata

    Examples
    --------
    >>> parser = LibraryParser()
    >>> compounds = parser.parse(Path('library_design.csv'))
    >>> len(compounds)
    1000
    """

    def parse(self, file_path: Path) -> List[Compound]:
        """
        Parse library design file into compounds.

        Parameters
        ----------
        file_path : Path
            Path to library design file (CSV or Excel)

        Returns
        -------
        List[Compound]
            List of compounds with sequences

        Raises
        ------
        FileNotFoundError
            If file doesn't exist
        ValueError
            If file format is invalid

        Examples
        --------
        >>> compounds = parser.parse(Path('library.csv'))
        >>> compounds[0].positional_block_sequence
        'Leu-Ala-...'
        """
        if not file_path.exists():
            raise FileNotFoundError(f"Library file not found: {file_path}")

        # Determine file type and read
        if file_path.suffix == '.csv':
            df = pd.read_csv(file_path)
        elif file_path.suffix in ['.xlsx', '.xls']:
            df = pd.read_excel(file_path)
        else:
            raise ValueError(f"Unsupported file type: {file_path.suffix}")

        if df.empty:
            raise ValueError(f"Library file is empty: {file_path}")

        # Parse compounds
        compounds = []
        for idx, row in df.iterrows():
            compound = self._parse_compound_row(row)
            if compound:
                compounds.append(compound)

        return compounds

    def _parse_compound_row(self, row: pd.Series) -> Compound:
        """
        Parse single compound from DataFrame row.

        Parameters
        ----------
        row : pd.Series
            Row from library design DataFrame

        Returns
        -------
        Compound
            Parsed compound (or None if invalid)
        """
        # Extract position columns
        pos_columns = [col for col in row.index if col.startswith('pos_')]

        if not pos_columns:
            return None

        # Build building blocks list (sorted by position/cycle)
        building_blocks_dict = {}
        for col in sorted(pos_columns):
            pos_num = int(col.split('_')[1])
            bb_code = row[col]

            # Handle NULL/empty positions
            if pd.isna(bb_code) or str(bb_code).upper() == 'NULL':
                building_block = BuildingBlock.create_null(cycle=pos_num)
            else:
                # Create building block
                building_block = BuildingBlock(
                    cycle=pos_num,
                    code=str(bb_code),
                    is_null=False
                )

            building_blocks_dict[pos_num] = building_block

        # Convert to sorted list (by cycle)
        building_blocks = [building_blocks_dict[pos] for pos in sorted(building_blocks_dict.keys())]

        # Create compound (Note: chromatogram would need to be added separately)
        # This parser only handles the building block sequence from library design
        # Chromatogram data comes from separate LC-MS files
        from ...domain.entities.chromatogram import Chromatogram
        import numpy as np

        # Placeholder chromatogram (would be replaced with actual data)
        placeholder_chrom = Chromatogram(
            time_points=np.array([0.0, 1.0]),
            counts=np.array([0.0, 0.0])
        )
        compound = Compound(building_blocks=building_blocks, chromatogram=placeholder_chrom)

        return compound

    def parse_with_metadata(
        self,
        file_path: Path
    ) -> List[Dict[str, Any]]:
        """
        Parse library with all metadata preserved.

        Parameters
        ----------
        file_path : Path
            Path to library design file

        Returns
        -------
        List[Dict[str, Any]]
            List of dicts with compound and metadata

        Examples
        --------
        >>> data = parser.parse_with_metadata(Path('library.csv'))
        >>> data[0]
        {
            'compound': <Compound>,
            'compound_id': 'cpd_001',
            'plate': 'Plate_1',
            'well': 'A01',
            ...
        }
        """
        if not file_path.exists():
            raise FileNotFoundError(f"Library file not found: {file_path}")

        # Read file
        if file_path.suffix == '.csv':
            df = pd.read_csv(file_path)
        elif file_path.suffix in ['.xlsx', '.xls']:
            df = pd.read_excel(file_path)
        else:
            raise ValueError(f"Unsupported file type: {file_path.suffix}")

        # Parse with metadata
        results = []
        for idx, row in df.iterrows():
            compound = self._parse_compound_row(row)
            if compound:
                # Convert row to dict
                metadata = row.to_dict()
                metadata['compound'] = compound
                results.append(metadata)

        return results
