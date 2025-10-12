"""
CSV parser for chromatogram data.

Parses CSV files containing chromatogram time-series data.
"""

import numpy as np
import pandas as pd
from pathlib import Path
from typing import Dict, List
from ...domain.entities.chromatogram import Chromatogram


class CSVParser:
    """
    Parse chromatogram data from CSV files.

    Expected CSV format:
    - Column 1: time (retention time values)
    - Remaining columns: signal values for each compound/variant
    - Header row with compound identifiers

    Examples
    --------
    >>> parser = CSVParser()
    >>> chromatograms = parser.parse(Path('data/chromatograms.csv'))
    >>> len(chromatograms)
    100
    """

    def parse(self, file_path: Path) -> Dict[str, Chromatogram]:
        """
        Parse CSV file into chromatograms.

        Parameters
        ----------
        file_path : Path
            Path to CSV file

        Returns
        -------
        Dict[str, Chromatogram]
            Map of compound_id -> Chromatogram

        Raises
        ------
        FileNotFoundError
            If file doesn't exist
        ValueError
            If CSV format is invalid

        Examples
        --------
        >>> chromatograms = parser.parse(Path('data.csv'))
        >>> chromatograms['compound_1'].time_points.shape
        (1000,)
        """
        if not file_path.exists():
            raise FileNotFoundError(f"CSV file not found: {file_path}")

        # Read CSV
        df = pd.read_csv(file_path)

        if df.empty:
            raise ValueError(f"CSV file is empty: {file_path}")

        # First column is time
        time_column = df.columns[0]
        time_points = df[time_column].values

        # Remaining columns are signals
        chromatograms = {}
        for col in df.columns[1:]:
            signal = df[col].values

            chrom = Chromatogram(
                time_points=time_points,
                counts=signal
            )
            chromatograms[col] = chrom

        return chromatograms

    def parse_with_variants(
        self,
        file_path: Path,
        variant_pattern: str = '_'
    ) -> Dict[str, List[Chromatogram]]:
        """
        Parse CSV with variant grouping.

        Parameters
        ----------
        file_path : Path
            Path to CSV file
        variant_pattern : str
            Separator between compound ID and variant ID

        Returns
        -------
        Dict[str, List[Chromatogram]]
            Map of compound_id -> list of variant chromatograms

        Notes
        -----
        Assumes column names like "compound1_v1", "compound1_v2", etc.
        Groups variants by splitting on variant_pattern.

        Examples
        --------
        >>> chroms_grouped = parser.parse_with_variants(Path('data.csv'))
        >>> chroms_grouped['compound1']  # List of variants
        [<Chromatogram>, <Chromatogram>, <Chromatogram>]
        """
        all_chroms = self.parse(file_path)

        # Group by compound ID (before variant pattern)
        grouped = {}
        for compound_id, chrom in all_chroms.items():
            # Extract base compound ID
            if variant_pattern in compound_id:
                base_id = compound_id.split(variant_pattern)[0]
            else:
                base_id = compound_id

            if base_id not in grouped:
                grouped[base_id] = []
            grouped[base_id].append(chrom)

        return grouped
