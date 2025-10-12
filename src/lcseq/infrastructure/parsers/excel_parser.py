"""
Excel parser for chromatogram data.

Parses Excel files containing chromatogram time-series data.
"""

import pandas as pd
from pathlib import Path
from typing import Dict, List, Optional
from ...domain.entities.chromatogram import Chromatogram


class ExcelParser:
    """
    Parse chromatogram data from Excel files.

    Expected Excel format:
    - First column: time (retention time values)
    - Remaining columns: signal values for each compound/variant
    - Header row with compound identifiers
    - Can handle multiple sheets

    Examples
    --------
    >>> parser = ExcelParser()
    >>> chromatograms = parser.parse(Path('data/chromatograms.xlsx'))
    >>> len(chromatograms)
    100
    """

    def parse(
        self,
        file_path: Path,
        sheet_name: Optional[str] = None
    ) -> Dict[str, Chromatogram]:
        """
        Parse Excel file into chromatograms.

        Parameters
        ----------
        file_path : Path
            Path to Excel file (.xlsx or .xls)
        sheet_name : Optional[str]
            Sheet name to parse (if None, uses first sheet)

        Returns
        -------
        Dict[str, Chromatogram]
            Map of compound_id -> Chromatogram

        Raises
        ------
        FileNotFoundError
            If file doesn't exist
        ValueError
            If Excel format is invalid

        Examples
        --------
        >>> chromatograms = parser.parse(Path('data.xlsx'), sheet_name='Sheet1')
        >>> chromatograms['compound_1'].counts.shape
        (1000,)
        """
        if not file_path.exists():
            raise FileNotFoundError(f"Excel file not found: {file_path}")

        # Read Excel
        if sheet_name:
            df = pd.read_excel(file_path, sheet_name=sheet_name)
        else:
            df = pd.read_excel(file_path)

        if df.empty:
            raise ValueError(f"Excel file/sheet is empty: {file_path}")

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

    def parse_all_sheets(
        self,
        file_path: Path
    ) -> Dict[str, Dict[str, Chromatogram]]:
        """
        Parse all sheets in Excel file.

        Parameters
        ----------
        file_path : Path
            Path to Excel file

        Returns
        -------
        Dict[str, Dict[str, Chromatogram]]
            Map of sheet_name -> (compound_id -> Chromatogram)

        Examples
        --------
        >>> all_data = parser.parse_all_sheets(Path('data.xlsx'))
        >>> all_data['Sheet1']['compound_1']
        <Chromatogram>
        """
        if not file_path.exists():
            raise FileNotFoundError(f"Excel file not found: {file_path}")

        # Get all sheet names
        excel_file = pd.ExcelFile(file_path)
        sheet_names = excel_file.sheet_names

        # Parse each sheet
        all_sheets = {}
        for sheet_name in sheet_names:
            chromatograms = self.parse(file_path, sheet_name=sheet_name)
            all_sheets[sheet_name] = chromatograms

        return all_sheets
