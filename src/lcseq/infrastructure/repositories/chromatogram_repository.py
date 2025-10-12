"""
Repository for loading and saving chromatogram data.

Implements persistence operations for chromatogram entities.
"""

from pathlib import Path
from typing import Dict, List, Optional
from ...domain.entities.chromatogram import Chromatogram
from ...domain.entities.compound import Compound
from ..parsers.csv_parser import CSVParser
from ..parsers.excel_parser import ExcelParser


class ChromatogramRepository:
    """
    Repository for chromatogram data persistence.

    Handles loading chromatograms from various file formats and
    associating them with compounds.

    Examples
    --------
    >>> repo = ChromatogramRepository()
    >>> chromatograms = repo.load_from_csv(Path('data.csv'))
    >>> len(chromatograms)
    100
    """

    def __init__(self):
        """Initialize repository with parsers."""
        self.csv_parser = CSVParser()
        self.excel_parser = ExcelParser()

    def load_from_csv(
        self,
        file_path: Path
    ) -> Dict[str, Chromatogram]:
        """
        Load chromatograms from CSV file.

        Parameters
        ----------
        file_path : Path
            Path to CSV file

        Returns
        -------
        Dict[str, Chromatogram]
            Map of compound_id -> Chromatogram

        Examples
        --------
        >>> chromatograms = repo.load_from_csv(Path('data.csv'))
        """
        return self.csv_parser.parse(file_path)

    def load_from_excel(
        self,
        file_path: Path,
        sheet_name: Optional[str] = None
    ) -> Dict[str, Chromatogram]:
        """
        Load chromatograms from Excel file.

        Parameters
        ----------
        file_path : Path
            Path to Excel file
        sheet_name : Optional[str]
            Sheet to load (if None, uses first sheet)

        Returns
        -------
        Dict[str, Chromatogram]
            Map of compound_id -> Chromatogram

        Examples
        --------
        >>> chromatograms = repo.load_from_excel(Path('data.xlsx'))
        """
        return self.excel_parser.parse(file_path, sheet_name)

    def load_and_associate(
        self,
        file_path: Path,
        compounds: List[Compound]
    ) -> List[Compound]:
        """
        Load chromatograms and associate with compounds.

        Parameters
        ----------
        file_path : Path
            Path to chromatogram file
        compounds : List[Compound]
            Compounds to associate with chromatograms

        Returns
        -------
        List[Compound]
            Compounds with associated chromatograms

        Notes
        -----
        Matches compounds to chromatograms by compound ID/sequence.

        Examples
        --------
        >>> compounds_with_data = repo.load_and_associate(
        ...     Path('chroms.csv'), compounds
        ... )
        """
        # Determine file type and load
        if file_path.suffix == '.csv':
            chromatograms = self.load_from_csv(file_path)
        elif file_path.suffix in ['.xlsx', '.xls']:
            chromatograms = self.load_from_excel(file_path)
        else:
            raise ValueError(f"Unsupported file type: {file_path.suffix}")

        # Associate with compounds
        for compound in compounds:
            compound_id = str(compound)
            if compound_id in chromatograms:
                compound.chromatogram = chromatograms[compound_id]

        return compounds
