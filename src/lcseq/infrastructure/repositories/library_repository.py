"""
Repository for library design data.

Implements persistence operations for library design and compound collections.
"""

from pathlib import Path
from typing import List
from ...domain.entities.compound import Compound
from ..parsers.library_parser import LibraryParser


class LibraryRepository:
    """
    Repository for library design persistence.

    Handles loading library designs from various file formats.

    Examples
    --------
    >>> repo = LibraryRepository()
    >>> compounds = repo.load(Path('library_design.csv'))
    >>> len(compounds)
    1000
    """

    def __init__(self):
        """Initialize repository with parser."""
        self.parser = LibraryParser()

    def load(self, file_path: Path) -> List[Compound]:
        """
        Load library design from file.

        Parameters
        ----------
        file_path : Path
            Path to library design file (CSV or Excel)

        Returns
        -------
        List[Compound]
            Parsed compounds from library

        Examples
        --------
        >>> compounds = repo.load(Path('library.csv'))
        """
        return self.parser.parse(file_path)

    def load_with_metadata(self, file_path: Path) -> list:
        """
        Load library with metadata preserved.

        Parameters
        ----------
        file_path : Path
            Path to library design file

        Returns
        -------
        list
            List of dicts with compounds and metadata

        Examples
        --------
        >>> data = repo.load_with_metadata(Path('library.csv'))
        >>> data[0]['compound']
        <Compound>
        >>> data[0]['plate']
        'Plate_1'
        """
        return self.parser.parse_with_metadata(file_path)
