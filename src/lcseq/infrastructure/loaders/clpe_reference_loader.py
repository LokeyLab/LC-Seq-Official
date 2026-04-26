"""
cLPE Reference Data Loader.

Loads reference data for cLPE (chromatographic Linear Peptide Equation)
validation from CSV files. This includes pre-computed AlogP and scaffold
grouping information.

Note: LogK is computed from OBSERVED retention times, not loaded from
reference data. The reference only provides structural properties.
"""

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import pandas as pd
import numpy as np

# Note: LogK is computed from observed RT, NOT loaded from reference data
# cLPE validation uses: AlogP (lipophilicity) + scaffold (grouping)


@dataclass
class CLPEReferenceData:
    """
    Container for cLPE reference data.

    Only contains pre-computed values needed for cLPE validation:
    - AlogP: calculated lipophilicity from structure
    - Scaffold group: stereochemistry grouping for regression

    LogK is computed from observed RT, NOT loaded from reference data.

    Attributes
    ----------
    alogp_map : Dict[str, float]
        Mapping from compound identifier to AlogP value
    scaffold_map : Dict[str, str]
        Mapping from compound identifier to scaffold group
    """
    alogp_map: Dict[str, float]
    scaffold_map: Dict[str, str]

    def get_compound_data(
        self,
        compound_id: str
    ) -> Tuple[Optional[float], Optional[str]]:
        """
        Get cLPE reference data for a compound.

        Parameters
        ----------
        compound_id : str
            Compound identifier (e.g., block sequence)

        Returns
        -------
        Tuple[Optional[float], Optional[str]]
            (alogp, scaffold_group)
        """
        return (
            self.alogp_map.get(compound_id),
            self.scaffold_map.get(compound_id)
        )


class CLPEReferenceLoader:
    """
    Loads cLPE reference data from CSV files.

    Only loads structural properties needed for cLPE validation:
    - AlogP: Calculated lipophilicity from structure
    - Scaffold group: Stereochemistry grouping for regression

    LogK is computed from OBSERVED retention times, not loaded from reference.

    The CSV file is expected to have columns:
    - Common_Name (N-->C): Compound name/sequence
    - AlogP: Calculated lipophilicity
    - All Stereochem (N-->C): Scaffold/stereochemistry grouping
    - BB1 Name, BB2 Name, BB3 Name: Building block codes

    Examples
    --------
    >>> loader = CLPEReferenceLoader()
    >>> ref_data = loader.load("test_data/raw_data.csv")
    >>> alogp, scaffold = ref_data.get_compound_data("Leu-Ala-Pro")
    """

    def load(self, csv_path: Path) -> CLPEReferenceData:
        """
        Load reference data from CSV file.

        Parameters
        ----------
        csv_path : Path
            Path to CSV file with cLPE reference data

        Returns
        -------
        CLPEReferenceData
            Loaded reference data (AlogP and scaffold only)
        """
        csv_path = Path(csv_path)
        if not csv_path.exists():
            raise FileNotFoundError(f"Reference data file not found: {csv_path}")

        df = pd.read_csv(csv_path)

        alogp_map: Dict[str, float] = {}
        scaffold_map: Dict[str, str] = {}

        for _, row in df.iterrows():
            # Get compound identifier
            compound_name = row.get("Common_Name (N-->C)")
            if pd.isna(compound_name):
                continue

            # Also create identifier from building blocks
            bb_names = []
            for bb_col in ["BB3 Name", "BB2 Name", "BB1 Name"]:  # N→C order
                bb_name = row.get(bb_col)
                if pd.notna(bb_name) and bb_name not in ["-", "null", "NULL", "AgxNull"]:
                    bb_names.append(str(bb_name).strip())

            # Use compound name as primary key, building block sequence as secondary
            keys = [compound_name]
            if bb_names:
                bb_sequence = "-".join(bb_names)
                keys.append(bb_sequence)

            # Extract only structural data (not RT or LogK)
            alogp = row.get("AlogP")
            scaffold = row.get("All Stereochem (N-->C)")

            # Store for each key
            for key in keys:
                if pd.notna(alogp) and not np.isnan(alogp):
                    alogp_map[key] = float(alogp)

                if pd.notna(scaffold):
                    scaffold_map[key] = str(scaffold).strip()

        return CLPEReferenceData(
            alogp_map=alogp_map,
            scaffold_map=scaffold_map
        )

    def load_and_apply(
        self,
        csv_path: Path,
        compounds: List
    ) -> CLPEReferenceData:
        """
        Load reference data and apply to compounds.

        Parameters
        ----------
        csv_path : Path
            Path to CSV file with cLPE reference data
        compounds : List[Compound]
            Compounds to apply reference data to

        Returns
        -------
        CLPEReferenceData
            Loaded reference data
        """
        ref_data = self.load(csv_path)

        for compound in compounds:
            # Try to match by compound_id first, then by block_support_sequence
            keys_to_try = []
            if compound.compound_id:
                keys_to_try.append(compound.compound_id)
            keys_to_try.append(compound.block_support_sequence)
            keys_to_try.append(compound.positional_block_sequence)

            for key in keys_to_try:
                alogp, scaffold = ref_data.get_compound_data(key)
                if alogp is not None or scaffold is not None:
                    compound.alogp = alogp
                    compound.scaffold_group = scaffold
                    break

        return ref_data
