"""HDF5 compound loader for LC-Seq data.

This module handles loading compound data from HDF5 files.
"""

from pathlib import Path
from typing import List
import h5py

from lcseq.domain.entities import BuildingBlock, Chromatogram, Compound


class HDF5CompoundLoader:
    """Loads compounds from HDF5 files.

    This loader reads chromatographic data stored in HDF5 format and constructs
    domain entities (Compounds with BuildingBlocks and Chromatograms).
    """

    def load_all(self, hdf5_path: Path) -> List[Compound]:
        """Load all compounds from HDF5 file.

        Args:
            hdf5_path: Path to HDF5 file containing compound data

        Returns:
            List of Compound objects with chromatographic data

        The HDF5 file is expected to have the following structure:
        - metadata/bb1_name: Array of building block names at position 1
        - metadata/bb2_name: Array of building block names at position 2
        - metadata/bb3_name: Array of building block names at position 3
        - chromatograms/time_points: Concatenated time point arrays
        - chromatograms/counts: Concatenated count arrays
        - chromatograms/lengths: Length of each chromatogram
        """
        print(f"Loading data from {hdf5_path}...")
        compounds = []

        with h5py.File(hdf5_path, "r") as f:
            # Read metadata
            bb1_names = f["metadata"]["bb1_name"][:]
            bb2_names = f["metadata"]["bb2_name"][:]
            bb3_names = f["metadata"]["bb3_name"][:]

            # Read chromatogram data
            all_time_points = f["chromatograms"]["time_points"][:]
            all_counts = f["chromatograms"]["counts"][:]
            lengths = f["chromatograms"]["lengths"][:]

            # Decode byte strings to UTF-8
            bb1_names = [
                n.decode("utf-8") if isinstance(n, bytes) else str(n)
                for n in bb1_names
            ]
            bb2_names = [
                n.decode("utf-8") if isinstance(n, bytes) else str(n)
                for n in bb2_names
            ]
            bb3_names = [
                n.decode("utf-8") if isinstance(n, bytes) else str(n)
                for n in bb3_names
            ]

            # Parse compounds
            total = len(bb1_names)
            offset = 0

            for i in range(total):
                # Extract chromatogram data for this compound
                length = lengths[i]
                time_points = all_time_points[offset : offset + length]
                counts = all_counts[offset : offset + length]
                offset += length

                # Create chromatogram entity
                chromatogram = Chromatogram(time_points=time_points, counts=counts)

                # Create building blocks (note: positions are 0-indexed internally)
                # bb3 is position 0, bb2 is position 1, bb1 is position 2
                building_blocks = [
                    BuildingBlock.from_code(0, bb3_names[i]),
                    BuildingBlock.from_code(1, bb2_names[i]),
                    BuildingBlock.from_code(2, bb1_names[i]),
                ]

                # Create compound entity
                compound = Compound(building_blocks, chromatogram)
                compounds.append(compound)

                # Progress tracking
                if (i + 1) % 5000 == 0:
                    print(f"  Loaded {i+1:,} compounds...")

        print(f"✓ Loaded {len(compounds):,} compounds")
        return compounds
