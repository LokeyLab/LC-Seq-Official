"""HDF5 compound loader for LC-Seq data.

This module handles loading and saving compound data from/to HDF5 files.
"""

import re
from pathlib import Path
from typing import List
import h5py
import numpy as np

from lcseq.domain.entities import BuildingBlock, Chromatogram, Compound


class HDF5CompoundLoader:
    """Loads compounds from HDF5 files.

    This loader reads chromatographic data stored in HDF5 format and constructs
    domain entities (Compounds with BuildingBlocks and Chromatograms).

    Supports libraries of any length (3-mer, 4-mer, 5-mer, etc.) by auto-detecting
    the number of positions from the metadata keys.
    """

    def load_all(self, hdf5_path: Path) -> List[Compound]:
        """Load all compounds from HDF5 file.

        Args:
            hdf5_path: Path to HDF5 file containing compound data

        Returns:
            List of Compound objects with chromatographic data

        The HDF5 file is expected to have the following structure:
        - metadata/bb1_name: Array of building block names at position 1 (N-terminus)
        - metadata/bb2_name: Array of building block names at position 2
        - metadata/bb3_name: Array of building block names at position 3
        - ... (additional positions for 4-mer, 5-mer, etc.)
        - metadata/bbN_name: Array of building block names at position N (C-terminus)
        - chromatograms/time_points: Concatenated time point arrays
        - chromatograms/counts: Concatenated count arrays
        - chromatograms/lengths: Length of each chromatogram

        Notes
        -----
        Building block positions in HDF5 are numbered 1 to N (N-terminus to C-terminus).
        Internally, compounds store blocks with cycle 0 at C-terminus.
        So bb[N] maps to cycle 0, bb[N-1] to cycle 1, ..., bb[1] to cycle N-1.
        """
        print(f"Loading data from {hdf5_path}...")
        compounds = []

        with h5py.File(hdf5_path, "r") as f:
            # Auto-detect number of positions from metadata keys
            # Keys follow pattern: bb1_name, bb2_name, bb3_name, ...
            bb_keys = sorted(
                [k for k in f["metadata"].keys() if re.match(r"bb\d+_name", k)],
                key=lambda k: int(re.search(r"\d+", k).group())
            )
            n_positions = len(bb_keys)

            if n_positions == 0:
                raise ValueError("No building block metadata found (expected bb1_name, bb2_name, ...)")

            print(f"  Detected {n_positions}-mer library ({n_positions} positions)")

            # Load all position metadata dynamically
            # bb_names_by_hdf5_position[0] = bb1_names, [1] = bb2_names, etc.
            bb_names_by_hdf5_position = []
            for key in bb_keys:
                names = f["metadata"][key][:]
                names = [
                    n.decode("utf-8") if isinstance(n, bytes) else str(n)
                    for n in names
                ]
                bb_names_by_hdf5_position.append(names)

            # Read chromatogram data
            all_time_points = f["chromatograms"]["time_points"][:]
            all_counts = f["chromatograms"]["counts"][:]
            lengths = f["chromatograms"]["lengths"][:]

            # Parse compounds
            total = len(bb_names_by_hdf5_position[0])
            offset = 0

            for i in range(total):
                # Extract chromatogram data for this compound
                length = lengths[i]
                time_points = all_time_points[offset : offset + length]
                counts = all_counts[offset : offset + length]
                offset += length

                # Create chromatogram entity
                chromatogram = Chromatogram(time_points=time_points, counts=counts)

                # Create building blocks dynamically
                # HDF5 positions: bb1 (N-terminus) ... bbN (C-terminus)
                # Internal cycles: cycle 0 (C-terminus) ... cycle N-1 (N-terminus)
                # So bbN maps to cycle 0, bbN-1 to cycle 1, ..., bb1 to cycle N-1
                building_blocks = []
                for cycle in range(n_positions):
                    # cycle 0 → bbN (index n_positions-1), cycle 1 → bbN-1, etc.
                    hdf5_position_index = n_positions - 1 - cycle
                    bb_name = bb_names_by_hdf5_position[hdf5_position_index][i]
                    building_blocks.append(BuildingBlock.from_code(cycle, bb_name))

                # Create compound entity
                compound = Compound(building_blocks, chromatogram)
                compounds.append(compound)

                # Progress tracking
                if (i + 1) % 5000 == 0:
                    print(f"  Loaded {i+1:,} compounds...")

        print(f"✓ Loaded {len(compounds):,} compounds ({n_positions}-mer)")
        return compounds

    def save_all(self, compounds: List[Compound], hdf5_path: Path) -> None:
        """Save compounds to HDF5 file.

        Args:
            compounds: List of Compound objects to save
            hdf5_path: Output path for HDF5 file

        Notes
        -----
        The HDF5 file will have the same structure as input files:
        - metadata/bb1_name: Array of building block names at position 1 (N-terminus)
        - metadata/bb2_name: Array of building block names at position 2
        - ... (additional positions)
        - metadata/bbN_name: Array of building block names at position N (C-terminus)
        - chromatograms/time_points: Concatenated time point arrays
        - chromatograms/counts: Concatenated count arrays
        - chromatograms/lengths: Length of each chromatogram
        """
        if not compounds:
            raise ValueError("No compounds to save")

        # Determine number of positions from first compound
        n_positions = len(compounds[0].building_blocks)

        # Build metadata arrays - convert cycle back to HDF5 position
        # Internal cycles: cycle 0 (C-terminus) ... cycle N-1 (N-terminus)
        # HDF5 positions: bb1 (N-terminus) ... bbN (C-terminus)
        # So cycle 0 → bbN, cycle 1 → bbN-1, ..., cycle N-1 → bb1
        bb_names_by_hdf5_position = [[] for _ in range(n_positions)]

        # Collect chromatogram data
        all_time_points = []
        all_counts = []
        lengths = []

        for compound in compounds:
            # Extract building block names
            for cycle in range(n_positions):
                hdf5_position_index = n_positions - 1 - cycle
                bb = compound.building_blocks[cycle]
                bb_names_by_hdf5_position[hdf5_position_index].append(bb.code)

            # Extract chromatogram data
            if compound.chromatogram is not None:
                time_points = compound.chromatogram.time_points
                counts = compound.chromatogram.counts
                all_time_points.extend(time_points)
                all_counts.extend(counts)
                lengths.append(len(time_points))
            else:
                lengths.append(0)

        # Write to HDF5
        print(f"Saving {len(compounds):,} compounds to {hdf5_path}...")

        with h5py.File(hdf5_path, "w") as f:
            # Create groups
            metadata_grp = f.create_group("metadata")
            chrom_grp = f.create_group("chromatograms")

            # Write building block metadata
            for i, bb_names in enumerate(bb_names_by_hdf5_position):
                key = f"bb{i + 1}_name"
                # Encode strings as bytes for HDF5
                encoded_names = [n.encode("utf-8") for n in bb_names]
                metadata_grp.create_dataset(key, data=encoded_names)

            # Write chromatogram data
            chrom_grp.create_dataset(
                "time_points",
                data=np.array(all_time_points, dtype=np.float64)
            )
            chrom_grp.create_dataset(
                "counts",
                data=np.array(all_counts, dtype=np.float64)
            )
            chrom_grp.create_dataset(
                "lengths",
                data=np.array(lengths, dtype=np.int64)
            )

        print(f"✓ Saved {len(compounds):,} compounds ({n_positions}-mer)")
