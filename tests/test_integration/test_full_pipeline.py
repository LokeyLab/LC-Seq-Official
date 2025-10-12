"""
End-to-end integration test for full analysis pipeline.

Tests the complete workflow from raw data to validated results.
"""

import pytest
import numpy as np
from pathlib import Path

from lcseq.domain.entities.building_block import BuildingBlock
from lcseq.domain.entities.compound import Compound
from lcseq.domain.entities.chromatogram import Chromatogram
from lcseq.application.dtos.analysis_request import AnalysisRequest
from lcseq.application.pipelines.full_analysis_pipeline import FullAnalysisPipeline


class TestFullAnalysisPipeline:
    """Test end-to-end analysis pipeline."""

    @pytest.fixture
    def mock_library(self):
        """Create mock library with synthetic compounds."""
        # Create building blocks
        leu = BuildingBlock(cycle=0, code='Leu', is_null=False)
        ala = BuildingBlock(cycle=1, code='Ala', is_null=False)
        val = BuildingBlock(cycle=2, code='Val', is_null=False)
        null = BuildingBlock.from_code(cycle=0, code='Null')

        # Create compounds at different levels
        compounds = []

        # Create null blocks for each position
        null0 = BuildingBlock.from_code(cycle=0, code='Null')
        null1 = BuildingBlock.from_code(cycle=1, code='Null')
        null2 = BuildingBlock.from_code(cycle=2, code='Null')

        # L0 (full null) - level 0
        l0 = Compound(
            building_blocks=[null0, null1, null2],
            chromatogram=None  # Will be added later
        )
        compounds.append(l0)

        # L1 compounds - level 1
        l1_a = Compound(
            building_blocks=[leu, null1, null2],
            chromatogram=None
        )
        compounds.append(l1_a)

        l1_b = Compound(
            building_blocks=[null0, ala, null2],
            chromatogram=None
        )
        compounds.append(l1_b)

        # L2 compound - level 2
        l2 = Compound(
            building_blocks=[leu, ala, null2],
            chromatogram=None
        )
        compounds.append(l2)

        # L3 compound - level 3 (maximal)
        l3 = Compound(
            building_blocks=[leu, ala, val],
            chromatogram=None
        )
        compounds.append(l3)

        return compounds

    @pytest.fixture
    def mock_chromatograms(self, mock_library):
        """Create synthetic chromatograms for test compounds."""
        # Create time axis
        time_points = np.linspace(0, 60, 1000)  # 0-60 minutes

        chromatograms = []

        for idx, compound in enumerate(mock_library):
            # Create synthetic signal with peaks
            # Base level: level 0 elutes earliest, level 3 elutes latest
            retention_time = 10 + compound.level * 5  # Spacing peaks

            # Generate Gaussian peak
            signal = np.zeros_like(time_points)
            peak_height = 100.0 / (compound.level + 1)  # Higher level = smaller peak
            peak_width = 1.0

            signal = peak_height * np.exp(-0.5 * ((time_points - retention_time) / peak_width) ** 2)

            # Add some baseline and noise
            signal += 5.0  # Baseline
            signal += np.random.normal(0, 0.5, size=signal.shape)  # Noise

            # Create chromatogram
            chrom = Chromatogram(
                time_points=time_points,
                counts=signal
            )
            chromatograms.append(chrom)

        return chromatograms

    def test_full_pipeline_execution(self, mock_library, mock_chromatograms):
        """
        Test complete pipeline execution.

        Verifies:
        1. Pipeline executes without errors
        2. Results are generated for all compounds
        3. Validation is performed
        4. Statistics are computed
        """
        # Attach chromatograms to compounds
        for compound, chrom in zip(mock_library, mock_chromatograms):
            compound.chromatogram = chrom

        # Create analysis request
        request = AnalysisRequest(
            data_path=Path('/dev/null'),  # Not used for this test
            output_path=Path('/tmp/test_results'),
            variant_mode='individual',
            hierarchy_mode='building_block',
            detection_params={'min_persistence': 0.05},
            validation_params={'retention_precision': 0.5},
            export_formats=['csv']
        )

        # Execute pipeline
        pipeline = FullAnalysisPipeline()
        response = pipeline.execute(mock_library, request)

        # Verify response structure
        assert response is not None
        assert response.request_id is not None
        assert len(response.compound_results) == len(mock_library)
        assert response.validation_summary is not None
        assert response.dataset_stats is not None

        # Verify some compounds were analyzed
        assert response.validation_summary.total_compounds == len(mock_library)

        # Verify no critical errors
        assert len(response.errors) == 0

        print(f"\\n✓ Pipeline executed successfully")
        print(f"  - Analyzed {len(mock_library)} compounds")
        print(f"  - Validation rate: {response.validation_summary.validation_rate:.1%}")
        print(f"  - Median purity: {response.validation_summary.median_purity:.3f}")

    def test_pipeline_with_empty_library(self):
        """Test pipeline handles empty library gracefully."""
        request = AnalysisRequest(
            data_path=Path('/dev/null'),
            output_path=Path('/tmp/test_results')
        )

        pipeline = FullAnalysisPipeline()
        response = pipeline.execute([], request)

        # Should return valid (empty) response, not crash
        assert response is not None
        assert len(response.compound_results) == 0
        assert response.validation_summary.total_compounds == 0

    def test_validation_services_integration(self, mock_library, mock_chromatograms):
        """
        Test that validation services are properly integrated.

        Verifies:
        1. Adaptive thresholds are computed
        2. Bayesian validation is applied
        3. Results contain purity, SNR, validation status
        """
        # Attach chromatograms
        for compound, chrom in zip(mock_library, mock_chromatograms):
            compound.chromatogram = chrom

        request = AnalysisRequest(
            data_path=Path('/dev/null'),
            output_path=Path('/tmp/test_results')
        )

        pipeline = FullAnalysisPipeline()
        response = pipeline.execute(mock_library, request)

        # Verify dataset stats include required keys
        required_stats = ['purity_p25', 'purity_p50', 'purity_p75', 'background']
        for key in required_stats:
            assert key in response.dataset_stats, f"Missing stat: {key}"

        # Verify compound results have validation info
        for result in response.compound_results:
            assert hasattr(result, 'validation_status')
            assert hasattr(result, 'purity')
            assert hasattr(result, 'snr')
            assert hasattr(result, 'purity_category')

        print(f"\\n✓ Validation integration verified")
        print(f"  - P25: {response.dataset_stats['purity_p25']:.3f}")
        print(f"  - P50: {response.dataset_stats['purity_p50']:.3f}")
        print(f"  - P75: {response.dataset_stats['purity_p75']:.3f}")


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
