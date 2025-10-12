"""
Tests for BayesianValidator service.

Tests Bayesian synthesis validation framework.
"""

import pytest
import numpy as np

from lcseq.domain.entities.building_block import BuildingBlock
from lcseq.domain.entities.chromatogram import Chromatogram
from lcseq.domain.entities.compound import Compound
from lcseq.domain.entities.peak import Peak, PeakType, ValidationStatus
from lcseq.domain.models.compound_hierarchy import CompoundHierarchy, HierarchyMode
from lcseq.domain.services.bayesian_validator import BayesianValidator
from lcseq.domain.services.purity_calculator import PurityCalculator
from lcseq.domain.services.snr_calculator import SNRCalculator


@pytest.fixture
def dataset_stats():
    """Create dataset statistics for testing."""
    return {
        'purity_p25': 0.5,
        'purity_p50': 0.7,
        'purity_p75': 0.85,
        'background': 10.0
    }


@pytest.fixture
def simple_hierarchy():
    """Create simple hierarchy for testing."""
    hierarchy = CompoundHierarchy(mode=HierarchyMode.BUILDING_BLOCK)

    time_points = np.linspace(0, 100, 200)
    counts = np.ones(200)
    chrom = Chromatogram(time_points=time_points, counts=counts)

    # L0
    l0 = Compound([BuildingBlock.from_code(0, "Null")], chrom)
    l0_peak = Peak(position=10.0, left_base=5, right_base=15, height=50, area=100, peak_type=PeakType.NULL)
    l0.detected_peaks = [l0_peak]
    l0.selected_peak = l0_peak
    l0.selected_peak.validation_status = ValidationStatus.VALIDATED

    # L1
    l1 = Compound([BuildingBlock.from_code(0, "Pro")], chrom)
    l1_peak = Peak(position=30.0, left_base=25, right_base=35, height=100, area=200, peak_type=PeakType.PUTATIVE_PRODUCT)
    l1.detected_peaks = [l1_peak]
    l1.selected_peak = l1_peak
    l1.selected_peak.validation_status = ValidationStatus.VALIDATED

    hierarchy.add_compound(l0)
    hierarchy.add_compound(l1)
    hierarchy.add_edge(l1, l0)

    return hierarchy


class TestBayesianValidator:
    """Test BayesianValidator service."""

    def test_validate_high_quality(self, simple_hierarchy, dataset_stats):
        """Test validation of high-quality compound."""
        validator = BayesianValidator()

        # Get L1 compound
        l1 = [c for c in simple_hierarchy.compounds if c.level == 1][0]

        # Add truncation peak to calculate purity
        trunc_peak = Peak(position=10.0, left_base=5, right_base=15, height=20, area=20, peak_type=PeakType.TRUNCATION)
        l1.detected_peaks.insert(0, trunc_peak)

        # High purity: product area / total area = 200 / 220 = 0.91 > p75
        status = validator.validate(l1, simple_hierarchy, dataset_stats)

        assert status == ValidationStatus.VALIDATED

    def test_validate_moderate_quality(self, simple_hierarchy, dataset_stats):
        """Test validation of moderate-quality compound."""
        validator = BayesianValidator()

        l1 = [c for c in simple_hierarchy.compounds if c.level == 1][0]

        # Add larger truncation peak for moderate purity
        trunc_peak = Peak(position=10.0, left_base=5, right_base=15, height=50, area=100, peak_type=PeakType.TRUNCATION)
        l1.detected_peaks.insert(0, trunc_peak)

        # Moderate purity: 200 / 300 = 0.67 (between p50 and p75)
        status = validator.validate(l1, simple_hierarchy, dataset_stats)

        assert status in [ValidationStatus.LIKELY_SUCCESS, ValidationStatus.UNCERTAIN]

    def test_validate_low_quality(self, dataset_stats):
        """Test validation of low-quality compound."""
        validator = BayesianValidator()

        hierarchy = CompoundHierarchy(mode=HierarchyMode.BUILDING_BLOCK)
        time_points = np.linspace(0, 100, 200)
        counts = np.ones(200)
        chrom = Chromatogram(time_points=time_points, counts=counts)

        compound = Compound([BuildingBlock.from_code(0, "Pro")], chrom)

        # Low purity: high truncation, low product
        trunc_peak = Peak(position=10.0, left_base=5, right_base=15, height=100, area=200, peak_type=PeakType.TRUNCATION)
        product_peak = Peak(position=30.0, left_base=25, right_base=35, height=30, area=50, peak_type=PeakType.PUTATIVE_PRODUCT)

        compound.detected_peaks = [trunc_peak, product_peak]
        compound.selected_peak = product_peak

        hierarchy.add_compound(compound)

        # Low purity: 50 / 250 = 0.20 < p25
        status = validator.validate(compound, hierarchy, dataset_stats)

        assert status in [ValidationStatus.LIKELY_FAILURE, ValidationStatus.FAILED]

    def test_validate_retention_violation(self, simple_hierarchy, dataset_stats):
        """Test validation with retention order violation."""
        validator = BayesianValidator()

        # Create descendant at higher retention time (violation)
        time_points = np.linspace(0, 100, 200)
        counts = np.ones(200)
        chrom = Chromatogram(time_points=time_points, counts=counts)

        l2 = Compound([BuildingBlock.from_code(0, "Pro"), BuildingBlock.from_code(1, "Leu")], chrom)
        l2_peak = Peak(position=25.0, left_base=20, right_base=30, height=100, area=200, peak_type=PeakType.PUTATIVE_PRODUCT)
        l2.detected_peaks = [l2_peak]
        l2.selected_peak = l2_peak

        # L1 is at 30.0, L2 at 25.0 → VIOLATION (ancestor should elute later)
        l1 = [c for c in simple_hierarchy.compounds if c.level == 1][0]

        simple_hierarchy.add_compound(l2)
        simple_hierarchy.add_edge(l2, l1)

        status = validator.validate(l2, simple_hierarchy, dataset_stats, retention_precision=1.0)

        # Should fail due to retention violation
        assert status == ValidationStatus.FAILED

    def test_validate_not_validated_no_peak(self, simple_hierarchy, dataset_stats):
        """Test that compound with no selected peak returns NOT_VALIDATED."""
        validator = BayesianValidator()

        time_points = np.linspace(0, 100, 200)
        counts = np.ones(200)
        chrom = Chromatogram(time_points=time_points, counts=counts)

        compound = Compound([BuildingBlock.from_code(0, "Val")], chrom)
        compound.detected_peaks = []
        compound.selected_peak = None

        simple_hierarchy.add_compound(compound)

        status = validator.validate(compound, simple_hierarchy, dataset_stats)

        assert status == ValidationStatus.NOT_VALIDATED


class TestBayesianValidatorHelperMethods:
    """Test helper methods."""

    def test_calculate_purity(self):
        """Test purity calculation."""
        time_points = np.linspace(0, 100, 200)
        counts = np.ones(200)
        chrom = Chromatogram(time_points=time_points, counts=counts)

        compound = Compound([BuildingBlock.from_code(0, "Pro")], chrom)

        # Add peaks: product=100, truncation=50, null=30 → purity = 100/180 ≈ 0.556
        compound.detected_peaks = [
            Peak(position=10.0, left_base=5, right_base=15, height=20, area=30, peak_type=PeakType.NULL),
            Peak(position=20.0, left_base=15, right_base=25, height=30, area=50, peak_type=PeakType.TRUNCATION),
            Peak(position=40.0, left_base=35, right_base=45, height=50, area=100, peak_type=PeakType.PUTATIVE_PRODUCT),
        ]

        purity = PurityCalculator.calculate(compound)

        assert purity == pytest.approx(100 / 180, rel=0.01)

    def test_calculate_purity_no_peaks(self):
        """Test purity calculation with no peaks."""
        time_points = np.linspace(0, 100, 200)
        counts = np.ones(200)
        chrom = Chromatogram(time_points=time_points, counts=counts)

        compound = Compound([BuildingBlock.from_code(0, "Pro")], chrom)
        compound.detected_peaks = []

        purity = PurityCalculator.calculate(compound)

        assert purity == 0.0

    def test_calculate_snr(self):
        """Test SNR calculation."""
        time_points = np.linspace(0, 100, 200)
        counts = np.ones(200)
        chrom = Chromatogram(time_points=time_points, counts=counts)

        compound = Compound([BuildingBlock.from_code(0, "Pro")], chrom)

        peak = Peak(position=40.0, left_base=35, right_base=45, height=100, area=200, peak_type=PeakType.PUTATIVE_PRODUCT)
        compound.selected_peak = peak

        snr = SNRCalculator.calculate(compound, background=10.0)

        assert snr == pytest.approx(10.0)

    def test_check_retention_order_valid(self):
        """Test retention order checking with valid order."""
        validator = BayesianValidator()

        hierarchy = CompoundHierarchy(mode=HierarchyMode.BUILDING_BLOCK)
        time_points = np.linspace(0, 100, 200)
        counts = np.ones(200)
        chrom = Chromatogram(time_points=time_points, counts=counts)

        # Descendant at 20, ancestor at 40 → valid
        desc = Compound([BuildingBlock.from_code(0, "Pro")], chrom)
        desc.selected_peak = Peak(position=20.0, left_base=15, right_base=25, height=50, area=100)

        anc = Compound([BuildingBlock.from_code(0, "Pro"), BuildingBlock.from_code(1, "Leu")], chrom)
        anc.selected_peak = Peak(position=40.0, left_base=35, right_base=45, height=60, area=120)

        hierarchy.add_compound(desc)
        hierarchy.add_compound(anc)
        hierarchy.add_edge(anc, desc)

        valid = validator._check_retention_order(anc, hierarchy, precision=1.0)

        assert valid is True

    def test_check_retention_order_invalid(self):
        """Test retention order checking with invalid order."""
        validator = BayesianValidator()

        hierarchy = CompoundHierarchy(mode=HierarchyMode.BUILDING_BLOCK)
        time_points = np.linspace(0, 100, 200)
        counts = np.ones(200)
        chrom = Chromatogram(time_points=time_points, counts=counts)

        # Descendant at 40, ancestor at 20 → INVALID
        desc = Compound([BuildingBlock.from_code(0, "Pro")], chrom)
        desc.selected_peak = Peak(position=40.0, left_base=35, right_base=45, height=50, area=100)

        anc = Compound([BuildingBlock.from_code(0, "Pro"), BuildingBlock.from_code(1, "Leu")], chrom)
        anc.selected_peak = Peak(position=20.0, left_base=15, right_base=25, height=60, area=120)

        hierarchy.add_compound(desc)
        hierarchy.add_compound(anc)
        hierarchy.add_edge(anc, desc)

        valid = validator._check_retention_order(anc, hierarchy, precision=1.0)

        assert valid is False

    def test_get_descendant_validation_fraction(self):
        """Test descendant validation fraction calculation."""
        validator = BayesianValidator()

        hierarchy = CompoundHierarchy(mode=HierarchyMode.BUILDING_BLOCK)
        time_points = np.linspace(0, 100, 200)
        counts = np.ones(200)
        chrom = Chromatogram(time_points=time_points, counts=counts)

        # Create compounds
        anc = Compound([BuildingBlock.from_code(0, "Pro"), BuildingBlock.from_code(1, "Leu")], chrom)

        desc1 = Compound([BuildingBlock.from_code(0, "Pro")], chrom)
        desc1.selected_peak = Peak(position=20.0, left_base=15, right_base=25, height=50, area=100)
        desc1.selected_peak.validation_status = ValidationStatus.VALIDATED

        desc2 = Compound([BuildingBlock.from_code(0, "Null"), BuildingBlock.from_code(1, "Leu")], chrom)
        desc2.selected_peak = Peak(position=25.0, left_base=20, right_base=30, height=40, area=80)
        desc2.selected_peak.validation_status = ValidationStatus.FAILED

        hierarchy.add_compound(anc)
        hierarchy.add_compound(desc1)
        hierarchy.add_compound(desc2)
        hierarchy.add_edge(anc, desc1)
        hierarchy.add_edge(anc, desc2)

        fraction = validator._get_descendant_validation_fraction(anc, hierarchy)

        # 1 validated out of 2 → 0.5
        assert fraction == pytest.approx(0.5)


class TestBayesianValidatorLikelihoodRatio:
    """Test likelihood ratio computation."""

    def test_compute_likelihood_ratio_high_purity(self):
        """Test likelihood ratio with high purity."""
        validator = BayesianValidator()

        lr = validator.compute_likelihood_ratio(
            purity=0.95,
            retention_order_valid=True,
            descendant_fraction=1.0
        )

        # Should be strongly positive (>10)
        assert lr > 10.0

    def test_compute_likelihood_ratio_low_purity(self):
        """Test likelihood ratio with low purity."""
        validator = BayesianValidator()

        lr = validator.compute_likelihood_ratio(
            purity=0.2,
            retention_order_valid=True,
            descendant_fraction=0.5
        )

        # Should be less positive or negative
        assert lr > 0  # Still positive due to retention order

    def test_compute_likelihood_ratio_retention_violation(self):
        """Test likelihood ratio with retention violation."""
        validator = BayesianValidator()

        # With retention violation, LR should be lower than with valid retention
        lr_violation = validator.compute_likelihood_ratio(
            purity=0.9,
            retention_order_valid=False,
            descendant_fraction=1.0
        )

        lr_valid = validator.compute_likelihood_ratio(
            purity=0.9,
            retention_order_valid=True,
            descendant_fraction=1.0
        )

        # Violation should reduce LR (but may not make it < 1 if other factors are strong)
        assert lr_violation < lr_valid
        # Should reduce by factor of ~0.05/0.95 ≈ 0.053
        assert lr_violation / lr_valid < 0.1


class TestBayesianValidatorEdgeCases:
    """Test edge cases."""

    def test_validate_zero_background(self, simple_hierarchy):
        """Test validation with zero background (edge case)."""
        validator = BayesianValidator()

        dataset_stats = {
            'purity_p25': 0.5,
            'purity_p50': 0.7,
            'purity_p75': 0.85,
            'background': 0.0  # Zero background
        }

        l1 = [c for c in simple_hierarchy.compounds if c.level == 1][0]

        # Should handle gracefully (SNR will be 0 or inf)
        status = validator.validate(l1, simple_hierarchy, dataset_stats)

        assert status in [ValidationStatus.FAILED, ValidationStatus.NOT_VALIDATED]

    def test_validate_perfect_purity(self, simple_hierarchy, dataset_stats):
        """Test validation with perfect purity (1.0)."""
        validator = BayesianValidator()

        l1 = [c for c in simple_hierarchy.compounds if c.level == 1][0]

        # Only product peak (purity = 1.0)
        product_peak = Peak(position=30.0, left_base=25, right_base=35, height=100, area=200, peak_type=PeakType.PUTATIVE_PRODUCT)
        l1.detected_peaks = [product_peak]
        l1.selected_peak = product_peak

        status = validator.validate(l1, simple_hierarchy, dataset_stats)

        assert status == ValidationStatus.VALIDATED

    def test_validate_all_impurities(self, dataset_stats):
        """Test validation with all impurities (no product)."""
        validator = BayesianValidator()

        hierarchy = CompoundHierarchy(mode=HierarchyMode.BUILDING_BLOCK)
        time_points = np.linspace(0, 100, 200)
        counts = np.ones(200)
        chrom = Chromatogram(time_points=time_points, counts=counts)

        compound = Compound([BuildingBlock.from_code(0, "Pro")], chrom)

        # Only truncation peaks (purity = 0)
        trunc_peak = Peak(position=10.0, left_base=5, right_base=15, height=50, area=100, peak_type=PeakType.TRUNCATION)
        compound.detected_peaks = [trunc_peak]
        compound.selected_peak = None  # No product

        hierarchy.add_compound(compound)

        status = validator.validate(compound, hierarchy, dataset_stats)

        assert status == ValidationStatus.NOT_VALIDATED  # No selected peak
