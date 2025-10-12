"""
Comprehensive tests for AnalysisConfiguration model.

Tests configuration parameters for LC-Seq analysis.
"""

import pytest
from lcseq.domain.models.analysis_configuration import (
    AnalysisConfiguration,
    AnalysisMode,
    HierarchyMode,
)


class TestAnalysisConfigurationCreation:
    """Test analysis configuration creation."""

    def test_create_default_configuration(self):
        """Test creating configuration with defaults."""
        config = AnalysisConfiguration()

        assert config.analysis_mode == AnalysisMode.INDIVIDUAL
        assert config.hierarchy_mode == HierarchyMode.BUILDING_BLOCK
        assert config.baseline_params["method"] == "asls"
        assert config.baseline_params["p"] == 0.001
        assert config.baseline_params["lambda"] == 1e7
        assert config.peak_detection_params["min_persistence"] == 0.05
        assert config.validation_params["purity_threshold"] == "auto"

    def test_create_with_custom_modes(self):
        """Test creating with custom analysis and hierarchy modes."""
        config = AnalysisConfiguration(
            analysis_mode=AnalysisMode.CONSENSUS,
            hierarchy_mode=HierarchyMode.MONOMER,
        )

        assert config.analysis_mode == AnalysisMode.CONSENSUS
        assert config.hierarchy_mode == HierarchyMode.MONOMER

    def test_create_with_custom_baseline_params(self):
        """Test creating with custom baseline parameters."""
        config = AnalysisConfiguration(
            baseline_params={"method": "asls", "p": 0.01, "lambda": 1e6}
        )

        assert config.baseline_params["p"] == 0.01
        assert config.baseline_params["lambda"] == 1e6

    def test_create_with_custom_peak_detection_params(self):
        """Test creating with custom peak detection parameters."""
        config = AnalysisConfiguration(
            peak_detection_params={
                "min_persistence": 0.1,
                "boundary_method": "valley_or_5pct",
            }
        )

        assert config.peak_detection_params["min_persistence"] == 0.1

    def test_create_with_custom_validation_params(self):
        """Test creating with custom validation parameters."""
        config = AnalysisConfiguration(
            validation_params={"purity_threshold": 0.7, "snr_threshold": 10.0}
        )

        assert config.validation_params["purity_threshold"] == 0.7
        assert config.validation_params["snr_threshold"] == 10.0


class TestDefaultFactory:
    """Test default() factory method."""

    def test_default_factory(self):
        """Test default factory creates standard configuration."""
        config = AnalysisConfiguration.default()

        assert config.analysis_mode == AnalysisMode.INDIVIDUAL
        assert config.hierarchy_mode == HierarchyMode.BUILDING_BLOCK
        assert config.baseline_params["p"] == 0.001
        assert config.baseline_params["lambda"] == 1e7

    def test_default_creates_independent_instances(self):
        """Test that default creates independent instances."""
        config1 = AnalysisConfiguration.default()
        config2 = AnalysisConfiguration.default()

        # Modify one
        config1.baseline_params["p"] = 0.01

        # Should not affect the other
        assert config2.baseline_params["p"] == 0.001


class TestConsensusFactory:
    """Test for_consensus_mode() factory method."""

    def test_consensus_mode_factory(self):
        """Test consensus mode factory."""
        config = AnalysisConfiguration.for_consensus_mode()

        assert config.analysis_mode == AnalysisMode.CONSENSUS
        assert config.validation_params["correlation_threshold"] == 0.8

    def test_consensus_mode_custom_threshold(self):
        """Test consensus mode with custom correlation threshold."""
        config = AnalysisConfiguration.for_consensus_mode(correlation_threshold=0.85)

        assert config.validation_params["correlation_threshold"] == 0.85

    def test_consensus_mode_has_all_defaults(self):
        """Test consensus mode includes all default parameters."""
        config = AnalysisConfiguration.for_consensus_mode()

        assert config.baseline_params["p"] == 0.001
        assert config.peak_detection_params["min_persistence"] == 0.05


class TestMonomerModeFactory:
    """Test for_monomer_mode() factory method."""

    def test_monomer_mode_factory(self):
        """Test monomer mode factory."""
        config = AnalysisConfiguration.for_monomer_mode()

        assert config.hierarchy_mode == HierarchyMode.MONOMER
        assert config.analysis_mode == AnalysisMode.INDIVIDUAL

    def test_monomer_mode_has_all_defaults(self):
        """Test monomer mode includes all default parameters."""
        config = AnalysisConfiguration.for_monomer_mode()

        assert config.baseline_params["p"] == 0.001
        assert config.peak_detection_params["min_persistence"] == 0.05


class TestParameterValidation:
    """Test parameter validation."""

    def test_invalid_baseline_method_raises_error(self):
        """Test invalid baseline method raises error."""
        with pytest.raises(ValueError, match="Invalid baseline method"):
            AnalysisConfiguration(baseline_params={"method": "invalid"})

    def test_invalid_asls_p_below_zero_raises_error(self):
        """Test AsLS p below 0 raises error."""
        with pytest.raises(ValueError, match="AsLS parameter p must be in"):
            AnalysisConfiguration(baseline_params={"p": -0.1})

    def test_invalid_asls_p_above_one_raises_error(self):
        """Test AsLS p above 1 raises error."""
        with pytest.raises(ValueError, match="AsLS parameter p must be in"):
            AnalysisConfiguration(baseline_params={"p": 1.5})

    def test_invalid_asls_lambda_negative_raises_error(self):
        """Test negative AsLS lambda raises error."""
        with pytest.raises(ValueError, match="AsLS parameter lambda must be positive"):
            AnalysisConfiguration(baseline_params={"lambda": -1e7})

    def test_invalid_asls_lambda_zero_raises_error(self):
        """Test zero AsLS lambda raises error."""
        with pytest.raises(ValueError, match="AsLS parameter lambda must be positive"):
            AnalysisConfiguration(baseline_params={"lambda": 0})

    def test_invalid_persistence_below_zero_raises_error(self):
        """Test min_persistence below 0 raises error."""
        with pytest.raises(ValueError, match="min_persistence must be in"):
            AnalysisConfiguration(peak_detection_params={"min_persistence": -0.1})

    def test_invalid_persistence_above_one_raises_error(self):
        """Test min_persistence above 1 raises error."""
        with pytest.raises(ValueError, match="min_persistence must be in"):
            AnalysisConfiguration(peak_detection_params={"min_persistence": 1.5})

    def test_invalid_correlation_threshold_below_zero(self):
        """Test correlation threshold below 0 raises error."""
        with pytest.raises(ValueError, match="correlation_threshold must be in"):
            AnalysisConfiguration(
                analysis_mode=AnalysisMode.CONSENSUS,
                validation_params={"correlation_threshold": -0.1},
            )

    def test_invalid_correlation_threshold_above_one(self):
        """Test correlation threshold above 1 raises error."""
        with pytest.raises(ValueError, match="correlation_threshold must be in"):
            AnalysisConfiguration(
                analysis_mode=AnalysisMode.CONSENSUS,
                validation_params={"correlation_threshold": 1.5},
            )


class TestConsensusCorrelationThreshold:
    """Test consensus mode correlation threshold handling."""

    def test_consensus_mode_adds_correlation_threshold(self):
        """Test consensus mode automatically adds correlation threshold."""
        config = AnalysisConfiguration(analysis_mode=AnalysisMode.CONSENSUS)

        assert "correlation_threshold" in config.validation_params
        assert config.validation_params["correlation_threshold"] == 0.8

    def test_consensus_mode_preserves_custom_correlation_threshold(self):
        """Test consensus mode preserves custom correlation threshold."""
        config = AnalysisConfiguration(
            analysis_mode=AnalysisMode.CONSENSUS,
            validation_params={"correlation_threshold": 0.9},
        )

        assert config.validation_params["correlation_threshold"] == 0.9

    def test_individual_mode_no_correlation_threshold(self):
        """Test individual mode doesn't require correlation threshold."""
        config = AnalysisConfiguration(analysis_mode=AnalysisMode.INDIVIDUAL)

        # Should not have correlation threshold in individual mode
        assert "correlation_threshold" not in config.validation_params


class TestCopyMethod:
    """Test configuration copying."""

    def test_copy_creates_independent_instance(self):
        """Test copy creates independent instance."""
        config1 = AnalysisConfiguration.default()
        config2 = config1.copy()

        # Modify original
        config1.baseline_params["p"] = 0.01

        # Copy should be unchanged
        assert config2.baseline_params["p"] == 0.001

    def test_copy_preserves_all_attributes(self):
        """Test copy preserves all configuration attributes."""
        config1 = AnalysisConfiguration(
            analysis_mode=AnalysisMode.CONSENSUS,
            hierarchy_mode=HierarchyMode.MONOMER,
            baseline_params={"p": 0.01, "lambda": 1e6},
            peak_detection_params={"min_persistence": 0.1},
            validation_params={"purity_threshold": 0.7},
            classification_params={"custom": "value"},
        )

        config2 = config1.copy()

        assert config2.analysis_mode == AnalysisMode.CONSENSUS
        assert config2.hierarchy_mode == HierarchyMode.MONOMER
        assert config2.baseline_params["p"] == 0.01
        assert config2.peak_detection_params["min_persistence"] == 0.1
        assert config2.validation_params["purity_threshold"] == 0.7
        assert config2.classification_params["custom"] == "value"

    def test_copy_creates_deep_copy_of_dicts(self):
        """Test copy creates deep copies of parameter dicts."""
        config1 = AnalysisConfiguration.default()
        config2 = config1.copy()

        # Modify dict in original
        config1.baseline_params["new_key"] = "new_value"

        # Should not affect copy
        assert "new_key" not in config2.baseline_params


class TestStringRepresentations:
    """Test string representations."""

    def test_repr_includes_key_info(self):
        """Test repr includes key configuration info."""
        config = AnalysisConfiguration(
            analysis_mode=AnalysisMode.CONSENSUS,
            hierarchy_mode=HierarchyMode.MONOMER,
        )

        repr_str = repr(config)
        assert "AnalysisConfiguration" in repr_str
        assert "consensus" in repr_str
        assert "monomer" in repr_str
        assert "asls" in repr_str

    def test_repr_shows_persistence(self):
        """Test repr shows min_persistence."""
        config = AnalysisConfiguration(
            peak_detection_params={"min_persistence": 0.1}
        )

        repr_str = repr(config)
        assert "min_persistence=0.1" in repr_str


class TestParameterDefaults:
    """Test parameter default values."""

    def test_baseline_defaults(self):
        """Test baseline parameter defaults."""
        config = AnalysisConfiguration()

        assert config.baseline_params["method"] == "asls"
        assert config.baseline_params["p"] == 0.001
        assert config.baseline_params["lambda"] == 1e7

    def test_peak_detection_defaults(self):
        """Test peak detection parameter defaults."""
        config = AnalysisConfiguration()

        assert config.peak_detection_params["min_persistence"] == 0.05
        assert config.peak_detection_params["boundary_method"] == "valley_or_5pct"

    def test_validation_defaults(self):
        """Test validation parameter defaults."""
        config = AnalysisConfiguration()

        assert config.validation_params["purity_threshold"] == "auto"
        assert config.validation_params["snr_threshold"] == "auto"

    def test_classification_params_empty_by_default(self):
        """Test classification params empty by default."""
        config = AnalysisConfiguration()

        assert config.classification_params == {}


class TestCustomParameterOverrides:
    """Test that custom parameters override defaults."""

    def test_custom_baseline_overrides_defaults(self):
        """Test custom baseline params override defaults."""
        config = AnalysisConfiguration(
            baseline_params={"method": "asls", "p": 0.005}
        )

        # Custom value
        assert config.baseline_params["p"] == 0.005
        # Still has other defaults
        assert config.baseline_params["lambda"] == 1e7

    def test_partial_override_preserves_unspecified_defaults(self):
        """Test partial override preserves other defaults."""
        config = AnalysisConfiguration(
            peak_detection_params={"min_persistence": 0.1}
        )

        assert config.peak_detection_params["min_persistence"] == 0.1
        assert config.peak_detection_params["boundary_method"] == "valley_or_5pct"
