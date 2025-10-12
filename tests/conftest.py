"""Pytest configuration and shared fixtures."""

import numpy as np
import pytest


@pytest.fixture
def simple_gaussian_chromatogram():
    """Generate a chromatogram with a single Gaussian peak.

    Returns:
        tuple: (time_points, counts) arrays
    """
    time = np.linspace(0, 100, 1000)
    # Gaussian peak centered at t=50, width=5, height=100
    signal = 100 * np.exp(-((time - 50) ** 2) / (2 * 5**2))
    # Add small baseline
    baseline = 10 * np.ones_like(time)
    # Add noise
    rng = np.random.default_rng(42)
    noise = rng.normal(0, 2, size=time.shape)
    counts = np.maximum(0, signal + baseline + noise).astype(np.int64)
    return time, counts


@pytest.fixture
def multi_peak_chromatogram():
    """Generate a chromatogram with multiple peaks.

    Returns:
        tuple: (time_points, counts) arrays with peaks at t=30, 50, 70
    """
    time = np.linspace(0, 100, 1000)
    # Three Gaussian peaks
    peak1 = 80 * np.exp(-((time - 30) ** 2) / (2 * 4**2))
    peak2 = 100 * np.exp(-((time - 50) ** 2) / (2 * 5**2))
    peak3 = 60 * np.exp(-((time - 70) ** 2) / (2 * 3**2))
    baseline = 10 * np.ones_like(time)
    rng = np.random.default_rng(42)
    noise = rng.normal(0, 2, size=time.shape)
    counts = np.maximum(0, peak1 + peak2 + peak3 + baseline + noise).astype(np.int64)
    return time, counts


@pytest.fixture
def flat_chromatogram():
    """Generate a flat chromatogram with no peaks.

    Returns:
        tuple: (time_points, counts) arrays
    """
    time = np.linspace(0, 100, 1000)
    baseline = 10 * np.ones_like(time)
    rng = np.random.default_rng(42)
    noise = rng.normal(0, 1, size=time.shape)
    counts = np.maximum(0, baseline + noise).astype(np.int64)
    return time, counts
