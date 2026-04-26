"""
cLPE (chromatographic Linear Peptide Equation) Validator.

Validates product peak selection by checking if the observed retention time
is consistent with the compound's lipophilicity (AlogP) based on a linear
regression model fitted per scaffold group.

Principle:
    LogK = slope × AlogP + intercept

Where LogK = log10((RT - t0) / t0) is the log capacity factor.

Compounds with the same scaffold should fall on a line. Outliers from this
relationship may indicate incorrect peak selection.

When an outlier is detected, the validator attempts to re-select a better
peak from the candidate peaks.
"""

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple
import numpy as np
from sklearn.linear_model import LinearRegression

from ..entities.compound import Compound
from ..entities.peak import Peak, PeakType


@dataclass
class CLPEModel:
    """Fitted cLPE model for a scaffold group."""
    scaffold: str
    slope: float
    intercept: float
    r_squared: float
    residual_std: float
    n_compounds: int


@dataclass
class CLPEResult:
    """Result of cLPE validation for a single compound."""
    observed_logk: float
    predicted_logk: float
    alogp: float
    residual: float
    z_score: float
    is_outlier: bool
    scaffold: str
    reselected_peak: Optional[Peak] = None


class CLPEValidator:
    """
    Domain service for cLPE-based peak validation and re-selection.

    The cLPE (chromatographic Linear Peptide Equation) principle states that
    compounds with the same scaffold should have retention times that correlate
    linearly with their calculated lipophilicity (AlogP).

    This service:
    1. Fits LogK ~ AlogP regression models per scaffold group
    2. Validates selected peaks against predicted LogK
    3. Attempts to re-select better peaks for outliers

    Parameters
    ----------
    t0 : float
        Column dead time in minutes (default: 1.0)
    outlier_threshold : float
        Number of standard deviations to flag as outlier (default: 2.5)
    min_group_size : int
        Minimum compounds per scaffold for model fitting (default: 5)

    Examples
    --------
    >>> validator = CLPEValidator(t0=1.0, outlier_threshold=2.5)
    >>> validator.fit_models(compounds, reference_data)
    >>> result, new_peak = validator.validate_and_reselect(compound)
    >>> if new_peak:
    ...     compound.selected_peak = new_peak
    """

    def __init__(
        self,
        t0: float = 1.0,
        outlier_threshold: float = 2.5,
        min_group_size: int = 5
    ):
        self.t0 = t0
        self.outlier_threshold = outlier_threshold
        self.min_group_size = min_group_size
        self.models: Dict[str, CLPEModel] = {}

    def compute_logk(self, rt: float) -> float:
        """
        Compute log capacity factor from retention time.

        Parameters
        ----------
        rt : float
            Retention time in minutes

        Returns
        -------
        float
            LogK = log10((RT - t0) / t0)
        """
        if rt <= self.t0:
            return float('-inf')
        return np.log10((rt - self.t0) / self.t0)

    def rt_from_logk(self, logk: float) -> float:
        """
        Compute retention time from log capacity factor.

        Parameters
        ----------
        logk : float
            Log capacity factor

        Returns
        -------
        float
            RT = t0 * (10^logk + 1)
        """
        return self.t0 * (10**logk + 1)

    def fit_models(
        self,
        compounds: List[Compound],
        alogp_map: Dict[str, float],
        scaffold_map: Dict[str, str],
        expected_logk_map: Optional[Dict[str, float]] = None
    ) -> Dict[str, CLPEModel]:
        """
        Fit cLPE regression models for each scaffold group.

        Parameters
        ----------
        compounds : List[Compound]
            Compounds with selected peaks
        alogp_map : Dict[str, float]
            Mapping from compound identifier to AlogP value
        scaffold_map : Dict[str, str]
            Mapping from compound identifier to scaffold group
        expected_logk_map : Dict[str, float], optional
            If provided, use these LogK values for fitting instead of
            computing from observed RT

        Returns
        -------
        Dict[str, CLPEModel]
            Fitted models keyed by scaffold
        """
        # Group data by scaffold
        scaffold_data: Dict[str, List[Tuple[float, float]]] = {}

        for compound in compounds:
            compound_id = self._get_compound_id(compound)

            if compound_id not in alogp_map:
                continue
            if compound_id not in scaffold_map:
                continue

            alogp = alogp_map[compound_id]
            scaffold = scaffold_map[compound_id]

            # Get LogK (from reference or observed)
            if expected_logk_map and compound_id in expected_logk_map:
                logk = expected_logk_map[compound_id]
            elif compound.selected_peak:
                logk = self.compute_logk(compound.selected_peak.position)
            else:
                continue

            if np.isinf(logk) or np.isnan(logk):
                continue
            if np.isnan(alogp):
                continue

            if scaffold not in scaffold_data:
                scaffold_data[scaffold] = []
            scaffold_data[scaffold].append((alogp, logk))

        # Fit model for each scaffold
        self.models = {}

        for scaffold, data in scaffold_data.items():
            if len(data) < self.min_group_size:
                continue

            X = np.array([d[0] for d in data]).reshape(-1, 1)
            y = np.array([d[1] for d in data])

            model = LinearRegression()
            model.fit(X, y)

            y_pred = model.predict(X)
            residuals = y - y_pred

            # R-squared
            ss_res = np.sum(residuals**2)
            ss_tot = np.sum((y - np.mean(y))**2)
            r_squared = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0.0

            self.models[scaffold] = CLPEModel(
                scaffold=scaffold,
                slope=model.coef_[0],
                intercept=model.intercept_,
                r_squared=r_squared,
                residual_std=np.std(residuals),
                n_compounds=len(data)
            )

        return self.models

    def validate_peak(
        self,
        compound: Compound,
        alogp: float,
        scaffold: str
    ) -> CLPEResult:
        """
        Validate if the selected peak is consistent with cLPE model.

        Parameters
        ----------
        compound : Compound
            Compound with selected peak
        alogp : float
            Compound's AlogP value
        scaffold : str
            Compound's scaffold group

        Returns
        -------
        CLPEResult
            Validation result including residual and outlier flag
        """
        if scaffold not in self.models:
            # No model for this scaffold, cannot validate
            return CLPEResult(
                observed_logk=float('nan'),
                predicted_logk=float('nan'),
                alogp=alogp,
                residual=float('nan'),
                z_score=float('nan'),
                is_outlier=False,
                scaffold=scaffold
            )

        model = self.models[scaffold]

        if not compound.selected_peak:
            return CLPEResult(
                observed_logk=float('nan'),
                predicted_logk=model.slope * alogp + model.intercept,
                alogp=alogp,
                residual=float('nan'),
                z_score=float('nan'),
                is_outlier=False,
                scaffold=scaffold
            )

        observed_logk = self.compute_logk(compound.selected_peak.position)
        predicted_logk = model.slope * alogp + model.intercept
        residual = observed_logk - predicted_logk

        # Z-score relative to model's residual standard deviation
        if model.residual_std > 0:
            z_score = residual / model.residual_std
        else:
            z_score = 0.0

        is_outlier = abs(z_score) > self.outlier_threshold

        return CLPEResult(
            observed_logk=observed_logk,
            predicted_logk=predicted_logk,
            alogp=alogp,
            residual=residual,
            z_score=z_score,
            is_outlier=is_outlier,
            scaffold=scaffold
        )

    def validate_and_reselect(
        self,
        compound: Compound,
        alogp: float,
        scaffold: str
    ) -> Tuple[CLPEResult, Optional[Peak]]:
        """
        Validate selected peak and attempt re-selection if outlier.

        If the selected peak is an outlier, tries to find a better peak
        from the candidate peaks (UNKNOWN or unassigned peaks).

        Parameters
        ----------
        compound : Compound
            Compound with detected peaks
        alogp : float
            Compound's AlogP value
        scaffold : str
            Compound's scaffold group

        Returns
        -------
        Tuple[CLPEResult, Optional[Peak]]
            (validation_result, new_peak) where new_peak is None if
            the original selection is acceptable or no better alternative found
        """
        result = self.validate_peak(compound, alogp, scaffold)

        if not result.is_outlier:
            return result, None

        # Try to find a better peak
        new_peak = self._find_best_clpe_peak(compound, alogp, scaffold)

        if new_peak and new_peak != compound.selected_peak:
            result.reselected_peak = new_peak
            return result, new_peak

        return result, None

    def _find_best_clpe_peak(
        self,
        compound: Compound,
        alogp: float,
        scaffold: str
    ) -> Optional[Peak]:
        """
        Find the peak with the best cLPE fit from candidates.

        Considers all detected peaks that could be the product peak.

        Parameters
        ----------
        compound : Compound
            Compound with detected peaks
        alogp : float
            Compound's AlogP value
        scaffold : str
            Compound's scaffold group

        Returns
        -------
        Optional[Peak]
            Best fitting peak, or None if no acceptable peak found
        """
        if scaffold not in self.models:
            return None

        model = self.models[scaffold]
        predicted_logk = model.slope * alogp + model.intercept

        # Candidate peaks: UNKNOWN or PUTATIVE_PRODUCT only
        # Exclude NULL (always at t0) and TRUNCATION (claimed by ancestors)
        # Re-selecting a TRUNCATION peak would be chemically inconsistent:
        # a full-length product cannot elute at the same time as its truncation
        candidates = [
            p for p in compound.detected_peaks
            if p.peak_type not in (PeakType.NULL, PeakType.TRUNCATION, PeakType.TRUNCATION_UNKNOWN)
            and p.area is not None and p.area > 0
            and p.is_accepted
        ]

        if not candidates:
            return None

        best_peak = None
        best_residual = float('inf')

        for peak in candidates:
            logk = self.compute_logk(peak.position)
            if np.isinf(logk):
                continue

            residual = abs(logk - predicted_logk)

            if residual < best_residual:
                best_residual = residual
                best_peak = peak

        # Only accept if within threshold
        if model.residual_std > 0:
            if best_residual <= self.outlier_threshold * model.residual_std:
                return best_peak

        return None

    def _get_compound_id(self, compound: Compound) -> str:
        """
        Get a unique identifier for a compound.

        Uses compound_id if available, otherwise constructs from
        building block sequence.
        """
        if compound.compound_id:
            return compound.compound_id

        # Fallback: construct from building blocks
        return compound.block_support_sequence

    def get_model_summary(self) -> str:
        """
        Get a summary of fitted models.

        Returns
        -------
        str
            Human-readable summary of all fitted models
        """
        lines = [f"cLPE Models ({len(self.models)} scaffolds):"]
        for scaffold, model in sorted(self.models.items()):
            lines.append(
                f"  {scaffold}: slope={model.slope:.3f}, "
                f"intercept={model.intercept:.3f}, "
                f"R²={model.r_squared:.3f}, "
                f"σ={model.residual_std:.3f}, "
                f"n={model.n_compounds}"
            )
        return "\n".join(lines)
