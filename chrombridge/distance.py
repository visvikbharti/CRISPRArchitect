"""
Chromatin 3D Distance Prediction and Donor Bridgeability Analysis
=================================================================

This module converts the abstract polymer physics predictions (mean-squared
distances in nm^2) into practical quantities that a molecular biologist
can act on:

    - **Mean 3D distance** between two loci (in nm and um)
    - **Uncertainty** (standard deviation, 5th/95th percentiles)
    - **Donor bridgeability**: Can a donor template of a given size
      physically span the distance between two loci?

Why does 3D distance matter for CRISPR?
---------------------------------------
When designing a multi-site editing experiment — for example, using a single
donor to correct mutations in both exon 40 and exon 60 of a large gene — the
donor template must physically reach both cut sites. If the two sites are
separated by 2 Mb on the genome, they are typically ~500-1500 nm apart in 3D.
A 3 kb cssDNA donor, on the other hand, forms a coil of only ~50 nm in
diameter. It is like trying to bridge two buildings with a piece of string.

The distribution of 3D distances
--------------------------------
The distance between two points on a Gaussian polymer chain is NOT a single
number — it is a probability distribution. This distribution arises because
the polymer is constantly changing shape due to thermal fluctuations (Brownian
motion). The chromatin fiber is jiggling, bending, and writhing, so the 3D
distance between two loci varies from cell to cell and from moment to moment.

For a Gaussian chain, the probability distribution of the end-to-end distance
R is a **Maxwell-Boltzmann-like** (or "Rayleigh") distribution in 3D:

    P(R) dR = (4*pi*R^2) * (3 / (2*pi*<R^2>))^(3/2) * exp(-3*R^2 / (2*<R^2>)) dR

This looks like a chi distribution with 3 degrees of freedom, scaled by
sqrt(<R^2>/3). Key properties:

    - Mean:    <R> = sqrt(8 / (3*pi)) * sqrt(<R^2>) ~ 0.921 * sqrt(<R^2>)
    - Median:  ~ 0.939 * sqrt(<R^2>)
    - Mode:    sqrt(2/3) * sqrt(<R^2>) ~ 0.816 * sqrt(<R^2>)
    - StdDev:  sqrt(1 - 8/(3*pi)) * sqrt(<R^2>) ~ 0.389 * sqrt(<R^2>)

The distribution is **right-skewed**: most of the time the loci are closer
than the RMS distance, but occasionally the polymer stretches out and they
end up farther apart.

References
----------
- Doi & Edwards, "The Theory of Polymer Dynamics", Oxford, 1988, Ch. 2.
- Rubinstein & Colby, "Polymer Physics", Oxford, 2003, Ch. 2.
- Bystricky et al., PNAS, 2004 (experimental 3D distances in yeast).
- Giorgetti et al., Cell, 2014 (distribution of chromatin distances).
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import List, Optional, Tuple

import sys
import os

_PARENT = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
if _PARENT not in sys.path:
    sys.path.insert(0, _PARENT)

from .polymer_model import (
    GaussianChainModel,
    WormLikeChainModel,
    ConfinedPolymerModel,
    donor_rms_size_nm,
    donor_radius_of_gyration_nm,
)


# =============================================================================
# Constants for the Rayleigh (Maxwell-Boltzmann) distribution of R
# =============================================================================

# For a 3D Gaussian chain, R = |end-to-end vector| follows a distribution
# where the key ratios relative to sqrt(<R^2>) are:
_MEAN_OVER_RMS = math.sqrt(8.0 / (3.0 * math.pi))    # ~0.9213
_MEDIAN_OVER_RMS = 0.9394                              # numerically computed
_STD_OVER_RMS = math.sqrt(1.0 - 8.0 / (3.0 * math.pi))  # ~0.3893

# Percentile multipliers (from the chi distribution with k=3):
# The CDF of R / sigma_component (where sigma_component = sqrt(<R^2>/3)) is
# the chi distribution with 3 degrees of freedom.
# P(R < x * sqrt(<R^2>)) = p  =>  we precompute x for common percentiles.
# These are computed from scipy.stats.chi(df=3).ppf(p) / sqrt(3).
_P5_OVER_RMS = 0.3935     # 5th percentile
_P95_OVER_RMS = 1.4289    # 95th percentile


# =============================================================================
# Data Classes
# =============================================================================

@dataclass
class DistancePrediction:
    """Predicted 3D distance between two genomic loci.

    All distances are derived from the polymer physics models in
    ``polymer_model.py``. The distribution assumes a 3D Gaussian chain
    (Maxwell-Boltzmann distribution for the scalar distance R).

    Attributes
    ----------
    genomic_distance_bp : int
        Input genomic separation in base pairs.
    chromatin_state : str
        Chromatin state used for the prediction.
    model_used : str
        Name of the polymer model used ("gaussian", "wlc", "confined").
    mean_3d_distance_nm : float
        Mean (expected value) of the 3D distance distribution, in nm.
    std_nm : float
        Standard deviation of the 3D distance distribution, in nm.
    median_nm : float
        Median of the 3D distance distribution, in nm.
    p5_nm : float
        5th percentile: the distance below which the loci are found only
        5% of the time. Represents a "close encounter" configuration.
    p95_nm : float
        95th percentile: the distance below which the loci are found 95%
        of the time. Represents a "stretched out" configuration.
    mean_3d_distance_um : float
        Mean 3D distance in micrometers (convenience conversion).
    rms_distance_nm : float
        Root-mean-square distance in nm (sqrt(<R^2>)).
    """

    genomic_distance_bp: int
    chromatin_state: str
    model_used: str
    mean_3d_distance_nm: float
    std_nm: float
    median_nm: float
    p5_nm: float
    p95_nm: float
    mean_3d_distance_um: float
    rms_distance_nm: float


@dataclass
class BridgeabilityResult:
    """Assessment of whether a donor template can physically bridge two loci.

    When we talk about "bridging", we mean: can a single donor molecule,
    floating freely inside the nucleus as a random coil, simultaneously
    interact with two genomic loci separated by some 3D distance?

    For this to work, the donor's coil diameter must be comparable to or
    larger than the inter-locus distance. In practice, this is almost
    never the case for loci separated by more than ~10-50 kb, because
    donor templates (1-5 kb) form tiny coils (~30-100 nm) while
    megabase-scale loci are separated by hundreds to thousands of nm.

    Attributes
    ----------
    feasible : bool
        True if bridgeability_ratio >= 0.5 (generous threshold).
        In practice, even a ratio of 0.5 is marginal.
    inter_locus_distance_nm : float
        Predicted mean 3D distance between the two loci, in nm.
    donor_coil_diameter_nm : float
        Estimated diameter of the donor coil (2 * Rg), in nm.
    bridgeability_ratio : float
        donor_coil_diameter / inter_locus_distance. Values:
        - > 1.0: donor is larger than the gap (bridging plausible)
        - 0.5 - 1.0: marginal (possible with favorable fluctuations)
        - 0.1 - 0.5: unlikely (requires rare conformations)
        - < 0.1: physically implausible
    explanation : str
        Detailed biological explanation of the result, written for a
        molecular biologist.
    """

    feasible: bool
    inter_locus_distance_nm: float
    donor_coil_diameter_nm: float
    bridgeability_ratio: float
    explanation: str


# =============================================================================
# Main Predictor Class
# =============================================================================

class ChromatinDistancePredictor:
    """Predict 3D chromatin distances and assess donor bridgeability.

    This is the main user-facing class for spatial analysis in ChromBridge.
    It wraps the polymer physics models and converts their output into
    biologically meaningful predictions with uncertainty estimates.

    Usage
    -----
    >>> predictor = ChromatinDistancePredictor()
    >>> result = predictor.predict_3d_distance(2_000_000, "euchromatin")
    >>> print(f"Mean distance: {result.mean_3d_distance_nm:.0f} nm "
    ...       f"({result.mean_3d_distance_um:.2f} um)")
    >>> print(f"Range (5th-95th): {result.p5_nm:.0f} - {result.p95_nm:.0f} nm")
    """

    def __init__(self) -> None:
        """Initialize with all three polymer models."""
        self._gaussian = GaussianChainModel()
        self._wlc = WormLikeChainModel()
        self._confined = ConfinedPolymerModel()

    def _get_model(self, model: str):
        """Select the polymer model by name.

        Parameters
        ----------
        model : str
            "gaussian", "wlc", or "confined".

        Returns
        -------
        One of GaussianChainModel, WormLikeChainModel, ConfinedPolymerModel.
        """
        models = {
            "gaussian": self._gaussian,
            "wlc": self._wlc,
            "confined": self._confined,
        }
        if model not in models:
            raise ValueError(
                f"Unknown model '{model}'. Choose from: {list(models.keys())}"
            )
        return models[model]

    def predict_3d_distance(
        self,
        genomic_distance_bp: int,
        chromatin_state: str = "euchromatin",
        model: str = "confined",
    ) -> DistancePrediction:
        """Predict the 3D spatial distance between two loci on the same chromosome.

        This is the primary method for getting a distance prediction. It
        returns not just a single number, but the full statistical picture:
        mean, median, standard deviation, and percentiles.

        **Why is there a distribution?**
        The chromatin fiber is constantly in thermal motion (jiggling). At
        any given moment, the 3D distance between two loci varies from cell
        to cell in a population, and even within a single cell over time.
        The polymer model predicts the *statistical* distribution of this
        distance, not a fixed value.

        Parameters
        ----------
        genomic_distance_bp : int
            Genomic separation between the two loci in base pairs. Must be
            on the same chromosome. For different chromosomes, use the
            translocation module instead.
        chromatin_state : str
            "euchromatin" (active, open chromatin — default), "heterochromatin"
            (silent, compact chromatin), or "intermediate".
        model : str
            Which polymer model to use:
            - "confined" (default): Most realistic. WLC + territory confinement.
              Recommended for most applications.
            - "wlc": Worm-Like Chain without confinement. Good for distances
              < 5 Mb but over-predicts at larger separations.
            - "gaussian": Simplest. Good for quick estimates but ignores
              stiffness at short distances and confinement at large distances.

        Returns
        -------
        DistancePrediction
            Dataclass with mean, std, median, p5, p95 distances in nm,
            plus convenience conversion to um.

        Example
        -------
        >>> predictor = ChromatinDistancePredictor()
        >>> # Predict distance for 2 Mb separation (e.g., exon 40 to exon 60)
        >>> pred = predictor.predict_3d_distance(2_000_000)
        >>> print(f"Mean: {pred.mean_3d_distance_nm:.0f} nm = "
        ...       f"{pred.mean_3d_distance_um:.2f} um")
        >>> print(f"5th-95th percentile: {pred.p5_nm:.0f} - {pred.p95_nm:.0f} nm")
        """
        polymer_model = self._get_model(model)
        r2 = polymer_model.mean_squared_distance(
            genomic_distance_bp, chromatin_state
        )
        rms = math.sqrt(r2)

        # Compute distribution statistics from the 3D Gaussian chain
        # distance distribution (Maxwell-Boltzmann / chi with k=3).
        mean_r = _MEAN_OVER_RMS * rms
        median_r = _MEDIAN_OVER_RMS * rms
        std_r = _STD_OVER_RMS * rms
        p5_r = _P5_OVER_RMS * rms
        p95_r = _P95_OVER_RMS * rms

        return DistancePrediction(
            genomic_distance_bp=int(genomic_distance_bp),
            chromatin_state=chromatin_state,
            model_used=model,
            mean_3d_distance_nm=mean_r,
            std_nm=std_r,
            median_nm=median_r,
            p5_nm=p5_r,
            p95_nm=p95_r,
            mean_3d_distance_um=mean_r / 1000.0,
            rms_distance_nm=rms,
        )

    def can_donor_bridge(
        self,
        genomic_distance_bp: int,
        donor_size_bp: int,
        donor_type: str = "cssDNA",
        chromatin_state: str = "euchromatin",
    ) -> BridgeabilityResult:
        """Assess whether a donor template can physically bridge two genomic loci.

        This answers the critical question: "If I want to use a single donor
        template to mediate editing at two cut sites that are X bp apart on
        the genome, is this physically possible?"

        **The short answer for most cases: NO.**

        Here is the intuition. A 3 kb cssDNA donor, when floating freely
        in the nucleus, forms a random coil with a radius of gyration of
        about 20 nm and a diameter of about 40 nm. This is its "reach".
        Two loci separated by 2 Mb on the genome are typically ~500-1500 nm
        apart in 3D space. The donor would need to stretch across a distance
        30-80 times its own size. That is like trying to bridge a swimming
        pool with a golf ball.

        The donor's coil diameter is computed from the Gaussian chain model
        for the appropriate DNA type (ssDNA or dsDNA), and the inter-locus
        distance comes from the confined polymer model.

        Parameters
        ----------
        genomic_distance_bp : int
            Genomic separation between the two target loci in bp.
        donor_size_bp : int
            Size of the donor template in bp (for dsDNA) or nt (for ssDNA).
        donor_type : str
            "dsDNA", "ssDNA", "cssDNA", or "lssDNA".
        chromatin_state : str
            Chromatin state at the target loci.

        Returns
        -------
        BridgeabilityResult
            Contains: feasibility boolean, distances, ratio, and a detailed
            biological explanation.

        Example
        -------
        >>> predictor = ChromatinDistancePredictor()
        >>> result = predictor.can_donor_bridge(
        ...     genomic_distance_bp=2_000_000,  # 2 Mb between exon 40 and exon 60
        ...     donor_size_bp=3000,              # 3 kb cssDNA donor
        ...     donor_type="cssDNA"
        ... )
        >>> print(result.explanation)
        """
        # Calculate inter-locus 3D distance
        distance_pred = self.predict_3d_distance(
            genomic_distance_bp, chromatin_state, model="confined"
        )
        inter_locus_nm = distance_pred.mean_3d_distance_nm

        # Calculate donor coil size
        # Diameter ~ 2 * Rg (radius of gyration)
        donor_rg = donor_radius_of_gyration_nm(donor_size_bp, donor_type)
        donor_diameter = 2.0 * donor_rg
        donor_rms = donor_rms_size_nm(donor_size_bp, donor_type)

        # Bridgeability ratio
        ratio = donor_diameter / inter_locus_nm if inter_locus_nm > 0 else float("inf")
        feasible = ratio >= 0.5

        # Build explanation
        explanation = _build_bridgeability_explanation(
            genomic_distance_bp=genomic_distance_bp,
            donor_size_bp=donor_size_bp,
            donor_type=donor_type,
            inter_locus_nm=inter_locus_nm,
            inter_locus_um=inter_locus_nm / 1000.0,
            donor_rms_nm=donor_rms,
            donor_rg_nm=donor_rg,
            donor_diameter_nm=donor_diameter,
            ratio=ratio,
            feasible=feasible,
            distance_pred=distance_pred,
        )

        return BridgeabilityResult(
            feasible=feasible,
            inter_locus_distance_nm=inter_locus_nm,
            donor_coil_diameter_nm=donor_diameter,
            bridgeability_ratio=ratio,
            explanation=explanation,
        )

    def compare_distances_table(
        self,
        distances_bp: Optional[List[int]] = None,
        chromatin_state: str = "euchromatin",
        donor_sizes_bp: Optional[List[Tuple[int, str]]] = None,
    ) -> str:
        """Generate a formatted comparison table of 3D distances and donor sizes.

        This produces a human-readable table showing how 3D distance scales
        with genomic distance, alongside donor coil sizes for context. It
        helps the user quickly see the scale mismatch between inter-locus
        distances and donor template sizes.

        Parameters
        ----------
        distances_bp : list of int, optional
            Genomic distances to tabulate. Defaults to a range from 1 kb
            to 100 Mb covering typical experimental scenarios.
        chromatin_state : str
            Chromatin state for predictions.
        donor_sizes_bp : list of (int, str) tuples, optional
            Donor templates to include for comparison. Each tuple is
            (size_bp, type_str). Defaults to common donor types.

        Returns
        -------
        str
            Formatted ASCII table.

        Example
        -------
        >>> predictor = ChromatinDistancePredictor()
        >>> print(predictor.compare_distances_table())
        """
        if distances_bp is None:
            distances_bp = [
                1_000,
                10_000,
                100_000,
                500_000,
                1_000_000,
                2_000_000,
                5_000_000,
                10_000_000,
                50_000_000,
                100_000_000,
            ]

        if donor_sizes_bp is None:
            donor_sizes_bp = [
                (200, "ssDNA"),
                (1000, "cssDNA"),
                (3000, "cssDNA"),
                (5000, "dsDNA"),
            ]

        # Header
        lines = []
        lines.append("=" * 95)
        lines.append(
            "3D CHROMATIN DISTANCE vs. GENOMIC DISTANCE"
        )
        lines.append(
            f"  Chromatin state: {chromatin_state} | Model: confined polymer"
        )
        lines.append("=" * 95)
        lines.append("")

        # Table header
        lines.append(
            f"{'Genomic dist':>14s}  {'Mean 3D (nm)':>13s}  "
            f"{'Std (nm)':>10s}  {'5th% (nm)':>10s}  "
            f"{'95th% (nm)':>11s}  {'Mean (um)':>10s}"
        )
        lines.append("-" * 95)

        for dist_bp in distances_bp:
            pred = self.predict_3d_distance(
                dist_bp, chromatin_state, model="confined"
            )
            # Format genomic distance nicely
            if dist_bp >= 1_000_000:
                dist_label = f"{dist_bp / 1_000_000:.0f} Mb"
            elif dist_bp >= 1_000:
                dist_label = f"{dist_bp / 1_000:.0f} kb"
            else:
                dist_label = f"{dist_bp} bp"

            lines.append(
                f"{dist_label:>14s}  {pred.mean_3d_distance_nm:>13.1f}  "
                f"{pred.std_nm:>10.1f}  {pred.p5_nm:>10.1f}  "
                f"{pred.p95_nm:>11.1f}  {pred.mean_3d_distance_um:>10.3f}"
            )

        # Donor sizes for comparison
        lines.append("")
        lines.append("-" * 95)
        lines.append(
            "DONOR TEMPLATE COIL SIZES (for comparison)"
        )
        lines.append("-" * 95)
        lines.append(
            f"{'Donor':>20s}  {'RMS size (nm)':>14s}  "
            f"{'Rg (nm)':>9s}  {'Diameter (nm)':>14s}"
        )
        lines.append("-" * 95)

        for size_bp, dtype in donor_sizes_bp:
            rms = donor_rms_size_nm(size_bp, dtype)
            rg = donor_radius_of_gyration_nm(size_bp, dtype)
            diameter = 2.0 * rg

            if size_bp >= 1_000:
                size_label = f"{size_bp / 1_000:.0f} kb {dtype}"
            else:
                size_label = f"{size_bp} bp {dtype}"

            lines.append(
                f"{size_label:>20s}  {rms:>14.1f}  "
                f"{rg:>9.1f}  {diameter:>14.1f}"
            )

        lines.append("")
        lines.append("=" * 95)
        lines.append(
            "NOTE: Donor coil diameters are typically 10-100x smaller than"
        )
        lines.append(
            "inter-locus 3D distances for Mb-scale separations. A single"
        )
        lines.append(
            "donor template CANNOT physically bridge two loci separated by"
        )
        lines.append(
            ">100 kb. Each locus requires its own donor and its own DSB."
        )
        lines.append("=" * 95)

        return "\n".join(lines)


# =============================================================================
# Private Helper: Build the bridgeability explanation string
# =============================================================================

def _build_bridgeability_explanation(
    genomic_distance_bp: int,
    donor_size_bp: int,
    donor_type: str,
    inter_locus_nm: float,
    inter_locus_um: float,
    donor_rms_nm: float,
    donor_rg_nm: float,
    donor_diameter_nm: float,
    ratio: float,
    feasible: bool,
    distance_pred: DistancePrediction,
) -> str:
    """Build a detailed, biology-oriented explanation of bridgeability.

    This is designed for molecular biologists, not physicists. It uses
    analogies and avoids jargon where possible while remaining precise.
    """
    # Format distances for display
    if genomic_distance_bp >= 1_000_000:
        dist_label = f"{genomic_distance_bp / 1_000_000:.1f} Mb"
    elif genomic_distance_bp >= 1_000:
        dist_label = f"{genomic_distance_bp / 1_000:.1f} kb"
    else:
        dist_label = f"{genomic_distance_bp} bp"

    if donor_size_bp >= 1_000:
        donor_label = f"{donor_size_bp / 1_000:.1f} kb"
    else:
        donor_label = f"{donor_size_bp} bp"

    parts = []

    parts.append(
        f"DONOR BRIDGEABILITY ANALYSIS"
    )
    parts.append(f"{'=' * 60}")
    parts.append(f"")
    parts.append(f"Target separation: {dist_label} (genomic / 1D)")
    parts.append(f"Donor template:    {donor_label} {donor_type}")
    parts.append(f"")

    # Inter-locus distance
    parts.append(f"--- Inter-locus 3D Distance ---")
    parts.append(
        f"  Mean 3D distance:   {inter_locus_nm:.0f} nm ({inter_locus_um:.2f} um)"
    )
    parts.append(
        f"  5th percentile:     {distance_pred.p5_nm:.0f} nm  "
        f"(closest 5% of configurations)"
    )
    parts.append(
        f"  95th percentile:    {distance_pred.p95_nm:.0f} nm  "
        f"(most extended configurations)"
    )
    parts.append(f"")

    # Donor size
    parts.append(f"--- Donor Coil Size ---")
    parts.append(
        f"  RMS end-to-end:     {donor_rms_nm:.1f} nm"
    )
    parts.append(
        f"  Radius of gyration: {donor_rg_nm:.1f} nm"
    )
    parts.append(
        f"  Effective diameter:  {donor_diameter_nm:.1f} nm  "
        f"(2 x radius of gyration)"
    )
    parts.append(f"")

    # Ratio and verdict
    parts.append(f"--- Bridgeability Assessment ---")
    parts.append(
        f"  Bridgeability ratio (donor diameter / inter-locus distance): "
        f"{ratio:.4f}"
    )

    if ratio >= 1.0:
        parts.append(
            f"  VERDICT: FEASIBLE - The donor coil is larger than the inter-locus"
        )
        parts.append(
            f"  distance. Bridging is physically plausible."
        )
    elif ratio >= 0.5:
        parts.append(
            f"  VERDICT: MARGINAL - The donor coil is roughly half the"
        )
        parts.append(
            f"  inter-locus distance. Bridging is possible but requires"
        )
        parts.append(
            f"  favorable thermal fluctuations."
        )
    elif ratio >= 0.1:
        parts.append(
            f"  VERDICT: UNLIKELY - The donor coil is {1/ratio:.0f}x smaller"
        )
        parts.append(
            f"  than the inter-locus distance. Bridging would require extremely"
        )
        parts.append(
            f"  rare conformations."
        )
    else:
        parts.append(
            f"  VERDICT: PHYSICALLY IMPOSSIBLE - The donor coil is {1/ratio:.0f}x"
        )
        parts.append(
            f"  smaller than the inter-locus distance."
        )

    parts.append(f"")

    # Biological explanation
    parts.append(f"--- Biological Interpretation ---")

    if not feasible:
        parts.append(
            f"  The {donor_label} {donor_type} donor forms a tiny random coil"
        )
        parts.append(
            f"  (diameter ~{donor_diameter_nm:.0f} nm) that is dwarfed by the"
        )
        parts.append(
            f"  ~{inter_locus_nm:.0f} nm gap between the two target loci."
        )
        parts.append(f"")
        parts.append(
            f"  Think of it this way: the donor is like a marble"
        )
        parts.append(
            f"  ({donor_diameter_nm:.0f} nm) trying to bridge a hallway"
        )
        parts.append(
            f"  ({inter_locus_nm:.0f} nm). No amount of clever design can"
        )
        parts.append(
            f"  make the marble span that distance."
        )
        parts.append(f"")
        parts.append(
            f"  RECOMMENDATION: For multi-site editing at these distances,"
        )
        parts.append(
            f"  you MUST use separate donor templates and separate DSBs at"
        )
        parts.append(
            f"  each target locus. A single donor cannot serve both sites."
        )
        parts.append(
            f"  Consider sequential editing (edit one site, clone, then edit"
        )
        parts.append(
            f"  the second site) to avoid simultaneous DSBs and the associated"
        )
        parts.append(
            f"  translocation risk."
        )
    else:
        parts.append(
            f"  The donor coil size is comparable to the inter-locus distance."
        )
        parts.append(
            f"  Single-donor bridging may be feasible, though homology search"
        )
        parts.append(
            f"  dynamics will still limit efficiency. Consider that the donor"
        )
        parts.append(
            f"  must find BOTH homology arms simultaneously, which compounds"
        )
        parts.append(
            f"  the search problem."
        )

    return "\n".join(parts)
