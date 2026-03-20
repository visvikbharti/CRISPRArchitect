"""
ChromBridge: 3D Chromatin Distance & Translocation Risk Module
==============================================================

ChromBridge connects the 1D world of genomic coordinates with the 3D reality
of chromatin organization inside the nucleus. This matters for CRISPR editing
because:

    1. **Donor bridging**: When you want a single donor template to mediate
       editing at two distant sites (e.g., exon 40 and exon 60), the donor
       must physically span the 3D distance between those loci. ChromBridge
       calculates that distance and shows whether bridging is feasible.

    2. **Translocation risk**: When you create two simultaneous DSBs, broken
       chromosome ends can be mis-joined, creating translocations, deletions,
       or inversions. The probability depends on the 3D spatial proximity of
       the breaks. ChromBridge quantifies this risk.

    3. **TAD awareness**: Topologically Associating Domains (TADs) partition
       chromosomes into self-interacting neighborhoods. Two loci within the
       same TAD interact far more frequently than loci in different TADs,
       which affects both homology search efficiency and translocation risk.

The physics under the hood:
    Chromosomes are long polymer chains (chromatin fibers). We model their
    3D spatial statistics using polymer physics — the same mathematics used
    to describe rubber, spaghetti, and any other flexible chain. The key
    insight is that the 3D distance between two points on a polymer grows
    much more slowly than their 1D separation along the chain, because the
    chain folds back on itself many times.

Modules
-------
- ``polymer_model``: Gaussian chain, worm-like chain, and confined polymer
  models for predicting mean-squared end-to-end distances of chromatin.
- ``distance``: Converts polymer predictions into practical 3D distance
  estimates with uncertainty, and evaluates donor bridgeability.
- ``translocation``: Translocation, deletion, and inversion risk scoring
  when two DSBs are created simultaneously.
- ``tad_analysis``: Estimates TAD membership and boundary effects on
  inter-locus contact frequency.

Quick Start
-----------
>>> from crisprarchitect.chrombridge import ChromatinDistancePredictor
>>> predictor = ChromatinDistancePredictor()
>>> result = predictor.predict_3d_distance(2_000_000, chromatin_state="euchromatin")
>>> print(f"Mean 3D distance: {result.mean_3d_distance_nm:.0f} nm "
...       f"({result.mean_3d_distance_um:.2f} um)")

>>> from crisprarchitect.chrombridge import TranslocationRiskPredictor
>>> risk_pred = TranslocationRiskPredictor()
>>> risk = risk_pred.estimate_risk(50_000_000, 52_000_000, same_chromosome=True)
>>> print(risk_pred.risk_summary(risk))
"""

# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

from .polymer_model import (
    GaussianChainModel,
    WormLikeChainModel,
    ConfinedPolymerModel,
)

from .distance import (
    ChromatinDistancePredictor,
    DistancePrediction,
    BridgeabilityResult,
)

from .translocation import (
    TranslocationRiskPredictor,
    TranslocationRisk,
    SafetyReport,
)

from .tad_analysis import (
    TADAnalyzer,
    TADInfo,
)

__all__ = [
    # Polymer models
    "GaussianChainModel",
    "WormLikeChainModel",
    "ConfinedPolymerModel",
    # Distance prediction
    "ChromatinDistancePredictor",
    "DistancePrediction",
    "BridgeabilityResult",
    # Translocation risk
    "TranslocationRiskPredictor",
    "TranslocationRisk",
    "SafetyReport",
    # TAD analysis
    "TADAnalyzer",
    "TADInfo",
]
