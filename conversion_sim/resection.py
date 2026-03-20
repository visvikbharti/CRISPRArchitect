"""
resection.py — Simulation of DNA End Resection at a CRISPR-Induced DSB
=======================================================================

Background biology
------------------
When CRISPR-Cas nuclease creates a double-strand break (DSB), the broken
DNA ends must be processed before any repair pathway can operate.  The
processing step relevant to HDR is called **5'→3' end resection**: endo-
and exo-nucleases chew back the 5'-terminated strand on each side of the
break, exposing long 3' single-stranded DNA (ssDNA) tails.  These 3' tails
are the substrate for RAD51 filament formation and subsequent homology
search.

Resection proceeds in two kinetically distinct phases:

1. **Short-range (initiation) resection — MRN/CtIP**
   The MRE11–RAD50–NBS1 (MRN) complex, together with the CtIP co-factor,
   performs the initial incision ~50–300 bp from the break end.  MRE11 has
   3'→5' exonuclease activity, so the actual mechanism is an *endonucleolytic*
   nick internal to the 5' strand, followed by 3'→5' digestion *back toward*
   the break.  The net result is a short 3' overhang of ~100–300 bp.
   (Symington, Annu. Rev. Genet. 2011; Cejka, Annu. Rev. Genet. 2015)

2. **Long-range (processive) resection — EXO1 or BLM-DNA2**
   After the short-range step clears the break end, one of two partially
   redundant enzyme systems extends the resection further:
   - **EXO1**: a 5'→3' exonuclease that processively degrades the 5' strand.
   - **BLM helicase + DNA2 nuclease**: BLM unwinds dsDNA ahead of the
     replication fork-like structure, and DNA2 cleaves the displaced 5'
     strand.
   Long-range resection can extend 1–5 kb or more, with a heavy right tail
   (occasional tracts of 10 kb+).  This distribution is well-modelled by a
   log-normal: most tracts are moderate, but a non-negligible fraction are
   very long.
   (Zhou et al., Mol. Cell 2014; Zhu et al., Cell 2008)

Effect of staggered cuts
-------------------------
Blunt-cutting nucleases (e.g., SpCas9) produce flush DSB ends where the
5' and 3' strands terminate at the same position.  In contrast, staggered-
cutting nucleases (e.g., Cas12a, enFnCas9, vCas9) produce a DSB with a
pre-existing **5' overhang** — the 5' strand extends beyond the 3' strand
by several nucleotides.

This 5' overhang has two consequences for resection:

* The MRN/CtIP short-range step is effectively **pre-done** on the overhang
  side: the 3' strand is already recessed by the overhang length.  We model
  this as adding the stagger length to the total resection on that side.

* The exposed 5' single-stranded flap is a preferred substrate for EXO1.
  Studies show that 5' overhangs stimulate EXO1 loading and processivity
  in vitro (Cannavo et al., Mol. Cell 2013).  We model this as a **20 %
  boost** to long-range resection on the overhang side.

The opposite side of a staggered break (the side with a 3' overhang/
recessed 5' end) is resected normally.

Modelling approach
------------------
We draw *n* independent samples for each side of the break:

    total_resection = short_range + long_range  [+ stagger adjustments]

* ``short_range`` ~ Normal(μ=200, σ=80), clipped to [50, 500].
* ``long_range``  ~ LogNormal(μ_log, σ_log), fitted so median ≈ 1500 bp,
  with σ_log = 0.8 to produce a heavy right tail up to ~15 kb.

All computations are fully vectorised with NumPy (no per-simulation loops).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Literal

import numpy as np

# ---------------------------------------------------------------------------
# Import biological constants from the shared constants module.
# ---------------------------------------------------------------------------
import sys, os
_PACKAGE_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if _PACKAGE_ROOT not in sys.path:
    sys.path.insert(0, _PACKAGE_ROOT)

from utils.constants import (
    SHORT_RESECTION_MEAN_BP,
    SHORT_RESECTION_STD_BP,
    LONG_RESECTION_MIN_BP,
    LONG_RESECTION_MAX_BP,
)


# ---------------------------------------------------------------------------
# Module-level constants derived from the literature
# ---------------------------------------------------------------------------

# We fit a LogNormal so that the *median* of the long-range resection
# distribution is ~1500 bp.  For a LogNormal(mu, sigma), the median is
# exp(mu).  We choose sigma_log = 0.8 to produce a fat right tail (some
# cells will have 5–15 kb resection tracts).
#
# median = exp(mu_log) => mu_log = ln(1500) ≈ 7.313
# Source: parameterised to match tract-length observations from
#         Zhou et al., Mol Cell 2014 and Symington 2011 review.
_LONG_RESECTION_MEDIAN_BP: float = 1500.0
_LONG_RESECTION_MU_LOG: float = float(np.log(_LONG_RESECTION_MEDIAN_BP))
_LONG_RESECTION_SIGMA_LOG: float = 0.8  # unitless; controls tail weight

# Stagger-induced boost to long-range resection (fractional).
# A 5' overhang stimulates EXO1 loading, increasing processive resection
# by approximately 20 %.
# Source: inferred from Cannavo et al., Mol Cell 2013 (in vitro EXO1
#         stimulation on 5'-flap substrates) and Chauhan et al., PNAS 2023
#         (enhanced HDR with vCas9 staggered cuts).
_STAGGER_LONG_RANGE_BOOST: float = 0.20  # 20 % increase


@dataclass
class ResectionResult:
    """Container for the output of a resection simulation.

    Attributes
    ----------
    left_bp : np.ndarray, shape (n_simulations,)
        Resection length (bp) on the **left** (PAM-proximal / 5' overhang)
        side of the DSB for each simulation.  "Left" and "right" are
        arbitrary labels; what matters biologically is that both sides are
        resected independently.
    right_bp : np.ndarray, shape (n_simulations,)
        Resection length (bp) on the **right** (PAM-distal) side.
    cut_type : str
        The cut geometry that was simulated ("blunt" or "staggered_5prime").
    stagger_bp : int
        The length of the 5' overhang (0 for blunt cuts).
    """

    left_bp: np.ndarray
    right_bp: np.ndarray
    cut_type: str
    stagger_bp: int


class ResectionSimulator:
    """Monte Carlo simulator of 5'→3' DNA end resection at a DSB.

    This class models the two-phase resection process that generates 3'
    single-stranded DNA overhangs after a CRISPR-induced double-strand
    break.  The overhangs are the obligate substrate for RAD51 filament
    formation (see ``filament.py``) and therefore for HDR.

    The simulator supports two cut geometries:

    * **blunt** — both strands cleaved at the same position (e.g., SpCas9).
      Both sides of the break undergo identical, independent resection.

    * **staggered_5prime** — the non-target strand is cleaved further from
      the PAM than the target strand, producing a 5' overhang of
      ``stagger_bp`` nucleotides (e.g., Cas12a with ~5 nt overhang,
      enFnCas9 with ~3 nt, vCas9 with ~6 nt).  On the overhang side the
      pre-existing single-stranded flap adds directly to the resection
      length and stimulates EXO1 processivity (+20 %).

    Parameters
    ----------
    seed : int or None
        Random seed for reproducibility.  Pass ``None`` for a random seed.

    Examples
    --------
    >>> sim = ResectionSimulator(seed=42)
    >>> result = sim.simulate(
    ...     cut_type="staggered_5prime",
    ...     stagger_bp=5,
    ...     n_simulations=10_000,
    ... )
    >>> print(f"Median left resection: {np.median(result.left_bp):.0f} bp")
    """

    def __init__(self, seed=None) -> None:
        self._rng = np.random.default_rng(seed)

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def simulate(
        self,
        cut_type: Literal["blunt", "staggered_5prime"] = "blunt",
        stagger_bp: int = 0,
        n_simulations: int = 10_000,
    ) -> ResectionResult:
        """Run *n_simulations* independent resection events.

        Each simulation draws short-range and long-range resection lengths
        for **both** sides of the DSB.  The two sides are treated as
        independent random variables because, biologically, MRN/CtIP and
        EXO1 operate on each broken end separately.

        Parameters
        ----------
        cut_type : {"blunt", "staggered_5prime"}
            Geometry of the DSB produced by the nuclease.

            * ``"blunt"`` — flush ends, e.g. SpCas9.
            * ``"staggered_5prime"`` — 5' overhang on one side, e.g. Cas12a,
              enFnCas9, vCas9.

        stagger_bp : int
            Length of the 5' overhang in base pairs.  Ignored when
            ``cut_type="blunt"``.  Typical values: 3 (enFnCas9), 5 (Cas12a),
            6 (vCas9).

        n_simulations : int
            Number of independent Monte Carlo trials.

        Returns
        -------
        ResectionResult
            A dataclass containing arrays of resection lengths for both
            sides and metadata about the cut geometry.

        Notes
        -----
        All arithmetic is vectorised.  10 000 simulations complete in < 1 ms
        on a modern laptop.
        """
        # Validate inputs
        if cut_type not in ("blunt", "staggered_5prime"):
            raise ValueError(
                f"cut_type must be 'blunt' or 'staggered_5prime', got {cut_type!r}"
            )
        if cut_type == "staggered_5prime" and stagger_bp <= 0:
            raise ValueError(
                "stagger_bp must be > 0 for staggered_5prime cuts; "
                f"got {stagger_bp}"
            )
        n = n_simulations

        # ------------------------------------------------------------------
        # STEP 1: Short-range resection (MRN/CtIP)
        # ------------------------------------------------------------------
        # Draw from a normal distribution centred at ~200 bp (the typical
        # length of the initial ssDNA tract created by MRN/CtIP).
        # Clip to biologically plausible bounds:
        #   - minimum 50 bp: below this the break end is not sufficiently
        #     processed for long-range nuclease loading.
        #   - maximum 500 bp: CtIP-dependent resection rarely exceeds this.
        # We draw independently for both sides of the break.
        short_left = self._draw_short_range(n)   # shape (n,)
        short_right = self._draw_short_range(n)

        # ------------------------------------------------------------------
        # STEP 2: Long-range resection (EXO1 or BLM-DNA2)
        # ------------------------------------------------------------------
        # Drawn from a log-normal distribution, producing a right-skewed
        # distribution: most cells have moderate resection (1–3 kb) but a
        # tail extends to 10 kb+.  Log-normal is a natural model for a
        # *multiplicative* stochastic process (nuclease encountering various
        # chromatin barriers of variable strength).
        long_left = self._draw_long_range(n)
        long_right = self._draw_long_range(n)

        # ------------------------------------------------------------------
        # STEP 3: Apply stagger effects (if applicable)
        # ------------------------------------------------------------------
        if cut_type == "staggered_5prime":
            # The 5' overhang exists on one side of the break.  By convention
            # we assign it to the LEFT side.  (In reality the overhang side
            # depends on the orientation of the guide RNA relative to the
            # locus, but for tract-length statistics the label is arbitrary.)

            # 3a. The overhang length is *added* to the left-side resection.
            #     Rationale: the 3' end on the overhang side is already
            #     recessed by ``stagger_bp`` nucleotides.  This is
            #     equivalent to a head start for resection.
            #     (The 5' strand protrudes, so the 3' strand is recessed.)
            short_left = short_left + stagger_bp

            # 3b. EXO1 loads more efficiently on a 5' single-stranded flap.
            #     We boost the long-range component on the overhang side by
            #     20 % (see module docstring for source).
            long_left = long_left * (1.0 + _STAGGER_LONG_RANGE_BOOST)

        # ------------------------------------------------------------------
        # STEP 4: Total resection = short + long
        # ------------------------------------------------------------------
        total_left = short_left + long_left    # shape (n,)
        total_right = short_right + long_right

        # Ensure non-negative (should always hold, but defensive)
        total_left = np.maximum(total_left, 0.0)
        total_right = np.maximum(total_right, 0.0)

        return ResectionResult(
            left_bp=total_left,
            right_bp=total_right,
            cut_type=cut_type,
            stagger_bp=stagger_bp if cut_type == "staggered_5prime" else 0,
        )

    # ------------------------------------------------------------------
    # Private helpers
    # ------------------------------------------------------------------

    def _draw_short_range(self, n: int) -> np.ndarray:
        """Sample short-range resection lengths from a clipped Normal.

        Biology: MRN/CtIP makes an initial ssDNA overhang of ~100–300 bp.
        We model this as Normal(μ=200, σ=80) clipped to [50, 500].

        Parameters
        ----------
        n : int
            Number of samples.

        Returns
        -------
        np.ndarray, shape (n,), dtype float64
            Short-range resection lengths in bp.
        """
        samples = self._rng.normal(
            loc=SHORT_RESECTION_MEAN_BP,   # 200 bp (Symington 2011)
            scale=SHORT_RESECTION_STD_BP,  # 80 bp
            size=n,
        )
        # Clip to biologically plausible range.
        # 50 bp: minimum processing needed for EXO1/BLM-DNA2 to load.
        # 500 bp: upper bound of CtIP-mediated resection in vivo.
        return np.clip(samples, 50.0, 500.0)

    def _draw_long_range(self, n: int) -> np.ndarray:
        """Sample long-range resection lengths from a LogNormal.

        Biology: EXO1 or BLM-DNA2 processively resect the 5' strand for
        1–5+ kb.  The distribution is right-skewed because processivity is
        limited by stochastic encounters with nucleosomes, cohesin, and
        other chromatin obstacles.  A log-normal naturally arises when the
        logarithm of the processivity distance is normally distributed —
        i.e., multiplicative noise from sequential barriers.

        We parameterise the log-normal so that:
            median = exp(μ_log) = 1 500 bp
            σ_log  = 0.8  (heavy right tail)

        This gives:
            mean  ≈ 1 500 * exp(0.8² / 2) ≈ 2 100 bp
            95th percentile ≈ exp(μ_log + 1.645 * σ_log) ≈ 5 600 bp

        Samples are clipped to [LONG_RESECTION_MIN_BP, LONG_RESECTION_MAX_BP]
        from constants.py (300 bp to 10 000 bp) to exclude physiologically
        implausible extremes.

        Parameters
        ----------
        n : int
            Number of samples.

        Returns
        -------
        np.ndarray, shape (n,), dtype float64
            Long-range resection lengths in bp.
        """
        # numpy's lognormal takes the mean and std of the *underlying*
        # normal distribution (i.e., the log of the variable).
        samples = self._rng.lognormal(
            mean=_LONG_RESECTION_MU_LOG,      # ln(1500) ≈ 7.31
            sigma=_LONG_RESECTION_SIGMA_LOG,   # 0.8
            size=n,
        )
        return np.clip(samples, LONG_RESECTION_MIN_BP, LONG_RESECTION_MAX_BP)
