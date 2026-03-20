"""
synthesis.py — Simulation of Template-Directed DNA Synthesis During HDR
=======================================================================

Background biology
------------------
After the RAD51 filament successfully invades a homologous donor template
(forming a displacement loop, or D-loop), the 3' end of the invading
strand serves as a primer for DNA synthesis.  **DNA Polymerase delta**
(Pol δ), together with its processivity clamp PCNA and the clamp loader
RFC, extends the invading 3' end by copying the donor template sequence
into the broken chromosome.

The length of this synthesis tract is the **gene conversion tract length**
— the stretch of DNA in the repaired chromosome that now carries
information copied from the donor.  This is the single most important
parameter for HDR-based genome editing, because it determines how far
from the cut site donor-encoded edits (SNPs, insertions, tags, etc.)
can be reliably incorporated.

Two competing HDR sub-pathways determine tract length
-----------------------------------------------------

**DSBR (double-strand break repair / double Holliday-junction pathway)**
    Both broken ends invade the donor, forming two Holliday junctions.
    Resolution can produce crossovers.  This pathway is suppressed in
    mitotic mammalian cells by the BLM helicase "dissolvasome."

**SDSA (Synthesis-Dependent Strand Annealing)** — the DOMINANT mitotic pathway
    Only one end invades.  After a limited amount of DNA synthesis, the
    D-loop is *collapsed* (disrupted): the newly synthesized strand is
    displaced from the donor template and anneals back to the second
    resected end on the other side of the break.  This pathway:
    - Produces only gene conversions (no crossovers)
    - Has a characteristic *geometric / exponential-like* tract-length
      distribution, because at each base pair of synthesis there is a
      fixed probability that the D-loop collapses and synthesis terminates.
    (Nassif et al., Mol. Cell. Biol. 1994; McMahill et al., PLoS Genet. 2007)

We model synthesis as an SDSA process because:
1. SDSA dominates in mitotically dividing mammalian cells.
2. The geometric tract-length distribution matches experimental observations
   (Elliott et al., Mol. Cell Biol. 1998; Kan et al., Mol. Cell 2017).

Modelling approach
------------------
At each base pair of synthesis, there is a probability *p* that the D-loop
collapses (SDSA displacement event).  The tract length is therefore
distributed as a **Geometric random variable** with success probability *p*.

    p = SDSA_DISPLACEMENT_PROB_PER_BP  (from constants.py, default 0.002)

This gives:
    mean tract length = 1/p = 500 bp
    median tract length = ln(2)/p ≈ 347 bp
    P(tract ≥ 1000 bp) = (1-p)^1000 ≈ 13.5 %

Adjustments for donor topology and cut geometry
-------------------------------------------------
* **Circular ssDNA (cssDNA) donors**: The circular topology may stabilise
  the D-loop by preventing end-fraying of the donor.  Topological
  constraint means the displaced strand has no free end to "escape."
  We model this as a 20 % reduction in displacement probability, which
  increases mean tract length by ~25 %.
  (Iyer et al., CRISPR J. 2022; inferred from improved HDR rates)

* **Staggered cuts**: More extensive resection (from the 5' overhang boost)
  produces longer RAD51 filaments, which stabilise the D-loop.  A more
  stable D-loop is less likely to collapse at each step.  We model this
  as a 15 % reduction in displacement probability.
  (Chauhan et al., PNAS 2023; inferred from vCas9-mediated HDR improvement)

Implementation note
-------------------
Rather than simulating one bp at a time (slow), we use NumPy's geometric
distribution to draw the number of "failures" (bp synthesised) before the
first "success" (D-loop collapse).  This is O(1) per simulation, enabling
millions of samples per second.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal, Union

import numpy as np

# ---------------------------------------------------------------------------
# Import biological constants
# ---------------------------------------------------------------------------
import sys, os
_PACKAGE_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if _PACKAGE_ROOT not in sys.path:
    sys.path.insert(0, _PACKAGE_ROOT)

from utils.constants import (
    SDSA_DISPLACEMENT_PROB_PER_BP,
    CONVERSION_TRACT_MIN_BP,
    CONVERSION_TRACT_MAX_BP,
)


# ---------------------------------------------------------------------------
# Module-level constants
# ---------------------------------------------------------------------------

# Fractional reduction in D-loop collapse probability for circular donors.
# Circular topology prevents donor-strand fraying, making the D-loop more
# kinetically stable.
# Source: inferred from improved HDR efficiency with cssDNA donors;
#         Iyer et al., CRISPR J. 2022.
_CIRCULAR_DLOOP_STABILITY_BOOST: float = 0.20  # 20 % reduction in p

# Fractional reduction in D-loop collapse probability for staggered cuts.
# Longer resection -> longer RAD51 filament -> more base-pairs of
# heteroduplex -> greater D-loop thermodynamic stability.
# Source: inferred from vCas9 (staggered cut) HDR enhancement;
#         Chauhan et al., PNAS 2023.
_STAGGER_DLOOP_STABILITY_BOOST: float = 0.15  # 15 % reduction in p


@dataclass
class SynthesisResult:
    """Container for synthesis simulation output.

    Attributes
    ----------
    tract_lengths_bp : np.ndarray, shape (n_simulations,)
        The gene conversion tract length (bp) for each simulation.  This is
        the number of base pairs of donor sequence incorporated into the
        repaired chromosome, measured from the point of strand invasion
        (approximately the cut site for a centrally placed homology arm).
    effective_displacement_prob : float
        The per-bp D-loop collapse probability used for this simulation,
        after adjustments for donor topology and cut geometry.
    donor_topology : str
        The donor template topology that was simulated.
    is_staggered : bool
        Whether a staggered cut was modelled.
    """

    tract_lengths_bp: np.ndarray
    effective_displacement_prob: float
    donor_topology: str
    is_staggered: bool


class SynthesisSimulator:
    """Monte Carlo simulator of template-directed DNA synthesis during SDSA.

    After a RAD51-coated 3' ssDNA tail invades the donor template, DNA
    Polymerase delta synthesises new DNA by copying the donor sequence.
    The length of this synthesis determines the **gene conversion tract
    length** — the portion of the repaired chromosome that carries donor-
    derived sequence.

    In the SDSA pathway (dominant in mitotic mammalian cells), the D-loop
    is inherently unstable.  At every base pair of synthesis, there is a
    probability *p* that helicases (e.g., RTEL1, BLM) dismantle the D-loop,
    terminating synthesis.  This creates a geometric distribution of tract
    lengths.

    Factors that modulate D-loop stability (and hence tract length):
    - **Donor topology**: circular ssDNA stabilises the D-loop (−20 % *p*)
    - **Cut stagger**: staggered cuts produce more resection and a longer
      RAD51 filament, stabilising the D-loop (−15 % *p*)
    - (Future extensions could add: Pol δ processivity factors, chromatin
      context, RAD51 mutants, etc.)

    Parameters
    ----------
    seed : int or None
        Random seed for reproducibility.

    Examples
    --------
    >>> synth = SynthesisSimulator(seed=42)
    >>> result = synth.simulate_synthesis(
    ...     donor_topology="circular_ssDNA",
    ...     is_staggered=True,
    ...     n_simulations=10_000,
    ... )
    >>> print(f"Median tract: {np.median(result.tract_lengths_bp):.0f} bp")
    """

    def __init__(self, seed=None) -> None:
        self._rng = np.random.default_rng(seed)

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def simulate_synthesis(
        self,
        donor_topology: Literal[
            "circular_ssDNA", "linear_ssDNA", "linear_dsDNA", "plasmid"
        ] = "linear_dsDNA",
        is_staggered: bool = False,
        n_simulations: int = 10_000,
    ) -> SynthesisResult:
        """Simulate SDSA-mediated DNA synthesis for *n_simulations* cells.

        Each simulation independently draws a tract length from a geometric
        distribution whose parameter *p* (probability of D-loop collapse
        per bp synthesised) is adjusted for the experimental configuration.

        The geometric distribution naturally models the SDSA process:
        - At bp 1: probability *p* of stopping, (1−p) of continuing.
        - At bp 2: probability *p* of stopping, ...
        - Tract length L has P(L = k) = (1−p)^(k−1) * p.
        - Mean = 1/p, Variance = (1−p)/p².

        Parameters
        ----------
        donor_topology : str
            Shape and chemistry of the donor template:

            * ``"circular_ssDNA"`` — a circularised single-stranded DNA
              donor (cssDNA).  The circular topology prevents end-fraying,
              producing a more stable D-loop.  This reduces the per-bp
              displacement probability by 20 %.

            * ``"linear_ssDNA"`` — a long single-stranded oligonucleotide
              (lssDNA) or long ssDNA donor.  No topology-mediated
              stabilisation.

            * ``"linear_dsDNA"`` — a double-stranded DNA donor (PCR product
              or restriction fragment).  Standard displacement probability.

            * ``"plasmid"`` — a circular dsDNA plasmid donor.  Treated
              identically to linear_dsDNA for D-loop stability (the
              circular dsDNA topology does not stabilise the D-loop in
              the same way as circular ssDNA, because strand invasion
              must unwind both strands).

        is_staggered : bool
            Whether the DSB was produced by a staggered-cutting nuclease.
            If True, the per-bp displacement probability is reduced by 15 %
            (more resection → longer filament → stabler D-loop).

        n_simulations : int
            Number of independent Monte Carlo trials.

        Returns
        -------
        SynthesisResult
            Dataclass containing the array of tract lengths and metadata.

        Notes
        -----
        The ``numpy.random.Generator.geometric`` function returns values
        in {1, 2, 3, ...}, which conveniently maps to "the D-loop collapsed
        after synthesising k bp."
        """
        # ------------------------------------------------------------------
        # Step 1: Compute the effective displacement probability
        # ------------------------------------------------------------------
        # Start with the baseline from constants.py:
        #   p = 0.002 per bp => mean tract = 1/0.002 = 500 bp
        p = SDSA_DISPLACEMENT_PROB_PER_BP

        # Circular ssDNA stabilises the D-loop.
        # Rationale: In a circular donor, the displaced non-template strand
        # has no free 5' or 3' end.  This topological constraint opposes
        # spontaneous D-loop reversal, because rewinding the donor requires
        # superhelical stress rather than simple strand reannealing from a
        # free end.  Net effect: D-loop persists longer → longer tracts.
        if donor_topology == "circular_ssDNA":
            p = p * (1.0 - _CIRCULAR_DLOOP_STABILITY_BOOST)
            # p goes from 0.002 → 0.0016, mean tract 500 → 625 bp

        # Staggered cuts produce longer resection and hence longer RAD51
        # filaments.  A longer filament means more base-pairs of
        # heteroduplex DNA in the D-loop, increasing its thermodynamic
        # stability (more hydrogen bonds to break for reversal).
        if is_staggered:
            p = p * (1.0 - _STAGGER_DLOOP_STABILITY_BOOST)
            # Further reduction: e.g., 0.0016 → 0.00136 for cssDNA +
            # stagger, giving mean tract ~735 bp

        # Ensure p is in a valid range (must be > 0 and <= 1).
        p = np.clip(p, 1e-6, 1.0)

        # ------------------------------------------------------------------
        # Step 2: Draw tract lengths from a geometric distribution
        # ------------------------------------------------------------------
        # numpy.geometric(p, size) returns values in {1, 2, 3, ...}
        # where the PMF is P(X=k) = (1-p)^(k-1) * p.
        # k represents "the D-loop collapsed after synthesising k bp."
        tract_lengths = self._rng.geometric(p, size=n_simulations)

        # ------------------------------------------------------------------
        # Step 3: Clip to biologically plausible bounds
        # ------------------------------------------------------------------
        # Minimum: CONVERSION_TRACT_MIN_BP (50 bp from constants.py).
        #   Very short tracts (< 50 bp) are unlikely to produce detectable
        #   gene conversion because mismatch repair may erase the
        #   heteroduplex before the next round of replication.
        #   (Mitchel et al., Genetics 2010)
        # Maximum: CONVERSION_TRACT_MAX_BP (5000 bp from constants.py).
        #   Extremely long tracts are rare because RTEL1 and BLM helicases
        #   actively disassemble D-loops.  Tracts beyond 5 kb are
        #   essentially never observed in mammalian mitotic cells.
        #   (Kan et al., Mol Cell 2017)
        tract_lengths = np.clip(
            tract_lengths,
            CONVERSION_TRACT_MIN_BP,
            CONVERSION_TRACT_MAX_BP,
        ).astype(np.float64)

        return SynthesisResult(
            tract_lengths_bp=tract_lengths,
            effective_displacement_prob=float(p),
            donor_topology=donor_topology,
            is_staggered=is_staggered,
        )
