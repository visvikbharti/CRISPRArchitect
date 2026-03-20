"""
filament.py — RAD51 Nucleoprotein Filament Formation on Resected ssDNA
======================================================================

Background biology
------------------
After 5'→3' end resection exposes 3' single-stranded DNA (ssDNA) tails
on both sides of a DSB, these tails must be decorated with RAD51 recombinase
before homology search and strand invasion can occur.  The assembly process
is tightly regulated:

1. **RPA coating (immediate)**
   Replication Protein A (RPA) — the major eukaryotic ssDNA-binding protein
   — binds the exposed ssDNA almost instantly.  RPA removes secondary
   structures (hairpins) and protects the ssDNA from nucleolytic degradation.
   In our simulation we treat RPA coating as instantaneous and complete; it
   is not a rate-limiting step.
   (Wold, Annu. Rev. Biochem. 1997)

2. **BRCA2-mediated RAD51 loading**
   The tumour suppressor BRCA2 physically delivers RAD51 monomers to the
   RPA-coated ssDNA.  BRCA2 has multiple BRC repeats, each of which binds
   one RAD51 monomer, and a C-terminal DNA-binding domain that positions
   RAD51 on the ssDNA.  This step displaces RPA.
   (Jensen et al., Nature 2010; Thorslund et al., Nat. Struct. Mol. Biol. 2010)

3. **RAD51 filament growth**
   Once a small nucleus of RAD51 monomers (approximately 5–8) has assembled,
   the filament can grow cooperatively in both directions along the ssDNA.
   Each RAD51 monomer covers **3 nucleotides** of ssDNA and stretches the
   backbone by ~50 % relative to B-form DNA.  The resulting right-handed
   helical filament has ~6 monomers per turn.
   (Ogawa et al., Science 1993; Yu et al., Mol. Cell Biol. 2001)

4. **Stochastic coverage**
   Filament assembly is *not* deterministic.  Some regions of ssDNA remain
   uncoated (covered by RPA rather than RAD51), producing gaps.  Filament
   coverage depends on BRCA2 availability, RAD51 concentration, and
   competition with RPA.  Single-molecule studies show that typical coverage
   is 70–95 % of the available ssDNA under physiological conditions.
   (Modesti et al., Nature 2007; Hilario et al., PNAS 2009)

5. **Minimum filament for strand invasion**
   A *minimum contiguous* RAD51 filament of approximately **5 monomers**
   (~15 nt) is required for a kinetically stable D-loop (strand-invasion
   intermediate).  Shorter filaments can sample homologous sequences but
   dissociate before productive strand exchange.
   (Qi et al., Cell 2015 — 8-nt microhomology sampling; Renkawitz et al.,
   Mol. Cell 2014)

Modelling approach
------------------
Given a resection length *L* (bp), the filament model returns:

* **filament_length_nt**: the number of nucleotides coated by RAD51.
  Drawn as ``coverage_fraction * L``, where ``coverage_fraction`` ~
  Beta(α=8, β=2) rescaled to the range [0.70, 0.95].  The Beta(8, 2)
  shape concentrates most of the probability mass near the upper end
  (mean ≈ 0.80 within the [0.70, 0.95] window), reflecting the fact
  that BRCA2 is an efficient loader under normal conditions.

* **filament_complete**: boolean flag indicating whether the filament
  meets the minimum-length requirement for strand invasion
  (≥ ``MIN_HOMOLOGY_FOR_INVASION_BP`` nucleotides, default 15 nt).

Both scalar and vectorised (array) inputs are supported.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Union

import numpy as np

# ---------------------------------------------------------------------------
# Import biological constants
# ---------------------------------------------------------------------------
import sys, os
_PACKAGE_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if _PACKAGE_ROOT not in sys.path:
    sys.path.insert(0, _PACKAGE_ROOT)

from utils.constants import (
    RAD51_FOOTPRINT_NT,
    RAD51_NUCLEATION_MIN_MONOMERS,
    MIN_HOMOLOGY_FOR_INVASION_BP,
)


# ---------------------------------------------------------------------------
# Module-level constants
# ---------------------------------------------------------------------------

# Parameters for the Beta distribution used to draw filament coverage.
# Beta(alpha=8, beta=2) has mean 0.8 and is left-skewed (most values
# near the top of the range), matching single-molecule observations that
# RAD51 filaments are usually nearly contiguous.
_COVERAGE_ALPHA: float = 8.0
_COVERAGE_BETA: float = 2.0

# Biological limits on filament coverage fraction.
# Minimum 70 %: even under sub-optimal conditions, BRCA2 loads enough
#   RAD51 to coat the majority of the ssDNA.
# Maximum 95 %: small RPA-bound gaps always remain; 100 % coating is
#   essentially unobserved in vivo.
# Source: Modesti et al., Nature 2007 (EM visualisation of filaments)
_COVERAGE_MIN: float = 0.70
_COVERAGE_MAX: float = 0.95


@dataclass
class FilamentResult:
    """Output of filament formation simulation.

    Attributes
    ----------
    filament_length_nt : np.ndarray or float
        The number of nucleotides coated by the RAD51 filament.  This is
        the fraction of the resected ssDNA that is competent for homology
        search.
    filament_complete : np.ndarray (bool) or bool
        True if the filament is long enough (≥ MIN_HOMOLOGY_FOR_INVASION_BP)
        to support stable strand invasion.
    coverage_fraction : np.ndarray or float
        The fraction of resected ssDNA coated by RAD51 (between
        _COVERAGE_MIN and _COVERAGE_MAX).
    n_monomers : np.ndarray or float
        The number of RAD51 monomers in the filament.  Each monomer covers
        RAD51_FOOTPRINT_NT (= 3) nucleotides.
    """

    filament_length_nt: Union[np.ndarray, float]
    filament_complete: Union[np.ndarray, bool]
    coverage_fraction: Union[np.ndarray, float]
    n_monomers: Union[np.ndarray, float]


class FilamentModel:
    """Stochastic model of RAD51 nucleoprotein filament formation.

    After DNA end resection, the exposed 3' ssDNA tail must be coated with
    RAD51 recombinase to form a *presynaptic filament*.  This filament is
    the molecular machine that performs the homology search: it samples
    dsDNA sequences throughout the nucleus until a sufficiently homologous
    donor template is found, then catalyses strand invasion (D-loop
    formation).

    **Why is filament coverage stochastic?**
    In reality, filament assembly is a kinetic competition between:
      - BRCA2-mediated RAD51 loading (promotes filament)
      - RPA re-binding (inhibits filament, fills gaps)
      - RAD51 ATP hydrolysis (causes filament disassembly from the rear)
    The net result is that each cell ends up with a different coverage
    fraction, which we sample from a Beta distribution.

    Parameters
    ----------
    seed : int or None
        Random seed for the internal NumPy generator.

    Examples
    --------
    >>> model = FilamentModel(seed=42)
    >>> result = model.form_filament(resection_length_bp=np.array([1800, 900]))
    >>> result.filament_complete
    array([ True,  True])
    """

    def __init__(self, seed=None) -> None:
        self._rng = np.random.default_rng(seed)

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def form_filament(
        self,
        resection_length_bp: Union[np.ndarray, float],
    ) -> FilamentResult:
        """Simulate RAD51 filament assembly on resected ssDNA.

        This method takes the total resection length (in bp, which equals
        the ssDNA length in nt for one strand) and returns the length of
        RAD51 filament that forms, plus a flag indicating whether the
        filament is long enough for strand invasion.

        The biological sequence of events being modelled:

        1. RPA rapidly coats the entire ssDNA tail (instantaneous here).
        2. BRCA2 delivers RAD51 monomers, displacing RPA.
        3. RAD51 monomers assemble cooperatively into a helical filament.
        4. The filament covers a stochastic fraction (70–95 %) of the
           available ssDNA.
        5. If the contiguous filament is shorter than ~15 nt (5 monomers),
           the filament is too unstable for productive strand invasion.

        Parameters
        ----------
        resection_length_bp : array-like or float
            Length of the 3' ssDNA overhang produced by resection, in bp
            (= nt for ssDNA).  Can be a single value or an array of values
            from a batch of Monte Carlo simulations.

        Returns
        -------
        FilamentResult
            Dataclass with filament_length_nt, filament_complete,
            coverage_fraction, and n_monomers.

        Notes
        -----
        We treat *resection_length_bp* as equivalent to ssDNA length in
        nucleotides, because each resected base pair exposes one nucleotide
        of the 3' strand.  (The complementary 5' strand has been degraded.)
        """
        # Ensure we are working with a NumPy array for vectorised math.
        resection = np.asarray(resection_length_bp, dtype=np.float64)
        scalar_input = resection.ndim == 0
        resection = np.atleast_1d(resection)

        n = resection.shape[0]

        # ------------------------------------------------------------------
        # Draw coverage fractions from a Beta distribution
        # ------------------------------------------------------------------
        # Raw Beta(8, 2) produces values in [0, 1] with mean 0.8.
        # We rescale to the biologically meaningful range [0.70, 0.95].
        raw_beta = self._rng.beta(_COVERAGE_ALPHA, _COVERAGE_BETA, size=n)

        # Linear rescaling: raw in [0,1] -> coverage in [_MIN, _MAX]
        coverage = _COVERAGE_MIN + raw_beta * (_COVERAGE_MAX - _COVERAGE_MIN)

        # ------------------------------------------------------------------
        # Compute filament length
        # ------------------------------------------------------------------
        # filament_length_nt = coverage * resection_length (in nt)
        # Since resection removes one strand, the exposed ssDNA length in
        # nucleotides equals the resection length in base pairs.
        filament_nt = coverage * resection  # shape (n,)

        # ------------------------------------------------------------------
        # Compute number of RAD51 monomers
        # ------------------------------------------------------------------
        # Each RAD51 monomer binds and stretches 3 nt of ssDNA.
        # (Ogawa et al., Science 1993)
        n_monomers = filament_nt / RAD51_FOOTPRINT_NT  # 3 nt per monomer

        # ------------------------------------------------------------------
        # Check minimum filament length for strand invasion
        # ------------------------------------------------------------------
        # Stable D-loop formation requires a contiguous RAD51 filament of
        # at least ~5 monomers (~15 nt).  Below this threshold the
        # filament-DNA complex is kinetically unstable: it can transiently
        # sample homologous sequences but cannot sustain the open D-loop
        # long enough for DNA synthesis to initiate.
        # (Qi et al., Cell 2015)
        filament_ok = filament_nt >= MIN_HOMOLOGY_FOR_INVASION_BP  # bool array

        # ------------------------------------------------------------------
        # Return results (unwrap scalar if input was scalar)
        # ------------------------------------------------------------------
        if scalar_input:
            return FilamentResult(
                filament_length_nt=float(filament_nt[0]),
                filament_complete=bool(filament_ok[0]),
                coverage_fraction=float(coverage[0]),
                n_monomers=float(n_monomers[0]),
            )
        return FilamentResult(
            filament_length_nt=filament_nt,
            filament_complete=filament_ok,
            coverage_fraction=coverage,
            n_monomers=n_monomers,
        )
