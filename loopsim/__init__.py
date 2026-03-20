"""
LoopSim -- Cohesin Loop Extrusion Simulator for HDR Homology Search
====================================================================

This module simulates **cohesin-mediated loop extrusion** at DNA double-strand
breaks (DSBs) and predicts the **homology search domain** available to the
RAD51 recombinase filament during homology-directed repair (HDR).

Biological Background
---------------------
When CRISPR-Cas9 (or any nuclease) creates a DSB in a mammalian cell, the
broken chromosome must be repaired. One high-fidelity pathway is **HDR**,
which requires a homologous template. But how does the cell *find* that
template in a nucleus containing 6 billion base pairs?

The landmark study by **Marin-Gonzalez et al., "Cohesin drives chromatin
scanning during the RAD51-mediated homology search," Science, 2025** revealed
the mechanism:

1. **DSB triggers cohesin loading.** After Cas9 cuts, the ATM kinase
   phosphorylates H2AX across megabases of flanking chromatin (creating the
   gamma-H2AX signal). This recruits the NIPBL cohesin-loader complex, which
   deposits **cohesin** (a ring-shaped SMC complex: SMC1/SMC3/RAD21/SA1-2)
   directly at the break site. (See also: Arnould et al., Nature, 2021.)

2. **Cohesin performs loop extrusion.** Once loaded, cohesin acts as a
   molecular motor that reels chromatin through itself. Picture a carabiner
   clipped at the DSB: the DNA fiber is pulled through the ring from both
   sides simultaneously at ~1 kb/sec. This creates an expanding chromatin
   loop anchored at the DSB.

3. **CTCF sites act as roadblocks.** The protein CTCF binds specific DNA
   motifs in an oriented manner (the CTCF motif is asymmetric). Cohesin
   stalls when it encounters a **convergent** CTCF site — one whose
   orientation faces the approaching cohesin motor. This is why mammalian
   genomes are organized into Topologically Associating Domains (TADs)
   bounded by convergent CTCF pairs.

4. **The extruded loop defines the homology search domain.** As chromatin
   is reeled through the DSB anchor point, the RAD51 filament (formed on
   the resected ssDNA at the break) can scan the passing chromatin for
   sequence complementarity. The region of chromatin that passes through
   the anchor during extrusion is the "search domain."

5. **Key insight for genome editing:** For endogenous donors (sister
   chromatid, ectopic homology on the same chromosome), the donor must
   lie within the search domain to be found efficiently. For **exogenous**
   donors (such as cssDNA delivered by electroporation), the donor is free
   in the nucleoplasm and finds the DSB by 3D diffusion — loop extrusion
   is not the rate-limiting step.

What This Module Simulates
--------------------------
- **ChromatinFiber**: A 1D lattice model of the chromatin fiber, with
  each bead representing 1 kb of chromatin. Beads carry CTCF binding
  information (position, orientation, occupancy probability) and
  chromatin state (open euchromatin vs. closed heterochromatin).

- **CohesinExtruder**: A stochastic simulation engine that models the
  loop extrusion process. Cohesin is loaded at the DSB and extrudes
  bidirectionally until stalled by CTCF, until reaching a boundary,
  or until a maximum time elapses.

- **HomologySearchPredictor**: Translates extrusion simulation results
  into predictions about which genomic regions are scanned by RAD51
  and the probability that a cis-located donor template falls within
  the search domain.

- **LoopSimulator**: The main user-facing interface that chains
  together fiber setup, cohesin simulation, and result analysis.
  Supports single-DSB and dual-DSB (paired nickase / dual-cut)
  scenarios.

- **Visualization**: Publication-quality plots including extrusion
  kymographs (space-time diagrams), scan probability tracks, and
  comparative dual-DSB domain plots.

Quick Start
-----------
>>> from crisprarchitect.loopsim import LoopSimulator
>>> sim = LoopSimulator(region_start_bp=0, region_end_bp=5_000_000)
>>> sim.setup_fiber(ctcf_spacing=120_000, ctcf_strength=0.8)
>>> sim.add_dsb(2_500_000)
>>> results = sim.run(n_simulations=1000)
>>> sim.plot_search_domain()
>>> sim.summary()

Dependencies
------------
- numpy (numerical simulation, vectorized operations)
- matplotlib (visualization)

References
----------
- Marin-Gonzalez A, et al. "Cohesin drives chromatin scanning during the
  RAD51-mediated homology search." Science, 2025.
- Arnould C, et al. "Loop extrusion as a mechanism for formation of DNA
  damage repair foci." Nature, 2021.
- Davidson IF, et al. "DNA loop extrusion by human cohesin." Science, 2019.
- Fudenberg G, et al. "Formation of chromosomal domains by loop extrusion."
  Cell Reports, 2016.
"""

from .chromatin_fiber import ChromatinFiber, CTCFSite
from .cohesin_extruder import CohesinExtruder, ExtrusionResult
from .homology_search import HomologySearchPredictor, SearchDomainResult
from .simulator import LoopSimulator, LoopSimResults
from .visualize import (
    plot_extrusion_kymograph,
    plot_scan_probability_track,
    plot_dual_dsb_comparison,
)

__all__ = [
    "ChromatinFiber",
    "CTCFSite",
    "CohesinExtruder",
    "ExtrusionResult",
    "HomologySearchPredictor",
    "SearchDomainResult",
    "LoopSimulator",
    "LoopSimResults",
    "plot_extrusion_kymograph",
    "plot_scan_probability_track",
    "plot_dual_dsb_comparison",
]

__version__ = "0.1.0"
