"""
models.py — Pre-Calibrated Parameter Sets for Common Experimental Configurations
=================================================================================

This module provides ready-to-use parameter dictionaries for the
:class:`~conversion_sim.simulator.ConversionSimulator`.  Each dictionary
encodes a specific experimental configuration (nuclease + donor + cell type)
that has been calibrated against published HDR data.

Using a preset
--------------
>>> from conversion_sim import ConversionSimulator
>>> from conversion_sim.models import PRESET_CONFIGS
>>>
>>> # Look up a preset
>>> params = PRESET_CONFIGS["default_enfncas9_cssdna_ipsc"]
>>>
>>> # Feed it directly into the simulator
>>> sim = ConversionSimulator(**params)
>>> results = sim.run()
>>> sim.summary()

Available presets
-----------------

``default_spcas9_cssdna_ipsc``
    Standard SpCas9 (blunt cut) with a circular single-stranded DNA (cssDNA)
    donor template in human induced pluripotent stem cells (iPSCs).

    This is the "workhorse" configuration for large knock-in experiments in
    iPSCs.  SpCas9 produces a predominantly blunt DSB at position -3 from
    the PAM (Jinek et al., Science 2012).  The cssDNA donor is circularised
    single-stranded DNA, which is resistant to intracellular exonucleases
    and achieves ~3x higher HDR than linear dsDNA (Iyer et al., CRISPR J.
    2022).  Homology arms of 300 nt each are optimal for cssDNA.

    iPSCs are challenging for HDR because:
    - Only ~35 % of cells are in S/G2 (HDR-competent phase)
    - Strong p53-mediated apoptotic response to DSBs
    - Baseline HDR efficiency is low (~8 %)
    But cssDNA donors partially compensate with their enhanced stability.

``default_enfncas9_cssdna_ipsc``
    Engineered *Francisella novicida* Cas9 (enFnCas9) paired with cssDNA
    in iPSCs.

    enFnCas9 is a rationally engineered variant from the Bhatt/Chakraborty
    lab with two key advantages over SpCas9:
    1. **Staggered cut**: produces a ~3 bp 5' overhang instead of a blunt
       end.  This pre-existing overhang accelerates 5'->3' resection
       (giving RAD51 a head start) and stabilises the D-loop during
       synthesis.
    2. **Higher specificity**: single-nucleobase mismatch discrimination
       reduces off-target editing.
    3. **Broadened PAM (NRG)**: more target-site flexibility.

    The combination of staggered cut + cssDNA donor is predicted to give
    the highest HDR rates in iPSCs among the presets defined here.
    (Nature Communications 2024; Chakraborty lab)

``default_spcas9_ssODN_hek293t``
    Standard SpCas9 with a short single-stranded oligodeoxynucleotide
    (ssODN) donor in HEK293T cells.

    This is the simplest and most common HDR configuration for small edits
    (point mutations, small tags) in easy-to-transfect cell lines.

    HEK293T cells are the "lab rat" of genome editing because:
    - High transfection efficiency
    - ~55 % cells in S/G2
    - p53 pathway is disrupted (SV40 large T antigen) -> less DSB-induced
      apoptosis
    - Baseline HDR ~25 %

    ssODN donors are short (typically 100–200 nt total, with 40–60 nt arms)
    and are thought to be incorporated partly through an SDSA-like mechanism
    and partly through single-strand template repair (SSTR), a distinct
    pathway involving FANCM and RAD51-independent annealing
    (Richardson et al., Nat. Biotechnol. 2016).

    NOTE: The ssODN pathway may not be purely SDSA.  The tract-length
    distribution from this preset should be interpreted cautiously for
    ssODNs, as the SSTR component is not explicitly modelled.

Parameter dictionary keys
-------------------------
Each preset is a dictionary with keys matching the constructor parameters
of :class:`ConversionSimulator`:

    cut_type : str
    overhang_length : int
    donor_topology : str
    homology_arm_length : int
    cell_type : str
    n_simulations : int
"""

from __future__ import annotations

from typing import Any, Dict

# ==========================================================================
# Preset configurations
# ==========================================================================

PRESET_CONFIGS: Dict[str, Dict[str, Any]] = {

    # ----------------------------------------------------------------------
    # PRESET 1: SpCas9 + cssDNA in iPSCs
    # ----------------------------------------------------------------------
    # Nuclease: SpCas9 — blunt cut, 0 bp overhang
    #   Source: Jinek et al., Science 2012
    # Donor: circular ssDNA (cssDNA), 300 bp homology arms
    #   Source: Iyer et al., CRISPR J. 2022
    # Cell type: iPSC — 8 % baseline HDR, 35 % in S/G2, p53-active
    #   Source: aggregated from multiple iPSC editing studies
    "default_spcas9_cssdna_ipsc": {
        "cut_type": "blunt",
        "overhang_length": 0,
        "donor_topology": "circular_ssDNA",
        "homology_arm_length": 300,       # Optimal for cssDNA (Iyer 2022)
        "cell_type": "iPSC",
        "n_simulations": 10_000,
    },

    # ----------------------------------------------------------------------
    # PRESET 2: enFnCas9 + cssDNA in iPSCs
    # ----------------------------------------------------------------------
    # Nuclease: enFnCas9 — staggered 5' cut, ~3 bp overhang
    #   Source: Nature Communications 2024 (Chakraborty lab)
    #   enFnCas9 is an engineered Francisella novicida Cas9 that produces
    #   staggered cuts, has broadened PAM (NRG), and shows improved HDR
    #   knock-in efficiency compared to SpCas9.
    # Donor: circular ssDNA (cssDNA), 300 bp homology arms
    #   Source: Iyer et al., CRISPR J. 2022
    # Cell type: iPSC
    #
    # The staggered cut provides several advantages:
    # - Pre-existing 5' overhang accelerates resection on one side
    # - Longer resection → longer RAD51 filament → more stable D-loop
    # - D-loop stability → longer gene conversion tracts
    # - Higher overall HDR rate (stagger_mult = 1 + 0.15 * 3 = 1.45)
    "default_enfncas9_cssdna_ipsc": {
        "cut_type": "staggered_5prime",
        "overhang_length": 3,             # ~3 bp 5' overhang (estimated)
        "donor_topology": "circular_ssDNA",
        "homology_arm_length": 300,
        "cell_type": "iPSC",
        "n_simulations": 10_000,
    },

    # ----------------------------------------------------------------------
    # PRESET 3: SpCas9 + ssODN in HEK293T
    # ----------------------------------------------------------------------
    # Nuclease: SpCas9 — blunt cut
    #   Source: Jinek et al., Science 2012
    # Donor: linear ssDNA (ssODN), short homology arms (60 bp each)
    #   Source: Richardson et al., Nat. Biotechnol. 2016
    #   ssODNs typically have 30–60 nt homology on each side of the edit.
    #   The optimal asymmetric design uses a longer PAM-distal arm, but
    #   for the symmetric model we use 60 bp.
    # Cell type: HEK293T — 25 % baseline HDR, 55 % in S/G2, p53-inactive
    #
    # CAVEAT: ssODN-mediated editing may proceed partly through SSTR
    # (single-strand template repair), a RAD51-independent pathway
    # involving FANCM and FEN1.  The SDSA model used here captures the
    # HDR component but may underestimate total precise editing for
    # short edits near the cut site.
    # (Richardson et al., 2016; Gallagher & Bhatt, unpublished)
    "default_spcas9_ssODN_hek293t": {
        "cut_type": "blunt",
        "overhang_length": 0,
        "donor_topology": "linear_ssDNA",
        "homology_arm_length": 60,        # Typical ssODN arm length
        "cell_type": "HEK293T",
        "n_simulations": 10_000,
    },
}


def list_presets() -> None:
    """Print a summary of all available preset configurations.

    This is a convenience function for interactive exploration.  It prints
    the name and key parameters of each preset in the ``PRESET_CONFIGS``
    dictionary.

    Examples
    --------
    >>> from conversion_sim.models import list_presets
    >>> list_presets()
    """
    print("Available ConversionSim preset configurations:")
    print("=" * 60)
    for name, params in PRESET_CONFIGS.items():
        print(f"\n  {name}")
        print(f"  {'─' * len(name)}")
        for key, val in params.items():
            print(f"    {key:25s} = {val}")
    print()


def get_preset(name: str) -> Dict[str, Any]:
    """Retrieve a preset configuration by name.

    Parameters
    ----------
    name : str
        Name of the preset.  Must be a key in ``PRESET_CONFIGS``.

    Returns
    -------
    dict
        A *copy* of the preset dictionary (safe to modify).

    Raises
    ------
    KeyError
        If the preset name is not found.
    """
    if name not in PRESET_CONFIGS:
        available = ", ".join(PRESET_CONFIGS.keys())
        raise KeyError(
            f"Preset {name!r} not found.  Available presets: {available}"
        )
    return dict(PRESET_CONFIGS[name])  # shallow copy
