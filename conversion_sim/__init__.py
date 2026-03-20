"""
ConversionSim — Monte Carlo Simulator of HDR Gene Conversion Tracts
====================================================================

This sub-package models the entire molecular trajectory of Homology-Directed
Repair (HDR) after a CRISPR-induced double-strand break (DSB), from the
initial DNA end resection through RAD51 filament nucleation and strand
invasion, to template-directed DNA synthesis and D-loop collapse.

**Why a Monte Carlo approach?**
Every step of HDR is stochastic.  Resection lengths, filament coverage,
synthesis processivity, and D-loop stability all vary from cell to cell.
By sampling each step thousands of times we build probability distributions
that predict, for example, the chance that a donor-encoded SNP 800 bp from
the cut site will be incorporated into the genome.

Modules
-------
resection
    Simulates 5'→3' end resection by MRN/CtIP (short-range) and
    EXO1 / BLM-DNA2 (long-range), producing 3' ssDNA overhangs.
filament
    Models RAD51 nucleoprotein filament assembly on the resected ssDNA,
    including stochastic coverage and minimum-nucleus requirements.
synthesis
    Simulates Pol-δ–dependent DNA synthesis within the D-loop during
    Synthesis-Dependent Strand Annealing (SDSA), yielding gene conversion
    tract lengths.
simulator
    The main engine that chains resection → filament → invasion → synthesis
    for *n* Monte Carlo trials and exposes summary statistics and plots.
models
    Pre-calibrated parameter dictionaries for common experimental
    configurations (e.g., SpCas9 + cssDNA in iPSCs).

Quick start
-----------
>>> from conversion_sim import ConversionSimulator
>>> sim = ConversionSimulator(
...     cut_type="staggered_5prime",
...     overhang_length=3,
...     donor_topology="circular_ssDNA",
...     homology_arm_length=300,
...     cell_type="iPSC",
...     n_simulations=10_000,
... )
>>> results = sim.run()
>>> sim.summary()
"""

from .resection import ResectionSimulator
from .filament import FilamentModel
from .synthesis import SynthesisSimulator
from .simulator import ConversionSimulator
from .models import PRESET_CONFIGS

__all__ = [
    "ResectionSimulator",
    "FilamentModel",
    "SynthesisSimulator",
    "ConversionSimulator",
    "PRESET_CONFIGS",
]
