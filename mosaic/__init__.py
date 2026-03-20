"""
MOSAIC: Multi-locus Optimized Strategy for Allele-specific Integrated Correction
=================================================================================

MOSAIC is a module within the CRISPRArchitect toolkit that helps molecular
biologists decide the *best editing strategy* when a patient's cells carry
mutations at **two or more genomic loci** that must all be corrected in the
same cell clone.

Biological motivation
---------------------
Many monogenic diseases are caused by compound heterozygous mutations: two
different pathogenic variants in the same gene, one on each allele. Correcting
both mutations in patient-derived iPSCs requires careful planning because:

  1. Each double-strand break (DSB) introduced by a nuclease triggers a
     p53-dependent DNA-damage response that can select against correctly
     edited iPSC clones (Ihry et al., Nature Medicine, 2018).
  2. Introducing two DSBs simultaneously risks inter-locus chromosomal
     translocations, especially if the two sites are on the same chromosome.
  3. Different mutation types (transitions, transversions, indels) are
     amenable to different editing technologies (base editing, prime
     editing, HDR) with very different safety and efficiency profiles.

MOSAIC enumerates every feasible combination of editing technologies, scores
them for efficiency, safety, time, and cost, and produces a ranked
recommendation report.

Key classes
-----------
- ``GeneStructure``        – Represent and query a gene's exon/intron layout.
- ``Mutation``             – A single pathogenic variant.
- ``MutationClassifier``   – Determine which editors can fix a mutation.
- ``StrategyEnumerator``   – Enumerate all feasible multi-locus strategies.
- ``StrategyScorer``       – Score and rank strategies.
- ``StrategyReporter``     – Produce a human-readable decision report.
"""

from .gene_structure import GeneStructure, ExonInfo
from .mutation_classifier import Mutation, MutationClassifier
from .strategy_enumerator import EditingStrategy, StrategyEnumerator
from .scorer import StrategyScorer
from .reporter import StrategyReporter

__all__ = [
    "GeneStructure",
    "ExonInfo",
    "Mutation",
    "MutationClassifier",
    "EditingStrategy",
    "StrategyEnumerator",
    "StrategyScorer",
    "StrategyReporter",
]
