"""
CRISPRArchitect — Multi-Modality Genome Editing Strategy Optimizer
===================================================================

A transcript-aware computational platform for unified genome editing
strategy design. Evaluates base editing, prime editing, and HDR within
a single pipeline for correcting pathogenic variants in iPSCs.

Core Pipeline (v3)
------------------
core.pipeline    : End-to-end strategy pipeline (TOPSIS + sensitivity)
core.sequence    : Transcript mapping, variant normalization, HGVS parsing
core.feasibility : Multi-nuclease PAM scanning, BE/PE/HDR feasibility
core.mosaic      : Strategy generation with consequence-aware scoring

Simulation Modules
------------------
conversion_sim : Gene conversion tract simulator (Monte Carlo)
chrombridge    : 3D chromatin distance & translocation risk predictor
topopred       : cssDNA secondary structure analyzer
loopsim        : Cohesin loop extrusion simulator

Legacy Module
-------------
mosaic         : v1 strategy optimizer (retained for backward compatibility)

Quick Start
-----------
>>> from core.pipeline.strategy_stage import StrategyPipeline
>>> from core.models import GenomicVariantInput
>>> pipeline = StrategyPipeline(cell_type="iPSC", nuclease="enFnCas9")
"""

__version__ = "3.0.0"
__author__ = "Vishal Bharti"
