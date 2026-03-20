"""
CRISPRArchitect — Multi-Site Genome Editing Strategy Optimizer
==============================================================

A computational toolkit for planning complex genome editing experiments,
particularly multi-site mutation correction in large genes.

Modules
-------
conversion_sim : Gene conversion tract simulator (Monte Carlo)
mosaic         : Multi-locus editing strategy optimizer
topopred       : cssDNA secondary structure analyzer
chrombridge    : 3D chromatin distance & translocation risk predictor
loopsim        : Cohesin loop extrusion simulator

Quick Start
-----------
>>> from crisprarchitect.conversion_sim import ConversionSimulator
>>> from crisprarchitect.mosaic.gene_structure import GeneStructure
>>> from crisprarchitect.chrombridge import ChromatinDistancePredictor
>>> from crisprarchitect.topopred import DonorAnalyzer
>>> from crisprarchitect.loopsim import LoopSimulator
>>> from crisprarchitect.utils.ensembl import fetch_gene
"""

__version__ = "0.1.0"
__author__ = "CRISPRArchitect Team"
