"""
CRISPRArchitect v2 Core
========================

Transcript-aware, consequence-guided genome editing strategy design.

This package extends the v1 MOSAIC framework with:

- Transcript-aware variant normalization via Ensembl
- Reference allele validation
- Coding and splice-site consequence annotation
- Modality-specific feasibility engines (BE, PE, HDR)
- Consequence-aware multi-objective scoring
- Benchmark evaluation framework

Subpackages
-----------
sequence    Transcript mapping, variant normalization, coding annotation
feasibility PAM scanning, base/prime/HDR feasibility engines
mosaic      Strategy generation and consequence-aware scoring
pipeline    End-to-end orchestration
"""

__version__ = "2.0.0"
