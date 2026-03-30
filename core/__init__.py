"""
CRISPRArchitect Core (v3)
==========================

Transcript-aware, consequence-guided, multi-nuclease genome editing
strategy design with TOPSIS multi-criteria decision analysis.

v3 additions over v2:
- Multi-nuclease engine (SpCas9, enFnCas9, SpCas9-NG, SpRY, Cas12a)
- Multi-editor evaluation (ABE7.10, ABE8e, BE4max, CBE4, PE2, PE3)
- TOPSIS scoring with Monte Carlo sensitivity analysis
- CFD/MIT off-target specificity scoring
- HGVS parser with ClinVar batch ingestion

Subpackages
-----------
sequence    Transcript mapping, variant normalization, HGVS parsing
feasibility Multi-nuclease PAM scanning, BE/PE/HDR feasibility engines
mosaic      Strategy generation and consequence-aware scoring
pipeline    End-to-end orchestration with TOPSIS + sensitivity
"""

__version__ = "3.0.0"
