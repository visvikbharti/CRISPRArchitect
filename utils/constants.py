"""
Biological Constants for CRISPRArchitect
========================================

Every constant in this file is derived from published experimental data.
References are provided inline. These constants parameterize the simulation
models throughout the toolkit.

IMPORTANT: These are best-estimate values from the literature. Actual values
vary by cell type, locus, and experimental conditions. The simulation models
use distributions around these central values where appropriate.
"""

# =============================================================================
# DNA PHYSICAL PROPERTIES
# =============================================================================

# Rise per base pair in B-form DNA (nanometers)
# Source: Standard B-DNA geometry
DNA_RISE_PER_BP_NM = 0.34

# Persistence length of dsDNA (nanometers)
# Source: Hagerman, Ann Rev Biophys Biophys Chem, 1988
# dsDNA behaves as a stiff rod below this length scale
DSDNA_PERSISTENCE_LENGTH_NM = 50.0

# Persistence length of ssDNA (nanometers)
# Source: Murphy et al., Biophys J, 2004
# ssDNA is much more flexible than dsDNA
SSDNA_PERSISTENCE_LENGTH_NM = 1.5

# Rise per nucleotide in ssDNA (nanometers)
# Source: Approximate, varies with sequence
SSDNA_RISE_PER_NT_NM = 0.59

# Contour length per nucleotide in ssDNA (nanometers)
# Source: Murphy et al., Biophys J, 2004
SSDNA_CONTOUR_PER_NT_NM = 0.63


# =============================================================================
# END RESECTION PARAMETERS
# =============================================================================

# Short-range resection by MRN/CtIP
# Typical initial resection length (base pairs)
# Source: Symington, Ann Rev Genet, 2011; Cejka, Ann Rev Genet, 2015
# MRN/CtIP creates a short ssDNA overhang of ~100-300 bp
SHORT_RESECTION_MEAN_BP = 200
SHORT_RESECTION_STD_BP = 80

# Long-range resection by EXO1 or BLM-DNA2
# Source: Symington, 2011; Zhou et al., Mol Cell, 2014
# Long-range resection can extend 1-5 kb, sometimes further
LONG_RESECTION_MEAN_BP = 2000
LONG_RESECTION_STD_BP = 1000
LONG_RESECTION_MIN_BP = 300
LONG_RESECTION_MAX_BP = 10000

# Resection rate (nucleotides per second)
# Source: Zhu et al., Cell, 2008 (in vitro rates)
# EXO1: ~100-200 nt/s in vitro; slower in vivo due to chromatin
RESECTION_RATE_NT_PER_SEC = 50  # Conservative in vivo estimate


# =============================================================================
# RAD51 FILAMENT PARAMETERS
# =============================================================================

# RAD51 monomer binding footprint (nucleotides)
# Source: Ogawa et al., Science, 1993; Yu et al., Mol Cell Biol, 2001
# Each RAD51 monomer covers 3 nucleotides of ssDNA
RAD51_FOOTPRINT_NT = 3

# RAD51 filament nucleation minimum (monomers)
# Source: Estimated from in vitro studies
# A minimum nucleus of ~5-8 RAD51 monomers is needed to start filament growth
RAD51_NUCLEATION_MIN_MONOMERS = 5

# RAD51 filament growth rate (monomers per second)
# Source: Estimated from single-molecule studies
RAD51_GROWTH_RATE_PER_SEC = 10

# Minimum homology for stable strand invasion (base pairs)
# Source: Qi et al., Cell, 2015 (8-nt microhomology sampling)
# Stable invasion requires ~15-20 bp of continuous homology
MIN_HOMOLOGY_FOR_INVASION_BP = 15


# =============================================================================
# GENE CONVERSION / SYNTHESIS PARAMETERS
# =============================================================================

# Gene conversion tract length in mammalian cells
# Source: Elliott et al., Mol Cell Biol, 1998
# Most tracts are 200-2000 bp, with median ~500 bp
CONVERSION_TRACT_MEAN_BP = 500
CONVERSION_TRACT_STD_BP = 400
CONVERSION_TRACT_MIN_BP = 50
CONVERSION_TRACT_MAX_BP = 5000

# DNA polymerase processivity during HDR synthesis
# Source: Estimated from gene conversion tract distributions
# Pol delta extends from the invading 3' end
SYNTHESIS_PROCESSIVITY_MEAN_BP = 600
SYNTHESIS_PROCESSIVITY_STD_BP = 350

# SDSA displacement probability per bp synthesized
# Source: Derived from tract length distributions
# Probability of D-loop collapse increases with synthesis length
# This creates the exponential-like decay in tract length distribution
SDSA_DISPLACEMENT_PROB_PER_BP = 0.002  # ~0.2% per bp -> median ~350 bp


# =============================================================================
# CUT STRUCTURE PARAMETERS
# =============================================================================

# SpCas9 blunt cut: both strands cut at position -3 from PAM
# Source: Jinek et al., Science, 2012
SPAS9_CUT_POSITION = -3  # bp upstream of PAM on target strand

# SpCas9 stagger: non-target strand can be cut at -3 to -5
# Source: Shou et al., Cell Discovery, 2019
SPCAS9_STAGGER_RANGE = (0, 1)  # 0-1 bp stagger (mostly blunt)

# vCas9 staggered cut: 5' overhangs of >=6 bp
# Source: Chauhan et al., PNAS, 2023
VCAS9_STAGGER_RANGE = (4, 8)  # 4-8 bp 5' overhang

# enFnCas9: improved HDR suggests some staggering
# Source: Nature Communications, 2024 (inferred from improved HDR rates)
# Exact stagger data not published; using conservative estimate
ENFNCAS9_STAGGER_RANGE = (2, 5)  # Estimated 2-5 bp 5' overhang

# Cas12a/Cpf1: well-characterized 5' overhangs
# Source: Zetsche et al., Cell, 2015; Stella et al., Nature, 2017
CAS12A_STAGGER_RANGE = (4, 5)  # 4-5 bp 5' overhang


# =============================================================================
# STAGGER EFFECT ON HDR
# =============================================================================

# HDR enhancement factor from staggered cuts
# Source: Chauhan et al., PNAS, 2023
# vCas9 (6bp stagger) increased precise editing 1.4-2.8x (mean 1.9x)
# Relationship: longer overhang -> more resection -> more HDR
HDR_ENHANCEMENT_PER_BP_OVERHANG = 0.15  # ~15% increase per bp of overhang

# Baseline HDR fraction (blunt cut, no stagger)
# Source: Multiple studies; varies enormously by cell type
BASELINE_HDR_FRACTION_HEK293T = 0.25
BASELINE_HDR_FRACTION_IPSC = 0.08
BASELINE_HDR_FRACTION_K562 = 0.20


# =============================================================================
# DONOR TEMPLATE PARAMETERS
# =============================================================================

# Donor topology effectiveness multipliers (relative to linear dsDNA = 1.0)
# Source: Iyer et al., CRISPR Journal, 2022
DONOR_TOPOLOGY_MULTIPLIER = {
    "linear_dsDNA": 1.0,
    "plasmid_dsDNA": 0.8,
    "linear_ssDNA": 1.5,         # ~1.5x better than linear dsDNA
    "circular_ssDNA": 3.0,       # ~2x better than linear ssDNA (Iyer et al.)
    "AAV_ssDNA": 4.0,            # Best for large knock-ins in iPSCs
}

# Optimal homology arm lengths (bp)
# Source: Iyer et al., 2022 (300 nt optimal for cssDNA)
# Source: Richardson et al., 2016 (asymmetric design for ssODN)
OPTIMAL_HA_LENGTH_CSSDNA = 300   # bp per arm
OPTIMAL_HA_LENGTH_LSSDNA = 300
OPTIMAL_HA_LENGTH_DSDNA = 800
OPTIMAL_HA_LENGTH_AAV = 800

# cssDNA nuclease resistance factor
# Circular ssDNA is resistant to exonucleases, increasing effective half-life
# Source: Iyer et al., 2022 (circularization of lssDNA improved HDR)
CSSDNA_HALFLIFE_MULTIPLIER = 3.0  # ~3x longer intracellular half-life vs linear


# =============================================================================
# 3D GENOME / POLYMER PHYSICS CONSTANTS
# =============================================================================

# Hi-C contact frequency power-law exponent
# P(contact) ~ s^(-gamma), where s = genomic distance
# Source: Lieberman-Aiden et al., Science, 2009
HICCONTACT_POWER_LAW_GAMMA = 1.08

# Chromatin fiber compaction ratio (bp per nm)
# Source: Varies by chromatin state
# Euchromatin: ~10 bp/nm (30nm fiber)
# Heterochromatin: ~40 bp/nm (compacted)
CHROMATIN_COMPACTION_EUCHROMATIN = 10  # bp per nm
CHROMATIN_COMPACTION_HETEROCHROMATIN = 40

# Kuhn length of chromatin fiber (nm)
# Source: Bystricky et al., PNAS, 2004 (yeast); estimated for mammalian
CHROMATIN_KUHN_LENGTH_NM = 300  # ~30 kb per Kuhn segment

# Gaussian chain scaling exponent for 3D distance
# <R^2> = b^2 * N^(2*nu), where N = number of Kuhn segments
# Source: Polymer physics; nu = 0.5 for ideal chain, ~0.33 for confined
POLYMER_SCALING_EXPONENT = 0.5  # Ideal chain (approximate)

# Nuclear diameter (micrometers)
# Source: Typical mammalian cell nucleus
NUCLEAR_DIAMETER_UM = 10.0

# Translocation frequency scaling
# Translocation frequency ~ Hi-C contact frequency (approximately)
# Source: Chiarle et al., Cell, 2011; Zhang et al., Cell, 2012
# Normalized so that loci at 1 Mb have translocation prob ~ 1e-3 per DSB pair
TRANSLOCATION_BASELINE_PROB_1MB = 1e-3


# =============================================================================
# CELL TYPE PARAMETERS
# =============================================================================

CELL_TYPE_PARAMS = {
    "iPSC": {
        "hdr_base_efficiency": 0.08,      # 8% baseline HDR
        "cell_cycle_s_g2_fraction": 0.35,  # 35% of cells in S/G2
        "p53_active": True,                # Strong p53 response to DSBs
        "viability_single_dsb": 0.55,      # 55% viability after 1 DSB
        "viability_dual_dsb": 0.30,        # 30% viability after 2 DSBs
        "description": "Human induced pluripotent stem cells",
    },
    "HEK293T": {
        "hdr_base_efficiency": 0.25,
        "cell_cycle_s_g2_fraction": 0.55,
        "p53_active": False,               # p53 pathway often disrupted
        "viability_single_dsb": 0.85,
        "viability_dual_dsb": 0.70,
        "description": "Human embryonic kidney 293T cells",
    },
    "K562": {
        "hdr_base_efficiency": 0.20,
        "cell_cycle_s_g2_fraction": 0.50,
        "p53_active": False,
        "viability_single_dsb": 0.80,
        "viability_dual_dsb": 0.65,
        "description": "Human chronic myeloid leukemia cells",
    },
    "T_cell": {
        "hdr_base_efficiency": 0.15,
        "cell_cycle_s_g2_fraction": 0.45,  # After activation
        "p53_active": True,
        "viability_single_dsb": 0.70,
        "viability_dual_dsb": 0.50,
        "description": "Primary human T cells (activated)",
    },
    "HSC": {
        "hdr_base_efficiency": 0.10,
        "cell_cycle_s_g2_fraction": 0.25,  # Mostly quiescent
        "p53_active": True,
        "viability_single_dsb": 0.60,
        "viability_dual_dsb": 0.35,
        "description": "Hematopoietic stem cells",
    },
}


# =============================================================================
# NUCLEASE PARAMETERS
# =============================================================================

NUCLEASE_PARAMS = {
    "SpCas9": {
        "pam": "NGG",
        "cut_type": "blunt",
        "stagger_bp": 0,
        "hdr_multiplier": 1.0,
        "specificity": "moderate",
        "description": "Standard S. pyogenes Cas9",
    },
    "enFnCas9": {
        "pam": "NRG",  # Broadened PAM
        "cut_type": "staggered_5prime",
        "stagger_bp": 3,  # Estimated
        "hdr_multiplier": 1.5,  # Improved HDR knock-in
        "specificity": "high",  # Single-nucleobase specificity
        "description": "Engineered F. novicida Cas9 (Chakraborty lab)",
    },
    "Cas12a": {
        "pam": "TTTV",
        "cut_type": "staggered_5prime",
        "stagger_bp": 5,
        "hdr_multiplier": 1.4,
        "specificity": "high",
        "description": "Acidaminococcus sp. Cas12a (Cpf1)",
    },
    "vCas9": {
        "pam": "NGG",
        "cut_type": "staggered_5prime",
        "stagger_bp": 6,
        "hdr_multiplier": 1.9,  # Mean from Chauhan et al., 2023
        "specificity": "moderate",
        "description": "Staggered-cut SpCas9 variant (MIT)",
    },
}
