"""
Biological Constants for CRISPRArchitect
========================================

Parameter Provenance Guide
--------------------------
Each constant is annotated with one of three evidence levels:

  [MEASURED]  Value directly from a published measurement with citation.
  [DERIVED]   Value computed from published data via a stated procedure.
  [ASSUMED]   Modeling assumption with stated rationale; no direct measurement.
              These parameters should be explored via sensitivity analysis.

IMPORTANT: These are best-estimate values from the literature. Actual values
vary by cell type, locus, and experimental conditions. The simulation models
use distributions around these central values where appropriate. Parameters
marked [ASSUMED] are explored via the TOPSIS sensitivity analysis framework.
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
# [DERIVED] Mean of the 100-300 bp range reported in:
#   Symington, Ann Rev Genet, 2011 (review of MRN/CtIP endonuclease)
#   Cejka, Ann Rev Genet, 2015 (in vitro reconstitution, ~100 bp products)
#   Shibata et al., Mol Cell, 2014 (in vivo imaging, ~200 bp initial resection)
# The std=80 is set so that 95% of draws fall within [40, 360] bp,
# spanning the observed biological range.
SHORT_RESECTION_MEAN_BP = 200
SHORT_RESECTION_STD_BP = 80

# Long-range resection by EXO1 or BLM-DNA2
# [DERIVED] Range from:
#   Symington, 2011: long-range resection can extend "several kilobases"
#   Zhou et al., Mol Cell, 2014: in vivo resection tracks reaching 3-5 kb
#   Gravel et al., Genes Dev, 2008: EXO1 processivity ~1-5 kb in vitro
# The mean=2000 and LogNormal(median=1500, sigma=0.8) used in the simulation
# are calibrated so that the 5th percentile ≈ 400 bp and 95th ≈ 6 kb,
# matching the observed range. See resection.py for the actual distribution.
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

# Gene conversion tract length parameters for long-donor SDSA
# [ASSUMED] These are modeling parameters, NOT direct measurements.
#
# IMPORTANT: No published study has directly measured the distribution of
# SDSA tract lengths for exogenous long-donor (cssDNA/dsDNA) HDR in
# mammalian cells. The values below are derived from FUNCTIONAL evidence:
#
# What is known from the literature:
#   - Elliott et al., MCB 1998 (PMID 9418857): Measured chromosomal gene
#     conversion tracts from I-SceI-induced DSBs in mouse cells. Found 80%
#     of tracts were <=58 bp, max ~511 bp. NOTE: This used endogenous
#     chromosomal substrates, NOT exogenous donors. Not directly applicable.
#   - Kan et al., Genome Res 2017 (PMID 28356322): Measured ODN-mediated
#     editing tracts. Average ~20 bp. NOTE: This is SSTR pathway, not SDSA.
#     Not applicable to long-donor HDR.
#   - Stark lab, G3 2017 (PMID 28179392): Human cell SDSA assay requires
#     >=350 bp of new synthesis for product formation, confirming that SDSA
#     can produce tracts of several hundred bp in human cells.
#   - Functional evidence from HDR experiments: Successful knock-ins with
#     500-1000 bp homology arms (IDT recommendations; reviewed in
#     Banan 2020, Mol Ther Methods Clin Dev 18:583-596) imply that donor
#     sequence hundreds of bp from the cut site is routinely incorporated.
#
# Our model: p = 0.002 gives mean = 500 bp, median = 347 bp.
# This is consistent with the functional evidence (HDR works with 300-1000
# bp arms) and the Stark lab SDSA assay (~350 bp synthesis minimum).
# The value is explicitly a modeling assumption that should be explored
# via sensitivity analysis (p range: 0.001 to 0.005).
CONVERSION_TRACT_MEAN_BP = 500     # [ASSUMED] See derivation above
CONVERSION_TRACT_STD_BP = 400      # [ASSUMED]
CONVERSION_TRACT_MIN_BP = 50       # [ASSUMED] Below this, mismatch repair erases
CONVERSION_TRACT_MAX_BP = 5000     # [ASSUMED] Upper bound; tracts >5 kb unobserved

# DNA polymerase processivity during HDR synthesis
# [ASSUMED] Estimated from the tract length model above.
SYNTHESIS_PROCESSIVITY_MEAN_BP = 600
SYNTHESIS_PROCESSIVITY_STD_BP = 350

# SDSA displacement probability per bp synthesized
# [ASSUMED] This is the CENTRAL modeling parameter of ConversionSim.
#
# The geometric distribution P(tract >= d) = (1-p)^d gives:
#   - Mean = 1/p, Median = ln(2)/p
#
# With p = 0.002: mean = 500 bp, median = 347 bp, 90th pctile ~1150 bp.
#
# Evidence basis (functional, not direct measurement):
#   1. Stark lab SDSA assay (G3, 2017): SDSA produces >=350 bp of new
#      synthesis in human cells, confirming tracts reach this range.
#   2. Successful HDR with 300-1000 bp homology arms in mammalian cells
#      implies routine incorporation at distances of several hundred bp.
#   3. Helicase regulation: BLM and RTEL1 disrupt D-loops after "a few
#      hundred nucleotides" of synthesis (Gallagher & Haber, ACS Chem
#      Biol, 2018, PMC5835394), consistent with our model's mean of 500 bp.
#
# SENSITIVITY: This parameter should be explored across p = 0.001 to 0.005
# (mean tracts 200 to 1000 bp). The TOPSIS sensitivity analysis framework
# does not currently vary simulation parameters, but ConversionSim supports
# user-configurable p via direct instantiation.
#
# NOT calibrated to:
#   - Elliott et al., MCB 1998 (endogenous substrate, short tracts <58 bp)
#   - Kan et al., Genome Res 2017 (ODN/SSTR pathway, tracts ~20 bp)
#   - Paquet et al., Nature 2016 (ssODN/SSTR pathway, tracts ~10-50 bp)
SDSA_DISPLACEMENT_PROB_PER_BP = 0.002


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
# [ASSUMED] Exact stagger data NOT published as of March 2026.
# Rationale: enFnCas9 shows improved HDR knock-in rates (Acharya et al.,
# Nat Commun 2024, 15:5471), which is consistent with 5' overhang generation (by analogy
# with vCas9 and Cas12a). FnCas9 crystal structure (Hirano et al., Cell,
# 2016) suggests a non-target-strand cleavage site offset of 2-4 bp.
# We use stagger_bp=3 as a conservative midpoint estimate.
# WARNING: This parameter should be updated when direct biochemical
# characterization of enFnCas9 cut sites is published. Users can override
# via NUCLEASE_PARAMS["enFnCas9"]["stagger_bp"].
ENFNCAS9_STAGGER_RANGE = (2, 5)  # Estimated; see WARNING above

# Cas12a/Cpf1: well-characterized 5' overhangs
# Source: Zetsche et al., Cell, 2015; Stella et al., Nature, 2017
CAS12A_STAGGER_RANGE = (4, 5)  # 4-5 bp 5' overhang


# =============================================================================
# STAGGER EFFECT ON HDR
# =============================================================================

# HDR enhancement factor from staggered cuts
# [DERIVED] Single-point linear fit to Chauhan et al., PNAS, 2023.
# Data: vCas9 produces >=6 bp 5' overhangs and achieves 1.4-2.8x (mean
# 1.9x) improvement in precise editing over WT SpCas9 (blunt cut).
# Linear model: enhancement = 1 + k * overhang_bp, where k = (1.9 - 1)/6 = 0.15.
#
# LIMITATIONS of this estimate:
#   1. Fit to a single data point (vCas9 at 6 bp). No intermediate lengths tested.
#   2. Assumes linearity; the true relationship may saturate or be step-like.
#   3. Enhancement may partly reflect vCas9 protein differences, not just overhang.
#   4. Extrapolation to other nucleases (enFnCas9, Cas12a) is untested.
# This parameter is explored in ConversionSim sensitivity analysis.
HDR_ENHANCEMENT_PER_BP_OVERHANG = 0.15

# Baseline HDR fraction (blunt cut, no stagger)
# Source: Multiple studies; varies enormously by cell type
BASELINE_HDR_FRACTION_HEK293T = 0.25
BASELINE_HDR_FRACTION_IPSC = 0.08
BASELINE_HDR_FRACTION_K562 = 0.20


# =============================================================================
# DONOR TEMPLATE PARAMETERS
# =============================================================================

# Donor topology effectiveness multipliers (relative to linear dsDNA = 1.0)
#
# Evidence for each value:
#   linear_dsDNA = 1.0   [MEASURED] Baseline reference. Standard plasmid-free dsDNA.
#   plasmid_dsDNA = 0.8  [DERIVED]  Plasmid backbone may trigger innate immune
#                        sensing (Bhatt et al., ACS Chem Biol, 2015). ~20% reduction
#                        commonly observed vs. linear dsDNA.
#   linear_ssDNA = 1.5   [MEASURED] Richardson et al., Nat Biotechnol, 2016:
#                        asymmetric ssODN achieves ~60% higher knock-in than dsDNA.
#                        Conservative central estimate.
#   circular_ssDNA = 3.0 [DERIVED]  Iyer et al., CRISPR J, 2022: cssDNA achieved
#                        ~18% HDR vs ~9.5% for linear ssDNA (ratio 1.9x).
#                        Relative to linear dsDNA baseline: 1.5 * 1.9 ≈ 2.85,
#                        rounded to 3.0. The Xie et al. (Nat Biotechnol, 2025)
#                        GATALYST system reports up to 70% in iPSCs, but this
#                        includes protocol optimization beyond topology alone.
#   AAV_ssDNA = 4.0      [ASSUMED]  AAV6 HDR template delivery achieves 10-50%
#                        knock-in in iPSCs (Martin et al., Cell Stem Cell, 2019;
#                        Dever et al., Nature, 2016). The 4.0x estimate reflects
#                        the combined benefit of nuclear delivery + ssDNA template.
#                        Highly variable; 2-6x range depending on MOI and locus.
DONOR_TOPOLOGY_MULTIPLIER = {
    "linear_dsDNA": 1.0,
    "plasmid_dsDNA": 0.8,
    "linear_ssDNA": 1.5,
    "circular_ssDNA": 3.0,
    "AAV_ssDNA": 4.0,
}

# Optimal homology arm lengths (bp)
# [ASSUMED] These are practical recommendations, not optimized values.
# Iyer et al. (CRISPR J, 2022) used 35-60 nt arms for short inserts and
# ~1 kb arms for endogenous locus tagging, but did NOT systematically
# optimize arm length. The 300 bp default for cssDNA/lssDNA is based on
# general HDR design guidelines (reviewed in Banan, 2020; IDT protocols)
# suggesting 200-1000 bp arms for dsDNA and shorter for ssDNA donors.
# Richardson et al. (Nat Biotechnol, 2016) optimized ssODN at ~90 bp.
OPTIMAL_HA_LENGTH_CSSDNA = 300   # [ASSUMED] bp per arm; see note above
OPTIMAL_HA_LENGTH_LSSDNA = 300   # [ASSUMED]
OPTIMAL_HA_LENGTH_DSDNA = 800    # [ASSUMED] Standard for dsDNA donors
OPTIMAL_HA_LENGTH_AAV = 800      # [ASSUMED] AAV packaging allows ~800 bp/arm

# cssDNA nuclease resistance factor
# [ASSUMED] Circular ssDNA lacks free 3'/5' ends, making it resistant to
# cellular exonucleases (RecBCD, ExoI, ExoIII homologs). Iyer et al.
# (CRISPR J, 2022) showed cssDNA achieves ~2x higher HDR than linear ssDNA.
# This is partially attributable to increased intracellular persistence.
# The 3.0x half-life estimate is extrapolated from the HDR improvement,
# assuming a ~linear relationship between donor persistence and HDR yield
# in the donor-limiting regime. No direct half-life measurement was reported.
CSSDNA_HALFLIFE_MULTIPLIER = 3.0


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
# [ASSUMED] Bystricky et al., PNAS, 2004 measured ~200-300 nm in yeast.
# Mammalian chromatin is debated: estimates range 100-500 nm depending
# on measurement method (FISH: ~200 nm, Mateos-Langerak et al., PNAS,
# 2009; Hi-C polymer fits: ~300 nm, Fudenberg et al., Cell Rep, 2016).
# We use 300 nm (~30 kb per Kuhn segment in euchromatin at 10 bp/nm).
CHROMATIN_KUHN_LENGTH_NM = 300

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

# Cell-type parameters for HDR efficiency and DSB toxicity modeling.
#
# Evidence sources for iPSC parameters (primary focus of this tool):
#   hdr_base_efficiency: 5-15% range for RNP + ssODN in iPSCs
#     (Paquet et al., Nature 2016: 5-30% for point mutations;
#      median ~10% without enhancement). We use 0.08 as a conservative
#     estimate for unenhanced HDR with cssDNA donors.
#   cell_cycle_s_g2_fraction: iPSCs have an unusually short G1 phase
#     (~30-40% S/G2 per Becker et al., PNAS 2006; Ghule et al., MCB 2011).
#     We use 0.35 as a conservative estimate.
#   p53_active: iPSCs retain functional p53 that selects against DSB-bearing
#     cells (Ihry et al., Nat Med 2018; Haapaniemi et al., Nat Med 2018).
#   viability_single_dsb: Ihry et al., 2018 Fig 2: ~40-60% colony survival
#     after single Cas9 cut in iPSCs. We use 0.55.
#   viability_dual_dsb: [ASSUMED] Approximately viability^1.5 for two
#     independent DSBs, accounting for synergistic toxicity. 0.55^1.5 ≈ 0.41;
#     we use 0.30 to account for translocation-mediated additional lethality.
#
# For other cell types, [ASSUMED] values are order-of-magnitude estimates
# from the indicated literature. Users should substitute measured values
# for their specific experimental system.

CELL_TYPE_PARAMS = {
    "iPSC": {
        "hdr_base_efficiency": 0.08,      # [DERIVED] Paquet 2016, conservative
        "cell_cycle_s_g2_fraction": 0.35,  # [MEASURED] Becker 2006; Ghule 2011
        "p53_active": True,                # [MEASURED] Ihry 2018; Haapaniemi 2018
        "viability_single_dsb": 0.55,      # [MEASURED] Ihry 2018 Fig 2
        "viability_dual_dsb": 0.30,        # [ASSUMED] See derivation above
        "description": "Human induced pluripotent stem cells",
    },
    "HEK293T": {
        "hdr_base_efficiency": 0.25,       # [MEASURED] Ran et al., Nat Protoc 2013
        "cell_cycle_s_g2_fraction": 0.55,  # [MEASURED] Rapidly dividing line
        "p53_active": False,               # [MEASURED] SV40 LT inactivates p53
        "viability_single_dsb": 0.85,      # [ASSUMED] Robust line, high viability
        "viability_dual_dsb": 0.70,        # [ASSUMED] ~0.85^1.3
        "description": "Human embryonic kidney 293T cells",
    },
    "K562": {
        "hdr_base_efficiency": 0.20,       # [MEASURED] DeWitt et al., Nat Biotechnol 2016
        "cell_cycle_s_g2_fraction": 0.50,  # [MEASURED] Actively dividing CML line
        "p53_active": False,               # [MEASURED] Homozygous TP53 frameshift
        "viability_single_dsb": 0.80,      # [ASSUMED] Robust line
        "viability_dual_dsb": 0.65,        # [ASSUMED]
        "description": "Human chronic myeloid leukemia cells",
    },
    "T_cell": {
        "hdr_base_efficiency": 0.15,       # [MEASURED] Roth et al., Nature 2018
        "cell_cycle_s_g2_fraction": 0.45,  # [ASSUMED] After CD3/CD28 activation
        "p53_active": True,                # [ASSUMED] Primary cells, intact p53
        "viability_single_dsb": 0.70,      # [ASSUMED] Fragile primary cells
        "viability_dual_dsb": 0.50,        # [ASSUMED]
        "description": "Primary human T cells (activated)",
    },
    "HSC": {
        "hdr_base_efficiency": 0.10,       # [MEASURED] Dever et al., Nature 2016
        "cell_cycle_s_g2_fraction": 0.25,  # [MEASURED] Mostly quiescent
        "p53_active": True,                # [ASSUMED] Primary cells, intact p53
        "viability_single_dsb": 0.60,      # [ASSUMED]
        "viability_dual_dsb": 0.35,        # [ASSUMED]
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
        "reference": "Jinek et al., Science, 2012",
    },
    "enFnCas9": {
        "pam": "NRG",  # [MEASURED] Broadened PAM (R = A or G); Nat Commun 2024
        "cut_type": "staggered_5prime",  # [ASSUMED] See ENFNCAS9_STAGGER_RANGE
        "stagger_bp": 3,  # [ASSUMED] Not directly measured; see WARNING in stagger section
        "hdr_multiplier": 1.5,  # [ASSUMED] Improved HDR observed but not quantified as fold-change
        "specificity": "high",  # [MEASURED] Single-nucleobase specificity; Nat Commun 2024
        "description": "Engineered F. novicida Cas9 (Chakraborty lab)",
        "reference": "Hirano et al., Cell, 2016; Acharya et al., Nat Commun, 2024 (15:5471)",
    },
    "SpCas9-NG": {
        "pam": "NG",  # Relaxed PAM
        "cut_type": "blunt",
        "stagger_bp": 0,
        "hdr_multiplier": 1.0,
        "specificity": "low",  # Broader PAM reduces specificity
        "description": "SpCas9-NG with relaxed NG PAM requirement",
        "reference": "Nishimasu et al., Science, 2018",
    },
    "SpRY": {
        "pam": "NNN",  # Near-PAMless (NRN > NYN preference)
        "cut_type": "blunt",
        "stagger_bp": 0,
        "hdr_multiplier": 0.8,  # Slightly reduced activity vs SpCas9
        "specificity": "very_low",  # Near-PAMless increases off-target risk
        "description": "SpRY near-PAMless Cas9 variant",
        "reference": "Walton et al., Science, 2020",
    },
    "Cas12a": {
        "pam": "TTTV",
        "cut_type": "staggered_5prime",
        "stagger_bp": 5,
        "hdr_multiplier": 1.4,
        "specificity": "high",
        "description": "Acidaminococcus sp. Cas12a (Cpf1)",
        "reference": "Zetsche et al., Cell, 2015",
    },
    "vCas9": {
        "pam": "NGG",
        "cut_type": "staggered_5prime",
        "stagger_bp": 6,
        "hdr_multiplier": 1.9,  # Mean from Chauhan et al., 2023
        "specificity": "moderate",
        "description": "Staggered-cut SpCas9 variant (MIT)",
        "reference": "Chauhan et al., PNAS, 2023",
    },
}


# =============================================================================
# BASE EDITOR PROFILES
# =============================================================================
# Each editor profile defines:
#   - editor_type: "ABE" or "CBE"
#   - window_start, window_end: 1-indexed editing window in 20-mer protospacer
#   - compatible_nucleases: list of nucleases this editor has been paired with
#   - efficiency_class: relative activity level ("high", "moderate", "low")
#   - evidence_tier: "A" = experimentally characterized, "B" = extrapolated
#
# References are cited per-editor for traceability.

BASE_EDITOR_PROFILES = {
    # ── Adenine Base Editors ────────────────────────────────────────────
    "ABE7.10": {
        "editor_type": "ABE",
        "window_start": 4,
        "window_end": 7,
        "compatible_nucleases": ["SpCas9"],
        "efficiency_class": "moderate",
        "evidence_tier": "A",
        "description": "Original ABE with TadA-TadA* heterodimer",
        "reference": "Gaudelli et al., Nature, 2017",
    },
    "ABE8e": {
        "editor_type": "ABE",
        "window_start": 3,
        "window_end": 9,
        # Window note: Canonical SpCas9-ABE window is positions 4-8.
        # ABE8e is more processive than ABE7.10 (evolved TadA-8e monomer)
        # and shows extended activity at positions 3 and 9 (Richter et al.,
        # 2020, Fig 2). We use 3-9 as the EXTENDED window to maximize
        # rescue of borderline cases. The canonical 4-8 window captures
        # the highest-efficiency positions; 3 and 9 have reduced but
        # nonzero activity. This is a deliberate modeling choice to be
        # inclusive rather than conservative.
        "compatible_nucleases": ["SpCas9", "enFnCas9", "SpCas9-NG", "SpRY"],
        "efficiency_class": "high",
        "evidence_tier": "A",
        "description": "ABE8e with evolved TadA-8e monomer, extended window (3-9), higher activity",
        "reference": "Richter et al., Nat Biotechnol, 2020 (PMID 32433547)",
    },
    "ABE8e-SpCas9-NG": {
        "editor_type": "ABE",
        "window_start": 3,
        "window_end": 9,
        "compatible_nucleases": ["SpCas9-NG"],
        "efficiency_class": "moderate",  # Reduced vs SpCas9 due to NG PAM
        "evidence_tier": "B",  # Inferred from ABE8e + SpCas9-NG combination
        "description": "ABE8e fused with SpCas9-NG for NG PAM access",
        "reference": "Richter et al., 2020; Nishimasu et al., 2018",
    },
    "ABE8e-SpRY": {
        "editor_type": "ABE",
        "window_start": 3,
        "window_end": 9,
        "compatible_nucleases": ["SpRY"],
        "efficiency_class": "low",  # SpRY has reduced on-target activity
        "evidence_tier": "B",
        "description": "ABE8e fused with near-PAMless SpRY",
        "reference": "Walton et al., Science, 2020; Richter et al., 2020",
    },
    "ABE8e-enFnCas9": {
        "editor_type": "ABE",
        "window_start": 3,
        "window_end": 9,
        "compatible_nucleases": ["enFnCas9"],
        "efficiency_class": "moderate",
        "evidence_tier": "B",  # enFnCas9 + ABE8e not directly published
        "description": "ABE8e paired with enFnCas9 (NRG PAM)",
        "reference": "Richter et al., 2020; Acharya et al., Nat Commun, 2024 (15:5471)",
    },
    # ── Cytosine Base Editors ───────────────────────────────────────────
    "BE4max": {
        "editor_type": "CBE",
        "window_start": 4,
        "window_end": 8,
        "compatible_nucleases": ["SpCas9"],
        "efficiency_class": "high",
        "evidence_tier": "A",
        "description": "Optimized CBE with APOBEC1 + UGI",
        "reference": "Koblan et al., Nat Biotechnol, 2018",
    },
    "BE4max-SpCas9-NG": {
        "editor_type": "CBE",
        "window_start": 4,
        "window_end": 8,
        "compatible_nucleases": ["SpCas9-NG"],
        "efficiency_class": "moderate",
        "evidence_tier": "B",
        "description": "BE4max with SpCas9-NG for relaxed PAM",
        "reference": "Koblan et al., 2018; Nishimasu et al., 2018",
    },
    "BE4max-SpRY": {
        "editor_type": "CBE",
        "window_start": 4,
        "window_end": 8,
        "compatible_nucleases": ["SpRY"],
        "efficiency_class": "low",
        "evidence_tier": "B",
        "description": "BE4max with near-PAMless SpRY",
        "reference": "Walton et al., Science, 2020; Koblan et al., 2018",
    },
    "BE4max-enFnCas9": {
        "editor_type": "CBE",
        "window_start": 4,
        "window_end": 8,
        "compatible_nucleases": ["enFnCas9"],
        "efficiency_class": "moderate",
        "evidence_tier": "B",
        "description": "BE4max paired with enFnCas9 (NRG PAM)",
        "reference": "Koblan et al., 2018; Acharya et al., Nat Commun, 2024 (15:5471)",
    },
    # ── Legacy aliases (backward compatibility) ─────────────────────────
    "ABE": {
        "editor_type": "ABE",
        "window_start": 4,
        "window_end": 7,
        "compatible_nucleases": ["SpCas9"],
        "efficiency_class": "moderate",
        "evidence_tier": "A",
        "description": "Alias for ABE7.10 (legacy)",
        "reference": "Gaudelli et al., Nature, 2017",
    },
    "CBE": {
        "editor_type": "CBE",
        "window_start": 4,
        "window_end": 8,
        "compatible_nucleases": ["SpCas9"],
        "efficiency_class": "moderate",
        "evidence_tier": "A",
        "description": "Alias for BE4max (legacy)",
        "reference": "Komor et al., Nature, 2016",
    },
}


# =============================================================================
# EDITOR-NUCLEASE COMPATIBILITY MATRIX
# =============================================================================
# Which editors work with which nucleases, and with what efficiency modifier.
# Efficiency modifier: 1.0 = full activity, 0.7 = reduced, etc.
# These are literature-informed estimates (Tier B evidence).

EDITOR_NUCLEASE_EFFICIENCY = {
    # (editor_name, nuclease) -> relative_efficiency_modifier
    ("ABE7.10", "SpCas9"): 1.0,
    ("ABE8e", "SpCas9"): 1.0,
    ("ABE8e", "SpCas9-NG"): 0.7,    # Reduced due to NG PAM binding
    ("ABE8e", "SpRY"): 0.5,          # SpRY has lower on-target activity
    ("ABE8e", "enFnCas9"): 0.8,      # Inferred; not directly published
    ("BE4max", "SpCas9"): 1.0,
    ("BE4max", "SpCas9-NG"): 0.7,
    ("BE4max", "SpRY"): 0.5,
    ("BE4max", "enFnCas9"): 0.8,
}
