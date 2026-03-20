"""
Editing Strategy Enumerator for MOSAIC
=======================================

Given a gene structure, a set of mutations, and experimental constraints
(cell type, nuclease), this module enumerates every *feasible* multi-locus
editing strategy and annotates each one with predicted efficiency, safety,
and practical parameters.

Why enumerate all strategies?
-----------------------------
When a patient's iPSCs carry two pathogenic mutations in the same gene,
there are many possible ways to correct both:

  * Use one big donor template spanning both sites (if they are close)
  * Use two separate donors in one transfection (simultaneous dual HDR)
  * Correct one mutation, expand clones, then correct the second
    (sequential HDR)
  * Use base editing for one or both (if the mutation types allow it)
  * Use prime editing for one or both
  * A hybrid: base-edit one site, HDR the other
  * Delete the affected exon(s) (if the protein can tolerate it)

Each approach has different efficiency, safety, time, and cost
characteristics.  MOSAIC enumerates all *biologically feasible*
combinations so the user can make an informed decision.

Key biological constraints
--------------------------
1. **Donor size limit**: cssDNA donors are practical up to ~5-10 kb.
   AAV donors are limited to ~4.7 kb including ITRs.  A single donor
   can span two mutation sites only if the genomic distance (including
   introns) between them is within this limit.
   Source: Iyer et al., CRISPR Journal, 2022.

2. **Simultaneous DSBs and translocation risk**: Two DSBs on the same
   chromosome create a risk of interstitial deletion (loss of the
   intervening segment) or translocation.  The probability scales with
   3D proximity (Hi-C contact frequency).
   Source: Chiarle et al., Cell, 2011; Zhang et al., Cell, 2012.

3. **p53 toxicity in iPSCs**: Each DSB activates the p53 pathway,
   leading to apoptosis or G1 arrest.  iPSCs are particularly
   sensitive.  Clones that survive may have acquired p53 mutations,
   compromising their utility.
   Source: Ihry et al., Nature Medicine, 2018.

4. **Base/prime editing availability**: These DSB-free approaches are
   only applicable to certain mutation types (see mutation_classifier).

5. **Efficiency compounding**: When two independent editing events must
   both succeed in the same cell, the overall efficiency is the product
   of the individual efficiencies (assuming independence).  For
   simultaneous approaches: P(both) = P(site1) * P(site2).  For
   sequential approaches: each round is independent, but you need to
   find a correctly edited clone before proceeding to the next round.
"""

from __future__ import annotations

import sys
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

import numpy as np

# Import from sibling modules
from .gene_structure import GeneStructure
from .mutation_classifier import Mutation, MutationClassifier

# Import constants from the toolkit-wide constants module.
# We use a relative import that works whether this is run as part of the
# crisprarchitect package or from the mosaic subpackage.
try:
    from ..utils.constants import (
        CELL_TYPE_PARAMS,
        NUCLEASE_PARAMS,
        DONOR_TOPOLOGY_MULTIPLIER,
        TRANSLOCATION_BASELINE_PROB_1MB,
        HICCONTACT_POWER_LAW_GAMMA,
        BASELINE_HDR_FRACTION_IPSC,
    )
except ImportError:
    # Fallback for direct execution or testing outside package context.
    # These values mirror the constants.py file.
    CELL_TYPE_PARAMS = {
        "iPSC": {
            "hdr_base_efficiency": 0.08,
            "cell_cycle_s_g2_fraction": 0.35,
            "p53_active": True,
            "viability_single_dsb": 0.55,
            "viability_dual_dsb": 0.30,
        },
        "HEK293T": {
            "hdr_base_efficiency": 0.25,
            "cell_cycle_s_g2_fraction": 0.55,
            "p53_active": False,
            "viability_single_dsb": 0.85,
            "viability_dual_dsb": 0.70,
        },
    }
    NUCLEASE_PARAMS = {
        "SpCas9": {"hdr_multiplier": 1.0, "stagger_bp": 0},
        "vCas9": {"hdr_multiplier": 1.9, "stagger_bp": 6},
    }
    DONOR_TOPOLOGY_MULTIPLIER = {
        "circular_ssDNA": 3.0,
        "linear_ssDNA": 1.5,
        "linear_dsDNA": 1.0,
        "AAV_ssDNA": 4.0,
        "plasmid_dsDNA": 0.8,
    }
    TRANSLOCATION_BASELINE_PROB_1MB = 1e-3
    HICCONTACT_POWER_LAW_GAMMA = 1.08
    BASELINE_HDR_FRACTION_IPSC = 0.08


# ===========================================================================
# EditingStrategy dataclass
# ===========================================================================

@dataclass
class EditingStrategy:
    """A fully specified multi-locus editing strategy.

    Each instance represents one *complete plan* for correcting all
    target mutations in a patient's cells.  The ``StrategyEnumerator``
    creates these; the ``StrategyScorer`` ranks them.

    Attributes
    ----------
    name : str
        Short machine-readable strategy identifier (e.g.,
        ``"SEQUENTIAL_HDR"``).
    description : str
        Human-readable description of the strategy.
    mutations_addressed : list of Mutation
        The mutations that this strategy corrects.
    requires_dsb : bool
        Whether any editing step requires a double-strand break.
        ``False`` for pure base-editing or prime-editing strategies.
    num_dsbs : int
        Total number of DSBs across all editing rounds.  0 for DSB-free
        approaches; 1 for single-site HDR or hybrid; 2 for dual HDR.
    num_editing_rounds : int
        Number of separate transfection/electroporation rounds.
        1 for simultaneous approaches; 2 for sequential.
    donor_templates_needed : int
        Number of distinct donor DNA molecules to synthesize.  0 for
        base/prime editing; 1 for single-template HDR; 2 for dual HDR.
    donor_type : str
        Type of donor template (e.g., ``"cssDNA"``, ``"ssODN"``,
        ``"none"``).
    requires_selection : bool
        Whether a selection marker (e.g., puromycin cassette) is
        recommended to enrich for correctly edited clones.
    estimated_efficiency : float
        Predicted fraction of transfected cells that will have ALL
        target mutations corrected on the desired alleles.  This is the
        "bottom line" number: how many clones must you screen to find
        one perfect clone?
    safety_score : float
        0-1 score where 1 = safest.  Accounts for DSB count, p53
        activation, translocation risk, and large-deletion risk.
    time_weeks : int
        Estimated wall-clock time from start of editing to validated
        clone, including clonal expansion and genotyping.
    screening_burden : int
        Estimated number of clones to screen to find one correctly
        edited clone.  Computed as ``ceil(1 / estimated_efficiency)``.
    feasibility_warnings : list of str
        Any biological or practical warnings (e.g., translocation risk,
        donor size exceeded).
    """

    name: str
    description: str
    mutations_addressed: List[Mutation] = field(default_factory=list)
    requires_dsb: bool = False
    num_dsbs: int = 0
    num_editing_rounds: int = 1
    donor_templates_needed: int = 0
    donor_type: str = "none"
    requires_selection: bool = False
    estimated_efficiency: float = 0.0
    safety_score: float = 1.0
    time_weeks: int = 4
    screening_burden: int = 1
    feasibility_warnings: List[str] = field(default_factory=list)


# ===========================================================================
# StrategyEnumerator
# ===========================================================================

class StrategyEnumerator:
    """Enumerate all feasible editing strategies for a multi-locus correction.

    This class is the core decision engine of MOSAIC.  Given the gene
    structure, the patient's mutations, the cell type, and the available
    nuclease, it generates a list of ``EditingStrategy`` objects
    representing every biologically plausible approach to correcting all
    mutations.

    The enumerator does NOT rank strategies — that is the job of
    ``StrategyScorer``.  The enumerator only asks "is this strategy
    physically and biologically feasible?" and, if so, populates it
    with estimated parameters.

    Parameters and sources
    ----------------------
    - Base editing efficiency: 30-60% for CBE, 40-70% for ABE in
      HEK293T; ~50% lower in iPSCs.
      Source: Komor et al., Nature, 2016; Gaudelli et al., Nature, 2017;
      Koblan et al., Nature Biotechnology, 2018.
    - Prime editing efficiency: 5-50%, median ~20% in HEK293T; ~30-50%
      lower in iPSCs.
      Source: Anzalone et al., Nature, 2019; Chen et al., Cell, 2021.
    - HDR efficiency: from CELL_TYPE_PARAMS (constants.py), multiplied
      by nuclease and donor-topology factors.
    - Time per editing round: ~4 weeks (transfection + recovery +
      clonal expansion + genotyping).  This is a practical estimate for
      iPSC workflows.
      Source: Standard iPSC gene-editing protocols (Ran et al., Nature
      Protocols, 2013; Bruntraeger et al., Methods, 2019).
    """

    def __init__(self) -> None:
        self._classifier = MutationClassifier()

    # -----------------------------------------------------------------
    # Main enumeration
    # -----------------------------------------------------------------

    def enumerate_strategies(
        self,
        gene: GeneStructure,
        mutations: List[Mutation],
        cell_type: str = "iPSC",
        nuclease: str = "SpCas9",
    ) -> List[EditingStrategy]:
        """Enumerate all feasible strategies for correcting the given mutations.

        Parameters
        ----------
        gene : GeneStructure
            The gene model containing the mutations.
        mutations : list of Mutation
            The pathogenic mutations to correct (typically 2).
        cell_type : str
            Cell type being edited.  Must be a key in
            ``CELL_TYPE_PARAMS`` (e.g., ``"iPSC"``, ``"HEK293T"``).
        nuclease : str
            Nuclease to use for HDR-based approaches.  Must be a key in
            ``NUCLEASE_PARAMS``.

        Returns
        -------
        list of EditingStrategy
            All feasible strategies, unsorted.

        Raises
        ------
        ValueError
            If ``cell_type`` or ``nuclease`` is not recognized.
        """
        if cell_type not in CELL_TYPE_PARAMS:
            raise ValueError(
                f"Unknown cell type '{cell_type}'. "
                f"Available: {list(CELL_TYPE_PARAMS.keys())}"
            )
        if nuclease not in NUCLEASE_PARAMS:
            raise ValueError(
                f"Unknown nuclease '{nuclease}'. "
                f"Available: {list(NUCLEASE_PARAMS.keys())}"
            )

        ct = CELL_TYPE_PARAMS[cell_type]
        nuc = NUCLEASE_PARAMS[nuclease]

        strategies: List[EditingStrategy] = []

        # --- Classify every mutation ---
        for m in mutations:
            if not m.mutation_type:
                m.mutation_type = self._classifier.classify(
                    m.ref_allele, m.alt_allele
                )

        be_amenable = [self._classifier.base_editing_amenable(m) for m in mutations]
        pe_amenable = [self._classifier.prime_editing_amenable(m) for m in mutations]
        all_be = all(be_amenable)
        all_pe = all(pe_amenable)
        any_be = any(be_amenable)
        any_pe = any(pe_amenable)

        # --- Compute key distances ---
        if len(mutations) >= 2:
            exon_a = mutations[0].exon_number
            exon_b = mutations[1].exon_number
            try:
                genomic_dist = gene.genomic_distance(exon_a, exon_b)
            except ValueError:
                genomic_dist = float("inf")
        else:
            genomic_dist = 0

        # --- HDR efficiency calculation ---
        # Base HDR rate for this cell type, modified by nuclease and donor.
        # Formula:
        #   eff = hdr_base * nuclease_multiplier * donor_multiplier * viability
        # Source: constants.py values derived from multiple publications.
        hdr_base = ct["hdr_base_efficiency"]
        nuc_mult = nuc["hdr_multiplier"]
        donor_mult_cssdna = DONOR_TOPOLOGY_MULTIPLIER.get("circular_ssDNA", 3.0)
        viability_1dsb = ct["viability_single_dsb"]
        viability_2dsb = ct["viability_dual_dsb"]

        single_site_hdr_eff = min(
            hdr_base * nuc_mult * donor_mult_cssdna * viability_1dsb,
            0.95,  # Cap at 95% — no editing is perfect
        )

        # --- Base editing efficiency ---
        # ABE: ~40-70% in HEK293T, ~20-35% in iPSCs
        # CBE: ~30-60% in HEK293T, ~15-30% in iPSCs
        # Source: Koblan et al., Nat Biotech, 2018; Zuo et al., Nat Methods, 2020
        # We use a conservative central estimate scaled by cell type.
        base_edit_eff_hek = 0.45  # Central estimate for HEK293T
        # iPSCs have lower transfection efficiency and different chromatin;
        # scale by the ratio of S/G2 fractions (base editing is cell-cycle
        # independent, but delivery efficiency correlates with growth rate).
        cell_scale = ct["cell_cycle_s_g2_fraction"] / 0.55  # Normalised to HEK293T
        base_edit_eff = min(base_edit_eff_hek * cell_scale, 0.70)

        # --- Prime editing efficiency ---
        # PE2/PE3: ~5-50% in HEK293T, typically ~20% median
        # In iPSCs: typically 30-50% lower than HEK293T
        # Source: Anzalone et al., Nature, 2019; Chen et al., Cell, 2021
        prime_edit_eff_hek = 0.20  # Conservative central estimate
        prime_edit_eff = min(prime_edit_eff_hek * cell_scale, 0.40)

        # --- Time constants ---
        # One editing round in iPSCs: ~4 weeks
        #   1 week: electroporation + recovery
        #   1 week: clonal seeding + expansion
        #   2 weeks: genotyping + validation
        # Source: Standard iPSC protocols (Bruntraeger et al., Methods, 2019)
        weeks_per_round = 4

        # --- Translocation risk ---
        # For two simultaneous DSBs on the same chromosome separated by
        # `d` bp, translocation probability approximated as:
        #   P_trans ~ P_baseline * (d / 1e6)^(-gamma)
        # But for same-chromosome loci, the more relevant risk is
        # interstitial deletion (loss of segment between the two cuts).
        # Source: Chiarle et al., Cell, 2011
        if genomic_dist > 0 and genomic_dist < float("inf"):
            # Probability scales inversely with distance on same chromosome.
            # For nearby loci (<1 Mb), the interstitial deletion risk is
            # significant (~1-10% of DSB outcomes).
            # Source: Kraft et al., Cell Stem Cell, 2015
            trans_prob = min(
                TRANSLOCATION_BASELINE_PROB_1MB
                * (genomic_dist / 1e6) ** (-HICCONTACT_POWER_LAW_GAMMA),
                0.15,  # Cap at 15%
            )
        else:
            trans_prob = 0.0

        # =================================================================
        # Strategy (a): SINGLE_TEMPLATE_HDR
        # =================================================================
        # Only feasible if both mutations can be spanned by a single donor.
        # The donor must include the genomic sequence (with introns)
        # between the two mutation sites, plus homology arms on each end.
        #
        # Practical limit for cssDNA: ~5-10 kb total insert.
        # We use 5 kb as the default maximum (conservative).
        # Source: Iyer et al., CRISPR Journal, 2022.

        MAX_SINGLE_DONOR_BP = 5_000  # Conservative cssDNA limit

        if len(mutations) >= 2 and genomic_dist <= MAX_SINGLE_DONOR_BP:
            # A single cssDNA donor can span both sites!
            # Only one DSB is needed (between or near the two mutations).
            eff = single_site_hdr_eff  # Only one DSB event
            screening = max(1, int(np.ceil(1.0 / max(eff, 1e-6))))

            strategies.append(EditingStrategy(
                name="SINGLE_TEMPLATE_HDR",
                description=(
                    "A single cssDNA donor template spans both mutation "
                    f"sites (genomic distance: {genomic_dist:,} bp). "
                    "Only one DSB is needed. The donor encodes both "
                    "corrections plus homology arms. This is the most "
                    "efficient HDR approach when feasible."
                ),
                mutations_addressed=list(mutations),
                requires_dsb=True,
                num_dsbs=1,
                num_editing_rounds=1,
                donor_templates_needed=1,
                donor_type="cssDNA",
                requires_selection=False,
                estimated_efficiency=eff,
                safety_score=0.5,  # One DSB: moderate safety
                time_weeks=weeks_per_round,
                screening_burden=screening,
                feasibility_warnings=[
                    f"Genomic distance ({genomic_dist:,} bp) must fit "
                    f"within cssDNA synthesis limit (~{MAX_SINGLE_DONOR_BP:,} bp)."
                ],
            ))

        # =================================================================
        # Strategy (b): DUAL_TEMPLATE_SIMULTANEOUS_HDR
        # =================================================================
        # Two separate donors + two sgRNAs delivered in one transfection.
        # Both DSBs happen simultaneously.
        #
        # Efficiency: product of individual HDR efficiencies.
        # Key risk: two simultaneous DSBs on the same chromosome can
        # cause interstitial deletion or translocation.

        if len(mutations) >= 2:
            # Each site corrected independently via HDR
            eff_per_site = single_site_hdr_eff
            # Joint efficiency = product (independent events)
            eff = eff_per_site ** 2
            screening = max(1, int(np.ceil(1.0 / max(eff, 1e-6))))

            warnings = []
            if genomic_dist < float("inf"):
                warnings.append(
                    f"Two simultaneous DSBs on {gene.chromosome} separated "
                    f"by {genomic_dist:,} bp. Estimated interstitial "
                    f"deletion / translocation risk: {trans_prob:.1%}."
                )
            if ct.get("p53_active", False):
                warnings.append(
                    "Two simultaneous DSBs in iPSCs will strongly activate "
                    "p53, causing significant cell death and potential "
                    "selection for p53-mutant clones "
                    "(Ihry et al., Nature Medicine, 2018)."
                )

            # Safety: penalised for 2 DSBs + translocation risk
            safety = max(0.1, 0.5 - 0.2 * (trans_prob / 0.05))
            if ct.get("p53_active", False):
                safety *= 0.6  # Extra penalty for p53-active cells

            strategies.append(EditingStrategy(
                name="DUAL_TEMPLATE_SIMULTANEOUS_HDR",
                description=(
                    "Two cssDNA donors and two sgRNAs delivered "
                    "simultaneously in one electroporation. Both sites "
                    "are cut and repaired by HDR in the same cell cycle. "
                    f"Per-site HDR efficiency: {eff_per_site:.1%}. "
                    f"Joint efficiency: {eff:.2%}."
                ),
                mutations_addressed=list(mutations),
                requires_dsb=True,
                num_dsbs=2,
                num_editing_rounds=1,
                donor_templates_needed=2,
                donor_type="cssDNA",
                requires_selection=True,
                estimated_efficiency=eff,
                safety_score=round(safety, 3),
                time_weeks=weeks_per_round,
                screening_burden=screening,
                feasibility_warnings=warnings,
            ))

        # =================================================================
        # Strategy (c): SEQUENTIAL_HDR
        # =================================================================
        # Correct one mutation, expand and genotype clones, then correct
        # the second mutation in a validated clone.
        #
        # Pros: avoids simultaneous DSBs -> no translocation risk
        # Cons: takes twice as long; two rounds of clonal expansion

        if len(mutations) >= 2:
            eff_per_round = single_site_hdr_eff
            # For sequential editing, we need to succeed in round 1,
            # then succeed again in round 2.  The overall probability of
            # getting a doubly-corrected clone is:
            #   P = P(round1) * P(round2)
            # But we screen after each round, so the practical efficiency
            # is determined by the worse of the two rounds.
            eff = eff_per_round  # Per-round screening efficiency
            screening_per_round = max(
                1, int(np.ceil(1.0 / max(eff_per_round, 1e-6)))
            )

            # Safety: only one DSB at a time, no translocation risk
            safety = 0.6  # Single DSB: moderate-good safety
            if ct.get("p53_active", False):
                safety = 0.5  # p53 activation still a concern per round

            strategies.append(EditingStrategy(
                name="SEQUENTIAL_HDR",
                description=(
                    "Correct one mutation per editing round. Round 1: "
                    f"correct mutation in exon {mutations[0].exon_number}. "
                    "Screen, validate, expand a correct clone. Round 2: "
                    f"correct mutation in exon {mutations[1].exon_number}. "
                    "Avoids simultaneous DSBs and translocation risk. "
                    f"Efficiency per round: {eff_per_round:.1%}."
                ),
                mutations_addressed=list(mutations),
                requires_dsb=True,
                num_dsbs=1,  # Per round (total 2 but never simultaneous)
                num_editing_rounds=2,
                donor_templates_needed=2,
                donor_type="cssDNA",
                requires_selection=False,
                estimated_efficiency=eff,
                safety_score=round(safety, 3),
                time_weeks=2 * weeks_per_round,
                screening_burden=screening_per_round,
                feasibility_warnings=[
                    "Total time is doubled due to two sequential rounds. "
                    "Each round includes clonal expansion and genotyping."
                ],
            ))

        # =================================================================
        # Strategy (d): DUAL_BASE_EDITING
        # =================================================================
        # If BOTH mutations are base-editable, we can correct both
        # without any DSBs.  This is the safest approach.
        #
        # Both base edits can be delivered simultaneously (two sgRNAs +
        # one base editor) or sequentially.  Simultaneous is preferred
        # because base editing does not cause DSBs, so there is no
        # translocation risk.
        #
        # Source: Koblan et al., Nat Biotech, 2018 (multiplexed base editing)

        if all_be and len(mutations) >= 2:
            # Simultaneous dual base editing
            eff = base_edit_eff ** 2  # Product of independent events
            screening = max(1, int(np.ceil(1.0 / max(eff, 1e-6))))

            strategies.append(EditingStrategy(
                name="DUAL_BASE_EDITING",
                description=(
                    "Both mutations are transitions amenable to base "
                    "editing. Deliver two sgRNAs with a single base "
                    "editor (ABE or CBE as appropriate) in one "
                    "electroporation. NO double-strand breaks are "
                    "introduced, making this the safest approach for "
                    f"iPSCs. Per-site efficiency: ~{base_edit_eff:.0%}. "
                    f"Joint efficiency: {eff:.1%}."
                ),
                mutations_addressed=list(mutations),
                requires_dsb=False,
                num_dsbs=0,
                num_editing_rounds=1,
                donor_templates_needed=0,
                donor_type="none",
                requires_selection=False,
                estimated_efficiency=eff,
                safety_score=0.95,  # Highest safety: no DSBs at all
                time_weeks=weeks_per_round,
                screening_burden=screening,
                feasibility_warnings=[
                    "Requires suitable PAM sites within the base-editing "
                    "window (positions 4-8 of the protospacer) at both "
                    "mutation sites. Verify computationally before "
                    "proceeding.",
                    "If both mutations require the SAME type of base editor "
                    "(both ABE or both CBE), one editor protein suffices. "
                    "If they require DIFFERENT editors (one ABE, one CBE), "
                    "you will need to deliver both editors or edit "
                    "sequentially.",
                ],
            ))

        # =================================================================
        # Strategy (e): HYBRID_BASE_EDIT_PLUS_HDR
        # =================================================================
        # If one mutation is base-editable but the other is not, use
        # base editing for the amenable site and HDR for the other.
        # This results in only ONE DSB total (at the HDR site).

        if any_be and not all_be and len(mutations) >= 2:
            # Figure out which mutation gets base editing and which gets HDR
            be_idx = next(i for i, v in enumerate(be_amenable) if v)
            hdr_idx = next(i for i, v in enumerate(be_amenable) if not v)
            be_mut = mutations[be_idx]
            hdr_mut = mutations[hdr_idx]

            # These can be done simultaneously or sequentially.
            # Simultaneous: one base editor + one nuclease + one donor
            # in one electroporation.  Only one DSB (at the HDR site).
            eff = base_edit_eff * single_site_hdr_eff
            screening = max(1, int(np.ceil(1.0 / max(eff, 1e-6))))

            safety = 0.6  # One DSB: moderate
            if ct.get("p53_active", False):
                safety = 0.5

            strategies.append(EditingStrategy(
                name="HYBRID_BASE_EDIT_PLUS_HDR",
                description=(
                    f"Base-edit the mutation in exon {be_mut.exon_number} "
                    f"({be_mut.mutation_type}) using "
                    f"{'ABE' if be_mut.mutation_type == 'transition_AG' else 'CBE'}, "
                    f"and correct the mutation in exon {hdr_mut.exon_number} "
                    f"({hdr_mut.mutation_type}) via HDR with a cssDNA donor. "
                    "Only one DSB is introduced (at the HDR site). "
                    "Good compromise between efficiency and safety."
                ),
                mutations_addressed=list(mutations),
                requires_dsb=True,
                num_dsbs=1,
                num_editing_rounds=1,
                donor_templates_needed=1,
                donor_type="cssDNA",
                requires_selection=False,
                estimated_efficiency=eff,
                safety_score=round(safety, 3),
                time_weeks=weeks_per_round,
                screening_burden=screening,
                feasibility_warnings=[
                    "Requires co-delivery of a base editor AND a nuclease "
                    "(+ donor). This is technically feasible via "
                    "electroporation of RNPs + mRNA, but increases the "
                    "complexity of the transfection.",
                    "Consider sequential approach (base edit first, HDR "
                    "second) to avoid any p53-related complications "
                    "during the base-editing step.",
                ],
            ))

        # =================================================================
        # Strategy (f): PRIME_EDITING (dual)
        # =================================================================
        # If both mutations are prime-editable, we can correct both
        # without DSBs.  Prime editing is less efficient than base
        # editing but can handle a wider range of mutations.

        if all_pe and len(mutations) >= 2:
            eff = prime_edit_eff ** 2
            screening = max(1, int(np.ceil(1.0 / max(eff, 1e-6))))

            strategies.append(EditingStrategy(
                name="DUAL_PRIME_EDITING",
                description=(
                    "Both mutations can be corrected by prime editing "
                    "(PE3 or PE5). Two pegRNAs + prime editor delivered "
                    "simultaneously. NO double-strand breaks. Lower "
                    "efficiency than base editing but applicable to a "
                    "wider range of mutation types (substitutions, small "
                    f"indels). Per-site efficiency: ~{prime_edit_eff:.0%}. "
                    f"Joint efficiency: {eff:.1%}."
                ),
                mutations_addressed=list(mutations),
                requires_dsb=False,
                num_dsbs=0,
                num_editing_rounds=1,
                donor_templates_needed=0,
                donor_type="none",
                requires_selection=False,
                estimated_efficiency=eff,
                safety_score=0.90,  # Very safe: no DSBs
                time_weeks=weeks_per_round,
                screening_burden=screening,
                feasibility_warnings=[
                    "Prime editing efficiency is highly locus-dependent. "
                    "The estimates here are based on median published "
                    "values. Validate at your specific target sites.",
                    "Consider PE5 (with MLH1dn) for improved efficiency "
                    "at the cost of slightly higher indel rates "
                    "(Chen et al., Cell, 2021).",
                ],
            ))

        # =================================================================
        # Strategy (g): EXON_DELETION
        # =================================================================
        # If the affected exons can be deleted without destroying protein
        # function (e.g., in-frame exon skipping), we can use two
        # sgRNAs flanking the exon(s) to create a deletion via NHEJ.
        # No donor template is needed.
        #
        # This is only applicable to a small subset of genes/exons where
        # the protein can tolerate the loss.  MOSAIC flags this as an
        # option but warns that functional validation is required.

        if len(mutations) >= 2:
            exon_nums = sorted(set(m.exon_number for m in mutations))
            total_exonic_to_delete = gene.exonic_distance(
                min(exon_nums), max(exon_nums)
            )
            # Exon deletion is only practical if:
            #   1. The total exonic sequence lost is small (<= ~1000 bp,
            #      ~330 aa) — larger deletions are unlikely to preserve
            #      protein function.
            #   2. The deletion maintains the reading frame (must be
            #      divisible by 3). We check modulo 3 here.
            in_frame = (total_exonic_to_delete % 3 == 0)

            eff = 0.10  # NHEJ-based deletion efficiency estimate
            screening = max(1, int(np.ceil(1.0 / max(eff, 1e-6))))

            warnings = [
                "EXON DELETION REQUIRES FUNCTIONAL VALIDATION. The protein "
                "must tolerate loss of the deleted exon(s). This must be "
                "confirmed experimentally (e.g., via minigene assay or "
                "patient phenotype data).",
            ]
            if not in_frame:
                warnings.append(
                    f"WARNING: Deleting exons {min(exon_nums)}-"
                    f"{max(exon_nums)} removes {total_exonic_to_delete} bp "
                    f"of exonic sequence, which is NOT divisible by 3. "
                    f"This will cause a frameshift and likely a "
                    f"non-functional protein."
                )
            if total_exonic_to_delete > 1000:
                warnings.append(
                    f"Large exonic deletion ({total_exonic_to_delete} bp = "
                    f"~{total_exonic_to_delete // 3} amino acids). "
                    f"Protein function is unlikely to be preserved."
                )

            strategies.append(EditingStrategy(
                name="EXON_DELETION",
                description=(
                    f"Delete exon(s) {min(exon_nums)}-{max(exon_nums)} "
                    f"using two sgRNAs flanking the target region. "
                    f"NHEJ-mediated deletion removes "
                    f"{total_exonic_to_delete} bp of exonic sequence. "
                    f"No donor template needed. "
                    f"Reading frame {'PRESERVED' if in_frame else 'DISRUPTED'}."
                ),
                mutations_addressed=list(mutations),
                requires_dsb=True,
                num_dsbs=2,
                num_editing_rounds=1,
                donor_templates_needed=0,
                donor_type="none",
                requires_selection=False,
                estimated_efficiency=eff,
                # Safety score: 2 DSBs but no donor complexity
                safety_score=0.3 if ct.get("p53_active") else 0.5,
                time_weeks=weeks_per_round,
                screening_burden=screening,
                feasibility_warnings=warnings,
            ))

        # =================================================================
        # Strategy (h): SEQUENTIAL_BASE_AND_PRIME
        # =================================================================
        # If one mutation is base-editable and the other is prime-editable
        # (but not base-editable), we can do both without any DSBs by
        # combining the two DSB-free approaches in sequence.
        #
        # Why sequential rather than simultaneous?  Because the base
        # editor and prime editor are different proteins — delivering
        # both simultaneously increases complexity and may cause
        # interference (e.g., the nCas9-RT of PE and the nCas9-deaminase
        # of BE competing for the same Cas9 scaffold).

        # Check: at least one mutation is BE-amenable, at least one is
        # PE-amenable (but not BE-amenable).
        if len(mutations) >= 2:
            be_only = [
                i for i in range(len(mutations))
                if be_amenable[i]
            ]
            pe_only_not_be = [
                i for i in range(len(mutations))
                if pe_amenable[i] and not be_amenable[i]
            ]

            if be_only and pe_only_not_be:
                # We have at least one mutation for BE and one for PE
                be_mut = mutations[be_only[0]]
                pe_mut = mutations[pe_only_not_be[0]]

                # Sequential: two rounds, each with its own efficiency
                eff_be = base_edit_eff
                eff_pe = prime_edit_eff
                screening_be = max(1, int(np.ceil(1.0 / max(eff_be, 1e-6))))
                screening_pe = max(1, int(np.ceil(1.0 / max(eff_pe, 1e-6))))

                strategies.append(EditingStrategy(
                    name="SEQUENTIAL_BASE_AND_PRIME",
                    description=(
                        "Two sequential DSB-free editing rounds. "
                        f"Round 1: Base-edit the mutation in exon "
                        f"{be_mut.exon_number} ({be_mut.mutation_type}) "
                        f"using {'ABE' if be_mut.mutation_type == 'transition_AG' else 'CBE'}. "
                        f"Round 2: Prime-edit the mutation in exon "
                        f"{pe_mut.exon_number} ({pe_mut.mutation_type}). "
                        "NO DSBs in either round. Safest sequential "
                        "approach, but requires more time."
                    ),
                    mutations_addressed=list(mutations),
                    requires_dsb=False,
                    num_dsbs=0,
                    num_editing_rounds=2,
                    donor_templates_needed=0,
                    donor_type="none",
                    requires_selection=False,
                    # Per-round efficiency is the bottleneck
                    estimated_efficiency=min(eff_be, eff_pe),
                    safety_score=0.92,  # Very safe: no DSBs
                    time_weeks=2 * weeks_per_round,
                    screening_burden=max(screening_be, screening_pe),
                    feasibility_warnings=[
                        "Requires two separate transfections with different "
                        "editor proteins. Plan for 8+ weeks total.",
                        "Base editing round should be performed first (higher "
                        "efficiency, less screening).",
                    ],
                ))

        # =================================================================
        # Strategy: HYBRID_PRIME_EDIT_PLUS_HDR
        # =================================================================
        # If one mutation is prime-editable but not base-editable, and
        # the other requires HDR, we can prime-edit one and HDR the other.

        if len(mutations) >= 2:
            pe_not_be = [
                i for i in range(len(mutations))
                if pe_amenable[i] and not be_amenable[i]
            ]
            needs_hdr = [
                i for i in range(len(mutations))
                if not pe_amenable[i] and not be_amenable[i]
            ]

            if pe_not_be and needs_hdr:
                pe_mut = mutations[pe_not_be[0]]
                hdr_mut = mutations[needs_hdr[0]]

                # Sequential recommended: prime edit first (no DSB),
                # then HDR second.
                eff = min(prime_edit_eff, single_site_hdr_eff)
                screening = max(
                    1,
                    int(np.ceil(1.0 / max(prime_edit_eff, 1e-6))),
                    int(np.ceil(1.0 / max(single_site_hdr_eff, 1e-6))),
                )

                safety = 0.55  # One DSB in one of the two rounds
                if ct.get("p53_active", False):
                    safety = 0.45

                strategies.append(EditingStrategy(
                    name="HYBRID_PRIME_EDIT_PLUS_HDR",
                    description=(
                        f"Prime-edit the mutation in exon "
                        f"{pe_mut.exon_number} (round 1, no DSB), then "
                        f"HDR-correct the mutation in exon "
                        f"{hdr_mut.exon_number} (round 2, one DSB). "
                        f"Only one DSB total across both rounds."
                    ),
                    mutations_addressed=list(mutations),
                    requires_dsb=True,
                    num_dsbs=1,
                    num_editing_rounds=2,
                    donor_templates_needed=1,
                    donor_type="cssDNA",
                    requires_selection=False,
                    estimated_efficiency=eff,
                    safety_score=round(safety, 3),
                    time_weeks=2 * weeks_per_round,
                    screening_burden=screening,
                    feasibility_warnings=[
                        "Perform prime editing first to avoid p53 "
                        "activation from the HDR DSB affecting the "
                        "prime-editing round.",
                    ],
                ))

        return strategies

    def summarize_strategies(
        self, strategies: List[EditingStrategy]
    ) -> str:
        """Return a concise text summary of all enumerated strategies.

        Parameters
        ----------
        strategies : list of EditingStrategy
            Strategies to summarize.

        Returns
        -------
        str
            Formatted summary string.
        """
        lines = [
            f"Enumerated {len(strategies)} feasible strategies:",
            f"{'=' * 60}",
        ]
        for i, s in enumerate(strategies, 1):
            dsb_label = f"{s.num_dsbs} DSB(s)" if s.requires_dsb else "No DSBs"
            lines.append(
                f"  {i}. {s.name:40s} | {dsb_label:10s} | "
                f"Eff: {s.estimated_efficiency:.2%} | "
                f"Safety: {s.safety_score:.2f} | "
                f"{s.time_weeks} wk"
            )
        return "\n".join(lines)


# =========================================================================
# Quick self-test
# =========================================================================

if __name__ == "__main__":
    from .gene_structure import build_demo_gene
    from .mutation_classifier import Mutation

    gene = build_demo_gene()
    m1 = Mutation(exon_number=40, position=5_000_000, ref_allele="G",
                  alt_allele="A", name="c.5000G>A")
    m2 = Mutation(exon_number=60, position=8_000_000, ref_allele="CT",
                  alt_allele="-", name="c.8000delCT")

    enumerator = StrategyEnumerator()
    strategies = enumerator.enumerate_strategies(
        gene, [m1, m2], cell_type="iPSC", nuclease="SpCas9"
    )
    print(enumerator.summarize_strategies(strategies))
