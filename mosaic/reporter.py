"""
Strategy Report Generator for MOSAIC
======================================

This module produces a comprehensive, human-readable report that
summarises the gene structure, mutation analysis, and ranked editing
strategies.  The report is designed to be read by a molecular biologist
who needs to make an informed decision about which multi-locus editing
approach to pursue.

The report includes:
  1. Gene structure summary (exon/intron architecture)
  2. Mutation analysis (classification, editor compatibility)
  3. Ranked strategies with scores, pros, and cons
  4. Recommended strategy with detailed rationale
  5. Experimental design suggestions
  6. Safety warnings and quality control recommendations
  7. Literature references for key claims
"""

from __future__ import annotations

import textwrap
from datetime import datetime
from typing import Dict, List, Optional

import numpy as np

from .gene_structure import GeneStructure
from .mutation_classifier import Mutation, MutationClassifier
from .scorer import ScoredStrategy, StrategyScorer
from .strategy_enumerator import EditingStrategy, StrategyEnumerator


class StrategyReporter:
    """Generate formatted text reports for MOSAIC analysis.

    The reporter is stateless — it receives all data as arguments and
    produces a formatted string.  This makes it easy to integrate into
    notebooks, scripts, or web applications.

    Usage
    -----
    >>> reporter = StrategyReporter()
    >>> report = reporter.generate_report(
    ...     gene=gene,
    ...     mutations=[m1, m2],
    ...     ranked_strategies=ranked,
    ...     cell_type="iPSC",
    ...     nuclease="SpCas9",
    ... )
    >>> print(report)
    """

    def __init__(self, line_width: int = 78) -> None:
        """Initialise the reporter.

        Parameters
        ----------
        line_width : int
            Maximum line width for text wrapping (default 78).
        """
        self.line_width = line_width
        self._classifier = MutationClassifier()

    # -----------------------------------------------------------------
    # Main report generation
    # -----------------------------------------------------------------

    def generate_report(
        self,
        gene: GeneStructure,
        mutations: List[Mutation],
        ranked_strategies: List[ScoredStrategy],
        cell_type: str = "iPSC",
        nuclease: str = "SpCas9",
    ) -> str:
        """Generate the full MOSAIC decision report.

        Parameters
        ----------
        gene : GeneStructure
            Gene model.
        mutations : list of Mutation
            Patient mutations.
        ranked_strategies : list of ScoredStrategy
            Strategies scored and ranked by ``StrategyScorer``.
        cell_type : str
            Cell type being edited.
        nuclease : str
            Nuclease platform.

        Returns
        -------
        str
            Complete formatted report.
        """
        sections = [
            self._header(cell_type, nuclease),
            self._gene_section(gene, mutations),
            self._mutation_section(mutations),
            self._distance_section(gene, mutations),
            self._ranking_section(ranked_strategies),
            self._recommendation_section(ranked_strategies, cell_type),
            self._experimental_design_section(
                ranked_strategies, gene, mutations, cell_type, nuclease
            ),
            self._warnings_section(ranked_strategies, cell_type),
            self._references_section(),
            self._footer(),
        ]
        return "\n".join(sections)

    # -----------------------------------------------------------------
    # Individual sections
    # -----------------------------------------------------------------

    def _header(self, cell_type: str, nuclease: str) -> str:
        """Report header with title and metadata."""
        now = datetime.now().strftime("%Y-%m-%d %H:%M")
        border = "=" * self.line_width
        return "\n".join([
            border,
            "MOSAIC: Multi-locus Optimized Strategy for Allele-specific "
            "Integrated Correction",
            "CRISPRArchitect Toolkit — Strategy Decision Report",
            border,
            f"Generated:   {now}",
            f"Cell type:   {cell_type}",
            f"Nuclease:    {nuclease}",
            border,
            "",
        ])

    def _gene_section(
        self, gene: GeneStructure, mutations: List[Mutation]
    ) -> str:
        """Gene structure summary section."""
        lines = [
            "SECTION 1: GENE STRUCTURE",
            "-" * 40,
            gene.summary(),
            "",
        ]
        return "\n".join(lines)

    def _mutation_section(self, mutations: List[Mutation]) -> str:
        """Mutation analysis section."""
        lines = [
            "SECTION 2: MUTATION ANALYSIS",
            "-" * 40,
        ]

        for i, m in enumerate(mutations, 1):
            lines.append(f"Mutation {i}: {m.name}")
            lines.append(self._classifier.describe_mutation(m))
            lines.append("")

        # Summary table
        lines.append("Correction strategy compatibility matrix:")
        lines.append(
            f"  {'Mutation':<30} {'Base Edit':>10} {'Prime Edit':>12} "
            f"{'HDR Required':>13}"
        )
        lines.append("  " + "-" * 67)
        for m in mutations:
            be = "Yes" if self._classifier.base_editing_amenable(m) else "No"
            pe = "Yes" if self._classifier.prime_editing_amenable(m) else "No"
            hdr = "Yes" if self._classifier.hdr_required(m) else "No"
            lines.append(
                f"  {m.name:<30} {be:>10} {pe:>12} {hdr:>13}"
            )
        lines.append("")
        return "\n".join(lines)

    def _distance_section(
        self, gene: GeneStructure, mutations: List[Mutation]
    ) -> str:
        """Genomic distance analysis between mutation sites."""
        lines = [
            "SECTION 3: INTER-SITE DISTANCE ANALYSIS",
            "-" * 40,
        ]

        if len(mutations) < 2:
            lines.append("Only one mutation — distance analysis not applicable.")
            lines.append("")
            return "\n".join(lines)

        exon_a = mutations[0].exon_number
        exon_b = mutations[1].exon_number

        try:
            gd = gene.genomic_distance(exon_a, exon_b)
            ed = gene.exonic_distance(exon_a, exon_b)
            intd = gene.intronic_distance(exon_a, exon_b)
            can_span = gene.can_span_with_single_donor(exon_a, exon_b)
            can_span_10k = gene.can_span_with_single_donor(
                exon_a, exon_b, max_donor_size=10_000
            )
        except ValueError as e:
            lines.append(f"ERROR computing distances: {e}")
            lines.append("")
            return "\n".join(lines)

        lines.extend([
            f"Mutation sites: exon {exon_a} and exon {exon_b}",
            f"Exons between:  {abs(exon_b - exon_a) - 1} intervening exons",
            f"",
            f"Genomic distance (including introns): {gd:>12,} bp "
            f"({gd / 1000:.1f} kb)",
            f"Exonic distance (coding only):        {ed:>12,} bp "
            f"({ed / 1000:.1f} kb)",
            f"Intronic distance (non-coding):       {intd:>12,} bp "
            f"({intd / 1000:.1f} kb)",
            f"Intronic fraction:                    "
            f"{intd / gd * 100 if gd > 0 else 0:.1f}%",
            f"",
            f"Single cssDNA donor (5 kb limit):     "
            f"{'FEASIBLE' if can_span else 'NOT FEASIBLE'}",
            f"Single cssDNA donor (10 kb limit):    "
            f"{'FEASIBLE' if can_span_10k else 'NOT FEASIBLE'}",
        ])

        if not can_span:
            lines.extend([
                f"",
                f"  -> The two mutation sites are too far apart on the ",
                f"     chromosome for a single donor template. The introns ",
                f"     between exon {exon_a} and exon {exon_b} inflate the ",
                f"     genomic distance to {gd:,} bp, far beyond the ",
                f"     practical limit for cssDNA (~5-10 kb).",
                f"  -> Each site must be addressed independently.",
            ])

        lines.append("")
        return "\n".join(lines)

    def _ranking_section(self, ranked: List[ScoredStrategy]) -> str:
        """Strategy ranking table."""
        lines = [
            "SECTION 4: STRATEGY RANKING",
            "-" * 40,
            f"{'Rank':<5} {'Strategy':<38} {'Overall':>8} "
            f"{'Safety':>7} {'Effic.':>7} {'Time':>5} {'Cost':>5}",
            "-" * self.line_width,
        ]

        for ss in ranked:
            s = ss.strategy
            lines.append(
                f"{ss.rank:<5} {s.name:<38} "
                f"{ss.overall_score:>8.3f} "
                f"{ss.safety_score:>7.3f} "
                f"{ss.efficiency_score:>7.3f} "
                f"{ss.time_score:>5.2f} "
                f"{ss.cost_score:>5.2f}"
            )

        lines.append("")

        # Detailed breakdown for each strategy
        lines.append("DETAILED STRATEGY DESCRIPTIONS:")
        lines.append("")

        for ss in ranked:
            s = ss.strategy
            lines.append(f"  #{ss.rank}: {s.name}")
            # Wrap description text
            wrapped = textwrap.fill(
                s.description,
                width=self.line_width - 6,
                initial_indent="      ",
                subsequent_indent="      ",
            )
            lines.append(wrapped)
            lines.append(
                f"      DSBs: {s.num_dsbs} | "
                f"Rounds: {s.num_editing_rounds} | "
                f"Donors: {s.donor_templates_needed} "
                f"({s.donor_type}) | "
                f"Screening: ~{s.screening_burden} clones"
            )
            if s.feasibility_warnings:
                lines.append("      Warnings:")
                for w in s.feasibility_warnings:
                    w_wrapped = textwrap.fill(
                        w,
                        width=self.line_width - 10,
                        initial_indent="        - ",
                        subsequent_indent="          ",
                    )
                    lines.append(w_wrapped)
            lines.append("")

        return "\n".join(lines)

    def _recommendation_section(
        self, ranked: List[ScoredStrategy], cell_type: str
    ) -> str:
        """Recommended strategy with detailed rationale."""
        lines = [
            "SECTION 5: RECOMMENDATION",
            "-" * 40,
        ]

        if not ranked:
            lines.append("No feasible strategies were identified.")
            lines.append("")
            return "\n".join(lines)

        best = ranked[0]
        s = best.strategy

        lines.extend([
            f"RECOMMENDED STRATEGY: {s.name}",
            f"Overall score: {best.overall_score:.3f} "
            f"(Safety: {best.safety_score:.3f}, "
            f"Efficiency: {best.efficiency_score:.3f}, "
            f"Time: {best.time_score:.2f}, "
            f"Cost: {best.cost_score:.2f})",
            "",
        ])

        # Rationale
        lines.append("Rationale:")

        if s.num_dsbs == 0:
            rationale = (
                "This strategy introduces NO double-strand breaks, making it "
                "the safest option for iPSCs. DSB-free approaches avoid the "
                "p53-mediated toxicity that causes massive cell death and "
                "potential selection for p53-mutant clones in iPSCs "
                "(Ihry et al., Nature Medicine, 2018). They also eliminate "
                "the risk of large deletions at the cut site (Kosicki et al., "
                "Nature Biotechnology, 2018) and chromosomal translocations."
            )
        elif s.num_dsbs == 1:
            rationale = (
                "This strategy introduces only ONE double-strand break, "
                "which is a reasonable safety profile for iPSC editing. "
                "While a single DSB does activate p53 and carries some risk "
                "of large deletions at the cut site, these risks are "
                "manageable with proper quality control (karyotyping, "
                "p53 sequencing of validated clones)."
            )
        else:
            rationale = (
                "This strategy requires multiple DSBs. While it may be the "
                "most efficient approach in terms of timeline, extra care "
                "must be taken to validate clones for karyotypic stability, "
                "absence of translocations, and intact p53 pathway."
            )

        wrapped = textwrap.fill(
            rationale,
            width=self.line_width - 2,
            initial_indent="  ",
            subsequent_indent="  ",
        )
        lines.append(wrapped)
        lines.append("")

        # Alternative strategies worth considering
        if len(ranked) > 1:
            second = ranked[1]
            score_diff = best.overall_score - second.overall_score
            if score_diff < 0.10:
                lines.append("NOTE: The second-ranked strategy is close:")
                lines.append(
                    f"  #{second.rank}: {second.strategy.name} "
                    f"(score: {second.overall_score:.3f}, "
                    f"delta: {score_diff:.3f})"
                )
                lines.append(
                    "  Consider both options based on your lab's specific "
                    "expertise and available reagents."
                )
                lines.append("")

        return "\n".join(lines)

    def _experimental_design_section(
        self,
        ranked: List[ScoredStrategy],
        gene: GeneStructure,
        mutations: List[Mutation],
        cell_type: str,
        nuclease: str,
    ) -> str:
        """Experimental design suggestions for the recommended strategy."""
        lines = [
            "SECTION 6: EXPERIMENTAL DESIGN SUGGESTIONS",
            "-" * 40,
        ]

        if not ranked:
            lines.append("No strategies to design for.")
            lines.append("")
            return "\n".join(lines)

        best = ranked[0]
        s = best.strategy

        # --- sgRNA design ---
        lines.append("sgRNA / Guide Design:")
        if "BASE" in s.name:
            lines.extend([
                "  - For base editing: the target base must fall within the",
                "    editing window (positions 4-8 of the protospacer, counting",
                "    from the PAM-distal end). Use tools like BE-Designer",
                "    (http://www.rgenome.net/be-designer/) to identify suitable",
                "    guide RNAs.",
                "  - Check for bystander edits: other C or A bases within the",
                "    editing window may also be deaminated.",
                "    Source: Komor et al., Nature, 2016; Gaudelli et al., ",
                "    Nature, 2017.",
            ])
        if "PRIME" in s.name:
            lines.extend([
                "  - For prime editing: design pegRNAs with the PrimeDesign",
                "    tool (https://primedesign.pinellolab.partners.org/).",
                "  - Optimal PBS length: 10-16 nt; RTT length depends on edit",
                "    (typically 10-30 nt for point mutations).",
                "  - Consider PE3 nicking guide design for improved efficiency.",
                "    Source: Anzalone et al., Nature, 2019.",
            ])
        if "HDR" in s.name:
            lines.extend([
                f"  - For HDR with {nuclease}: design sgRNAs that cut as close",
                "    as possible to the mutation site (ideally within 10 bp).",
                "    HDR efficiency drops sharply with distance from the cut",
                "    to the edit (Paquet et al., Nature, 2016).",
                f"  - Check for off-targets using Cas-OFFinder or CRISPOR.",
            ])
        lines.append("")

        # --- Donor design ---
        if s.donor_templates_needed > 0:
            lines.append(f"Donor Template Design ({s.donor_type}):")
            if s.donor_type == "cssDNA":
                lines.extend([
                    "  - cssDNA (circular single-stranded DNA) provides ~3x",
                    "    higher HDR efficiency than linear dsDNA donors",
                    "    (Iyer et al., CRISPR Journal, 2022).",
                    "  - Optimal homology arm length: 300 nt per arm",
                    "    (Iyer et al., 2022).",
                    "  - The donor must encode the GENOMIC sequence (including",
                    "    introns if the edit spans exon-intron boundaries),",
                    "    NOT the cDNA. HDR uses chromosomal homology.",
                    "  - Include a silent PAM-blocking mutation in the donor to",
                    "    prevent re-cutting of the corrected allele by the",
                    "    nuclease. This is critical for maintaining the HDR",
                    "    product.",
                    "    Source: Paquet et al., Nature, 2016.",
                    f"  - Number of donors to synthesize: {s.donor_templates_needed}",
                ])
            lines.append("")

        # --- Screening strategy ---
        lines.append("Screening Strategy:")
        lines.extend([
            f"  - Estimated clones to screen: ~{s.screening_burden}",
            f"    (based on predicted efficiency of {s.estimated_efficiency:.1%})",
            "  - Genotyping method: PCR + Sanger sequencing across each",
            "    mutation site, on both alleles.",
            "  - CRITICAL: Also sequence the non-edited allele to check for",
            "    NHEJ-induced indels (if HDR was used). Allele-specific PCR",
            "    or amplicon deep sequencing is recommended.",
        ])
        if s.num_dsbs > 0:
            lines.extend([
                "  - Karyotype validated clones (G-banding or SNP array) to",
                "    detect chromosomal rearrangements.",
                "    Source: International Stem Cell Initiative, Nat Biotech, 2011.",
            ])
        if s.num_dsbs >= 2 and s.num_editing_rounds == 1:
            lines.extend([
                "  - ADDITIONAL: For simultaneous dual-DSB strategies, perform",
                "    junction PCR to rule out interstitial deletion between",
                "    the two cut sites. Design primers flanking the expected",
                "    deletion breakpoints.",
            ])
        lines.append("")

        # --- Timeline ---
        lines.append("Estimated Timeline:")
        if s.num_editing_rounds == 1:
            lines.extend([
                f"  Week 1:   Electroporation + recovery",
                f"  Week 2:   Single-cell plating + clonal expansion",
                f"  Week 3-4: Colony picking + genotyping",
                f"  Week 5-6: Validation (Sanger confirmation, karyotyping)",
                f"  TOTAL:    ~{s.time_weeks + 2} weeks to validated clone",
            ])
        elif s.num_editing_rounds == 2:
            lines.extend([
                f"  Round 1 (Weeks 1-4):",
                f"    Week 1:   Electroporation + recovery",
                f"    Week 2:   Single-cell plating",
                f"    Week 3-4: Genotyping + expansion of correct clone",
                f"  Round 2 (Weeks 5-8):",
                f"    Week 5:   Second electroporation",
                f"    Week 6:   Single-cell plating",
                f"    Week 7-8: Genotyping + expansion",
                f"  Validation (Weeks 9-10):",
                f"    Sanger confirmation, karyotyping, mycoplasma test",
                f"  TOTAL:    ~{s.time_weeks + 2} weeks to validated clone",
            ])
        lines.append("")

        return "\n".join(lines)

    def _warnings_section(
        self, ranked: List[ScoredStrategy], cell_type: str
    ) -> str:
        """Safety warnings and quality control recommendations."""
        lines = [
            "SECTION 7: SAFETY WARNINGS AND QUALITY CONTROL",
            "-" * 40,
        ]

        has_dsb = any(ss.strategy.requires_dsb for ss in ranked[:3])
        has_dual_dsb = any(
            ss.strategy.num_dsbs >= 2 and ss.strategy.num_editing_rounds == 1
            for ss in ranked[:3]
        )

        if cell_type == "iPSC":
            lines.extend([
                "iPSC-SPECIFIC WARNINGS:",
                "",
                "1. p53 SELECTION RISK",
                "   Double-strand breaks activate the p53 pathway in iPSCs,",
                "   causing significant cell death. Surviving clones may have",
                "   acquired inactivating TP53 mutations under this selective",
                "   pressure. TP53-mutant iPSCs are genomically unstable and",
                "   UNSUITABLE for clinical applications or reliable disease",
                "   modelling.",
                "   Action: Sequence TP53 exons 4-9 in all validated clones.",
                "   Source: Ihry et al., Nature Medicine, 2018.",
                "",
            ])

        if has_dsb:
            lines.extend([
                "2. LARGE DELETION RISK",
                "   Even a single Cas9-induced DSB can cause large deletions",
                "   (>1 kb, sometimes >10 kb) at the cut site. These are",
                "   often UNDETECTABLE by standard short-range PCR genotyping.",
                "   Action: Use long-range PCR (>5 kb amplicons) or Southern",
                "   blot to confirm structural integrity around the edit site.",
                "   Source: Kosicki et al., Nature Biotechnology, 2018.",
                "",
            ])

        if has_dual_dsb:
            lines.extend([
                "3. TRANSLOCATION / INTERSTITIAL DELETION RISK",
                "   Two simultaneous DSBs on the same chromosome can cause",
                "   loss of the intervening genomic segment (interstitial",
                "   deletion) or chromosomal translocation / inversion.",
                "   These events produce cells with apparently correct",
                "   genotypes at each edit site but catastrophic structural",
                "   damage elsewhere.",
                "   Action: Perform G-banding karyotyping or SNP array on",
                "   all candidate clones. Consider FISH with probes spanning",
                "   the inter-site region.",
                "   Source: Leibowitz et al., Nature Genetics, 2021.",
                "",
            ])

        lines.extend([
            "RECOMMENDED QUALITY CONTROL PIPELINE:",
            "  1. PCR + Sanger sequencing at all edit sites (both alleles)",
            "  2. Long-range PCR or Southern blot to rule out large deletions",
            "  3. TP53 sequencing (exons 4-9)",
            "  4. G-banding karyotype (at least 20 metaphases)",
            "  5. Mycoplasma testing",
            "  6. STR profiling (to confirm cell line identity)",
            "  7. Pluripotency marker verification (SSEA-4, TRA-1-60, OCT4)",
            "",
        ])

        return "\n".join(lines)

    def _references_section(self) -> str:
        """Key literature references cited in the report."""
        lines = [
            "SECTION 8: KEY REFERENCES",
            "-" * 40,
            "",
            "Editing technologies:",
            "  - Komor AC et al. Programmable editing of a target base in",
            "    genomic DNA without double-stranded DNA cleavage.",
            "    Nature. 2016;533(7603):420-424. [CBE]",
            "",
            "  - Gaudelli NM et al. Programmable base editing of A*T to G*C",
            "    in genomic DNA without DNA cleavage.",
            "    Nature. 2017;551(7681):464-471. [ABE]",
            "",
            "  - Anzalone AV et al. Search-and-replace genome editing without",
            "    double-strand breaks or donor DNA.",
            "    Nature. 2019;576(7785):149-157. [Prime editing]",
            "",
            "  - Chen PJ et al. Enhanced prime editing systems by manipulating",
            "    cellular determinants of editing outcomes.",
            "    Cell. 2021;184(22):5635-5652. [PE4/PE5]",
            "",
            "Donor templates:",
            "  - Iyer S et al. Efficient homology-directed repair with circular",
            "    single-stranded DNA donors.",
            "    CRISPR Journal. 2022;5(5):685-701. [cssDNA donors]",
            "",
            "  - Richardson CD et al. Enhancing homology-directed genome editing",
            "    by catalytically active and inactive CRISPR-Cas9 using",
            "    asymmetric donor DNA.",
            "    Nature Biotechnology. 2016;34(3):339-344. [Asymmetric donors]",
            "",
            "Safety and toxicity:",
            "  - Ihry RJ et al. p53 inhibits CRISPR-Cas9 engineering in human",
            "    pluripotent stem cells.",
            "    Nature Medicine. 2018;24(7):939-946. [p53 selection in iPSCs]",
            "",
            "  - Kosicki M et al. Repair of double-strand breaks induced by",
            "    CRISPR-Cas9 leads to large deletions and complex",
            "    rearrangements.",
            "    Nature Biotechnology. 2018;36(8):765-771. [Large deletions]",
            "",
            "  - Leibowitz ML et al. Chromothripsis as an on-target consequence",
            "    of CRISPR-Cas9 genome editing.",
            "    Nature Genetics. 2021;53(6):895-905. [Chromothripsis]",
            "",
            "Nuclease variants:",
            "  - Chauhan VP et al. Staggered-cut Cas9 variants generate",
            "    5-prime overhangs that enhance homology-directed repair.",
            "    PNAS. 2023. [vCas9]",
            "",
            "3D genome and translocation risk:",
            "  - Lieberman-Aiden E et al. Comprehensive mapping of long-range",
            "    interactions reveals folding principles of the human genome.",
            "    Science. 2009;326(5950):289-293. [Hi-C contact frequency]",
            "",
            "  - Zhang Y et al. Spatial organization of the mouse genome and",
            "    its role in recurrent chromosomal translocations.",
            "    Cell. 2012;148(5):908-921. [Translocation-contact correlation]",
            "",
        ]
        return "\n".join(lines)

    def _footer(self) -> str:
        """Report footer."""
        border = "=" * self.line_width
        return "\n".join([
            border,
            "End of MOSAIC Report",
            "This report is generated by CRISPRArchitect and is intended as",
            "a decision-support tool.  All predictions are based on published",
            "literature values and simplified models.  Experimental validation",
            "is always required.",
            border,
        ])


# =========================================================================
# COMPLETE EXAMPLE: 60-exon gene, mutations in exon 40 and exon 60
# =========================================================================
#
# This example demonstrates the full MOSAIC workflow for the user's
# actual use case: a large multi-exon gene with compound heterozygous
# mutations in widely separated exons.

if __name__ == "__main__":

    # =====================================================================
    # STEP 1: Build the gene structure
    # =====================================================================
    # We use the demo 60-exon gene (inspired by COL7A1, USH2A, RYR1).
    # In a real workflow, you would provide actual exon coordinates from
    # UCSC or Ensembl.

    from .gene_structure import build_demo_gene

    print("Building gene structure...")
    gene = build_demo_gene()
    print(gene.summary())
    print()

    # =====================================================================
    # STEP 2: Define the patient's mutations
    # =====================================================================
    # Scenario: compound heterozygous — one mutation per allele.
    #
    # Mutation 1: G>A missense in exon 40
    #   - This is a transition: the patient has A where the reference has G.
    #   - Correction needed: A -> G (ABE can do this).
    #   - Alternatively, on the antisense strand: T -> C (also ABE).
    #
    # Mutation 2: 2-bp deletion (CT deleted) in exon 60
    #   - This causes a frameshift.
    #   - Not amenable to base editing (it is a deletion, not a substitution).
    #   - Amenable to prime editing (2 bp insertion to restore the deleted CT).
    #   - Amenable to HDR (any mutation can be corrected by HDR).

    m1 = Mutation(
        exon_number=40,
        position=5_000_000,       # Genomic coordinate (example)
        ref_allele="G",
        alt_allele="A",
        name="c.5000G>A (p.Gly1667Asp) [exon 40]",
    )

    m2 = Mutation(
        exon_number=60,
        position=8_000_000,       # Genomic coordinate (example)
        ref_allele="CT",
        alt_allele="-",
        name="c.8000delCT (p.Leu2667fs) [exon 60]",
    )

    mutations = [m1, m2]

    print("Mutations defined:")
    classifier = MutationClassifier()
    for m in mutations:
        print(classifier.describe_mutation(m))
        print()

    # =====================================================================
    # STEP 3: Enumerate feasible strategies
    # =====================================================================
    # The enumerator will check:
    #   - Genomic distance between exon 40 and exon 60 (likely very large,
    #     >100 kb, so single-template HDR is infeasible).
    #   - Mutation classifications (m1 = transition_AG -> ABE-amenable;
    #     m2 = deletion -> prime-editing-amenable but not base-editable).
    #   - Cell type constraints (iPSC: p53-active, low viability with DSBs).

    print("Enumerating strategies...")
    enumerator = StrategyEnumerator()
    strategies = enumerator.enumerate_strategies(
        gene=gene,
        mutations=mutations,
        cell_type="iPSC",
        nuclease="SpCas9",
    )
    print(enumerator.summarize_strategies(strategies))
    print()

    # =====================================================================
    # STEP 4: Score and rank strategies
    # =====================================================================
    # The scorer weights safety highest (0.40) because we are working
    # with iPSCs.

    print("Scoring and ranking...")
    scorer = StrategyScorer(
        weight_efficiency=0.30,
        weight_safety=0.40,
        weight_time=0.15,
        weight_cost=0.15,
    )
    ranked = scorer.rank_strategies(strategies, cell_type="iPSC")
    print(scorer.format_ranking(ranked))
    print()

    # =====================================================================
    # STEP 5: Generate the full report
    # =====================================================================
    reporter = StrategyReporter()
    report = reporter.generate_report(
        gene=gene,
        mutations=mutations,
        ranked_strategies=ranked,
        cell_type="iPSC",
        nuclease="SpCas9",
    )
    print(report)
