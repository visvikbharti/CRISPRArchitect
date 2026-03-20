"""
Gene Structure Representation for MOSAIC
=========================================

This module models the exon-intron architecture of a gene so that MOSAIC
can calculate genomic distances, determine whether two mutation sites can
be spanned by a single donor template, and estimate translocation risk
from inter-locus DSBs.

Why gene structure matters for editing strategy
------------------------------------------------
When two mutations sit in different exons of the same gene, the *genomic*
distance between them (including introns) determines:

  * Whether a single homology-directed repair (HDR) donor template can
    span both sites.  Circular single-stranded DNA (cssDNA) donors have
    a practical upper limit of ~5-10 kb total insert + homology arms
    (Iyer et al., CRISPR Journal, 2022).  AAV donors are limited to
    ~4.7 kb total.

  * The translocation risk if two DSBs are introduced simultaneously.
    Translocation frequency scales roughly with 3D contact probability,
    which in turn scales as a power law of genomic separation:
        P(translocation) ~ d^(-gamma),  gamma ~ 1.08
    (Lieberman-Aiden et al., Science, 2009; Zhang et al., Cell, 2012).

  * The feasibility of exon-deletion strategies.  If two mutations lie
    in adjacent small exons, it may be simpler to delete the exons and
    rely on exon skipping (if the reading frame is preserved).

Biology primer on gene structure
---------------------------------
Human protein-coding genes consist of:
  - **Exons**: the portions of the gene that are retained in the mature
    mRNA after splicing.  Median human exon size is ~150 bp, though
    first and last exons tend to be larger.
    Source: Sakharkar et al., In Silico Biology, 2004.
  - **Introns**: intervening sequences that are removed by the
    spliceosome.  Human intron sizes vary enormously, from <100 bp to
    >500 kb.  Median ~1.5 kb, but the distribution is very right-skewed.
    Source: Sakharkar et al., 2004; Lander et al., Nature, 2001.

For large multi-exon genes (e.g., 60+ exons like COL7A1, USH2A, RYR1),
the total genomic span can be hundreds of kilobases, even though the
coding sequence (CDS) is only 5-15 kb.  This means that two mutations
that are only a few kb apart in the mRNA may be separated by hundreds
of kb on the chromosome — far too large for a single donor.
"""

from __future__ import annotations

import sys
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

import numpy as np


# ---------------------------------------------------------------------------
# ExonInfo dataclass
# ---------------------------------------------------------------------------

@dataclass
class ExonInfo:
    """A single exon within a gene.

    Attributes
    ----------
    exon_number : int
        1-based exon number (exon 1 is the first exon in the mRNA).
    start : int
        Genomic start coordinate (0-based, inclusive).  For genes on the
        minus strand the ``start`` is still the lower genomic coordinate
        (i.e., coordinates are always stored in ascending order on the
        reference assembly, following BED convention).
    end : int
        Genomic end coordinate (0-based, exclusive — BED convention).
    size : int
        Length of the exon in base pairs.  Automatically computed as
        ``end - start`` if not provided.
    sequence : Optional[str]
        The nucleotide sequence of the exon (5' -> 3' on the sense
        strand).  ``None`` if not yet fetched.
    """

    exon_number: int
    start: int
    end: int
    size: int = 0
    sequence: Optional[str] = None

    def __post_init__(self) -> None:
        """Auto-compute size from coordinates if the caller left it at 0."""
        if self.size == 0:
            self.size = self.end - self.start


# ---------------------------------------------------------------------------
# GeneStructure class
# ---------------------------------------------------------------------------

class GeneStructure:
    """Represent the exon/intron architecture of a gene.

    This is the foundation for all genomic-distance calculations in MOSAIC.
    You can build a ``GeneStructure`` either from manually provided exon
    coordinates (``from_manual``) or — in future releases — from an Ensembl
    gene ID.

    Parameters
    ----------
    gene_name : str
        HGNC gene symbol (e.g. ``"COL7A1"``).
    chromosome : str
        Chromosome name (e.g. ``"chr3"``).
    strand : str
        ``"+"`` or ``"-"``.
    exons : List[ExonInfo]
        Sorted list of ``ExonInfo`` objects ordered by exon number.
        Coordinates must be in ascending genomic order.

    Notes
    -----
    Introns are derived automatically from the gaps between consecutive
    exons.  Intron *i* is defined as the region between the end of exon
    *i* and the start of exon *i+1*.
    """

    def __init__(
        self,
        gene_name: str,
        chromosome: str,
        strand: str,
        exons: List[ExonInfo],
    ) -> None:
        self.gene_name: str = gene_name
        self.chromosome: str = chromosome
        self.strand: str = strand

        # Sort exons by genomic start coordinate to guarantee order.
        self.exons: List[ExonInfo] = sorted(exons, key=lambda e: e.start)

        # Build a quick lookup: exon_number -> index in self.exons
        self._exon_index: Dict[int, int] = {
            e.exon_number: i for i, e in enumerate(self.exons)
        }

        # Derive introns from inter-exon gaps.
        self.introns: List[Dict[str, int]] = self._derive_introns()

    # -----------------------------------------------------------------
    # Construction helpers
    # -----------------------------------------------------------------

    @classmethod
    def from_manual(
        cls,
        gene_name: str,
        exons_list: List[Dict[str, int]],
        chromosome: str = "chrN",
        strand: str = "+",
    ) -> "GeneStructure":
        """Create a ``GeneStructure`` from manually provided exon coords.

        This is the simplest way to build a gene model when you already
        know the exon boundaries (e.g., from UCSC Genome Browser or
        Ensembl).

        Parameters
        ----------
        gene_name : str
            HGNC gene symbol.
        exons_list : list of dict
            Each dict must have keys ``"start"``, ``"end"``, and
            ``"number"`` (1-based exon number).  Coordinates are
            0-based, half-open (BED convention).
        chromosome : str, optional
            Chromosome name (default ``"chrN"`` placeholder).
        strand : str, optional
            ``"+"`` or ``"-"`` (default ``"+"``).

        Returns
        -------
        GeneStructure
            Fully constructed gene model with derived introns.

        Example
        -------
        >>> gs = GeneStructure.from_manual(
        ...     "MYH7",
        ...     [
        ...         {"number": 1, "start": 1000, "end": 1200},
        ...         {"number": 2, "start": 3000, "end": 3180},
        ...     ],
        ... )
        """
        exon_infos: List[ExonInfo] = []
        for e in exons_list:
            exon_infos.append(
                ExonInfo(
                    exon_number=e["number"],
                    start=e["start"],
                    end=e["end"],
                )
            )
        return cls(
            gene_name=gene_name,
            chromosome=chromosome,
            strand=strand,
            exons=exon_infos,
        )

    # -----------------------------------------------------------------
    # Intron derivation
    # -----------------------------------------------------------------

    def _derive_introns(self) -> List[Dict[str, int]]:
        """Derive introns from the gaps between consecutive exons.

        An intron is defined as the genomic interval between the end of
        exon *i* and the start of exon *i+1*.  This matches the standard
        GT-AG splice-site convention: the intron begins immediately after
        the donor site at the end of the upstream exon and ends just
        before the acceptor site at the start of the downstream exon.

        Returns
        -------
        list of dict
            Each dict has keys ``"intron_number"`` (1-based, intron 1 is
            between exon 1 and exon 2), ``"start"``, ``"end"``, and
            ``"size"``.
        """
        introns: List[Dict[str, int]] = []
        for i in range(len(self.exons) - 1):
            intron_start = self.exons[i].end
            intron_end = self.exons[i + 1].start
            introns.append(
                {
                    "intron_number": i + 1,
                    "start": intron_start,
                    "end": intron_end,
                    "size": intron_end - intron_start,
                }
            )
        return introns

    # -----------------------------------------------------------------
    # Distance calculations
    # -----------------------------------------------------------------

    def _validate_exon_numbers(self, exon_a: int, exon_b: int) -> Tuple[int, int]:
        """Return (low, high) exon numbers after validation.

        Raises ``ValueError`` if either exon number is not present in
        the gene model.
        """
        if exon_a not in self._exon_index:
            raise ValueError(
                f"Exon {exon_a} not found in gene {self.gene_name}. "
                f"Available exons: {sorted(self._exon_index.keys())[:5]}..."
            )
        if exon_b not in self._exon_index:
            raise ValueError(
                f"Exon {exon_b} not found in gene {self.gene_name}. "
                f"Available exons: {sorted(self._exon_index.keys())[:5]}..."
            )
        lo = min(exon_a, exon_b)
        hi = max(exon_a, exon_b)
        return lo, hi

    def genomic_distance(self, exon_a: int, exon_b: int) -> int:
        """Total genomic distance between the starts of two exons.

        This is the distance on the chromosome, *including* all introns
        and exons that lie between the two target exons.  It determines
        whether a single donor template can physically span both sites,
        and it feeds into the translocation-risk calculation.

        Parameters
        ----------
        exon_a, exon_b : int
            1-based exon numbers.

        Returns
        -------
        int
            Genomic distance in base pairs.

        Biology
        -------
        For a 60-exon gene where exon 40 and exon 60 are separated by
        20 intervening introns, the genomic distance can easily be
        200-500 kb — far beyond the capacity of any single donor
        template.  In contrast, two adjacent exons may be separated by
        only 1-2 kb, easily spanned by cssDNA.
        """
        lo, hi = self._validate_exon_numbers(exon_a, exon_b)
        idx_lo = self._exon_index[lo]
        idx_hi = self._exon_index[hi]
        exon_lo = self.exons[idx_lo]
        exon_hi = self.exons[idx_hi]
        return exon_hi.end - exon_lo.start

    def exonic_distance(self, exon_a: int, exon_b: int) -> int:
        """Total exonic (coding) sequence between two exons, inclusive.

        This sums the sizes of all exons from ``exon_a`` through
        ``exon_b`` inclusive.  It tells you how much coding sequence
        separates the two mutation sites, which is relevant for:

          * Designing cDNA-based donors (which lack introns).
          * Estimating the impact of an exon-deletion strategy on the
            protein product (how many amino acids would be lost).

        Parameters
        ----------
        exon_a, exon_b : int
            1-based exon numbers.

        Returns
        -------
        int
            Total exonic base pairs.
        """
        lo, hi = self._validate_exon_numbers(exon_a, exon_b)
        idx_lo = self._exon_index[lo]
        idx_hi = self._exon_index[hi]
        total = 0
        for i in range(idx_lo, idx_hi + 1):
            total += self.exons[i].size
        return total

    def intronic_distance(self, exon_a: int, exon_b: int) -> int:
        """Total intronic sequence between two exons.

        This is ``genomic_distance - exonic_distance`` and represents
        the non-coding sequence that inflates the genomic span.  In many
        human genes, introns account for >95% of the genomic distance
        between distant exons.

        Parameters
        ----------
        exon_a, exon_b : int
            1-based exon numbers.

        Returns
        -------
        int
            Total intronic base pairs between the two exons.
        """
        return self.genomic_distance(exon_a, exon_b) - self.exonic_distance(exon_a, exon_b)

    def can_span_with_single_donor(
        self,
        exon_a: int,
        exon_b: int,
        max_donor_size: int = 10_000,
    ) -> bool:
        """Check if both mutation sites can be covered by one donor template.

        A single donor can span two sites only if the entire genomic
        region between (and including) the two target exons is smaller
        than the maximum donor size.  For cssDNA, practical limits are
        ~5-10 kb of insert sequence (Iyer et al., CRISPR Journal, 2022).
        For AAV, the packaging limit is ~4.7 kb.

        Parameters
        ----------
        exon_a, exon_b : int
            1-based exon numbers where the mutations reside.
        max_donor_size : int, optional
            Maximum donor insert size in bp (default 10 000 for cssDNA).
            This should include the homology arms.  Typical values:
              - cssDNA: 5000-10000 bp
              - AAV: 4700 bp
              - plasmid: essentially unlimited, but HDR efficiency drops
                for large donors.

        Returns
        -------
        bool
            ``True`` if the genomic distance is within the donor limit.

        Biology
        -------
        The key insight is that HDR donors must contain the **genomic**
        sequence (with introns), not just the cDNA.  You cannot provide
        a cDNA-only donor because the cell's HDR machinery recombines
        based on homology to the *chromosomal* DNA, which includes
        introns.  (Exception: some specialized gene-trap or exon-trap
        strategies, which MOSAIC does not currently model.)
        """
        gd = self.genomic_distance(exon_a, exon_b)
        return gd <= max_donor_size

    # -----------------------------------------------------------------
    # Summary
    # -----------------------------------------------------------------

    def summary(self) -> str:
        """Return a human-readable summary of the gene structure.

        The summary includes total genomic span, number of exons, total
        coding sequence length, largest and smallest exons and introns,
        and some derived statistics useful for editing strategy design.

        Returns
        -------
        str
            Multi-line formatted summary string.
        """
        total_genomic = self.exons[-1].end - self.exons[0].start
        total_exonic = sum(e.size for e in self.exons)
        total_intronic = total_genomic - total_exonic
        exon_sizes = [e.size for e in self.exons]
        intron_sizes = [i["size"] for i in self.introns] if self.introns else [0]

        lines = [
            f"Gene Structure Summary: {self.gene_name}",
            f"{'=' * 50}",
            f"Chromosome:          {self.chromosome} ({self.strand} strand)",
            f"Genomic span:        {total_genomic:,} bp "
            f"({total_genomic / 1000:.1f} kb)",
            f"Number of exons:     {len(self.exons)}",
            f"Number of introns:   {len(self.introns)}",
            f"Total exonic seq:    {total_exonic:,} bp "
            f"({total_exonic / 1000:.1f} kb)",
            f"Total intronic seq:  {total_intronic:,} bp "
            f"({total_intronic / 1000:.1f} kb)",
            f"Intronic fraction:   {total_intronic / total_genomic * 100:.1f}%",
            f"",
            f"Exon sizes:   min={min(exon_sizes)} bp, "
            f"median={int(np.median(exon_sizes))} bp, "
            f"max={max(exon_sizes)} bp",
            f"Intron sizes: min={min(intron_sizes):,} bp, "
            f"median={int(np.median(intron_sizes)):,} bp, "
            f"max={max(intron_sizes):,} bp",
        ]
        return "\n".join(lines)

    def __repr__(self) -> str:
        return (
            f"GeneStructure(gene_name={self.gene_name!r}, "
            f"exons={len(self.exons)}, "
            f"span={self.exons[-1].end - self.exons[0].start:,} bp)"
        )


# =========================================================================
# DEMO GENE: a realistic 60-exon gene for testing
# =========================================================================

def build_demo_gene() -> GeneStructure:
    """Build a demo 60-exon gene with realistic human exon/intron sizes.

    This gene model is inspired by large multi-exon genes such as:
      - **COL7A1** (118 exons, mutations cause dystrophic epidermolysis
        bullosa)
      - **USH2A** (72 exons, mutations cause Usher syndrome)
      - **RYR1** (106 exons, mutations cause malignant hyperthermia)

    Exon and intron sizes are drawn from distributions matching human
    genome averages:
      - Exon sizes: log-normal with median ~150 bp
        Source: Sakharkar et al., In Silico Biology, 2004.
      - Intron sizes: log-normal with median ~2 kb, but with occasional
        very large introns (>100 kb).  The first intron is typically
        the largest.
        Source: Sakharkar et al., 2004; Lander et al., Nature, 2001.

    The gene is placed on chr3 starting at position 1,000,000
    (arbitrary).  Strand is "+".

    Returns
    -------
    GeneStructure
        A fully built gene model with 60 exons.

    Notes
    -----
    We use a fixed random seed (42) so that the demo gene is identical
    every time this function is called.  This ensures reproducible
    testing.
    """
    rng = np.random.RandomState(42)

    num_exons = 60
    gene_start = 1_000_000  # Arbitrary genomic start on chr3

    # --- Generate exon sizes ---
    # Median human exon size ~ 150 bp (Sakharkar et al., 2004).
    # We model this as a log-normal distribution:
    #   ln(size) ~ Normal(mu, sigma)  where exp(mu) = median
    #   mu = ln(150), sigma chosen to give a range of ~50-500 bp
    exon_mu = np.log(150)
    exon_sigma = 0.45  # Gives a spread roughly consistent with Sakharkar
    exon_sizes = rng.lognormal(mean=exon_mu, sigma=exon_sigma, size=num_exons)
    exon_sizes = np.clip(exon_sizes, 50, 800).astype(int)

    # First and last exons tend to be larger (UTR-containing).
    # Source: General genomics knowledge; first exon often contains 5'UTR.
    exon_sizes[0] = int(exon_sizes[0] * 2.0)   # First exon larger
    exon_sizes[-1] = int(exon_sizes[-1] * 2.5)  # Last exon larger (3'UTR)

    # --- Generate intron sizes ---
    # Human introns have a very wide size distribution.  Median ~ 1.5-2 kb
    # but 10-20% of introns can be >20 kb, and the first intron is often
    # the largest (median ~7 kb; Bradnam & Korf, PLoS Comp Bio, 2008).
    # We model with a bimodal strategy:
    #   - Most introns: log-normal, median ~2 kb
    #   - A few large introns: log-normal, median ~50 kb
    num_introns = num_exons - 1
    intron_sizes = np.zeros(num_introns, dtype=int)

    for i in range(num_introns):
        if i == 0:
            # First intron is typically the largest in human genes.
            # Source: Bradnam & Korf, PLoS Comp Bio, 2008
            intron_sizes[i] = int(rng.lognormal(mean=np.log(80_000), sigma=0.4))
        elif rng.random() < 0.15:
            # ~15% of introns are large (>20 kb)
            intron_sizes[i] = int(rng.lognormal(mean=np.log(50_000), sigma=0.6))
        else:
            # Most introns are 0.5-5 kb
            intron_sizes[i] = int(rng.lognormal(mean=np.log(2_000), sigma=0.7))
        # Enforce a minimum intron size (human minimum ~60 bp for functional
        # splicing; the U2-type spliceosome needs ~80 nt minimum).
        # Source: Wieringa et al., EMBO J, 1984
        intron_sizes[i] = max(intron_sizes[i], 80)

    # --- Build exon coordinates ---
    exons_list: List[Dict[str, int]] = []
    current_pos = gene_start

    for i in range(num_exons):
        exon_start = current_pos
        exon_end = exon_start + int(exon_sizes[i])
        exons_list.append(
            {"number": i + 1, "start": exon_start, "end": exon_end}
        )
        current_pos = exon_end
        # Add intron after every exon except the last one.
        if i < num_introns:
            current_pos += int(intron_sizes[i])

    return GeneStructure.from_manual(
        gene_name="DEMO60",
        exons_list=exons_list,
        chromosome="chr3",
        strand="+",
    )


# =========================================================================
# Quick self-test when run as a script
# =========================================================================

if __name__ == "__main__":
    gene = build_demo_gene()
    print(gene.summary())
    print()
    print(f"Genomic distance exon 40 -> exon 60:  "
          f"{gene.genomic_distance(40, 60):,} bp")
    print(f"Exonic distance exon 40 -> exon 60:   "
          f"{gene.exonic_distance(40, 60):,} bp")
    print(f"Intronic distance exon 40 -> exon 60: "
          f"{gene.intronic_distance(40, 60):,} bp")
    print(f"Can single donor span exon 40-60?     "
          f"{gene.can_span_with_single_donor(40, 60)}")
    print(f"Can single donor span exon 59-60?     "
          f"{gene.can_span_with_single_donor(59, 60)}")
