"""
Coding Annotator — Variant consequence and HGVS annotation
===========================================================

Given a variant, its transcript coordinate, and the CDS sequence, this
module determines the functional consequence (synonymous, missense,
nonsense, frameshift, etc.) and generates HGVS-like c. and p. notation.

Biology context
---------------
The consequence of a coding variant depends on:
  1. Whether it changes the encoded amino acid (synonymous vs. missense)
  2. Whether it introduces a premature stop codon (nonsense)
  3. Whether it shifts the reading frame (frameshift for indels not
     divisible by 3)
  4. Whether it is near a splice junction (splice_donor, splice_acceptor,
     or splice_region)

Splice site definitions (ACMG/AMP standards):
  - splice_donor:   first 2 intronic bases at the 3' end of an exon
                     (GT of the donor site), or the last 2 exonic bases
  - splice_acceptor: last 2 intronic bases at the 5' end of an exon
                     (AG of the acceptor site), or the first 2 exonic bases
  - splice_region:  within 3-8 bp of an exon boundary (either side)

For simplicity, this module uses exonic distances:
  - <= 2 bp from exon boundary = potential splice donor/acceptor
  - 3-8 bp from exon boundary = splice region

Three-letter amino acid codes are used in p. notation per HGVS guidelines.

References
----------
- Richards et al., Genetics in Medicine, 2015 (ACMG/AMP variant guidelines)
- den Dunnen et al., Human Mutation, 2016 (HGVS nomenclature v15.11)
- McLaren et al., Genome Biology, 2016 (Ensembl VEP consequence hierarchy)
"""

from __future__ import annotations

from typing import Optional

from core.models import (
    CodingAnnotation,
    ConsequenceType,
    TranscriptCoordinate,
)

# v1 genetic code
try:
    from utils.sequence import GENETIC_CODE
except ImportError:
    try:
        from crisprarchitect.utils.sequence import GENETIC_CODE
    except ImportError:
        GENETIC_CODE = {}


# Three-letter amino acid codes for HGVS p. notation
# Reference: IUPAC-IUB Commission on Biochemical Nomenclature, 1968
AA_THREE_LETTER = {
    'A': 'Ala', 'R': 'Arg', 'N': 'Asn', 'D': 'Asp', 'C': 'Cys',
    'E': 'Glu', 'Q': 'Gln', 'G': 'Gly', 'H': 'His', 'I': 'Ile',
    'L': 'Leu', 'K': 'Lys', 'M': 'Met', 'F': 'Phe', 'P': 'Pro',
    'S': 'Ser', 'T': 'Thr', 'W': 'Trp', 'Y': 'Tyr', 'V': 'Val',
    '*': 'Ter',  # stop codon
}

# Splice distance thresholds (bp from exon boundary)
# ACMG PVS1 decision tree; Tayoun et al., Human Mutation, 2018
SPLICE_DONOR_ACCEPTOR_THRESHOLD = 2  # within 2 bp = canonical splice site
SPLICE_REGION_THRESHOLD = 8          # within 3-8 bp = splice region


class CodingAnnotator:
    """Annotates variants with coding consequence and HGVS notation.

    Takes a variant's transcript coordinate, the CDS sequence, and the
    variant alleles, and produces a CodingAnnotation with:
      - ConsequenceType (synonymous, missense, nonsense, etc.)
      - HGVS c. notation (e.g., "c.910C>T")
      - HGVS p. notation (e.g., "p.Arg304Ter")
      - Splice proximity assessment
      - Reference and alternate codons/amino acids

    Parameters
    ----------
    cds_sequence : str
        The full CDS sequence in transcript orientation (starts with ATG).
        Used for codon lookup and reading-frame analysis.

    Examples
    --------
    >>> annotator = CodingAnnotator("ATGCCCAAG...")
    >>> result = annotator.annotate(coord, "C", "T")
    >>> print(result.consequence, result.hgvs_c, result.hgvs_p)
    """

    def __init__(self, cds_sequence: str = ""):
        self.cds = cds_sequence.upper() if cds_sequence else ""

    def annotate(
        self,
        coord: TranscriptCoordinate,
        ref_allele: str,
        alt_allele: str,
    ) -> CodingAnnotation:
        """Annotate a variant with coding consequence.

        Parameters
        ----------
        coord : TranscriptCoordinate
            The transcript coordinate from TranscriptMapper.
        ref_allele : str
            Reference allele in transcript orientation.
        alt_allele : str
            Alternate allele in transcript orientation.

        Returns
        -------
        CodingAnnotation
            Full annotation including consequence type and HGVS strings.
        """
        ref = ref_allele.upper()
        alt = alt_allele.upper()

        # Check splice proximity first (applies to all variants)
        splice_prox = self._assess_splice_proximity(coord)

        # Handle non-coding / intronic positions
        if not coord.in_cds or coord.cds_position < 1:
            return self._non_coding_annotation(coord, ref, alt, splice_prox)

        # Determine if this is an SNV or indel
        is_insertion = (ref == "-" or ref == "")
        is_deletion = (alt == "-" or alt == "")
        is_snv = (
            not is_insertion
            and not is_deletion
            and len(ref) == 1
            and len(alt) == 1
        )
        is_mnv = (
            not is_insertion
            and not is_deletion
            and len(ref) == len(alt)
            and len(ref) > 1
        )
        is_complex_indel = (
            not is_snv
            and not is_mnv
            and not is_insertion
            and not is_deletion
        )

        # Check for splice-site override (within 2 bp of exon boundary)
        if splice_prox is not None:
            min_dist = min(coord.distance_to_exon_start, coord.distance_to_exon_end)
            if min_dist <= SPLICE_DONOR_ACCEPTOR_THRESHOLD:
                return self._splice_annotation(coord, ref, alt, splice_prox)

        if is_snv:
            return self._annotate_snv(coord, ref, alt, splice_prox)
        elif is_insertion:
            return self._annotate_insertion(coord, ref, alt, splice_prox)
        elif is_deletion:
            return self._annotate_deletion(coord, ref, alt, splice_prox)
        elif is_mnv or is_complex_indel:
            return self._annotate_complex(coord, ref, alt, splice_prox)
        else:
            return CodingAnnotation(
                consequence=ConsequenceType.UNKNOWN,
                message=f"Unable to classify variant: {ref}>{alt}",
                splice_proximity=splice_prox,
            )

    # ----- SNV annotation ---------------------------------------------------

    def _annotate_snv(
        self,
        coord: TranscriptCoordinate,
        ref: str,
        alt: str,
        splice_prox: Optional[str],
    ) -> CodingAnnotation:
        """Annotate a single nucleotide variant."""
        cds_pos = coord.cds_position  # 1-based
        codon_idx = coord.codon_index  # 1-based
        codon_pos = coord.codon_position  # 1, 2, or 3

        # Retrieve reference codon from CDS
        codon_start = (codon_idx - 1) * 3  # 0-based in CDS string
        if self.cds and codon_start + 3 <= len(self.cds):
            ref_codon = self.cds[codon_start:codon_start + 3]
        else:
            ref_codon = coord.reference_codon

        if not ref_codon or len(ref_codon) != 3:
            return CodingAnnotation(
                consequence=ConsequenceType.UNKNOWN,
                hgvs_c=f"c.{cds_pos}{ref}>{alt}",
                message="Could not determine reference codon.",
                splice_proximity=splice_prox,
            )

        # Build alternate codon
        codon_list = list(ref_codon)
        codon_list[codon_pos - 1] = alt  # codon_pos is 1-based
        alt_codon = ''.join(codon_list)

        ref_aa = GENETIC_CODE.get(ref_codon, "?")
        alt_aa = GENETIC_CODE.get(alt_codon, "?")

        # Determine consequence
        if ref_aa == alt_aa:
            consequence = ConsequenceType.SYNONYMOUS
        elif alt_aa == "*":
            consequence = ConsequenceType.NONSENSE
        else:
            consequence = ConsequenceType.MISSENSE

        # Override with splice consequence if very close to boundary
        if splice_prox and consequence == ConsequenceType.SYNONYMOUS:
            min_dist = min(
                coord.distance_to_exon_start, coord.distance_to_exon_end
            )
            if min_dist <= SPLICE_DONOR_ACCEPTOR_THRESHOLD:
                # Even synonymous changes at splice sites are significant
                consequence = (
                    ConsequenceType.SPLICE_DONOR
                    if coord.distance_to_exon_end <= SPLICE_DONOR_ACCEPTOR_THRESHOLD
                    else ConsequenceType.SPLICE_ACCEPTOR
                )

        # Build HGVS notation
        hgvs_c = f"c.{cds_pos}{ref}>{alt}"

        ref_aa_3 = AA_THREE_LETTER.get(ref_aa, ref_aa)
        alt_aa_3 = AA_THREE_LETTER.get(alt_aa, alt_aa)

        if consequence == ConsequenceType.SYNONYMOUS:
            hgvs_p = f"p.{ref_aa_3}{codon_idx}="
        elif consequence == ConsequenceType.NONSENSE:
            hgvs_p = f"p.{ref_aa_3}{codon_idx}{alt_aa_3}"
        elif consequence == ConsequenceType.MISSENSE:
            hgvs_p = f"p.{ref_aa_3}{codon_idx}{alt_aa_3}"
        else:
            hgvs_p = ""

        message_parts = [
            f"{consequence.value} variant",
            f"codon {codon_idx}: {ref_codon} ({ref_aa}) -> {alt_codon} ({alt_aa})",
        ]
        if splice_prox:
            message_parts.append(f"splice proximity: {splice_prox}")

        return CodingAnnotation(
            consequence=consequence,
            hgvs_c=hgvs_c,
            hgvs_p=hgvs_p,
            reference_codon=ref_codon,
            alternate_codon=alt_codon,
            reference_aa=ref_aa,
            alternate_aa=alt_aa,
            codon_index=codon_idx,
            codon_position=codon_pos,
            splice_proximity=splice_prox,
            message="; ".join(message_parts),
        )

    # ----- Indel annotation -------------------------------------------------

    def _annotate_insertion(
        self,
        coord: TranscriptCoordinate,
        ref: str,
        alt: str,
        splice_prox: Optional[str],
    ) -> CodingAnnotation:
        """Annotate an insertion variant."""
        inserted_bases = alt if alt != "-" else ""
        ins_len = len(inserted_bases)
        cds_pos = coord.cds_position

        if ins_len % 3 == 0:
            consequence = ConsequenceType.INFRAME_INSERTION
            msg = f"In-frame insertion of {ins_len} bp ({ins_len // 3} codons)"
        else:
            consequence = ConsequenceType.FRAMESHIFT
            msg = f"Frameshift insertion of {ins_len} bp"

        hgvs_c = f"c.{cds_pos}_{cds_pos + 1}ins{inserted_bases}"

        return CodingAnnotation(
            consequence=consequence,
            hgvs_c=hgvs_c,
            hgvs_p="p.?",  # frameshift protein notation is complex
            codon_index=coord.codon_index,
            codon_position=coord.codon_position,
            splice_proximity=splice_prox,
            message=msg,
        )

    def _annotate_deletion(
        self,
        coord: TranscriptCoordinate,
        ref: str,
        alt: str,
        splice_prox: Optional[str],
    ) -> CodingAnnotation:
        """Annotate a deletion variant."""
        del_len = len(ref) if ref != "-" else 0
        cds_pos = coord.cds_position

        if del_len % 3 == 0:
            consequence = ConsequenceType.INFRAME_DELETION
            msg = f"In-frame deletion of {del_len} bp ({del_len // 3} codons)"
        else:
            consequence = ConsequenceType.FRAMESHIFT
            msg = f"Frameshift deletion of {del_len} bp"

        if del_len == 1:
            hgvs_c = f"c.{cds_pos}del"
        elif del_len > 1:
            hgvs_c = f"c.{cds_pos}_{cds_pos + del_len - 1}del"
        else:
            hgvs_c = f"c.{cds_pos}del"

        return CodingAnnotation(
            consequence=consequence,
            hgvs_c=hgvs_c,
            hgvs_p="p.?",
            codon_index=coord.codon_index,
            codon_position=coord.codon_position,
            splice_proximity=splice_prox,
            message=msg,
        )

    def _annotate_complex(
        self,
        coord: TranscriptCoordinate,
        ref: str,
        alt: str,
        splice_prox: Optional[str],
    ) -> CodingAnnotation:
        """Annotate a complex (MNV or delins) variant."""
        ref_len = len(ref) if ref != "-" else 0
        alt_len = len(alt) if alt != "-" else 0
        net_change = alt_len - ref_len

        if net_change % 3 == 0:
            consequence = ConsequenceType.INFRAME_INSERTION if net_change > 0 \
                else ConsequenceType.INFRAME_DELETION if net_change < 0 \
                else ConsequenceType.MISSENSE
            msg = f"Complex variant: {ref}>{alt} (net {net_change:+d} bp, in-frame)"
        else:
            consequence = ConsequenceType.FRAMESHIFT
            msg = f"Complex variant: {ref}>{alt} (net {net_change:+d} bp, frameshift)"

        cds_pos = coord.cds_position
        hgvs_c = f"c.{cds_pos}_{cds_pos + ref_len - 1}delins{alt}"

        return CodingAnnotation(
            consequence=consequence,
            hgvs_c=hgvs_c,
            hgvs_p="p.?",
            codon_index=coord.codon_index,
            codon_position=coord.codon_position,
            splice_proximity=splice_prox,
            message=msg,
        )

    # ----- Splice / non-coding helpers --------------------------------------

    def _assess_splice_proximity(
        self, coord: TranscriptCoordinate
    ) -> Optional[str]:
        """Assess whether a position is near a splice junction.

        Returns
        -------
        str or None
            "near_exon_start" or "near_exon_end" if within the splice
            region threshold, otherwise None.
        """
        if coord.exon_number < 1:
            # Intronic — check distance to nearest boundary
            min_dist = min(
                coord.distance_to_exon_start,
                coord.distance_to_exon_end,
            )
            if min_dist <= SPLICE_REGION_THRESHOLD:
                return "intronic_splice_region"
            return None

        # Exonic
        dist_start = coord.distance_to_exon_start
        dist_end = coord.distance_to_exon_end

        if dist_start <= SPLICE_REGION_THRESHOLD:
            return "near_exon_start"
        if dist_end <= SPLICE_REGION_THRESHOLD:
            return "near_exon_end"

        return None

    def _splice_annotation(
        self,
        coord: TranscriptCoordinate,
        ref: str,
        alt: str,
        splice_prox: str,
    ) -> CodingAnnotation:
        """Build annotation for a splice-site variant."""
        dist_start = coord.distance_to_exon_start
        dist_end = coord.distance_to_exon_end
        min_dist = min(dist_start, dist_end)

        if min_dist <= SPLICE_DONOR_ACCEPTOR_THRESHOLD:
            if dist_end <= dist_start:
                consequence = ConsequenceType.SPLICE_DONOR
            else:
                consequence = ConsequenceType.SPLICE_ACCEPTOR
        else:
            consequence = ConsequenceType.SPLICE_REGION

        cds_pos = coord.cds_position
        hgvs_c = f"c.{cds_pos}{ref}>{alt}" if cds_pos > 0 else ""

        return CodingAnnotation(
            consequence=consequence,
            hgvs_c=hgvs_c,
            splice_proximity=splice_prox,
            codon_index=coord.codon_index,
            codon_position=coord.codon_position,
            message=(
                f"{consequence.value}: {min_dist} bp from exon boundary "
                f"(exon {coord.exon_number})"
            ),
        )

    def _non_coding_annotation(
        self,
        coord: TranscriptCoordinate,
        ref: str,
        alt: str,
        splice_prox: Optional[str],
    ) -> CodingAnnotation:
        """Build annotation for non-coding (intronic/UTR) variants."""
        if coord.exon_number == -1:
            # Intronic
            if splice_prox:
                min_dist = min(
                    coord.distance_to_exon_start,
                    coord.distance_to_exon_end,
                )
                if min_dist <= SPLICE_DONOR_ACCEPTOR_THRESHOLD:
                    consequence = ConsequenceType.SPLICE_DONOR
                elif min_dist <= SPLICE_REGION_THRESHOLD:
                    consequence = ConsequenceType.SPLICE_REGION
                else:
                    consequence = ConsequenceType.INTRONIC
            else:
                consequence = ConsequenceType.INTRONIC

            return CodingAnnotation(
                consequence=consequence,
                splice_proximity=splice_prox,
                message=(
                    f"Intronic variant at genomic position "
                    f"{coord.genomic_position}"
                ),
            )

        # Exonic but non-CDS (UTR)
        return CodingAnnotation(
            consequence=ConsequenceType.NON_CODING,
            splice_proximity=splice_prox,
            message=(
                f"Non-coding exonic variant in exon {coord.exon_number}"
            ),
        )


# ---------------------------------------------------------------------------
# Self-test
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    print("=" * 60)
    print("CodingAnnotator — self-test")
    print("=" * 60)

    from core.models import TranscriptCoordinate as TC

    # Build a test CDS: 12 codons = 36 bp
    # ATG GCC AAG TAT TGC CAA GGT GAC TGA xxx xxx
    # M   A   K   Y   C   Q   G   D   *
    cds = "ATGGCCAAGTATGCCAAGGTGACTGAAAAAAAAAA"
    # Note: position 25-27 = TGA = stop codon (*)

    annotator = CodingAnnotator(cds)

    # Test 1: synonymous variant
    # Codon 2 (GCC, Ala), position 3 in codon = C
    # GCC -> GCT = Ala -> Ala (synonymous)
    coord1 = TC(
        genomic_position=100,
        exon_number=1,
        transcript_position=6,
        cds_position=6,
        codon_index=2,
        codon_position=3,
        reference_codon="GCC",
        reference_aa="A",
        distance_to_exon_start=5,
        distance_to_exon_end=94,
    )
    r1 = annotator.annotate(coord1, "C", "T")
    assert r1.consequence == ConsequenceType.SYNONYMOUS, \
        f"Expected SYNONYMOUS, got {r1.consequence}"
    assert r1.hgvs_c == "c.6C>T"
    assert "=" in r1.hgvs_p  # synonymous notation
    print(f"[PASS] Synonymous: {r1.hgvs_c} {r1.hgvs_p} ({r1.message})")

    # Test 2: missense variant
    # Codon 3 (AAG, Lys), position 1 in codon = A
    # AAG -> GAG = Lys -> Glu (missense)
    coord2 = TC(
        genomic_position=200,
        exon_number=1,
        transcript_position=7,
        cds_position=7,
        codon_index=3,
        codon_position=1,
        reference_codon="AAG",
        reference_aa="K",
        distance_to_exon_start=6,
        distance_to_exon_end=93,
    )
    r2 = annotator.annotate(coord2, "A", "G")
    assert r2.consequence == ConsequenceType.MISSENSE, \
        f"Expected MISSENSE, got {r2.consequence}"
    assert r2.alternate_aa == "E"  # Glu
    print(f"[PASS] Missense: {r2.hgvs_c} {r2.hgvs_p} ({r2.message})")

    # Test 3: nonsense variant
    # Codon 3 (AAG, Lys), position 2 in codon = A
    # AAG -> TAG = Lys -> Stop (nonsense)
    coord3 = TC(
        genomic_position=201,
        exon_number=1,
        transcript_position=8,
        cds_position=8,
        codon_index=3,
        codon_position=2,
        reference_codon="AAG",
        reference_aa="K",
        distance_to_exon_start=7,
        distance_to_exon_end=92,
    )
    r3 = annotator.annotate(coord3, "A", "T")
    assert r3.consequence == ConsequenceType.NONSENSE, \
        f"Expected NONSENSE, got {r3.consequence}"
    assert r3.alternate_aa == "*"
    assert "Ter" in r3.hgvs_p
    print(f"[PASS] Nonsense: {r3.hgvs_c} {r3.hgvs_p} ({r3.message})")

    # Test 4: frameshift insertion
    coord4 = TC(
        genomic_position=300,
        exon_number=2,
        transcript_position=10,
        cds_position=10,
        codon_index=4,
        codon_position=1,
        reference_codon="TAT",
        reference_aa="Y",
        distance_to_exon_start=9,
        distance_to_exon_end=40,
    )
    r4 = annotator.annotate(coord4, "-", "A")
    assert r4.consequence == ConsequenceType.FRAMESHIFT
    print(f"[PASS] Frameshift insertion: {r4.hgvs_c} ({r4.message})")

    # Test 5: in-frame deletion (3 bp)
    r5 = annotator.annotate(coord4, "TAT", "-")
    assert r5.consequence == ConsequenceType.INFRAME_DELETION
    print(f"[PASS] In-frame deletion: {r5.hgvs_c} ({r5.message})")

    # Test 6: splice proximity
    coord6 = TC(
        genomic_position=400,
        exon_number=3,
        transcript_position=20,
        cds_position=20,
        codon_index=7,
        codon_position=2,
        reference_codon="GGT",
        reference_aa="G",
        distance_to_exon_start=1,  # 1 bp from exon start
        distance_to_exon_end=50,
    )
    r6 = annotator.annotate(coord6, "G", "A")
    assert r6.splice_proximity is not None
    assert "near_exon_start" in r6.splice_proximity
    print(f"[PASS] Splice proximity: {r6.splice_proximity} ({r6.message})")

    # Test 7: intronic variant
    coord7 = TC(
        genomic_position=500,
        exon_number=-1,
        transcript_position=-1,
        cds_position=-1,
        codon_index=-1,
        codon_position=-1,
        reference_codon="",
        reference_aa="",
        distance_to_exon_start=50,
        distance_to_exon_end=50,
        in_cds=False,
    )
    r7 = annotator.annotate(coord7, "A", "G")
    assert r7.consequence == ConsequenceType.INTRONIC
    print(f"[PASS] Intronic: {r7.message}")

    # Test 8: in-frame insertion (3 bp)
    r8 = annotator.annotate(coord4, "-", "ACG")
    assert r8.consequence == ConsequenceType.INFRAME_INSERTION
    print(f"[PASS] In-frame insertion: {r8.hgvs_c} ({r8.message})")

    print("\nAll CodingAnnotator tests passed.")
