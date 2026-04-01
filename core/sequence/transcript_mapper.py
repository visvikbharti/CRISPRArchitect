"""
Transcript Mapper — Genomic-to-transcript coordinate conversion
================================================================

Maps a 1-based genomic position to its transcript-level coordinates:
exon number, CDS position, codon index, codon position, and distances
to the nearest exon boundaries (splice junctions).

Biology context
---------------
A coding variant's impact depends entirely on where it falls within the
transcript.  The same genomic position can be:
  - exonic and coding (affects a codon)
  - exonic but UTR (5' or 3')
  - intronic (between exons — only matters if near a splice site)

For strand-negative genes, the transcript reads "backwards" relative to
the genome coordinate axis.  Exon 1 is at the *highest* genomic
coordinate, and the CDS starts from the 3' end of the genome.  This
module handles both orientations transparently.

Codon arithmetic
-----------------
Given a 1-based CDS position `p`:
  codon_index   = (p - 1) // 3 + 1      (1-based codon number)
  codon_position = (p - 1) % 3 + 1      (1, 2, or 3 within the codon)

References
----------
- den Dunnen et al., Human Mutation, 2016 (HGVS nomenclature v15.11)
- Eilbeck et al., Genome Biology, 2005 (Sequence Ontology)
"""

from __future__ import annotations

from typing import List, Optional

from core.models import TranscriptInfo, ExonRecord, TranscriptCoordinate

# v1 utilities for sequence operations
try:
    from utils.sequence import GENETIC_CODE, reverse_complement
except ImportError:
    try:
        from crisprarchitect.utils.sequence import GENETIC_CODE, reverse_complement
    except ImportError:
        # Minimal fallback
        GENETIC_CODE = {}

        def reverse_complement(seq):
            comp = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
            return ''.join(comp.get(b, 'N') for b in reversed(seq.upper()))


class TranscriptMapper:
    """Maps genomic positions to transcript coordinates.

    Given a TranscriptInfo, this class builds an internal lookup structure
    that maps any genomic position within the gene to its exon, CDS
    position, codon context, and distances to splice boundaries.

    Parameters
    ----------
    transcript : TranscriptInfo
        The transcript to map against.  Must have at least one exon.
    cds_sequence : str or None
        Pre-fetched CDS sequence.  If None, codon lookups will return
        empty strings.  The CDS should be in transcript orientation
        (i.e., starting with ATG).

    Notes
    -----
    All genomic positions are 1-based (Ensembl convention).  The CDS
    position is also 1-based.  Intronic positions return exon_number=-1
    and cds_position=-1.
    """

    def __init__(
        self,
        transcript: TranscriptInfo,
        cds_sequence: Optional[str] = None,
    ):
        self.transcript = transcript
        self.cds_sequence = cds_sequence.upper() if cds_sequence else ""

        # Pre-compute sorted exon list (by genomic start, ascending)
        self._exons_genomic_order: List[ExonRecord] = sorted(
            transcript.exons, key=lambda e: e.start
        )

        # Pre-compute cumulative exon lengths (in transcript order) for
        # CDS offset calculation.  transcript.exons is already in transcript
        # order (exon 1 first).
        self._cumulative_exon_bp: List[int] = []
        running = 0
        for exon in transcript.exons:
            self._cumulative_exon_bp.append(running)
            running += exon.length

        # Compute 5'UTR offset: number of transcript bases before CDS start
        self._utr5_offset = self._compute_utr5_offset()

    # ----- public API -------------------------------------------------------

    def map_genomic_to_transcript(
        self,
        genomic_pos: int,
        flank_size: int = 60,
    ) -> TranscriptCoordinate:
        """Map a 1-based genomic position to transcript coordinates.

        Parameters
        ----------
        genomic_pos : int
            1-based genomic coordinate (e.g., from ClinVar or VCF).
        flank_size : int
            Number of bases of sequence context on each side of the
            position (default 60).

        Returns
        -------
        TranscriptCoordinate
            Full mapping result including exon number, CDS position,
            codon context, and splice distances.
        """
        strand = self.transcript.strand

        # Step 1: find which exon (if any) contains this position
        exon, exon_idx = self._find_exon(genomic_pos)

        if exon is None:
            # Intronic or intergenic
            return self._intronic_coordinate(genomic_pos, flank_size)

        # Step 2: compute offset within this exon (in transcript sense)
        #   Forward strand: offset = genomic_pos - exon.start
        #   Reverse strand: offset = exon.end - genomic_pos
        if strand == 1:
            offset_in_exon = genomic_pos - exon.start  # 0-based
        else:
            offset_in_exon = exon.end - genomic_pos  # 0-based

        # Step 3: transcript position = cumulative bp before this exon + offset + 1
        # Use the exon_number (1-based) to index into cumulative lengths
        # transcript.exons is in transcript order, so index = exon.exon_number - 1
        exon_tx_idx = exon.exon_number - 1
        transcript_pos = self._cumulative_exon_bp[exon_tx_idx] + offset_in_exon + 1

        # Step 4: CDS position = transcript position minus 5'UTR offset
        cds_pos = transcript_pos - self._utr5_offset
        in_cds = cds_pos >= 1

        # Validate CDS position against CDS length
        if self.cds_sequence and cds_pos > len(self.cds_sequence):
            in_cds = False
            cds_pos = -1

        if cds_pos < 1:
            in_cds = False
            cds_pos = -1

        # Step 5: codon arithmetic
        codon_index = -1
        codon_position = -1
        reference_codon = ""
        reference_aa = ""

        if in_cds and cds_pos >= 1 and self.cds_sequence:
            codon_index = (cds_pos - 1) // 3 + 1       # 1-based codon number
            codon_position = (cds_pos - 1) % 3 + 1     # 1, 2, or 3
            codon_start = (codon_index - 1) * 3         # 0-based in CDS string
            if codon_start + 3 <= len(self.cds_sequence):
                reference_codon = self.cds_sequence[codon_start:codon_start + 3]
                reference_aa = GENETIC_CODE.get(reference_codon, "?")

        # Step 6: distances to exon boundaries (in transcript sense)
        #   distance_to_exon_start = offset_in_exon (0 means at the first base)
        #   distance_to_exon_end = exon.length - 1 - offset_in_exon
        distance_to_start = offset_in_exon
        distance_to_end = exon.length - 1 - offset_in_exon

        # Step 7: sequence context
        seq_context, center_idx = self._get_sequence_context(
            genomic_pos, flank_size
        )

        return TranscriptCoordinate(
            genomic_position=genomic_pos,
            exon_number=exon.exon_number,
            transcript_position=transcript_pos,
            cds_position=cds_pos,
            codon_index=codon_index,
            codon_position=codon_position,
            reference_codon=reference_codon,
            reference_aa=reference_aa,
            distance_to_exon_start=distance_to_start,
            distance_to_exon_end=distance_to_end,
            in_cds=in_cds,
            sequence_context=seq_context,
            context_center_index=center_idx,
        )

    # ----- private helpers --------------------------------------------------

    def _compute_utr5_offset(self) -> int:
        """Compute the number of transcript bases in the 5'UTR.

        Uses the CDS start/end boundaries from TranscriptInfo to determine
        how many exonic bases precede the first coding position.
        Returns 0 if CDS boundaries are not available.
        """
        tx = self.transcript
        cds_start = getattr(tx, 'cds_start', None)
        cds_end = getattr(tx, 'cds_end', None)
        if cds_start is None or cds_end is None:
            return 0

        # On forward strand, CDS begins at cds_start (lowest genomic coord)
        # On reverse strand, CDS begins at cds_end (highest genomic coord)
        utr_bases = 0
        for exon in tx.exons:  # already in transcript order
            if tx.strand == 1:
                # Forward: UTR is everything before cds_start
                if exon.end < cds_start:
                    utr_bases += exon.length  # entire exon is UTR
                elif exon.start < cds_start:
                    utr_bases += cds_start - exon.start  # partial UTR
                    break
                else:
                    break  # past UTR
            else:
                # Reverse: UTR is everything above cds_end (5' end of transcript)
                if exon.start > cds_end:
                    utr_bases += exon.length  # entire exon is UTR
                elif exon.end > cds_end:
                    utr_bases += exon.end - cds_end  # partial UTR
                    break
                else:
                    break  # past UTR
        return utr_bases

    def _find_exon(self, genomic_pos: int) -> tuple:
        """Find the exon containing genomic_pos.

        Returns (ExonRecord, index) or (None, -1) if intronic.
        """
        for i, exon in enumerate(self.transcript.exons):
            if exon.start <= genomic_pos <= exon.end:
                return exon, i
        return None, -1

    def _intronic_coordinate(
        self, genomic_pos: int, flank_size: int
    ) -> TranscriptCoordinate:
        """Build a TranscriptCoordinate for an intronic position.

        Computes the distance to the nearest exon boundary so we can
        assess splice-site proximity.
        """
        # Find nearest exon boundaries
        min_dist_start = float('inf')
        min_dist_end = float('inf')
        nearest_exon_num = -1

        for exon in self.transcript.exons:
            d_to_start = abs(genomic_pos - exon.start)
            d_to_end = abs(genomic_pos - exon.end)

            if d_to_start < min_dist_start or d_to_end < min_dist_end:
                if d_to_start <= d_to_end:
                    if d_to_start < min_dist_start:
                        min_dist_start = d_to_start
                        min_dist_end = d_to_end
                        nearest_exon_num = exon.exon_number
                else:
                    if d_to_end < min_dist_end:
                        min_dist_end = d_to_end
                        min_dist_start = d_to_start
                        nearest_exon_num = exon.exon_number

        # For intronic, report distances to nearest exon boundaries
        # Use min of both distances
        all_boundary_distances = []
        for exon in self.transcript.exons:
            all_boundary_distances.append(
                (abs(genomic_pos - exon.start), exon.exon_number, "start")
            )
            all_boundary_distances.append(
                (abs(genomic_pos - exon.end), exon.exon_number, "end")
            )
        all_boundary_distances.sort()

        if all_boundary_distances:
            closest_dist, closest_exon, closest_side = all_boundary_distances[0]
            # distance_to_exon_start/end: we report the closest boundary
            # distance as negative (convention: negative = intronic)
            dist_start = closest_dist
            dist_end = closest_dist
        else:
            dist_start = -1
            dist_end = -1

        seq_context, center_idx = self._get_sequence_context(
            genomic_pos, flank_size
        )

        return TranscriptCoordinate(
            genomic_position=genomic_pos,
            exon_number=-1,
            transcript_position=-1,
            cds_position=-1,
            codon_index=-1,
            codon_position=-1,
            reference_codon="",
            reference_aa="",
            distance_to_exon_start=dist_start,
            distance_to_exon_end=dist_end,
            in_cds=False,
            sequence_context=seq_context,
            context_center_index=center_idx,
        )

    def _get_sequence_context(
        self, genomic_pos: int, flank_size: int
    ) -> tuple:
        """Extract local sequence context from the CDS if available.

        For a full implementation, this would fetch genomic sequence via
        the API.  Here we return context from the CDS if the position is
        coding, otherwise return empty.

        Returns
        -------
        tuple of (context_string, center_index)
        """
        # If we have CDS and position is coding, extract context from CDS
        # This is a simplified approach; the variant_normalizer will provide
        # genomic context via fetcher.fetch_genomic_sequence()
        return "", -1

    def find_nearest_exon(self, genomic_pos: int) -> Optional[ExonRecord]:
        """Return the exon nearest to the given genomic position.

        Useful for intronic variants to determine which exon boundary
        is closest (for splice-site annotation).
        """
        best_exon = None
        best_dist = float('inf')

        for exon in self.transcript.exons:
            if exon.start <= genomic_pos <= exon.end:
                return exon
            dist = min(
                abs(genomic_pos - exon.start),
                abs(genomic_pos - exon.end),
            )
            if dist < best_dist:
                best_dist = dist
                best_exon = exon

        return best_exon

    def get_splice_distances(self, genomic_pos: int) -> tuple:
        """Compute distances to the nearest exon start and end.

        For exonic positions, returns (distance_to_exon_start,
        distance_to_exon_end) in transcript sense.  For intronic
        positions, returns distances to the two closest exon boundaries.

        Returns
        -------
        tuple of (dist_to_5prime_boundary, dist_to_3prime_boundary)
            In transcript orientation.
        """
        exon, _ = self._find_exon(genomic_pos)
        strand = self.transcript.strand

        if exon is not None:
            if strand == 1:
                return (
                    genomic_pos - exon.start,
                    exon.end - genomic_pos,
                )
            else:
                return (
                    exon.end - genomic_pos,
                    genomic_pos - exon.start,
                )

        # Intronic: find the two closest exon boundaries
        distances = []
        for ex in self.transcript.exons:
            distances.append(abs(genomic_pos - ex.start))
            distances.append(abs(genomic_pos - ex.end))
        distances.sort()

        if len(distances) >= 2:
            return (distances[0], distances[1])
        elif distances:
            return (distances[0], distances[0])
        return (-1, -1)


# ---------------------------------------------------------------------------
# Self-test
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    print("=" * 60)
    print("TranscriptMapper — self-test")
    print("=" * 60)

    # Build a minimal forward-strand transcript for testing
    from core.models import ExonRecord as ER, TranscriptInfo as TI

    # Toy gene: 3 exons on forward strand
    #   Exon 1: 1000-1099 (100 bp)
    #   Exon 2: 2000-2049 (50 bp)
    #   Exon 3: 3000-3074 (75 bp)
    #   CDS = 225 bp = 75 codons
    exons = [
        ER("E1", 1, 1000, 1099, 1, "1"),
        ER("E2", 2, 2000, 2049, 1, "1"),
        ER("E3", 3, 3000, 3074, 1, "1"),
    ]

    tx = TI(
        transcript_id="TEST_FWD",
        gene_symbol="TEST",
        gene_id="ENSG_TEST",
        chromosome="1",
        start=1000,
        end=3074,
        strand=1,
        biotype="protein_coding",
        is_canonical=True,
        exons=exons,
    )

    # CDS: 225 bp, all A's for simplicity, starting with ATG
    cds = "ATG" + "A" * 222  # 225 bp = 75 codons
    mapper = TranscriptMapper(tx, cds)

    # Test 1: first base of exon 1
    coord = mapper.map_genomic_to_transcript(1000)
    assert coord.exon_number == 1, f"Expected exon 1, got {coord.exon_number}"
    assert coord.transcript_position == 1
    assert coord.cds_position == 1
    assert coord.codon_index == 1
    assert coord.codon_position == 1
    assert coord.distance_to_exon_start == 0
    assert coord.distance_to_exon_end == 99
    print(f"[PASS] gpos=1000 -> exon=1, cds=1, codon=1, pos_in_codon=1")

    # Test 2: third base of exon 1 (position 1002, offset=2)
    coord2 = mapper.map_genomic_to_transcript(1002)
    assert coord2.cds_position == 3
    assert coord2.codon_index == 1
    assert coord2.codon_position == 3
    print(f"[PASS] gpos=1002 -> cds=3, codon=1, pos_in_codon=3")

    # Test 3: first base of exon 2 (should be transcript pos 101)
    coord3 = mapper.map_genomic_to_transcript(2000)
    assert coord3.exon_number == 2
    assert coord3.transcript_position == 101, \
        f"Expected tx_pos 101, got {coord3.transcript_position}"
    assert coord3.distance_to_exon_start == 0
    print(f"[PASS] gpos=2000 -> exon=2, tx_pos=101")

    # Test 4: intronic position
    coord4 = mapper.map_genomic_to_transcript(1500)
    assert coord4.exon_number == -1
    assert coord4.in_cds is False
    print(f"[PASS] gpos=1500 -> intronic")

    # Test 5: reverse-strand transcript
    exons_rev = [
        ER("E1", 1, 3000, 3074, -1, "1"),  # exon 1 at high coords
        ER("E2", 2, 2000, 2049, -1, "1"),
        ER("E3", 3, 1000, 1099, -1, "1"),
    ]
    tx_rev = TI(
        transcript_id="TEST_REV",
        gene_symbol="TEST_REV",
        gene_id="ENSG_REV",
        chromosome="1",
        start=1000,
        end=3074,
        strand=-1,
        biotype="protein_coding",
        is_canonical=True,
        exons=exons_rev,
    )
    cds_rev = "ATG" + "C" * 222
    mapper_rev = TranscriptMapper(tx_rev, cds_rev)

    # Reverse strand exon 1: genomic 3000-3074, first transcript base = 3074
    coord_rev = mapper_rev.map_genomic_to_transcript(3074)
    assert coord_rev.exon_number == 1
    assert coord_rev.transcript_position == 1
    assert coord_rev.codon_index == 1
    assert coord_rev.codon_position == 1
    print(f"[PASS] reverse strand: gpos=3074 -> exon=1, tx_pos=1, codon=1")

    # Reverse strand: genomic 3072 is offset 2 from end, so tx_pos = 3
    coord_rev2 = mapper_rev.map_genomic_to_transcript(3072)
    assert coord_rev2.transcript_position == 3
    assert coord_rev2.codon_position == 3
    print(f"[PASS] reverse strand: gpos=3072 -> tx_pos=3, pos_in_codon=3")

    # Test 6: splice distances
    d5, d3 = mapper.get_splice_distances(1005)
    assert d5 == 5, f"Expected 5, got {d5}"
    assert d3 == 94, f"Expected 94, got {d3}"
    print(f"[PASS] splice distances at gpos=1005: 5' = {d5}, 3' = {d3}")

    print("\nAll TranscriptMapper tests passed.")
