"""
Reference Validator — Verify user-provided ref alleles against Ensembl
======================================================================

Before annotating a variant, we must confirm that the user's reference
allele matches what the genome actually has at that position.  Mismatches
indicate:
  - wrong coordinate (off-by-one, wrong chromosome)
  - wrong strand (user gave the complement)
  - wrong genome build (GRCh37 vs GRCh38)
  - a genuine error in the input data

For reverse-strand transcripts, the user typically provides the allele in
transcript (coding) orientation, but the genome stores the forward-strand
base.  This validator checks both orientations and reports which match.

Biology context
---------------
The reference genome (GRCh38) stores all sequences on the forward strand.
When a gene is on the reverse strand, its coding sequence is the reverse
complement of the genomic sequence.  ClinVar and HGVS notation report
variants in the transcript (coding) orientation, so for reverse-strand
genes the alleles are complemented relative to the genome.

Example:
  Gene NF1 is on chromosome 17 reverse strand.
  If the transcript shows "c.910C>T", the genomic forward strand has G>A.

References
----------
- den Dunnen et al., Human Mutation, 2016 (HGVS nomenclature)
- Dalgleish et al., Human Mutation, 2010 (reference sequence conventions)
"""

from __future__ import annotations

from typing import Optional

from core.models import ReferenceValidation

# v1 utilities
try:
    from utils.sequence import complement, reverse_complement
except ImportError:
    try:
        from crisprarchitect.utils.sequence import complement, reverse_complement
    except ImportError:
        def complement(base):
            comp = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
            return comp.get(base.upper(), 'N')

        def reverse_complement(seq):
            return ''.join(complement(b) for b in reversed(seq))


class ReferenceValidator:
    """Validates that a user-provided ref allele matches the genome.

    This class uses a fetcher to retrieve the actual genomic sequence at
    the variant position and compares it with the user-provided reference
    allele.  It handles reverse-strand genes by checking the complement.

    Parameters
    ----------
    fetcher : TranscriptFetcher
        An instance of TranscriptFetcher for sequence retrieval.

    Examples
    --------
    >>> from core.sequence.fetcher import TranscriptFetcher
    >>> fetcher = TranscriptFetcher()
    >>> validator = ReferenceValidator(fetcher)
    >>> result = validator.validate(
    ...     chromosome="7", position=117559593, ref_allele="G",
    ...     strand=1, ref_length=1
    ... )
    >>> print(result.is_valid, result.message)
    """

    def __init__(self, fetcher):
        """
        Parameters
        ----------
        fetcher : TranscriptFetcher
            Used to fetch genomic sequence from Ensembl.
        """
        self.fetcher = fetcher

    def validate(
        self,
        chromosome: str,
        position: int,
        ref_allele: str,
        strand: int = 1,
        ref_length: Optional[int] = None,
    ) -> ReferenceValidation:
        """Validate the reference allele at a genomic position.

        Parameters
        ----------
        chromosome : str
            Chromosome (e.g., "17", "X").
        position : int
            1-based genomic coordinate.
        ref_allele : str
            User-provided reference allele (in transcript orientation
            for reverse-strand genes, or genomic orientation for
            forward-strand genes).
        strand : int
            Transcript strand: 1 (forward) or -1 (reverse).
        ref_length : int or None
            Length of reference allele to fetch.  If None, uses
            len(ref_allele).  Set explicitly for insertions where
            ref_allele is "-".

        Returns
        -------
        ReferenceValidation
            Contains is_valid, expected vs. provided alleles on both
            genomic and transcript strands, and a human-readable message.
        """
        # Handle insertion (ref_allele is "-" or empty)
        if ref_allele in ("-", ""):
            return ReferenceValidation(
                is_valid=True,
                expected_ref_genomic="-",
                provided_ref_genomic="-",
                expected_ref_transcript="-",
                provided_ref_transcript="-",
                message="Insertion — no reference allele to validate.",
            )

        allele_len = ref_length if ref_length is not None else len(ref_allele)
        ref_allele_upper = ref_allele.upper()

        # Fetch the genomic sequence at this position (always forward strand)
        try:
            genomic_ref = self.fetcher.fetch_genomic_sequence(
                chromosome=chromosome,
                start=position,
                end=position + allele_len - 1,
                strand=1,  # always fetch forward strand
            )
        except Exception as e:
            return ReferenceValidation(
                is_valid=False,
                expected_ref_genomic="?",
                provided_ref_genomic=ref_allele_upper,
                expected_ref_transcript="?",
                provided_ref_transcript=ref_allele_upper,
                message=f"Failed to fetch reference sequence: {e}",
            )

        genomic_ref = genomic_ref.upper()

        # Compute transcript-strand reference
        if strand == -1:
            transcript_ref = reverse_complement(genomic_ref)
        else:
            transcript_ref = genomic_ref

        # The user's allele could be in genomic or transcript orientation.
        # Check both:
        #   1. User allele matches genomic forward strand directly
        #   2. User allele matches transcript strand (reverse-complemented
        #      for reverse-strand genes)
        matches_genomic = (ref_allele_upper == genomic_ref)
        matches_transcript = (ref_allele_upper == transcript_ref)

        is_valid = matches_genomic or matches_transcript

        # Determine what the user provided in genomic orientation
        if matches_transcript and not matches_genomic and strand == -1:
            # User gave transcript-orientation allele
            provided_genomic = reverse_complement(ref_allele_upper)
            provided_transcript = ref_allele_upper
        else:
            provided_genomic = ref_allele_upper
            if strand == -1:
                provided_transcript = reverse_complement(ref_allele_upper)
            else:
                provided_transcript = ref_allele_upper

        # Build message
        if is_valid:
            if matches_genomic and matches_transcript:
                msg = (
                    f"Reference allele '{ref_allele_upper}' matches genomic "
                    f"sequence at chr{chromosome}:{position}."
                )
            elif matches_genomic:
                msg = (
                    f"Reference allele '{ref_allele_upper}' matches genomic "
                    f"forward strand at chr{chromosome}:{position}."
                )
            else:
                msg = (
                    f"Reference allele '{ref_allele_upper}' matches transcript "
                    f"strand (reverse complement of genomic '{genomic_ref}') "
                    f"at chr{chromosome}:{position}."
                )
        else:
            msg = (
                f"MISMATCH: provided ref='{ref_allele_upper}' does not match "
                f"genomic='{genomic_ref}' or transcript='{transcript_ref}' "
                f"at chr{chromosome}:{position} (strand={strand}). "
                f"Check coordinates, strand, or genome build."
            )

        return ReferenceValidation(
            is_valid=is_valid,
            expected_ref_genomic=genomic_ref,
            provided_ref_genomic=provided_genomic,
            expected_ref_transcript=transcript_ref,
            provided_ref_transcript=provided_transcript,
            message=msg,
        )


# ---------------------------------------------------------------------------
# Self-test
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    print("=" * 60)
    print("ReferenceValidator — self-test")
    print("=" * 60)

    # Test 1: Mock fetcher for offline testing
    class MockFetcher:
        """Returns pre-defined sequences for testing."""

        def __init__(self):
            self._sequences = {
                # chr7:117559593 (CFTR, forward strand) = "G"
                ("7", 117559593, 117559593, 1): "G",
                # chr17:31232193 (NF1, reverse strand) genomic fwd = "G"
                # transcript strand (reverse complement) = "C"
                ("17", 31232193, 31232193, 1): "G",
                # Multi-base example
                ("1", 100, 102, 1): "ACG",
            }

        def fetch_genomic_sequence(self, chromosome, start, end, strand=1):
            key = (chromosome, start, end, 1)  # always forward
            if key in self._sequences:
                return self._sequences[key]
            raise Exception(f"Mock: no data for {key}")

    mock = MockFetcher()
    validator = ReferenceValidator(mock)

    # Test 1: forward strand, correct allele
    r1 = validator.validate("7", 117559593, "G", strand=1)
    assert r1.is_valid, f"Expected valid, got: {r1.message}"
    print(f"[PASS] Forward strand correct: {r1.message}")

    # Test 2: forward strand, wrong allele
    r2 = validator.validate("7", 117559593, "T", strand=1)
    assert not r2.is_valid, f"Expected invalid, got: {r2.message}"
    print(f"[PASS] Forward strand mismatch: {r2.message}")

    # Test 3: reverse strand, user provides transcript-orientation allele
    # NF1 genomic = G, transcript = C (reverse complement)
    r3 = validator.validate("17", 31232193, "C", strand=-1)
    assert r3.is_valid, f"Expected valid (transcript match), got: {r3.message}"
    print(f"[PASS] Reverse strand transcript match: {r3.message}")

    # Test 4: reverse strand, user provides genomic-orientation allele
    r4 = validator.validate("17", 31232193, "G", strand=-1)
    assert r4.is_valid, f"Expected valid (genomic match), got: {r4.message}"
    print(f"[PASS] Reverse strand genomic match: {r4.message}")

    # Test 5: insertion (ref = "-")
    r5 = validator.validate("1", 100, "-", strand=1)
    assert r5.is_valid
    print(f"[PASS] Insertion validation: {r5.message}")

    # Test 6: multi-base reference
    r6 = validator.validate("1", 100, "ACG", strand=1, ref_length=3)
    assert r6.is_valid
    print(f"[PASS] Multi-base ref: {r6.message}")

    # Test 7: fetch failure
    r7 = validator.validate("99", 1, "A", strand=1)
    assert not r7.is_valid
    assert "Failed" in r7.message
    print(f"[PASS] Fetch failure handled: {r7.message}")

    print("\nAll ReferenceValidator tests passed.")
