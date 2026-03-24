"""
Variant Normalizer — End-to-end variant annotation pipeline
=============================================================

Orchestrates the full variant normalization workflow:

    GenomicVariantInput
         |
         v
    TranscriptFetcher   -->  TranscriptInfo
         |
         v
    TranscriptMapper    -->  TranscriptCoordinate
         |
         v
    ReferenceValidator  -->  ReferenceValidation
         |
         v
    CodingAnnotator     -->  CodingAnnotation
         |
         v
    NormalizedVariant  (+ v1 Mutation bridge)

This is the single entry point that downstream modules (feasibility
assessors, strategy generators) call to get a fully annotated variant.

v1 compatibility
-----------------
The NormalizedVariant includes a v1_mutation field that holds a v1
Mutation dataclass for backward compatibility with v1's
StrategyEnumerator and StrategyScorer.

CRITICAL: v1 Mutation uses 0-based positions; v2 uses 1-based.  The
bridge conversion subtracts 1 explicitly.

References
----------
- den Dunnen et al., Human Mutation, 2016 (HGVS nomenclature)
- McLaren et al., Genome Biology, 2016 (VEP annotation pipeline)
"""

from __future__ import annotations

import logging
from typing import Optional

from core.models import (
    GenomicVariantInput,
    NormalizedVariant,
    TranscriptInfo,
    TranscriptCoordinate,
    CodingAnnotation,
    ReferenceValidation,
    ConsequenceType,
)
from core.sequence.fetcher import TranscriptFetcher
from core.sequence.transcript_mapper import TranscriptMapper
from core.sequence.reference_validator import ReferenceValidator
from core.sequence.coding_annotation import CodingAnnotator

# v1 bridge — import Mutation for backward compatibility
try:
    from mosaic.mutation_classifier import Mutation as V1Mutation
except ImportError:
    try:
        from crisprarchitect.mosaic.mutation_classifier import Mutation as V1Mutation
    except ImportError:
        V1Mutation = None

# v1 sequence utilities
try:
    from utils.sequence import reverse_complement
except ImportError:
    try:
        from crisprarchitect.utils.sequence import reverse_complement
    except ImportError:
        def reverse_complement(seq):
            comp = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
            return ''.join(comp.get(b, 'N') for b in reversed(seq.upper()))


logger = logging.getLogger(__name__)


class VariantNormalizer:
    """End-to-end variant annotation pipeline.

    Given a GenomicVariantInput (raw user input), this class:
    1. Fetches the transcript structure from Ensembl
    2. Maps the genomic position to transcript coordinates
    3. Validates the reference allele against the genome
    4. Annotates the coding consequence
    5. Creates a v1-compatible Mutation object for backward compat

    The result is a NormalizedVariant that contains everything needed
    for editing feasibility assessment.

    Parameters
    ----------
    fetcher : TranscriptFetcher or None
        Pre-configured fetcher.  If None, a new one is created.
    species : str
        Species name for Ensembl lookups.

    Examples
    --------
    >>> normalizer = VariantNormalizer()
    >>> variant_input = GenomicVariantInput(
    ...     chromosome="7", position=117559593,
    ...     ref_allele="G", alt_allele="A",
    ...     gene_symbol="CFTR"
    ... )
    >>> result = normalizer.normalize(variant_input)
    >>> print(result.coding.consequence, result.coding.hgvs_c)
    """

    def __init__(
        self,
        fetcher: Optional[TranscriptFetcher] = None,
        species: str = "homo_sapiens",
    ):
        self.fetcher = fetcher or TranscriptFetcher(species=species)
        self.validator = ReferenceValidator(self.fetcher)
        self._transcript_cache: dict = {}

    def normalize(
        self,
        variant: GenomicVariantInput,
        flank_size: int = 60,
    ) -> NormalizedVariant:
        """Normalize and annotate a genomic variant.

        This is the main public method.  It runs the full pipeline:
        fetch transcript -> map coordinates -> validate ref -> annotate
        coding consequence -> build v1 Mutation bridge.

        Parameters
        ----------
        variant : GenomicVariantInput
            Raw variant input with chromosome, position, ref/alt alleles,
            and gene symbol or transcript ID.
        flank_size : int
            Number of flanking bases to include in sequence context
            (default: 60, giving a 121-bp window).

        Returns
        -------
        NormalizedVariant
            Fully annotated variant with all fields populated.

        Raises
        ------
        ValueError
            If neither gene_symbol nor transcript_id is provided.
        """
        # Step 0: validate input
        if not variant.gene_symbol and not variant.transcript_id:
            raise ValueError(
                "GenomicVariantInput must have gene_symbol or transcript_id"
            )

        # Step 1: fetch transcript
        transcript = self._get_transcript(variant)

        # Step 2: fetch CDS sequence
        cds_seq = self._fetch_cds(transcript.transcript_id)

        # Step 3: map genomic position to transcript coordinates
        mapper = TranscriptMapper(transcript, cds_seq)
        coord = mapper.map_genomic_to_transcript(
            variant.position, flank_size=flank_size
        )

        # Step 4: validate reference allele
        ref_validation = self.validator.validate(
            chromosome=variant.chromosome,
            position=variant.position,
            ref_allele=variant.ref_allele,
            strand=transcript.strand,
        )

        # Step 5: determine alleles in transcript orientation
        ref_tx, alt_tx = self._alleles_in_transcript_orientation(
            variant.ref_allele, variant.alt_allele, transcript.strand
        )

        # Step 6: annotate coding consequence
        annotator = CodingAnnotator(cds_seq)
        coding = annotator.annotate(coord, ref_tx, alt_tx)

        # Step 7: fetch local genomic sequence context
        local_seq, edit_idx = self._fetch_local_sequence(
            variant, flank_size
        )

        # Step 8: build v1 Mutation bridge
        v1_mut = self._build_v1_mutation(variant, coord, coding)

        return NormalizedVariant(
            input=variant,
            transcript=transcript,
            transcript_coord=coord,
            coding=coding,
            ref_validation=ref_validation,
            v1_mutation=v1_mut,
            local_sequence=local_seq,
            local_seq_edit_index=edit_idx,
        )

    # ----- Transcript fetching ----------------------------------------------

    def _get_transcript(self, variant: GenomicVariantInput) -> TranscriptInfo:
        """Fetch and cache the transcript for a variant."""
        # Use transcript_id if provided, otherwise look up by gene symbol
        if variant.transcript_id:
            cache_key = variant.transcript_id
            if cache_key not in self._transcript_cache:
                self._transcript_cache[cache_key] = \
                    self.fetcher.fetch_by_transcript_id(variant.transcript_id)
            return self._transcript_cache[cache_key]

        cache_key = variant.gene_symbol.upper()
        if cache_key not in self._transcript_cache:
            self._transcript_cache[cache_key] = \
                self.fetcher.fetch_by_gene(variant.gene_symbol)
        return self._transcript_cache[cache_key]

    def _fetch_cds(self, transcript_id: str) -> str:
        """Fetch CDS sequence, returning empty string on failure."""
        try:
            return self.fetcher.fetch_cds_sequence(transcript_id)
        except Exception as e:
            logger.warning(
                "Could not fetch CDS for %s: %s", transcript_id, e
            )
            return ""

    # ----- Allele orientation -----------------------------------------------

    @staticmethod
    def _alleles_in_transcript_orientation(
        ref_allele: str,
        alt_allele: str,
        strand: int,
    ) -> tuple:
        """Convert alleles to transcript orientation if on reverse strand.

        For reverse-strand genes, alleles given in genomic orientation need
        to be reverse-complemented.  If alleles are already in transcript
        orientation (common for ClinVar / HGVS input), this is handled by
        the reference validator upstream.

        We assume the user provides alleles in transcript orientation (the
        HGVS convention), so for reverse-strand genes we do NOT flip them.
        The ReferenceValidator will catch mismatches.

        Parameters
        ----------
        ref_allele, alt_allele : str
            User-provided alleles.
        strand : int
            Transcript strand (1 or -1).

        Returns
        -------
        tuple of (ref_transcript, alt_transcript)
        """
        # Special cases: insertion/deletion markers
        ref = ref_allele.upper()
        alt = alt_allele.upper()

        if ref == "-" or alt == "-":
            return ref, alt

        # We assume alleles are in transcript orientation (HGVS convention).
        # If the user provided genomic-orientation alleles for a reverse-strand
        # gene, the reference validator will flag the mismatch and the user
        # can correct. We do not auto-flip to avoid silent errors.
        return ref, alt

    # ----- Local sequence context -------------------------------------------

    def _fetch_local_sequence(
        self, variant: GenomicVariantInput, flank_size: int
    ) -> tuple:
        """Fetch genomic sequence around the variant position.

        Returns
        -------
        tuple of (sequence_string, edit_index)
            The edit_index is the 0-based index of the variant position
            within the returned sequence.
        """
        try:
            start = max(1, variant.position - flank_size)
            end = variant.position + flank_size
            seq = self.fetcher.fetch_genomic_sequence(
                chromosome=variant.chromosome,
                start=start,
                end=end,
                strand=1,  # always fetch forward strand
            )
            edit_idx = variant.position - start  # 0-based index
            return seq, edit_idx
        except Exception as e:
            logger.warning(
                "Could not fetch local sequence for chr%s:%d: %s",
                variant.chromosome, variant.position, e
            )
            return "", -1

    # ----- v1 bridge --------------------------------------------------------

    @staticmethod
    def _build_v1_mutation(
        variant: GenomicVariantInput,
        coord: TranscriptCoordinate,
        coding: CodingAnnotation,
    ) -> object:
        """Create a v1 Mutation object for backward compatibility.

        v1 Mutation fields:
          - exon_number: int (1-based)
          - position: int (0-BASED — critical difference from v2!)
          - ref_allele: str
          - alt_allele: str
          - mutation_type: str (auto-classified)
          - name: str

        The position conversion is:
          v1_position = v2_genomic_position - 1
        """
        if V1Mutation is None:
            logger.debug(
                "v1 Mutation class not available; skipping bridge"
            )
            return None

        exon_num = coord.exon_number if coord.exon_number > 0 else 1

        # CRITICAL: v1 uses 0-based positions, v2 uses 1-based
        v1_position = variant.position - 1

        # Build a descriptive name
        name = variant.name or coding.hgvs_c or (
            f"chr{variant.chromosome}:{variant.position}"
            f"{variant.ref_allele}>{variant.alt_allele}"
        )

        try:
            v1_mut = V1Mutation(
                exon_number=exon_num,
                position=v1_position,
                ref_allele=variant.ref_allele,
                alt_allele=variant.alt_allele,
                name=name,
            )
            return v1_mut
        except Exception as e:
            logger.warning("Could not create v1 Mutation: %s", e)
            return None


# ---------------------------------------------------------------------------
# Self-test
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    print("=" * 60)
    print("VariantNormalizer — self-test")
    print("=" * 60)

    # ---- Offline tests with mock objects ----

    # Create a minimal mock fetcher that doesn't hit the network
    class MockTranscriptFetcher:
        """Provides canned responses for offline testing."""

        def __init__(self):
            from core.models import ExonRecord as ER, TranscriptInfo as TI

            self._transcript = TI(
                transcript_id="ENST00000003084",
                gene_symbol="CFTR",
                gene_id="ENSG00000001626",
                chromosome="7",
                start=117480025,
                end=117668665,
                strand=1,
                biotype="protein_coding",
                is_canonical=True,
                exons=[
                    ER("E1", 1, 117480025, 117480124, 1, "7"),   # 100 bp
                    ER("E2", 2, 117530895, 117531024, 1, "7"),   # 130 bp
                    ER("E3", 3, 117559500, 117559650, 1, "7"),   # 151 bp
                ],
            )
            self.species = "homo_sapiens"

        def fetch_by_gene(self, gene_symbol):
            if gene_symbol.upper() == "CFTR":
                return self._transcript
            raise Exception(f"Mock: unknown gene {gene_symbol}")

        def fetch_by_transcript_id(self, tid):
            return self._transcript

        def fetch_cds_sequence(self, tid):
            # 381 bp = 100 + 130 + 151
            return "ATG" + "GCC" * 126  # 381 bp

        def fetch_genomic_sequence(self, chromosome, start, end, strand=1):
            length = end - start + 1
            # Return a plausible sequence
            return "ACGT" * (length // 4 + 1)[:length] if length > 0 else ""

    mock_fetcher = MockTranscriptFetcher()
    normalizer = VariantNormalizer(fetcher=mock_fetcher)

    # Test 1: basic normalization
    variant_input = GenomicVariantInput(
        chromosome="7",
        position=117559593,
        ref_allele="G",
        alt_allele="A",
        gene_symbol="CFTR",
        name="test_variant",
    )

    result = normalizer.normalize(variant_input, flank_size=30)

    assert result.transcript.gene_symbol == "CFTR"
    assert result.transcript_coord.exon_number == 3, \
        f"Expected exon 3, got {result.transcript_coord.exon_number}"
    assert result.input is variant_input
    assert result.local_sequence != ""
    print(f"[PASS] Basic normalization: exon={result.transcript_coord.exon_number}, "
          f"cds_pos={result.transcript_coord.cds_position}")
    print(f"       Consequence: {result.coding.consequence.value}")
    print(f"       HGVS: {result.coding.hgvs_c}")
    print(f"       Local seq length: {len(result.local_sequence)}")

    # Test 2: v1 Mutation bridge
    if V1Mutation is not None:
        assert result.v1_mutation is not None
        # v1 uses 0-based positions
        assert result.v1_mutation.position == variant_input.position - 1, \
            (f"v1 position should be {variant_input.position - 1}, "
             f"got {result.v1_mutation.position}")
        print(f"[PASS] v1 bridge: position={result.v1_mutation.position} (0-based), "
              f"type={result.v1_mutation.mutation_type}")
    else:
        print("[SKIP] v1 Mutation not importable — bridge test skipped")

    # Test 3: intronic position
    variant_intronic = GenomicVariantInput(
        chromosome="7",
        position=117500000,  # between exon 1 and exon 2
        ref_allele="A",
        alt_allele="G",
        gene_symbol="CFTR",
    )
    result_int = normalizer.normalize(variant_intronic)
    assert result_int.transcript_coord.exon_number == -1
    assert not result_int.transcript_coord.in_cds
    assert result_int.coding.consequence in (
        ConsequenceType.INTRONIC,
        ConsequenceType.SPLICE_REGION,
        ConsequenceType.SPLICE_DONOR,
        ConsequenceType.NON_CODING,
    )
    print(f"[PASS] Intronic: consequence={result_int.coding.consequence.value}")

    # Test 4: insertion
    variant_ins = GenomicVariantInput(
        chromosome="7",
        position=117559593,
        ref_allele="-",
        alt_allele="ACG",
        gene_symbol="CFTR",
    )
    result_ins = normalizer.normalize(variant_ins)
    # Insertion ref is "-", so ref validation should pass trivially
    assert result_ins.ref_validation.is_valid
    print(f"[PASS] Insertion: consequence={result_ins.coding.consequence.value}")

    # Test 5: transcript_id lookup
    variant_by_tx = GenomicVariantInput(
        chromosome="7",
        position=117559593,
        ref_allele="G",
        alt_allele="A",
        transcript_id="ENST00000003084",
    )
    result_tx = normalizer.normalize(variant_by_tx)
    assert result_tx.transcript.transcript_id == "ENST00000003084"
    print(f"[PASS] Transcript ID lookup: {result_tx.transcript.transcript_id}")

    # Test 6: missing gene symbol and transcript_id should raise ValueError
    variant_bad = GenomicVariantInput(
        chromosome="7",
        position=100,
        ref_allele="A",
        alt_allele="G",
    )
    try:
        normalizer.normalize(variant_bad)
        assert False, "Should have raised ValueError"
    except ValueError as e:
        print(f"[PASS] Missing gene/transcript raises ValueError: {e}")

    # ---- Online test (optional, skipped if no internet) ----
    try:
        print("\n--- Optional online test ---")
        live_normalizer = VariantNormalizer()

        # CFTR p.Gly551Asp (G551D) — a well-known CF mutation
        # chr7:117559593 G>A (GRCh38)
        live_variant = GenomicVariantInput(
            chromosome="7",
            position=117559593,
            ref_allele="G",
            alt_allele="A",
            gene_symbol="CFTR",
            name="CFTR G551D",
        )
        live_result = live_normalizer.normalize(live_variant)
        print(f"  Gene: {live_result.transcript.gene_symbol}")
        print(f"  Transcript: {live_result.transcript.transcript_id}")
        print(f"  Exon: {live_result.transcript_coord.exon_number}")
        print(f"  CDS pos: {live_result.transcript_coord.cds_position}")
        print(f"  Consequence: {live_result.coding.consequence.value}")
        print(f"  HGVS c: {live_result.coding.hgvs_c}")
        print(f"  HGVS p: {live_result.coding.hgvs_p}")
        print(f"  Ref valid: {live_result.ref_validation.is_valid}")
        print(f"  Ref msg: {live_result.ref_validation.message}")
        if live_result.v1_mutation:
            print(f"  v1 Mutation: pos={live_result.v1_mutation.position} "
                  f"(0-based), type={live_result.v1_mutation.mutation_type}")
        print("[PASS] Online normalization succeeded")
    except Exception as e:
        print(f"[SKIP] Online test: {e}")

    print("\nAll VariantNormalizer tests passed.")
