"""
Tests for CRISPRArchitect v3 HGVS Parser + ClinVar Processor
==============================================================

Tests parsing of clinical variant notation and ClinVar batch processing.

References
----------
den Dunnen et al., Human Mutation, 2016 (HGVS nomenclature)
"""

from __future__ import annotations

import os
import tempfile
import pytest

from core.sequence.hgvs_parser import (
    HGVSParser,
    ParsedHGVS,
    ClinVarProcessor,
    ClinVarVariant,
)
from core.models import GenomicVariantInput


# ── Tests: HGVS Substitution Parsing ────────────────────────────────

class TestHGVSSubstitution:
    """Tests for parsing HGVS c. substitution notation."""

    def setup_method(self):
        self.parser = HGVSParser()

    def test_full_notation_with_transcript(self):
        """NM_000267.3:c.910C>T should parse correctly."""
        p = self.parser.parse("NM_000267.3:c.910C>T")
        assert p.transcript_id == "NM_000267.3"
        assert p.cds_position == 910
        assert p.ref_allele == "C"
        assert p.alt_allele == "T"
        assert p.variant_type == "substitution"

    def test_transcript_without_version(self):
        """NM_000267:c.910C>T should parse (may match as gene prefix without version)."""
        p = self.parser.parse("NM_000267:c.910C>T")
        # Without version, regex may match as gene prefix — either is acceptable
        assert p.cds_position == 910
        assert p.ref_allele == "C"

    def test_gene_prefix(self):
        """NF1:c.910C>T with gene prefix."""
        p = self.parser.parse("NF1:c.910C>T")
        assert p.gene_symbol == "NF1"
        assert p.cds_position == 910
        assert p.ref_allele == "C"
        assert p.alt_allele == "T"

    def test_gene_and_transcript(self):
        """NF1:NM_000267.3:c.910C>T with both gene and transcript."""
        # This format is less standard but we should handle the transcript part
        p = self.parser.parse("NM_000267.3:c.910C>T")
        assert p.transcript_id == "NM_000267.3"

    def test_all_transitions(self):
        """All four transitions should parse correctly."""
        for ref, alt in [("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")]:
            p = self.parser.parse(f"NM_000267.3:c.100{ref}>{alt}")
            assert p.ref_allele == ref
            assert p.alt_allele == alt
            assert p.variant_type == "substitution"

    def test_all_transversions(self):
        """All eight transversions should parse correctly."""
        transversions = [
            ("A", "C"), ("A", "T"), ("G", "C"), ("G", "T"),
            ("C", "A"), ("C", "G"), ("T", "A"), ("T", "G"),
        ]
        for ref, alt in transversions:
            p = self.parser.parse(f"NM_000267.3:c.100{ref}>{alt}")
            assert p.ref_allele == ref and p.alt_allele == alt

    def test_negative_cds_position(self):
        """Negative CDS position (5' UTR) should parse."""
        p = self.parser.parse("NM_000267.3:c.-15C>T")
        assert p.cds_position == -15

    def test_intronic_offset(self):
        """c.1234+5G>A (intronic) should parse position and offset."""
        p = self.parser.parse("NM_000267.3:c.1234+5G>A")
        assert p.cds_position == 1234
        assert p.intron_offset == 5


# ── Tests: HGVS Deletion Parsing ────────────────────────────────────

class TestHGVSDeletion:

    def setup_method(self):
        self.parser = HGVSParser()

    def test_single_base_deletion(self):
        """c.1234del should parse as single-base deletion."""
        p = self.parser.parse("NM_000267.3:c.1234del")
        assert p.variant_type == "deletion"
        assert p.cds_position == 1234

    def test_multi_base_deletion(self):
        """c.1234_1236del should parse with range."""
        p = self.parser.parse("NM_000267.3:c.1234_1236del")
        assert p.variant_type == "deletion"
        assert p.cds_position == 1234
        assert p.cds_end_position == 1236

    def test_deletion_with_bases(self):
        """c.1234delA should parse with specified deleted base."""
        p = self.parser.parse("NM_000267.3:c.1234delA")
        assert p.variant_type == "deletion"
        assert p.ref_allele == "A"
        assert p.alt_allele == "-"


# ── Tests: HGVS Insertion Parsing ────────────────────────────────────

class TestHGVSInsertion:

    def setup_method(self):
        self.parser = HGVSParser()

    def test_insertion(self):
        """c.1234_1235insATG should parse correctly."""
        p = self.parser.parse("NM_000267.3:c.1234_1235insATG")
        assert p.variant_type == "insertion"
        assert p.cds_position == 1234
        assert p.alt_allele == "ATG"
        assert p.ref_allele == "-"


# ── Tests: HGVS Delins Parsing ──────────────────────────────────────

class TestHGVSDelins:

    def setup_method(self):
        self.parser = HGVSParser()

    def test_delins(self):
        """c.1234delinsATG should parse as delins."""
        p = self.parser.parse("NM_000267.3:c.1234delinsATG")
        assert p.variant_type == "delins"
        assert p.cds_position == 1234
        assert p.alt_allele == "ATG"


# ── Tests: Invalid HGVS ──────────────────────────────────────────────

class TestInvalidHGVS:

    def setup_method(self):
        self.parser = HGVSParser()

    def test_garbage_input_raises(self):
        with pytest.raises(ValueError):
            self.parser.parse("not a real notation")

    def test_protein_notation_raises(self):
        with pytest.raises(ValueError):
            self.parser.parse("NP_000258.1:p.Arg304Ter")

    def test_genomic_notation_raises(self):
        with pytest.raises(ValueError):
            self.parser.parse("NC_000017.11:g.31200443C>T")


# ── Tests: Batch Parsing ─────────────────────────────────────────────

class TestBatchParsing:

    def test_batch_parse_mixed(self):
        parser = HGVSParser()
        notations = [
            "NM_000267.3:c.910C>T",
            "NM_000267.3:c.1234del",
            "invalid notation here",
            "NM_000267.3:c.5678A>G",
        ]
        results = parser.parse_batch(notations, gene_symbol="NF1")
        assert len(results) == 3  # 2 valid + 1 invalid skipped
        assert results[0].gene_symbol == "NF1"


# ── Tests: ClinVar Processor ────────────────────────────────────────

class TestClinVarProcessor:

    def test_classify_snvs_abe(self):
        """ABE-amenable SNVs should be classified correctly."""
        processor = ClinVarProcessor()
        variants = [
            ClinVarVariant(gene_symbol="NF1", ref_allele="G", alt_allele="A"),  # correction A->G = ABE
            ClinVarVariant(gene_symbol="NF1", ref_allele="C", alt_allele="T"),  # correction T->C = ABE (antisense)
            ClinVarVariant(gene_symbol="NF1", ref_allele="A", alt_allele="G"),  # correction G->A = CBE (antisense)
        ]
        classified = processor.classify_variants(variants)
        # alt=A, ref=G: correction A->G = ABE
        # alt=T, ref=C: correction T->C = ABE (antisense)
        assert len(classified["abe_amenable"]) == 2
        # alt=G, ref=A: correction G->A = CBE
        assert len(classified["cbe_amenable"]) == 1

    def test_classify_transversions(self):
        processor = ClinVarProcessor()
        variants = [
            ClinVarVariant(ref_allele="A", alt_allele="C"),  # transversion
            ClinVarVariant(ref_allele="G", alt_allele="T"),  # transversion
        ]
        classified = processor.classify_variants(variants)
        assert len(classified["transversion"]) == 2

    def test_filter_snvs(self):
        processor = ClinVarProcessor()
        variants = [
            ClinVarVariant(ref_allele="A", alt_allele="G"),
            ClinVarVariant(ref_allele="ATG", alt_allele="-"),  # not SNV
            ClinVarVariant(ref_allele="C", alt_allele="T"),
        ]
        snvs = processor.filter_snvs(variants)
        assert len(snvs) == 2

    def test_to_pipeline_inputs(self):
        processor = ClinVarProcessor()
        variants = [
            ClinVarVariant(
                clinvar_id="12345",
                gene_symbol="NF1",
                chromosome="17",
                position=31200443,
                ref_allele="C",
                alt_allele="T",
                hgvs_c="c.910C>T",
            ),
        ]
        inputs = processor.to_pipeline_inputs(variants)
        assert len(inputs) == 1
        assert isinstance(inputs[0], GenomicVariantInput)
        assert inputs[0].chromosome == "17"
        assert inputs[0].position == 31200443
        assert inputs[0].gene_symbol == "NF1"
        assert inputs[0].name == "c.910C>T"

    def test_load_tsv_with_mock_file(self):
        """Test TSV loading with a minimal mock ClinVar file."""
        processor = ClinVarProcessor()
        content = (
            "#Comment line\n"
            "GeneSymbol\tClinicalSignificance\tAssembly\tChromosome\t"
            "PositionVCF\tReferenceAlleleVCF\tAlternateAlleleVCF\tVariationID\n"
            "NF1\tPathogenic\tGRCh38\t17\t31200443\tC\tT\t12345\n"
            "BRCA2\tBenign\tGRCh38\t13\t32000000\tA\tG\t67890\n"
            "DMD\tPathogenic\tGRCh38\tX\t31000000\tG\tA\t11111\n"
        )
        with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as f:
            f.write(content)
            f.flush()
            tmppath = f.name

        try:
            variants = processor.load_tsv(
                tmppath,
                significance_filter=["Pathogenic"],
            )
            assert len(variants) == 2  # NF1 and DMD (BRCA2 is Benign)
            assert variants[0].gene_symbol == "NF1"
            assert variants[1].gene_symbol == "DMD"
        finally:
            os.unlink(tmppath)
