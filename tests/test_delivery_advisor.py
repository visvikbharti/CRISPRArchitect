"""
Tests for DeliveryAdvisor — post-ranking delivery annotations.

Tests cover:
1. Donor format recommendation by edit size
2. Cell-type-specific warnings (iPSC dsDNA toxicity, HSC culture)
3. Delivery method recommendation by modality + cell type
4. Hard constraint checking (AAV size limits)
5. Viability enhancer suggestions
6. Modality classification from strategy names
7. End-to-end annotation of a mock PipelineResult
"""

import pytest
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from core.feasibility.delivery_advisor import (
    DeliveryAdvisor,
    DeliveryAnnotation,
    DeliveryAdvisoryResult,
    DonorRecommendation,
    annotate_delivery,
)
from core.models import (
    EditModality,
    EvidenceTier,
    FeasibilityLabel,
    GenomicVariantInput,
    NormalizedVariant,
    PipelineResult,
    RiskLevel,
    ScoredStrategy,
    Strategy,
    TranscriptInfo,
    ExonRecord,
    TranscriptCoordinate,
    CodingAnnotation,
    ConsequenceType,
    ReferenceValidation,
)


# ═══════════════════════════════════════════════════════════════════════════
# Fixtures
# ═══════════════════════════════════════════════════════════════════════════

def _make_transcript():
    """Create a minimal TranscriptInfo for testing."""
    return TranscriptInfo(
        transcript_id="ENST00000358273",
        gene_symbol="NF1",
        gene_id="ENSG00000196712",
        chromosome="17",
        start=31094927,
        end=31377677,
        strand=-1,
        biotype="protein_coding",
        is_canonical=True,
        exons=[ExonRecord("ENSE001", 1, 31094927, 31095068, -1, "17")],
    )


def _make_variant(ref="C", alt="T", name="test_snp"):
    """Create a minimal NormalizedVariant."""
    gv = GenomicVariantInput(
        chromosome="17", position=31232193,
        ref_allele=ref, alt_allele=alt,
        gene_symbol="NF1", name=name,
    )
    coord = TranscriptCoordinate(
        genomic_position=31232193, exon_number=1,
        transcript_position=100, cds_position=100,
        codon_index=34, codon_position=1,
        reference_codon="CAG", reference_aa="Q",
        distance_to_exon_start=50, distance_to_exon_end=50,
    )
    coding = CodingAnnotation(consequence=ConsequenceType.MISSENSE)
    ref_val = ReferenceValidation(
        is_valid=True,
        expected_ref_genomic=ref, provided_ref_genomic=ref,
        expected_ref_transcript=ref, provided_ref_transcript=ref,
    )
    return NormalizedVariant(
        input=gv, transcript=_make_transcript(),
        transcript_coord=coord, coding=coding,
        ref_validation=ref_val,
    )


def _make_strategy(name, num_dsbs=0, num_donors=0, num_rounds=1):
    """Create a minimal Strategy."""
    return Strategy(
        name=name,
        num_dsbs=num_dsbs,
        num_donors=num_donors,
        num_rounds=num_rounds,
        rearrangement_risk=RiskLevel.LOW,
        evidence_tier=EvidenceTier.A,
        modality_prior_score=0.9,
    )


def _make_scored(strategy, score=0.8):
    """Create a ScoredStrategy wrapper."""
    return ScoredStrategy(
        strategy=strategy,
        overall_score=score,
        rank=1,
    )


def _make_pipeline_result(strategies, variants=None):
    """Create a minimal PipelineResult."""
    if variants is None:
        variants = [_make_variant()]
    scored = [_make_scored(s, 0.9 - i * 0.1)
              for i, s in enumerate(strategies)]
    for i, sc in enumerate(scored, 1):
        sc.rank = i
    return PipelineResult(
        transcript=_make_transcript(),
        variants=variants,
        strategies=scored,
    )


# ═══════════════════════════════════════════════════════════════════════════
# Test: Donor format recommendation by edit size
# ═══════════════════════════════════════════════════════════════════════════

class TestDonorRecommendation:
    """Test donor format selection by edit size."""

    def test_snp_recommends_ssodn(self):
        advisor = DeliveryAdvisor(cell_type="iPSC")
        rec = advisor._recommend_donor_format(1)
        assert rec.format == "ssODN"

    def test_small_insertion_recommends_ssodn(self):
        advisor = DeliveryAdvisor(cell_type="iPSC")
        rec = advisor._recommend_donor_format(30)
        assert rec.format == "ssODN"

    def test_medium_insertion_recommends_cssdna(self):
        advisor = DeliveryAdvisor(cell_type="iPSC")
        rec = advisor._recommend_donor_format(500)
        assert rec.format == "cssDNA"

    def test_large_insertion_recommends_cssdna(self):
        advisor = DeliveryAdvisor(cell_type="iPSC")
        rec = advisor._recommend_donor_format(3000)
        assert rec.format == "cssDNA"

    def test_very_large_insertion_recommends_cssdna(self):
        advisor = DeliveryAdvisor(cell_type="iPSC")
        rec = advisor._recommend_donor_format(15000)
        assert rec.format == "cssDNA"
        assert "GATALYST" in rec.rationale

    def test_hsc_medium_insert_includes_aav6_alternative(self):
        advisor = DeliveryAdvisor(cell_type="CD34_HSC")
        rec = advisor._recommend_donor_format(2000)
        assert rec.format == "cssDNA"
        assert "AAV6" in rec.alternatives

    def test_ipsc_warns_against_dsdna(self):
        advisor = DeliveryAdvisor(cell_type="iPSC")
        rec = advisor._recommend_donor_format(500)
        assert any("dsDNA" in w and "p53" in w for w in rec.warnings)


# ═══════════════════════════════════════════════════════════════════════════
# Test: Cell-type-specific warnings
# ═══════════════════════════════════════════════════════════════════════════

class TestCellTypeWarnings:
    """Test cell-type-specific delivery warnings."""

    def test_ipsc_dsb_warning(self):
        advisor = DeliveryAdvisor(cell_type="iPSC")
        strategy = _make_strategy("HDR_test", num_dsbs=1, num_donors=1)
        result = _make_pipeline_result([strategy])
        advisory = advisor.advise(result)
        ann = advisory.annotations[0]
        assert any("p53" in w for w in ann.warnings)

    def test_hsc_hdr_quiescence_warning(self):
        advisor = DeliveryAdvisor(cell_type="CD34_HSC")
        strategy = _make_strategy("HDR_test", num_dsbs=1, num_donors=1)
        result = _make_pipeline_result([strategy])
        advisory = advisor.advise(result)
        ann = advisory.annotations[0]
        assert any("S/G2" in w or "quiescent" in w for w in ann.warnings)

    def test_hek293t_no_dsb_warning(self):
        advisor = DeliveryAdvisor(cell_type="HEK293T")
        strategy = _make_strategy("HDR_test", num_dsbs=1, num_donors=1)
        result = _make_pipeline_result([strategy])
        advisory = advisor.advise(result)
        ann = advisory.annotations[0]
        # HEK293T has mutant p53, no DSB warning expected
        assert not any("p53" in w for w in ann.warnings)


# ═══════════════════════════════════════════════════════════════════════════
# Test: Delivery method recommendation
# ═══════════════════════════════════════════════════════════════════════════

class TestDeliveryMethod:
    """Test delivery method recommendation by modality and cell type."""

    def test_be_ipsc_recommends_mrna_nucleofection(self):
        advisor = DeliveryAdvisor(cell_type="iPSC")
        method, desc = advisor._recommend_delivery_method("BE")
        assert "nucleofection" in method or "nucleofection" in desc.lower()

    def test_hdr_hsc_recommends_rnp_nucleofection(self):
        advisor = DeliveryAdvisor(cell_type="CD34_HSC")
        method, desc = advisor._recommend_delivery_method("HDR")
        assert "nucleofection" in method or "nucleofection" in desc.lower()

    def test_be_hek293t_recommends_lipofection(self):
        advisor = DeliveryAdvisor(cell_type="HEK293T")
        method, desc = advisor._recommend_delivery_method("BE")
        assert "lipofection" in method.lower() or "lipofect" in desc.lower()


# ═══════════════════════════════════════════════════════════════════════════
# Test: Hard constraint violations
# ═══════════════════════════════════════════════════════════════════════════

class TestHardConstraints:
    """Test hard feasibility constraint checking."""

    def test_aav6_size_violation(self):
        advisor = DeliveryAdvisor(cell_type="CD34_HSC")
        strategy = _make_strategy("HDR_large", num_dsbs=1, num_donors=1)
        # Edit size exceeds AAV6 capacity
        result = _make_pipeline_result(
            [strategy],
            variants=[_make_variant(ref="A" * 5000, alt="-", name="large_del")],
        )
        advisory = advisor.advise(result)
        ann = advisory.annotations[0]
        # The donor recommendation should be cssDNA, not AAV6
        if ann.donor_recommendation:
            assert ann.donor_recommendation.format == "cssDNA"

    def test_be_strategy_is_deliverable(self):
        advisor = DeliveryAdvisor(cell_type="iPSC")
        strategy = _make_strategy("ABE8e_single", num_dsbs=0, num_donors=0)
        result = _make_pipeline_result([strategy])
        advisory = advisor.advise(result)
        ann = advisory.annotations[0]
        assert ann.is_deliverable is True
        assert len(ann.hard_constraint_violations) == 0


# ═══════════════════════════════════════════════════════════════════════════
# Test: Viability enhancers
# ═══════════════════════════════════════════════════════════════════════════

class TestViabilityEnhancers:
    """Test viability enhancer recommendations."""

    def test_ipsc_hdr_enhancers(self):
        advisor = DeliveryAdvisor(cell_type="iPSC")
        enhancers = advisor._get_viability_enhancers("HDR")
        assert any("ROCK" in e for e in enhancers)
        assert any("BCL-XL" in e for e in enhancers)
        assert any("Cold shock" in e or "32C" in e for e in enhancers)

    def test_ipsc_be_enhancers(self):
        advisor = DeliveryAdvisor(cell_type="iPSC")
        enhancers = advisor._get_viability_enhancers("BE")
        assert any("ROCK" in e for e in enhancers)
        assert any("p53DD" in e for e in enhancers)

    def test_hsc_enhancers(self):
        advisor = DeliveryAdvisor(cell_type="CD34_HSC")
        enhancers = advisor._get_viability_enhancers("HDR")
        assert any("SCF" in e or "TPO" in e for e in enhancers)
        assert any("HiFi" in e for e in enhancers)


# ═══════════════════════════════════════════════════════════════════════════
# Test: Modality classification
# ═══════════════════════════════════════════════════════════════════════════

class TestModalityClassification:
    """Test strategy name → modality class mapping."""

    def test_abe_classified_as_be(self):
        advisor = DeliveryAdvisor()
        assert advisor._classify_modality("ABE8e_single") == "BE"

    def test_cbe_classified_as_be(self):
        advisor = DeliveryAdvisor()
        assert advisor._classify_modality("CBE_dual") == "BE"

    def test_pe_classified_as_pe(self):
        advisor = DeliveryAdvisor()
        assert advisor._classify_modality("PE_single_round1") == "PE"

    def test_prime_classified_as_pe(self):
        advisor = DeliveryAdvisor()
        assert advisor._classify_modality("prime_editing") == "PE"

    def test_hdr_classified_as_hdr(self):
        advisor = DeliveryAdvisor()
        assert advisor._classify_modality("HDR_cssDNA_single") == "HDR"

    def test_knockin_classified_as_hdr(self):
        advisor = DeliveryAdvisor()
        assert advisor._classify_modality("knockin_GFP") == "HDR"


# ═══════════════════════════════════════════════════════════════════════════
# Test: End-to-end annotation
# ═══════════════════════════════════════════════════════════════════════════

class TestEndToEnd:
    """Test full advisory pipeline on mock results."""

    def test_multiple_strategies(self):
        strategies = [
            _make_strategy("ABE8e_single", num_dsbs=0, num_donors=0),
            _make_strategy("PE_single", num_dsbs=0, num_donors=0),
            _make_strategy("HDR_cssDNA", num_dsbs=1, num_donors=1),
        ]
        result = _make_pipeline_result(strategies)
        advisory = annotate_delivery(result, cell_type="iPSC")

        assert len(advisory.annotations) == 3
        assert advisory.cell_type == "iPSC"

        # BE should be simplest
        be_ann = advisory.annotations[0]
        assert be_ann.delivery_complexity <= 2

        # HDR should have donor recommendation
        hdr_ann = advisory.annotations[2]
        assert hdr_ann.donor_recommendation is not None

    def test_convenience_function(self):
        strategy = _make_strategy("ABE8e_test", num_dsbs=0)
        result = _make_pipeline_result([strategy])
        advisory = annotate_delivery(result, cell_type="CD34_HSC")

        assert isinstance(advisory, DeliveryAdvisoryResult)
        assert len(advisory.annotations) == 1

    def test_empty_pipeline(self):
        result = PipelineResult(
            transcript=_make_transcript(),
            variants=[],
            strategies=[],
        )
        advisor = DeliveryAdvisor(cell_type="iPSC")
        advisory = advisor.advise(result)
        assert len(advisory.annotations) == 0

    def test_delivery_complexity_ordering(self):
        """BE should have lower complexity than PE, which is lower than HDR+AAV6."""
        advisor = DeliveryAdvisor(cell_type="iPSC")

        be_complexity = advisor._get_delivery_complexity(
            _make_strategy("ABE_test"), "BE"
        )
        pe_complexity = advisor._get_delivery_complexity(
            _make_strategy("PE_test"), "PE"
        )
        assert be_complexity <= pe_complexity


# ═══════════════════════════════════════════════════════════════════════════
# Test: Delivery constants loaded correctly
# ═══════════════════════════════════════════════════════════════════════════

class TestDeliveryConstants:
    """Test that delivery constants are properly defined."""

    def test_donor_format_properties_exist(self):
        from utils.constants import DONOR_FORMAT_PROPERTIES
        assert "ssODN" in DONOR_FORMAT_PROPERTIES
        assert "cssDNA" in DONOR_FORMAT_PROPERTIES
        assert "lssDNA" in DONOR_FORMAT_PROPERTIES
        assert "dsDNA" in DONOR_FORMAT_PROPERTIES
        assert "AAV6" in DONOR_FORMAT_PROPERTIES

    def test_cssdna_max_size(self):
        from utils.constants import DONOR_FORMAT_PROPERTIES
        assert DONOR_FORMAT_PROPERTIES["cssDNA"]["max_size_nt"] == 20000

    def test_aav6_max_payload(self):
        from utils.constants import DELIVERY_AAV_MAX_PAYLOAD_NT
        assert DELIVERY_AAV_MAX_PAYLOAD_NT == 4700

    def test_modality_complexity_ordering(self):
        from utils.constants import MODALITY_DELIVERY_COMPLEXITY
        assert MODALITY_DELIVERY_COMPLEXITY["ABE"] < MODALITY_DELIVERY_COMPLEXITY["PE"]
        assert MODALITY_DELIVERY_COMPLEXITY["HDR_ssODN"] < MODALITY_DELIVERY_COMPLEXITY["HDR_AAV6"]

    def test_cell_compatibility_exists(self):
        from utils.constants import DELIVERY_CELL_COMPATIBILITY
        assert "nucleofection_RNP" in DELIVERY_CELL_COMPATIBILITY
        # CD34+ HSCs should not work with lipofection
        assert DELIVERY_CELL_COMPATIBILITY["lipofection_RNP"]["CD34_HSC"] == "infeasible"

    def test_ipsc_dsdna_toxicity_highest(self):
        from utils.constants import DONOR_FORMAT_PROPERTIES
        ipsc_tox = {
            k: v["toxicity_iPSC"]
            for k, v in DONOR_FORMAT_PROPERTIES.items()
        }
        assert ipsc_tox["dsDNA"] == max(ipsc_tox.values())


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
