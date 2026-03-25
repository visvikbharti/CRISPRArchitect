"""
Tests for CRISPRArchitect v3 Off-Target Scoring Engine
========================================================

Tests CFD scoring and guide specificity assessment.

References
----------
Doench et al., Nat Biotechnol, 2016 (CFD score)
Hsu et al., Nat Biotechnol, 2013 (MIT specificity)
"""

from __future__ import annotations

import pytest

from core.models import GuideCandidate
from core.feasibility.off_target import (
    cfd_score_single_mismatch,
    cfd_score_guide_vs_target,
    count_mismatches,
    OffTargetScorer,
    GuideSpecificity,
)


# ── Tests: CFD Single Mismatch ──────────────────────────────────────

class TestCFDSingleMismatch:

    def test_match_returns_1(self):
        """Matching bases should return 1.0 (no penalty)."""
        assert cfd_score_single_mismatch("A", "A", 10) == 1.0
        assert cfd_score_single_mismatch("G", "G", 5) == 1.0

    def test_mismatch_returns_less_than_1(self):
        """Mismatched bases should return < 1.0."""
        score = cfd_score_single_mismatch("A", "C", 10)
        assert 0.0 < score < 1.0

    def test_pam_distal_more_tolerant(self):
        """Mismatches at PAM-distal positions (1-5) should be more tolerated."""
        score_distal = cfd_score_single_mismatch("A", "C", 1)
        score_proximal = cfd_score_single_mismatch("A", "C", 20)
        assert score_distal > score_proximal

    def test_seed_region_less_tolerant(self):
        """Mismatches in seed region (positions 13-20) should be less tolerated."""
        score_nonseed = cfd_score_single_mismatch("A", "C", 5)
        score_seed = cfd_score_single_mismatch("A", "C", 18)
        assert score_nonseed > score_seed


# ── Tests: CFD Guide vs Target ──────────────────────────────────────

class TestCFDGuideVsTarget:

    def test_perfect_match_returns_1(self):
        """Identical guide and target should score 1.0."""
        guide = "ATCGATCGATCGATCGATCG"
        score = cfd_score_guide_vs_target(guide, guide)
        assert score == 1.0

    def test_one_mismatch_less_than_1(self):
        """Single mismatch should reduce score below 1.0."""
        guide = "ATCGATCGATCGATCGATCG"
        target = "ATCGATCGATCGATCGATCC"  # last base different
        score = cfd_score_guide_vs_target(guide, target)
        assert 0.0 < score < 1.0

    def test_many_mismatches_near_zero(self):
        """Many mismatches should produce very low score."""
        guide = "AAAAAAAAAAAAAAAAAAAA"
        target = "CCCCCCCCCCCCCCCCCCCC"
        score = cfd_score_guide_vs_target(guide, target)
        assert score < 0.01

    def test_wrong_length_returns_0(self):
        """Non-20-mer inputs should return 0."""
        assert cfd_score_guide_vs_target("ATCG", "ATCG") == 0.0


# ── Tests: Count Mismatches ──────────────────────────────────────────

class TestCountMismatches:

    def test_identical_sequences(self):
        assert count_mismatches("ATCG", "ATCG") == 0

    def test_all_different(self):
        assert count_mismatches("AAAA", "CCCC") == 4

    def test_one_mismatch(self):
        assert count_mismatches("ATCG", "ATCC") == 1

    def test_case_insensitive(self):
        assert count_mismatches("atcg", "ATCG") == 0


# ── Tests: Off-Target Scorer ────────────────────────────────────────

class TestOffTargetScorer:

    def test_scorer_initializes(self):
        scorer = OffTargetScorer(nuclease="SpCas9")
        assert scorer.evidence_tier == "A"

    def test_non_spcas9_is_tier_b(self):
        scorer = OffTargetScorer(nuclease="enFnCas9")
        assert scorer.evidence_tier == "B"

    def test_score_guide_returns_specificity(self):
        """Scoring a guide should return a GuideSpecificity object."""
        scorer = OffTargetScorer()
        guide = GuideCandidate(
            sequence_20mer="ATCGATCGATCGATCGATCG",
            pam_sequence="AGG",
            strand="+",
            cut_position=200,
            score=0.8,
        )
        # Build a sequence with the guide embedded
        seq = ("A" * 180) + "ATCGATCGATCGATCGATCG" + "AGG" + ("T" * 197)
        result = scorer.score_guide(guide, seq)
        assert isinstance(result, GuideSpecificity)
        assert 0 <= result.specificity_score <= 100

    def test_unique_guide_has_high_specificity(self):
        """A guide with no off-target matches should have high specificity."""
        scorer = OffTargetScorer()
        guide = GuideCandidate(
            sequence_20mer="GCTAGCTAGCTAGCTAGCTA",
            pam_sequence="GGG",
            strand="+",
            cut_position=200,
            score=0.8,
        )
        # Random sequence unlikely to have similar 20-mers
        import random
        random.seed(42)
        bases = "ATCG"
        seq = "".join(random.choice(bases) for _ in range(400))
        result = scorer.score_guide(guide, seq)
        # Should be high specificity (few or no off-targets)
        assert result.specificity_score > 50

    def test_score_guides_returns_sorted_list(self):
        """score_guides should return results sorted by specificity."""
        scorer = OffTargetScorer()
        guides = [
            GuideCandidate(sequence_20mer="ATCGATCGATCGATCGATCG", score=0.8),
            GuideCandidate(sequence_20mer="GCTAGCTAGCTAGCTAGCTA", score=0.7),
        ]
        seq = "A" * 400
        results = scorer.score_guides(guides, seq)
        assert len(results) == 2
        assert results[0].specificity_score >= results[1].specificity_score

    def test_off_target_counts(self):
        """Off-target counts should be reasonable."""
        scorer = OffTargetScorer(max_mismatches=3)
        guide = GuideCandidate(
            sequence_20mer="ATCGATCGATCGATCGATCG",
            pam_sequence="AGG",
            strand="+",
            cut_position=200,
            score=0.8,
        )
        # Sequence with the exact guide + PAM (should find at least 1 match)
        seq = ("G" * 180) + "ATCGATCGATCGATCGATCG" + "AGG" + ("C" * 197)
        result = scorer.score_guide(guide, seq)
        assert result.n_off_targets_0mm >= 1  # the on-target itself

    def test_empty_guide_returns_zero_specificity(self):
        """Short guide should return 0 specificity."""
        scorer = OffTargetScorer()
        guide = GuideCandidate(sequence_20mer="ATCG")  # too short
        result = scorer.score_guide(guide, "A" * 100)
        assert result.specificity_score == 0.0
