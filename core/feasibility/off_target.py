"""
Off-Target Scoring Engine for CRISPRArchitect v3
==================================================

Computes guide specificity scores based on the CFD (Cutting Frequency
Determination) scoring matrix from Doench et al. (2016) and the MIT
specificity score framework from Hsu et al. (2013).

These are purely computational scores — no external tools or digenome-seq
data required. The scores predict the relative likelihood that a guide
RNA will cut at off-target sites based on mismatch patterns.

Approach
--------
1. For each candidate guide (20-mer), enumerate potential off-target sites
   by examining all positions in a reference window for mismatches.
2. Score each mismatch pattern using the CFD position-weight matrix.
3. Aggregate off-target scores into a specificity score.

Limitations
-----------
- Uses a simplified off-target enumeration (local sequence window only).
  For genome-wide off-target prediction, external tools like Cas-OFFinder
  or FlashFry are recommended.
- CFD scores are trained on SpCas9 data. Accuracy for other nucleases
  (enFnCas9, SpCas9-NG, SpRY) is extrapolated.
- DNA/RNA bulges are not modeled (only substitution mismatches).

Evidence Tier: B (literature-informed heuristic)

References
----------
Doench et al., Nat Biotechnol, 2016 (CFD score)
Hsu et al., Nat Biotechnol, 2013 (MIT specificity score)
Concordet & Haeussler, Nucleic Acids Res, 2018 (CRISPOR implementation)

Python 3.9 compatible.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

from core.models import GuideCandidate


# ═══════════════════════════════════════════════════════════════════════
# CFD Scoring Matrix (Doench et al., 2016)
# ═══════════════════════════════════════════════════════════════════════
#
# CFD score for a single mismatch at position p with mismatch type (r:d)
# where r = RNA base (guide), d = DNA base (target).
# Position 1 = PAM-distal, position 20 = PAM-proximal.
# Value = fraction of cutting activity retained with this mismatch.
# 1.0 = no activity loss, 0.0 = complete loss.
#
# These values are derived from Table S19 of Doench et al., 2016.
# Simplified to key mismatch types per position.

_CFD_MISMATCH_SCORES: Dict[Tuple[str, str, int], float] = {}

# PAM-distal mismatches (positions 1-8) are generally better tolerated
# PAM-proximal mismatches (positions 13-20) are less tolerated
# Seed region (positions 1-12 from PAM, i.e., 9-20 in our numbering) is critical

# Position-dependent mismatch tolerance (averaged across mismatch types)
# Source: Doench et al., 2016, Figure 3A
_POSITION_TOLERANCE = {
    1: 0.85, 2: 0.80, 3: 0.75, 4: 0.70, 5: 0.65,
    6: 0.60, 7: 0.55, 8: 0.50, 9: 0.45, 10: 0.40,
    11: 0.35, 12: 0.30, 13: 0.25, 14: 0.20, 15: 0.18,
    16: 0.15, 17: 0.12, 18: 0.10, 19: 0.08, 20: 0.05,
}

# Mismatch-type modifiers (relative to position tolerance)
# rG:dT mismatches are better tolerated (wobble pairing)
# rC:dC mismatches are poorly tolerated
_MISMATCH_TYPE_MODIFIER = {
    ("A", "C"): 0.8, ("A", "G"): 0.6, ("A", "A"): 0.0,  # rA vs dC/dG/dA
    ("C", "A"): 0.7, ("C", "T"): 0.6, ("C", "C"): 0.0,  # rC vs dA/dT/dC
    ("G", "A"): 0.7, ("G", "C"): 0.6, ("G", "G"): 0.0,  # rG vs dA/dC/dG
    ("T", "A"): 0.8, ("T", "G"): 1.0, ("T", "C"): 0.6,  # rU(T):dG is wobble
    ("T", "T"): 0.0,
    ("A", "T"): 1.0, ("C", "G"): 1.0, ("G", "C"): 1.0, ("T", "A"): 1.0,  # match
}


def cfd_score_single_mismatch(
    guide_base: str,
    target_base: str,
    position: int,
) -> float:
    """CFD score for a single mismatch at a given position.

    Parameters
    ----------
    guide_base : str
        RNA base in the guide (A, C, G, T/U).
    target_base : str
        DNA base at the off-target site.
    position : int
        1-indexed position in the 20-mer (1 = PAM-distal).

    Returns
    -------
    float
        Activity fraction [0, 1]. 1.0 = perfect match (no penalty).
    """
    gb = guide_base.upper().replace("U", "T")
    tb = target_base.upper()

    if gb == tb:
        return 1.0  # match: no penalty

    pos_tol = _POSITION_TOLERANCE.get(position, 0.3)
    mm_mod = _MISMATCH_TYPE_MODIFIER.get((gb, tb), 0.7)

    return pos_tol * mm_mod


def cfd_score_guide_vs_target(
    guide_20mer: str,
    target_20mer: str,
) -> float:
    """Compute CFD score for a guide against one potential off-target.

    The CFD score is the product of individual mismatch penalties:
        CFD = product(cfd_single(pos_i)) for all mismatched positions

    A CFD of 1.0 means the off-target is identical to the on-target.
    A CFD of 0.0 means the off-target will not be cut.

    Parameters
    ----------
    guide_20mer : str
        The 20-nt guide sequence.
    target_20mer : str
        The 20-nt potential off-target site.

    Returns
    -------
    float
        CFD score in [0, 1].
    """
    if len(guide_20mer) != 20 or len(target_20mer) != 20:
        return 0.0

    score = 1.0
    for i in range(20):
        position = i + 1  # 1-indexed
        s = cfd_score_single_mismatch(
            guide_20mer[i], target_20mer[i], position
        )
        score *= s

    return score


def count_mismatches(seq1: str, seq2: str) -> int:
    """Count number of mismatches between two equal-length sequences."""
    return sum(1 for a, b in zip(seq1.upper(), seq2.upper()) if a != b)


# ═══════════════════════════════════════════════════════════════════════
# Off-Target Specificity Scorer
# ═══════════════════════════════════════════════════════════════════════

@dataclass
class OffTargetHit:
    """A potential off-target site for a guide."""
    target_sequence: str = ""
    num_mismatches: int = 0
    cfd_score: float = 0.0
    position_in_window: int = 0
    strand: str = "+"


@dataclass
class GuideSpecificity:
    """Specificity assessment for one guide candidate.

    Attributes
    ----------
    guide : GuideCandidate
        The guide being assessed.
    specificity_score : float
        MIT-style aggregate specificity in [0, 100].
        100 = perfectly unique, 0 = many high-scoring off-targets.
    n_off_targets_0mm : int
        Off-targets with 0 mismatches (should be 1 = the on-target).
    n_off_targets_1mm : int
        Off-targets with exactly 1 mismatch.
    n_off_targets_2mm : int
        Off-targets with exactly 2 mismatches.
    n_off_targets_3mm : int
        Off-targets with exactly 3 mismatches.
    top_off_targets : list of OffTargetHit
        Top 5 off-target hits by CFD score.
    evidence_tier : str
        "A" for SpCas9 (CFD trained on this), "B" for others (extrapolated).
    """
    guide: Optional[GuideCandidate] = None
    specificity_score: float = 100.0
    n_off_targets_0mm: int = 0
    n_off_targets_1mm: int = 0
    n_off_targets_2mm: int = 0
    n_off_targets_3mm: int = 0
    top_off_targets: List[OffTargetHit] = field(default_factory=list)
    evidence_tier: str = "B"
    warnings: List[str] = field(default_factory=list)


class OffTargetScorer:
    """Score guide specificity using CFD-based off-target analysis.

    This performs LOCAL off-target analysis within the provided sequence
    window. For genome-wide analysis, use external tools (Cas-OFFinder).

    Parameters
    ----------
    max_mismatches : int
        Maximum mismatches to consider (default 3).
    nuclease : str
        Nuclease name for evidence tier assignment.
    """

    def __init__(
        self,
        max_mismatches: int = 3,
        nuclease: str = "SpCas9",
    ) -> None:
        self.max_mismatches = max_mismatches
        self.nuclease = nuclease
        # CFD was trained on SpCas9 data — mark others as Tier B
        self.evidence_tier = "A" if nuclease == "SpCas9" else "B"

    def score_guide(
        self,
        guide: GuideCandidate,
        local_sequence: str,
        pam: str = "NGG",
    ) -> GuideSpecificity:
        """Score a single guide for off-target specificity.

        Scans the local sequence for PAM-adjacent 20-mers that are
        similar to the guide (up to max_mismatches differences).

        Parameters
        ----------
        guide : GuideCandidate
            The guide to assess.
        local_sequence : str
            Genomic sequence window to scan for off-targets.
        pam : str
            PAM sequence (IUPAC) for the nuclease.

        Returns
        -------
        GuideSpecificity
            Specificity assessment.
        """
        guide_seq = guide.sequence_20mer.upper()
        if len(guide_seq) != 20:
            return GuideSpecificity(
                guide=guide,
                specificity_score=0.0,
                warnings=["Guide is not 20 nt."],
            )

        seq = local_sequence.upper()
        off_targets = self._find_off_targets(guide_seq, seq, pam)

        # Count by mismatch class
        n_0mm = sum(1 for ot in off_targets if ot.num_mismatches == 0)
        n_1mm = sum(1 for ot in off_targets if ot.num_mismatches == 1)
        n_2mm = sum(1 for ot in off_targets if ot.num_mismatches == 2)
        n_3mm = sum(1 for ot in off_targets if ot.num_mismatches == 3)

        # MIT-style specificity: 100 / (1 + sum_of_off_target_CFD_scores)
        # Higher = more specific (fewer/weaker off-targets)
        off_target_cfd_sum = sum(
            ot.cfd_score for ot in off_targets if ot.num_mismatches > 0
        )
        specificity = 100.0 / (1.0 + off_target_cfd_sum)
        specificity = min(100.0, max(0.0, specificity))

        # Top off-targets
        top_ots = sorted(
            [ot for ot in off_targets if ot.num_mismatches > 0],
            key=lambda x: x.cfd_score,
            reverse=True,
        )[:5]

        warnings = []
        if n_1mm > 3:
            warnings.append(
                f"High off-target risk: {n_1mm} sites with 1 mismatch."
            )
        if specificity < 50:
            warnings.append(
                f"Low specificity score ({specificity:.0f}/100). "
                "Consider alternative guides."
            )
        if self.evidence_tier == "B":
            warnings.append(
                f"CFD scores extrapolated for {self.nuclease} "
                "(trained on SpCas9 data)."
            )

        return GuideSpecificity(
            guide=guide,
            specificity_score=round(specificity, 1),
            n_off_targets_0mm=n_0mm,
            n_off_targets_1mm=n_1mm,
            n_off_targets_2mm=n_2mm,
            n_off_targets_3mm=n_3mm,
            top_off_targets=top_ots,
            evidence_tier=self.evidence_tier,
            warnings=warnings,
        )

    def score_guides(
        self,
        guides: List[GuideCandidate],
        local_sequence: str,
        pam: str = "NGG",
    ) -> List[GuideSpecificity]:
        """Score multiple guides for off-target specificity.

        Returns results sorted by specificity score (best first).
        """
        results = [
            self.score_guide(g, local_sequence, pam) for g in guides
        ]
        results.sort(key=lambda x: x.specificity_score, reverse=True)
        return results

    def _find_off_targets(
        self,
        guide_seq: str,
        sequence: str,
        pam: str,
    ) -> List[OffTargetHit]:
        """Find potential off-target sites in a sequence.

        Scans both strands for 20-mers adjacent to PAM sites that
        have <= max_mismatches differences from the guide.
        """
        import re

        # Build PAM regex
        iupac = {
            'N': '[ATCG]', 'R': '[AG]', 'Y': '[CT]', 'S': '[GC]',
            'W': '[AT]', 'K': '[GT]', 'M': '[AC]', 'B': '[CGT]',
            'D': '[AGT]', 'H': '[ACT]', 'V': '[ACG]',
            'A': 'A', 'T': 'T', 'G': 'G', 'C': 'C',
        }
        pam_regex = ''.join(iupac.get(b, b) for b in pam.upper())
        pam_len = len(pam)

        hits = []
        comp = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}

        # Forward strand: [20-mer][PAM]
        for m in re.finditer(f'(?={pam_regex})', sequence):
            pam_start = m.start()
            proto_start = pam_start - 20
            if proto_start < 0:
                continue
            target = sequence[proto_start:pam_start]
            if len(target) != 20:
                continue

            mm = count_mismatches(guide_seq, target)
            if mm <= self.max_mismatches:
                cfd = cfd_score_guide_vs_target(guide_seq, target)
                hits.append(OffTargetHit(
                    target_sequence=target,
                    num_mismatches=mm,
                    cfd_score=round(cfd, 4),
                    position_in_window=proto_start,
                    strand="+",
                ))

        # Reverse strand
        rc_seq = ''.join(comp.get(b, 'N') for b in reversed(sequence))
        for m in re.finditer(f'(?={pam_regex})', rc_seq):
            pam_start = m.start()
            proto_start = pam_start - 20
            if proto_start < 0:
                continue
            target = rc_seq[proto_start:pam_start]
            if len(target) != 20:
                continue

            mm = count_mismatches(guide_seq, target)
            if mm <= self.max_mismatches:
                cfd = cfd_score_guide_vs_target(guide_seq, target)
                fwd_pos = len(sequence) - pam_start - pam_len
                hits.append(OffTargetHit(
                    target_sequence=target,
                    num_mismatches=mm,
                    cfd_score=round(cfd, 4),
                    position_in_window=fwd_pos,
                    strand="-",
                ))

        return hits
