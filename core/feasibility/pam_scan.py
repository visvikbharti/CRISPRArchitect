"""
Enhanced PAM Scanner for CRISPRArchitect v2
=============================================

Scans a local genomic sequence for PAM sites, extracts 20-mer protospacers,
computes cut positions, and returns ranked GuideCandidate objects.

Supports SpCas9 (NGG), enFnCas9 (NRG), Cas12a (TTTV), and vCas9 (NGG)
as first-class nucleases via NUCLEASE_PARAMS lookup.

Biological background
---------------------
- SpCas9 recognises a 5'-NGG-3' PAM 3' of the protospacer and cuts 3 bp
  upstream of the PAM on both strands (Jinek et al., Science, 2012).
- enFnCas9 recognises a broadened 5'-NRG-3' PAM (R = A or G) and generates
  staggered 5' overhangs of ~2-5 bp (Hirano et al., Cell, 2016;
  Chakraborty lab, Nat Commun, 2024).
- Guides with GC content <30% or >70% are rejected (Doench et al.,
  Nat Biotechnol, 2014; Hart et al., 2015).
- Protospacers with >=4 consecutive T's are flagged because the TTTT motif
  acts as a Pol III terminator and reduces sgRNA expression from U6
  promoters (Graf et al., 2019).

Python 3.9 compatible -- no slots=True, no match statements, no X | Y unions.
"""

from __future__ import annotations

import math
import sys
from typing import Dict, List, Optional

from core.models import GuideCandidate

try:
    from utils.sequence import find_pam_sites, gc_content, reverse_complement
    from utils.constants import NUCLEASE_PARAMS
except ImportError:
    from crisprarchitect.utils.sequence import find_pam_sites, gc_content, reverse_complement
    from crisprarchitect.utils.constants import NUCLEASE_PARAMS


# ─── Scoring constants ────────────────────────────────────────────────
# Maximum distance (bp) for which we still consider a guide useful
_MAX_DISTANCE_BP = 200
# Optimal GC for sgRNA activity (Doench et al., 2014)
_OPTIMAL_GC = 0.50
# Weight factors for composite score
_DISTANCE_WEIGHT = 0.60
_GC_WEIGHT = 0.25
_POLY_T_PENALTY = 0.15   # deducted if >=4 T run present

# Length of the protospacer (constant for Cas9 family)
_PROTOSPACER_LEN = 20


class EnhancedPAMScanner:
    """Identify and rank sgRNA guides around an edit site.

    Parameters
    ----------
    nuclease : str
        Key into ``NUCLEASE_PARAMS`` (default ``"SpCas9"``).

    Raises
    ------
    ValueError
        If the nuclease name is not found in ``NUCLEASE_PARAMS``.

    References
    ----------
    Jinek et al., Science, 2012 (SpCas9 PAM recognition and cut mechanism)
    Hirano et al., Cell, 2016 (FnCas9 structure and PAM specificity)
    Doench et al., Nat Biotechnol, 2014 (sgRNA design rules and GC optimality)
    """

    def __init__(self, nuclease: str = "SpCas9") -> None:
        if nuclease not in NUCLEASE_PARAMS:
            raise ValueError(
                f"Unknown nuclease '{nuclease}'. "
                f"Available: {list(NUCLEASE_PARAMS.keys())}"
            )
        self.nuclease = nuclease
        params = NUCLEASE_PARAMS[nuclease]
        self.pam: str = params["pam"]
        self.cut_type: str = params["cut_type"]
        self.stagger_bp: int = params["stagger_bp"]

    # ── public API ─────────────────────────────────────────────────────

    def scan(
        self,
        sequence: str,
        edit_position: int,
        window_bp: int = 200,
    ) -> List[GuideCandidate]:
        """Scan *sequence* for PAM sites and return ranked guides.

        Parameters
        ----------
        sequence : str
            Local genomic DNA sequence (5' to 3', sense strand).
        edit_position : int
            0-based index of the target edit within *sequence*.
        window_bp : int
            Only consider PAM hits whose predicted cut falls within
            ``edit_position +/- window_bp``.  Default 200 bp.

        Returns
        -------
        List[GuideCandidate]
            Guides sorted by composite score (best first).
        """
        pam_hits = find_pam_sites(sequence, pam=self.pam, strand="both")

        candidates: List[GuideCandidate] = []
        for hit in pam_hits:
            guide = self._build_candidate(sequence, hit, edit_position)
            if guide is None:
                continue
            # Distance filter
            if abs(guide.distance_to_edit) > window_bp:
                continue
            # GC filter (Doench et al., 2014 -- 30-70% range)
            if guide.gc_content < 0.30 or guide.gc_content > 0.70:
                continue
            # Poly-T filter (reject >= 4 consecutive T's as hard filter)
            if guide.has_poly_t:
                continue
            # Score
            guide.score = self._score(guide, window_bp)
            candidates.append(guide)

        candidates.sort(key=lambda g: g.score, reverse=True)
        return candidates

    # ── internals ──────────────────────────────────────────────────────

    def _build_candidate(
        self,
        sequence: str,
        hit: Dict,
        edit_position: int,
    ) -> Optional[GuideCandidate]:
        """Extract protospacer and compute cut position for one PAM hit.

        For 3'-PAM nucleases (Cas9 family):
          + strand:  [20-mer protospacer][PAM]  ->  cut is 3 bp upstream of PAM
          - strand:  [PAM][20-mer protospacer]  ->  protospacer is downstream
                     (reverse complement coordinates)

        Returns None if the protospacer cannot be fully extracted from the
        sequence (e.g., near edges).
        """
        pam_pos = hit["position"]     # 0-based position of PAM on fwd strand
        strand = hit["strand"]
        pam_seq = hit["pam_sequence"]
        pam_len = len(self.pam)
        seq_upper = sequence.upper()

        # --- extract 20-mer protospacer ---
        if strand == "+":
            # PAM is at pam_pos..pam_pos+pam_len on + strand
            # protospacer is 20 bp upstream of PAM on + strand
            proto_start = pam_pos - _PROTOSPACER_LEN
            proto_end = pam_pos
            if proto_start < 0:
                return None
            protospacer = seq_upper[proto_start:proto_end]
            # Cut position: 3 bp upstream of PAM on the protospacer strand
            # i.e., between positions (pam_pos - 3) and (pam_pos - 4)
            # We report the 0-based coordinate of the cut (upstream side)
            cut_pos = pam_pos - 3
        else:
            # '-' strand hit: find_pam_sites gives fwd_pos = len-match_start-pam_len
            # On the - strand the arrangement (in 5'->3' of - strand) is:
            #   [20-mer protospacer][PAM]
            # In forward-strand coordinates the PAM occupies pam_pos..pam_pos+pam_len
            # The protospacer (on the - strand) maps to fwd positions
            #   pam_pos+pam_len .. pam_pos+pam_len+20  (downstream on fwd)
            proto_start = pam_pos + pam_len
            proto_end = proto_start + _PROTOSPACER_LEN
            if proto_end > len(seq_upper):
                return None
            # The protospacer read 5'->3' on the - strand is the reverse complement
            protospacer = reverse_complement(seq_upper[proto_start:proto_end])
            # Cut on - strand: 3 bp upstream of PAM *on the protospacer strand*
            # In fwd coordinates that is pam_pos + pam_len + 3 - 1 = pam_pos + pam_len + 2
            cut_pos = pam_pos + pam_len + 2

        if len(protospacer) != _PROTOSPACER_LEN:
            return None

        gc = gc_content(protospacer)
        has_poly_t = "TTTT" in protospacer
        distance = cut_pos - edit_position  # signed

        return GuideCandidate(
            sequence_20mer=protospacer,
            pam_sequence=pam_seq,
            strand=strand,
            cut_position=cut_pos,
            distance_to_edit=distance,
            gc_content=round(gc, 4),
            position_in_window=-1,   # set later by BE engine
            has_poly_t=has_poly_t,
            score=0.0,
        )

    @staticmethod
    def _score(guide: GuideCandidate, window_bp: int) -> float:
        """Composite score: distance penalty + GC optimality bonus.

        Score is in [0, 1].  Higher is better.
        - Distance component: 1.0 when cut is at edit, decays linearly
          to 0 at *window_bp*.
        - GC component: 1.0 at 50% GC, falls off quadratically.
        """
        # Distance component (linear decay)
        abs_dist = abs(guide.distance_to_edit)
        dist_score = max(0.0, 1.0 - abs_dist / max(window_bp, 1))

        # GC component (quadratic penalty away from 0.5)
        gc_deviation = abs(guide.gc_content - _OPTIMAL_GC)
        gc_score = max(0.0, 1.0 - (gc_deviation / 0.20) ** 2)

        score = _DISTANCE_WEIGHT * dist_score + _GC_WEIGHT * gc_score
        # No poly-T guides reach here (hard-filtered), but keep the
        # penalty weight in the budget for future soft-filter mode.
        return round(score, 4)


# ═══════════════════════════════════════════════════════════════════════
# Self-test
# ═══════════════════════════════════════════════════════════════════════

if __name__ == "__main__":
    print("=" * 60)
    print("EnhancedPAMScanner  --  self-test")
    print("=" * 60)

    # Synthetic 400-bp sequence with a known NGG site near the middle
    # The edit is at position 200.
    # We embed an NGG PAM at position 203 (so protospacer is 183-203, cut at 200)
    import random
    random.seed(42)

    # Build a controlled sequence
    # Positions 0-179: random (avoiding accidental TTTT or extreme GC)
    def _safe_random_seq(length: int) -> str:
        bases = "ATCG"
        seq = []
        for _ in range(length):
            seq.append(random.choice(bases))
        # break poly-T runs
        result = "".join(seq)
        while "TTTT" in result:
            result = result.replace("TTTT", "ATCG", 1)
        return result

    left = _safe_random_seq(183)
    # 20-mer protospacer with good GC (10G/C out of 20 = 50%)
    proto = "ATCGATCGATCGATCGATCG"   # GC = 50%
    pam_ngg = "AGG"
    right = _safe_random_seq(400 - 183 - 20 - 3)
    synth_seq = left + proto + pam_ngg + right
    edit_pos = 200  # should be within 3 bp of cut

    print(f"Sequence length : {len(synth_seq)}")
    print(f"Edit position   : {edit_pos}")
    print(f"Embedded PAM at : {183 + 20} (NGG = {synth_seq[203:206]})")

    # SpCas9 scan
    scanner_sp = EnhancedPAMScanner("SpCas9")
    guides_sp = scanner_sp.scan(synth_seq, edit_pos, window_bp=200)
    print(f"\nSpCas9 guides found: {len(guides_sp)}")
    for i, g in enumerate(guides_sp[:5]):
        print(f"  #{i+1}  score={g.score:.3f}  dist={g.distance_to_edit:+d}"
              f"  GC={g.gc_content:.2f}  strand={g.strand}  PAM={g.pam_sequence}")

    # enFnCas9 scan (NRG PAM -- should find MORE sites than NGG)
    scanner_en = EnhancedPAMScanner("enFnCas9")
    guides_en = scanner_en.scan(synth_seq, edit_pos, window_bp=200)
    print(f"\nenFnCas9 guides found: {len(guides_en)}")
    for i, g in enumerate(guides_en[:5]):
        print(f"  #{i+1}  score={g.score:.3f}  dist={g.distance_to_edit:+d}"
              f"  GC={g.gc_content:.2f}  strand={g.strand}  PAM={g.pam_sequence}")

    # ── assertions ────────────────────────────────────────────────────
    assert len(guides_sp) > 0, "SpCas9 must find at least 1 guide"
    assert len(guides_en) >= len(guides_sp), (
        f"enFnCas9 (NRG) must find >= SpCas9 (NGG) sites, "
        f"got {len(guides_en)} vs {len(guides_sp)}"
    )
    # Best guide should be sorted first (highest score)
    if len(guides_sp) > 1:
        assert guides_sp[0].score >= guides_sp[1].score, "Guides must be sorted best-first"
    # No guide should have extreme GC
    for g in guides_sp:
        assert 0.30 <= g.gc_content <= 0.70, f"GC filter failed: {g.gc_content}"
    # No guide should have poly-T
    for g in guides_sp:
        assert not g.has_poly_t, "Poly-T filter failed"

    print("\nAll assertions passed.  PASS")
