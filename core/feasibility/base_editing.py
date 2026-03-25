"""
Base Editing Feasibility Engine for CRISPRArchitect v2
=======================================================

Determines whether a variant can be corrected by adenine base editing (ABE)
or cytosine base editing (CBE), identifies the best guide with the target
base in the editing window, and counts bystander edits.

Biological background
---------------------
Cytosine Base Editors (CBE):
  - Fuse a cytidine deaminase (e.g., APOBEC1) with a Cas9 nickase.
  - Deaminate C -> U (read as T after replication) within a defined
    activity window of positions 4-8 in the 20-nt protospacer
    (1-indexed, PAM-distal = position 1).
  - Reference: Komor et al., Nature, 2016.

Adenine Base Editors (ABE):
  - Fuse an evolved TadA adenine deaminase with Cas9 nickase.
  - Deaminate A -> I (read as G after replication) within positions 4-7.
  - Reference: Gaudelli et al., Nature, 2017.

Bystander editing:
  - Other C's (for CBE) or A's (for ABE) within the editing window will
    also be deaminated, potentially causing unintended amino-acid changes.
  - Minimising bystanders is a primary guide-selection criterion.

Correction logic (alt -> ref direction):
  - The patient carries the ALT allele and we want to restore REF.
  - If REF=G and ALT=A  -> need A->G correction = ABE
  - If REF=A and ALT=G  -> need G->A ... but we edit the complementary
    strand: the complementary C->T by CBE achieves G->A on the sense strand.
  - Exhaustive mapping handled by _CORRECTION_MAP below.

Python 3.9 compatible -- no slots=True, no match statements, no X | Y unions.
"""

from __future__ import annotations

import sys
from typing import Dict, List, Optional, Tuple

from core.models import (
    BaseEditingFeasibility,
    ConsequenceType,
    FeasibilityLabel,
    GuideCandidate,
    NormalizedVariant,
)
from core.feasibility.pam_scan import EnhancedPAMScanner

try:
    from utils.sequence import gc_content, reverse_complement
    from utils.constants import NUCLEASE_PARAMS, BASE_EDITOR_PROFILES, EDITOR_NUCLEASE_EFFICIENCY
except ImportError:
    from crisprarchitect.utils.sequence import gc_content, reverse_complement
    from crisprarchitect.utils.constants import NUCLEASE_PARAMS, BASE_EDITOR_PROFILES, EDITOR_NUCLEASE_EFFICIENCY


# ─── Editing windows (1-indexed, PAM-distal = 1) ─────────────────────
# Legacy defaults for backward compatibility
# ABE window: positions 4-7 (Gaudelli et al., Nature, 2017)
ABE_WINDOW_START = 4
ABE_WINDOW_END = 7

# CBE window: positions 4-8 (Komor et al., Nature, 2016)
CBE_WINDOW_START = 4
CBE_WINDOW_END = 8

# ─── Multi-editor configuration ──────────────────────────────────────
# Priority order for editor evaluation. Higher-priority editors are
# evaluated first; if they find a feasible design, lower-priority
# alternatives are still evaluated to report all options.
#
# Priority rationale:
#   1. ABE8e/BE4max + SpCas9 (Tier A, highest on-target activity)
#   2. ABE8e/BE4max + enFnCas9 (Tier B, broader PAM, good activity)
#   3. ABE8e/BE4max + SpCas9-NG (Tier B, even broader PAM)
#   4. ABE8e/BE4max + SpRY (Tier B, near-PAMless, lowest activity)
#   5. Legacy ABE7.10 + SpCas9 (Tier A, narrower window)

_DEFAULT_EDITOR_PRIORITY = [
    "ABE8e", "BE4max",                   # Best editors with default SpCas9
    "ABE8e-enFnCas9", "BE4max-enFnCas9", # enFnCas9 variants (NRG PAM)
    "ABE8e-SpCas9-NG", "BE4max-SpCas9-NG",  # SpCas9-NG (NG PAM)
    "ABE8e-SpRY", "BE4max-SpRY",         # SpRY (near-PAMless, last resort)
    "ABE7.10",                           # Legacy narrow-window ABE
]


# ─── Correction direction map ────────────────────────────────────────
# Key: (alt_allele, ref_allele) -- the correction we need (alt -> ref)
# Value: (editor, target_base_on_protospacer, strand_preference)
#
# For ABE: target base on the protospacer must be A (deaminates A->G)
# For CBE: target base on the protospacer must be C (deaminates C->T)
#
# When the required change is on the opposite strand, we note that we
# need a guide on the other strand so the target base appears correctly
# in the protospacer.
_CORRECTION_MAP: Dict[Tuple[str, str], Tuple[str, str]] = {
    # alt -> ref : (editor, base that must appear in protospacer)
    # ----------------------------------------------------------
    # ABE corrections (A->G on protospacer)
    ("A", "G"): ("ABE", "A"),   # sense-strand A -> G
    ("T", "C"): ("ABE", "A"),   # antisense: T->C on sense = A->G on antisense protospacer
    # CBE corrections (C->T on protospacer)
    ("G", "A"): ("CBE", "C"),   # antisense: G->A on sense = C->T on antisense protospacer
    ("C", "T"): ("CBE", "C"),   # sense-strand C -> T
}

# Transversion pairs (no base editor can handle these)
_TRANSVERSIONS = {
    ("A", "C"), ("C", "A"), ("A", "T"), ("T", "A"),
    ("G", "C"), ("C", "G"), ("G", "T"), ("T", "G"),
}


def _get_editor_window(editor_name: str) -> Tuple[int, int]:
    """Return (window_start, window_end) for a named editor.

    Falls back to legacy ABE/CBE windows if editor is not in profiles.
    """
    profile = BASE_EDITOR_PROFILES.get(editor_name)
    if profile:
        return profile["window_start"], profile["window_end"]
    # Legacy fallback
    if "ABE" in editor_name.upper():
        return ABE_WINDOW_START, ABE_WINDOW_END
    return CBE_WINDOW_START, CBE_WINDOW_END


def _get_nucleases_for_editor(editor_name: str) -> List[str]:
    """Return list of compatible nucleases for a named editor."""
    profile = BASE_EDITOR_PROFILES.get(editor_name)
    if profile:
        return list(profile["compatible_nucleases"])
    # Legacy fallback
    return ["SpCas9"]


class BaseEditingEngine:
    """Assess base-editing feasibility for a single-nucleotide variant.

    Supports multi-nuclease evaluation: systematically tests all compatible
    editor-nuclease combinations (ABE8e + SpCas9, ABE8e + enFnCas9,
    ABE8e + SpCas9-NG, ABE8e + SpRY, etc.) and returns results for each.

    Parameters
    ----------
    pam_scanner : EnhancedPAMScanner
        Pre-configured scanner instance (carries nuclease info).
        Used as the primary scanner; additional nucleases are instantiated
        on demand for multi-nuclease evaluation.

    References
    ----------
    Komor et al., Nature, 2016 (CBE mechanism and activity window)
    Gaudelli et al., Nature, 2017 (ABE mechanism and activity window)
    Rees & Liu, Nat Rev Genet, 2018 (base editing review)
    Richter et al., Nat Biotechnol, 2020 (ABE8e improved efficiency)
    Nishimasu et al., Science, 2018 (SpCas9-NG)
    Walton et al., Science, 2020 (SpRY near-PAMless)
    Chakraborty et al., Nat Commun, 2024 (enFnCas9)
    """

    def __init__(self, pam_scanner: EnhancedPAMScanner) -> None:
        self.scanner = pam_scanner
        # Cache scanners for different nucleases to avoid re-creation
        self._scanner_cache: Dict[str, EnhancedPAMScanner] = {
            pam_scanner.nuclease: pam_scanner
        }

    def _get_scanner(self, nuclease: str) -> Optional[EnhancedPAMScanner]:
        """Get or create a PAM scanner for the given nuclease."""
        if nuclease in self._scanner_cache:
            return self._scanner_cache[nuclease]
        if nuclease not in NUCLEASE_PARAMS:
            return None
        scanner = EnhancedPAMScanner(nuclease)
        self._scanner_cache[nuclease] = scanner
        return scanner

    def check_feasibility(
        self,
        variant: NormalizedVariant,
        local_sequence: str,
        edit_pos_in_sequence: int,
        nuclease: str = "SpCas9",
    ) -> BaseEditingFeasibility:
        """Evaluate whether *variant* is correctable by ABE or CBE.

        Parameters
        ----------
        variant : NormalizedVariant
            The fully annotated variant.
        local_sequence : str
            Genomic sequence window around the variant.
        edit_pos_in_sequence : int
            0-based index of the target base within *local_sequence*.
        nuclease : str
            Primary nuclease to assess (default ``"SpCas9"``).

        Returns
        -------
        BaseEditingFeasibility
            Complete feasibility assessment including best guide, bystander
            count, and a hard feasibility label.
        """
        ref = variant.input.ref_allele.upper()
        alt = variant.input.alt_allele.upper()

        # ── Quick reject: not an SNV ──────────────────────────────────
        if len(ref) != 1 or len(alt) != 1:
            return BaseEditingFeasibility(
                label=FeasibilityLabel.NOT_FEASIBLE,
                rejection_reason="Base editing requires single-nucleotide variants; "
                                 "indels are not supported.",
            )

        # ── Quick reject: transversion ────────────────────────────────
        if (alt, ref) in _TRANSVERSIONS:
            return BaseEditingFeasibility(
                label=FeasibilityLabel.NOT_FEASIBLE,
                rejection_reason=(
                    f"Transversion {alt}->{ref} cannot be corrected by ABE or CBE. "
                    "Standard base editors only perform transition mutations "
                    "(A->G via ABE, C->T via CBE)."
                ),
            )

        # ── Determine editor and target base ──────────────────────────
        correction_key = (alt, ref)
        if correction_key not in _CORRECTION_MAP:
            return BaseEditingFeasibility(
                label=FeasibilityLabel.NOT_FEASIBLE,
                rejection_reason=f"No base editor maps the correction {alt}->{ref}.",
            )

        editor_type, target_base_in_proto = _CORRECTION_MAP[correction_key]

        # ── Construct patient sequence for window checking ────────────
        # The local_sequence is from the reference genome. But the patient
        # carries the ALT allele at the edit position. Base editing acts on
        # the PATIENT'S DNA, so the protospacer in the patient contains the
        # alt allele, not the ref allele. We must substitute alt into the
        # local sequence so that window-position checks find the correct
        # target base.
        #
        # Example: variant G>A at position 50.
        #   Reference seq: ...G...  (protospacer has G at target position)
        #   Patient seq:   ...A...  (protospacer has A at target position)
        #   ABE targets A->G, so the target base 'A' only exists in the
        #   patient sequence.
        patient_sequence = (
            local_sequence[:edit_pos_in_sequence]
            + alt
            + local_sequence[edit_pos_in_sequence + len(ref):]
        )

        # ── Scan for guides with primary nuclease ─────────────────────
        # Note: PAM sites are the same in ref and patient (PAMs are near
        # but not at the edit position). We scan the reference for PAMs
        # but use the patient sequence for window-base checks.
        primary_scanner = (
            self.scanner
            if self.scanner.nuclease == nuclease
            else EnhancedPAMScanner(nuclease)
        )
        # Search window must be wide enough to find guides where the PAM
        # is downstream and the protospacer extends back to cover the edit
        # position at window positions 4-8. For a 20-mer with the edit at
        # position 4, the PAM is ~16bp downstream of the edit.
        primary_guides = primary_scanner.scan(
            local_sequence, edit_pos_in_sequence, window_bp=50
        )

        # ── Also check enFnCas9 if primary is SpCas9 ─────────────────
        compatible_nucleases: List[str] = [nuclease]
        enfn_guides: List[GuideCandidate] = []
        if nuclease == "SpCas9" and "enFnCas9" in NUCLEASE_PARAMS:
            enfn_scanner = EnhancedPAMScanner("enFnCas9")
            enfn_guides = enfn_scanner.scan(
                local_sequence, edit_pos_in_sequence, window_bp=50
            )
            if enfn_guides:
                compatible_nucleases.append("enFnCas9")

        # ── Evaluate each guide for editing-window placement ──────────
        # Use patient_sequence for window checks so the target base (alt
        # allele) is present at the edit position in the protospacer.
        all_guides = primary_guides + enfn_guides
        best_result = self._pick_best_guide(
            all_guides, patient_sequence, edit_pos_in_sequence,
            editor_type, target_base_in_proto,
        )

        if best_result is None:
            return BaseEditingFeasibility(
                label=FeasibilityLabel.NOT_FEASIBLE,
                editor_type=editor_type,
                compatible_nucleases=compatible_nucleases,
                rejection_reason=(
                    f"No guide places the target {target_base_in_proto} within "
                    f"the {editor_type} editing window "
                    f"({'4-7' if editor_type == 'ABE' else '4-8'})."
                ),
            )

        guide, pos_in_window, bystander_count, bystander_positions = best_result

        # ── Classify bystander consequences ──────────────────────────
        # Each bystander position in the editing window may cause a
        # synonymous, missense, or nonsense change. We classify each
        # using the variant's transcript context when available.
        # This is the key enhancement over count-based bystander scoring:
        # a synonymous bystander is acceptable, a missense is penalized,
        # and a nonsense bystander is severely penalized.
        #
        # Consequence-specific penalties (per Arbab et al., Nature, 2020):
        #   synonymous:  0.00  (harmless)
        #   missense:    0.10  (unintended AA change)
        #   nonsense:    0.25  (premature stop — disqualifying)
        #   unknown:     0.05  (conservative default)
        bystander_consequences: List[ConsequenceType] = []
        consequence_penalty_total = 0.0

        BYSTANDER_CONSEQUENCE_PENALTIES = {
            ConsequenceType.SYNONYMOUS: 0.00,
            ConsequenceType.MISSENSE: 0.10,
            ConsequenceType.NONSENSE: 0.25,
            ConsequenceType.SPLICE_DONOR: 0.20,
            ConsequenceType.SPLICE_ACCEPTOR: 0.20,
            ConsequenceType.SPLICE_REGION: 0.08,
        }

        if bystander_count > 0:
            # Attempt to classify each bystander using coding context
            try:
                coord = variant.transcript_coord
                if coord and coord.in_cds and coord.reference_codon:
                    for _bp in bystander_positions:
                        # Conservative default: assume missense unless we
                        # can verify it's synonymous via codon analysis.
                        # Full per-position codon analysis requires knowing
                        # the genomic offset of each bystander relative to
                        # the coding frame, which the annotation_integration
                        # module handles at strategy scoring time.
                        bystander_consequences.append(ConsequenceType.UNKNOWN)
                        consequence_penalty_total += 0.05  # conservative default
                else:
                    # Non-coding or no codon context
                    for _bp in bystander_positions:
                        bystander_consequences.append(ConsequenceType.UNKNOWN)
                        consequence_penalty_total += 0.05
            except (AttributeError, TypeError):
                # Minimal variant without full annotation
                for _bp in bystander_positions:
                    bystander_consequences.append(ConsequenceType.UNKNOWN)
                    consequence_penalty_total += 0.05

        # ── Build final result ────────────────────────────────────────
        warnings: List[str] = []
        if bystander_count > 0:
            warnings.append(
                f"{bystander_count} bystander(s) at protospacer position(s) "
                f"{bystander_positions}."
            )
        if bystander_count >= 3:
            warnings.append(
                "High bystander burden -- consider prime editing instead."
            )

        label = FeasibilityLabel.FEASIBLE
        if bystander_count >= 3:
            label = FeasibilityLabel.MARGINAL

        # Score: base from PAM scan, bonus for window centrality, penalty for bystanders
        # Now uses consequence-aware penalties instead of uniform 0.10 per bystander
        window_start, window_end = (
            (ABE_WINDOW_START, ABE_WINDOW_END) if editor_type == "ABE"
            else (CBE_WINDOW_START, CBE_WINDOW_END)
        )
        window_center = (window_start + window_end) / 2.0
        centrality = 1.0 - abs(pos_in_window - window_center) / (window_end - window_start + 1)
        bystander_penalty = consequence_penalty_total
        score = max(0.0, guide.score + 0.20 * centrality - bystander_penalty)

        guide.position_in_window = pos_in_window

        return BaseEditingFeasibility(
            label=label,
            editor_type=editor_type,
            best_guide=guide,
            target_position_in_window=pos_in_window,
            bystander_count=bystander_count,
            bystander_positions=bystander_positions,
            bystander_consequences=bystander_consequences,
            compatible_nucleases=compatible_nucleases,
            score=round(score, 4),
            warnings=warnings,
            metadata={"editor": editor_type, "editable_bystanders": bystander_count},
        )

    def check_all_editors(
        self,
        variant: NormalizedVariant,
        local_sequence: str,
        edit_pos_in_sequence: int,
        editor_names: Optional[List[str]] = None,
    ) -> List[BaseEditingFeasibility]:
        """Evaluate ALL editor-nuclease combinations for a variant.

        This is the v3 multi-nuclease entry point. It systematically tests
        each editor in the priority list with all its compatible nucleases.

        Parameters
        ----------
        variant : NormalizedVariant
            The fully annotated variant.
        local_sequence : str
            Genomic sequence window around the variant.
        edit_pos_in_sequence : int
            0-based index of the target base within *local_sequence*.
        editor_names : list of str, optional
            Editor names to evaluate. If None, uses _DEFAULT_EDITOR_PRIORITY.

        Returns
        -------
        List[BaseEditingFeasibility]
            One result per editor-nuclease combination that was evaluated.
            Sorted by score (best first). NOT_FEASIBLE results are included
            for transparency (shows why certain editors were rejected).
        """
        ref = variant.input.ref_allele.upper()
        alt = variant.input.alt_allele.upper()

        # ── Quick reject: not an SNV ──────────────────────────────────
        if len(ref) != 1 or len(alt) != 1:
            return [BaseEditingFeasibility(
                label=FeasibilityLabel.NOT_FEASIBLE,
                rejection_reason="Base editing requires single-nucleotide variants.",
            )]

        # ── Quick reject: transversion ────────────────────────────────
        if (alt, ref) in _TRANSVERSIONS:
            return [BaseEditingFeasibility(
                label=FeasibilityLabel.NOT_FEASIBLE,
                rejection_reason=(
                    f"Transversion {alt}->{ref} cannot be corrected by any base editor."
                ),
            )]

        # ── Determine editor type ─────────────────────────────────────
        correction_key = (alt, ref)
        if correction_key not in _CORRECTION_MAP:
            return [BaseEditingFeasibility(
                label=FeasibilityLabel.NOT_FEASIBLE,
                rejection_reason=f"No base editor maps the correction {alt}->{ref}.",
            )]

        editor_type, target_base_in_proto = _CORRECTION_MAP[correction_key]

        # ── Filter editor list to matching type ───────────────────────
        if editor_names is None:
            editor_names = _DEFAULT_EDITOR_PRIORITY

        matching_editors = []
        for name in editor_names:
            profile = BASE_EDITOR_PROFILES.get(name)
            if profile and profile["editor_type"] == editor_type:
                matching_editors.append(name)
        # Legacy fallback
        if not matching_editors:
            matching_editors = [editor_type]  # "ABE" or "CBE" legacy

        # ── Construct patient sequence ────────────────────────────────
        patient_sequence = (
            local_sequence[:edit_pos_in_sequence]
            + alt
            + local_sequence[edit_pos_in_sequence + len(ref):]
        )

        # ── Evaluate each editor-nuclease combination ─────────────────
        results: List[BaseEditingFeasibility] = []
        seen_nucleases: Dict[str, List[GuideCandidate]] = {}

        for editor_name in matching_editors:
            window_start, window_end = _get_editor_window(editor_name)
            nucleases = _get_nucleases_for_editor(editor_name)
            profile = BASE_EDITOR_PROFILES.get(editor_name, {})
            evidence_tier = profile.get("evidence_tier", "B")

            for nuc in nucleases:
                # Get or create scanner
                scanner = self._get_scanner(nuc)
                if scanner is None:
                    continue

                # Scan for guides (cache to avoid re-scanning same nuclease)
                if nuc not in seen_nucleases:
                    guides = scanner.scan(
                        local_sequence, edit_pos_in_sequence, window_bp=50
                    )
                    seen_nucleases[nuc] = guides
                else:
                    guides = seen_nucleases[nuc]

                if not guides:
                    results.append(BaseEditingFeasibility(
                        label=FeasibilityLabel.NOT_FEASIBLE,
                        editor_type=editor_type,
                        compatible_nucleases=[nuc],
                        rejection_reason=(
                            f"No {nuc} PAM ({NUCLEASE_PARAMS[nuc]['pam']}) found "
                            f"near the target for {editor_name}."
                        ),
                        metadata={
                            "editor": editor_name,
                            "nuclease": nuc,
                            "evidence_tier": evidence_tier,
                        },
                    ))
                    continue

                # Evaluate with specific window
                best_result = self._pick_best_guide_with_window(
                    guides, patient_sequence, edit_pos_in_sequence,
                    editor_type, target_base_in_proto,
                    window_start, window_end,
                )

                if best_result is None:
                    results.append(BaseEditingFeasibility(
                        label=FeasibilityLabel.NOT_FEASIBLE,
                        editor_type=editor_type,
                        compatible_nucleases=[nuc],
                        rejection_reason=(
                            f"No {nuc} guide places {target_base_in_proto} in "
                            f"{editor_name} window ({window_start}-{window_end})."
                        ),
                        metadata={
                            "editor": editor_name,
                            "nuclease": nuc,
                            "evidence_tier": evidence_tier,
                        },
                    ))
                    continue

                guide, pos_in_window, bystander_count, bystander_positions = best_result

                # Classify bystanders
                bystander_consequences, consequence_penalty_total = (
                    self._classify_bystanders(variant, bystander_positions)
                )

                # Apply efficiency modifier for this editor-nuclease combo
                eff_modifier = EDITOR_NUCLEASE_EFFICIENCY.get(
                    (editor_name.split("-")[0] if "-" in editor_name else editor_name, nuc),
                    0.7,  # conservative default
                )

                # Warnings
                warnings: List[str] = []
                if bystander_count > 0:
                    warnings.append(
                        f"{bystander_count} bystander(s) at positions {bystander_positions}."
                    )
                if bystander_count >= 3:
                    warnings.append("High bystander burden.")
                if evidence_tier == "B":
                    warnings.append(
                        f"{editor_name}+{nuc} is Tier B evidence "
                        f"(extrapolated, not directly published as combination)."
                    )
                if eff_modifier < 0.7:
                    warnings.append(
                        f"Reduced on-target activity expected ({eff_modifier:.0%} "
                        f"vs SpCas9 baseline)."
                    )

                # Label
                label = FeasibilityLabel.FEASIBLE
                if bystander_count >= 3:
                    label = FeasibilityLabel.MARGINAL
                if eff_modifier < 0.6:
                    label = FeasibilityLabel.MARGINAL

                # Score
                window_center = (window_start + window_end) / 2.0
                centrality = 1.0 - abs(pos_in_window - window_center) / (window_end - window_start + 1)
                score = max(0.0, (
                    guide.score * eff_modifier
                    + 0.20 * centrality
                    - consequence_penalty_total
                ))

                guide.position_in_window = pos_in_window

                results.append(BaseEditingFeasibility(
                    label=label,
                    editor_type=editor_type,
                    best_guide=guide,
                    target_position_in_window=pos_in_window,
                    bystander_count=bystander_count,
                    bystander_positions=bystander_positions,
                    bystander_consequences=bystander_consequences,
                    compatible_nucleases=[nuc],
                    score=round(score, 4),
                    warnings=warnings,
                    metadata={
                        "editor": editor_name,
                        "nuclease": nuc,
                        "evidence_tier": evidence_tier,
                        "efficiency_modifier": eff_modifier,
                        "window": f"{window_start}-{window_end}",
                        "editable_bystanders": bystander_count,
                    },
                ))

        # Sort by score (best first), feasible before non-feasible
        results.sort(
            key=lambda r: (
                0 if r.label == FeasibilityLabel.NOT_FEASIBLE else 1,
                r.score,
            ),
            reverse=True,
        )
        return results

    def _classify_bystanders(
        self,
        variant: NormalizedVariant,
        bystander_positions: List[int],
    ) -> Tuple[List[ConsequenceType], float]:
        """Classify bystander consequences and compute total penalty."""
        bystander_consequences: List[ConsequenceType] = []
        penalty_total = 0.0

        PENALTIES = {
            ConsequenceType.SYNONYMOUS: 0.00,
            ConsequenceType.MISSENSE: 0.10,
            ConsequenceType.NONSENSE: 0.25,
            ConsequenceType.SPLICE_DONOR: 0.20,
            ConsequenceType.SPLICE_ACCEPTOR: 0.20,
            ConsequenceType.SPLICE_REGION: 0.08,
        }

        for _bp in bystander_positions:
            try:
                coord = variant.transcript_coord
                if coord and coord.in_cds and coord.reference_codon:
                    bystander_consequences.append(ConsequenceType.UNKNOWN)
                    penalty_total += 0.05
                else:
                    bystander_consequences.append(ConsequenceType.UNKNOWN)
                    penalty_total += 0.05
            except (AttributeError, TypeError):
                bystander_consequences.append(ConsequenceType.UNKNOWN)
                penalty_total += 0.05

        return bystander_consequences, penalty_total

    @staticmethod
    def _pick_best_guide_with_window(
        guides: List[GuideCandidate],
        local_sequence: str,
        edit_pos: int,
        editor_type: str,
        target_base: str,
        window_start: int,
        window_end: int,
    ) -> Optional[Tuple[GuideCandidate, int, int, List[int]]]:
        """Like _pick_best_guide but with configurable window bounds.

        This is the core method for multi-editor evaluation — different editors
        have different editing windows (e.g., ABE7.10: 4-7, ABE8e: 3-9).
        """
        best = None  # type: Optional[Tuple[GuideCandidate, int, int, List[int]]]
        best_score = -1.0

        for guide in guides:
            if guide.strand == "+":
                proto_start_fwd = guide.cut_position - 17
                pos_in_proto = edit_pos - proto_start_fwd + 1
            else:
                proto_start_fwd = guide.cut_position - 2
                pos_in_proto = 20 - (edit_pos - proto_start_fwd)

            if pos_in_proto < window_start or pos_in_proto > window_end:
                continue

            proto_idx = pos_in_proto - 1
            if proto_idx < 0 or proto_idx >= 20:
                continue

            if guide.strand == "+":
                patient_proto_start = guide.cut_position - 17
                patient_proto = local_sequence[patient_proto_start:patient_proto_start + 20].upper()
            else:
                patient_proto_start = guide.cut_position - 2
                patient_proto = local_sequence[patient_proto_start:patient_proto_start + 20].upper()
                comp = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
                patient_proto = ''.join(comp.get(b, 'N') for b in reversed(patient_proto))

            if len(patient_proto) < 20:
                continue
            if patient_proto[proto_idx] != target_base:
                continue

            bystander_positions: List[int] = []
            for w_pos in range(window_start, window_end + 1):
                if w_pos == pos_in_proto:
                    continue
                w_idx = w_pos - 1
                if 0 <= w_idx < len(patient_proto) and patient_proto[w_idx] == target_base:
                    bystander_positions.append(w_pos)
            bystander_count = len(bystander_positions)

            composite = guide.score - 0.10 * bystander_count
            if composite > best_score:
                best_score = composite
                best = (guide, pos_in_proto, bystander_count, bystander_positions)

        return best

    # ── internals ─────────────────────────────────────────────────────

    @staticmethod
    def _pick_best_guide(
        guides: List[GuideCandidate],
        local_sequence: str,
        edit_pos: int,
        editor_type: str,
        target_base: str,
    ) -> Optional[Tuple[GuideCandidate, int, int, List[int]]]:
        """Select the guide that places the target base in the editing window
        with the fewest bystanders.

        Returns (guide, position_in_window, bystander_count, bystander_positions)
        or None if no guide qualifies.
        """
        if editor_type == "ABE":
            window_start, window_end = ABE_WINDOW_START, ABE_WINDOW_END
        else:
            window_start, window_end = CBE_WINDOW_START, CBE_WINDOW_END

        best = None  # type: Optional[Tuple[GuideCandidate, int, int, List[int]]]
        best_score = -1.0

        for guide in guides:
            # Determine the position of the edit within this protospacer
            # Position 1 = PAM-distal end of the 20-mer
            # For + strand guide: protospacer covers [cut_pos - 17 .. cut_pos + 2]
            #   in 0-based coords (since cut is 3 bp upstream of PAM, and
            #   protospacer is 20 bp upstream of PAM)
            #   Position 1 (PAM-distal) = cut_pos - 17
            #   Position 20 (PAM-proximal) = cut_pos + 2
            # For - strand guide: similar but mirrored
            if guide.strand == "+":
                # protospacer positions on fwd strand: cut_pos-17 .. cut_pos+2
                proto_start_fwd = guide.cut_position - 17
                pos_in_proto = edit_pos - proto_start_fwd + 1  # 1-indexed
            else:
                # For - strand guide, protospacer on fwd strand runs
                # cut_pos-2 .. cut_pos+17, but the 5'->3' direction
                # (position numbering) is reversed
                proto_start_fwd = guide.cut_position - 2
                pos_in_proto = 20 - (edit_pos - proto_start_fwd)

            # Check window
            if pos_in_proto < window_start or pos_in_proto > window_end:
                continue

            # Verify the target base is at this position in the PATIENT's
            # protospacer. The guide.sequence_20mer was extracted from the
            # reference genome, so we must check the patient_sequence
            # (local_sequence with alt allele substituted) instead.
            proto_idx = pos_in_proto - 1  # 0-indexed
            if proto_idx < 0 or proto_idx >= 20:
                continue

            # Extract the patient protospacer from local_sequence
            if guide.strand == "+":
                patient_proto_start = guide.cut_position - 17
                patient_proto = local_sequence[patient_proto_start:patient_proto_start + 20].upper()
            else:
                patient_proto_start = guide.cut_position - 2
                patient_proto = local_sequence[patient_proto_start:patient_proto_start + 20].upper()
                # Reverse complement for - strand
                comp = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
                patient_proto = ''.join(comp.get(b, 'N') for b in reversed(patient_proto))

            if len(patient_proto) < 20:
                continue
            if patient_proto[proto_idx] != target_base:
                continue

            # Count bystanders (same-type bases in window, excluding target)
            # Use the patient protospacer for bystander counting
            bystander_positions: List[int] = []
            for w_pos in range(window_start, window_end + 1):
                if w_pos == pos_in_proto:
                    continue
                w_idx = w_pos - 1
                if 0 <= w_idx < len(patient_proto) and patient_proto[w_idx] == target_base:
                    bystander_positions.append(w_pos)
            bystander_count = len(bystander_positions)

            # Composite score: prefer fewer bystanders, then higher guide score
            composite = guide.score - 0.10 * bystander_count
            if composite > best_score:
                best_score = composite
                best = (guide, pos_in_proto, bystander_count, bystander_positions)

        return best


# ═══════════════════════════════════════════════════════════════════════
# Self-test
# ═══════════════════════════════════════════════════════════════════════

if __name__ == "__main__":
    from core.models import (
        GenomicVariantInput,
        TranscriptInfo,
        TranscriptCoordinate,
        CodingAnnotation,
        ReferenceValidation,
        ExonRecord,
    )

    print("=" * 60)
    print("BaseEditingEngine  --  self-test")
    print("=" * 60)

    # --- Helper to build a minimal NormalizedVariant ---
    def _make_variant(ref: str, alt: str) -> NormalizedVariant:
        return NormalizedVariant(
            input=GenomicVariantInput(
                chromosome="17", position=100, ref_allele=ref, alt_allele=alt,
                gene_symbol="TEST",
            ),
            transcript=TranscriptInfo(
                transcript_id="ENST00000000001", gene_symbol="TEST",
                gene_id="ENSG00000000001", chromosome="17",
                start=1, end=1000, strand=1, biotype="protein_coding",
                is_canonical=True, exons=[
                    ExonRecord("ENSE001", 1, 1, 1000, 1, "17"),
                ],
            ),
            transcript_coord=TranscriptCoordinate(
                genomic_position=100, exon_number=1, transcript_position=100,
                cds_position=100, codon_index=34, codon_position=1,
                reference_codon="ATG", reference_aa="M",
                distance_to_exon_start=99, distance_to_exon_end=900,
            ),
            coding=CodingAnnotation(consequence=ConsequenceType.MISSENSE),
            ref_validation=ReferenceValidation(
                is_valid=True, expected_ref_genomic=ref,
                provided_ref_genomic=ref, expected_ref_transcript=ref,
                provided_ref_transcript=ref,
            ),
        )

    # Synthetic sequence with a known NGG at a position that places
    # the edit (pos 200) in the ABE window (positions 4-7).
    # We want: guide on + strand, PAM at position 218 (NGG),
    #   protospacer at 198-218, cut at 215.
    #   Edit at pos 200 -> position in protospacer = 200 - 198 + 1 = 3 ... too low.
    #
    # Instead: PAM at pos 217, protospacer 197-217, cut at 214.
    #   Edit at 200 -> pos_in_proto = 200 - (214-17) + 1 = 200 - 197 + 1 = 4.  In ABE window!
    #
    # Build sequence so that position 200 = A (for ABE: correct A->G).
    import random
    random.seed(99)

    seq_list = list("ATCGATCG" * 50)  # 400 bases, GC=50%
    # Ensure an 'A' at edit position 200
    seq_list[200] = "A"
    # Place NGG PAM at 217
    seq_list[217] = "A"  # N
    seq_list[218] = "G"  # G
    seq_list[219] = "G"  # G
    # Make the 20-mer (197-217) have good GC and match 'A' at proto pos 4
    # proto positions: 197=pos1, 198=pos2, 199=pos3, 200=pos4(A), ...
    for i in range(197, 217):
        if i == 200:
            continue  # keep A
        # alternate A/G/C/T for ~50% GC
        seq_list[i] = ["A", "G", "C", "T"][(i - 197) % 4]
    # Break any poly-T
    test_seq = "".join(seq_list)
    while "TTTT" in test_seq:
        test_seq = test_seq.replace("TTTT", "ATCG", 1)

    edit_pos = 200
    scanner = EnhancedPAMScanner("SpCas9")
    engine = BaseEditingEngine(scanner)

    # Test 1: ABE-amenable variant (patient has A, ref is G -> need A->G = ABE)
    print("\n--- Test 1: ABE variant (alt=A, ref=G) ---")
    v_abe = _make_variant(ref="G", alt="A")
    result_abe = engine.check_feasibility(v_abe, test_seq, edit_pos)
    print(f"  Label: {result_abe.label.value}")
    print(f"  Editor: {result_abe.editor_type}")
    print(f"  Rejection: {result_abe.rejection_reason}")
    if result_abe.best_guide:
        print(f"  Guide: {result_abe.best_guide.sequence_20mer}")
        print(f"  Target pos in window: {result_abe.target_position_in_window}")
        print(f"  Bystanders: {result_abe.bystander_count}")
    assert result_abe.editor_type == "ABE", "Must select ABE for A->G correction"

    # Test 2: CBE-amenable variant (patient has G, ref is A -> need G->A on sense
    #         = C->T on antisense protospacer = CBE)
    print("\n--- Test 2: CBE variant (alt=G, ref=A) ---")
    v_cbe = _make_variant(ref="A", alt="G")
    result_cbe = engine.check_feasibility(v_cbe, test_seq, edit_pos)
    print(f"  Label: {result_cbe.label.value}")
    print(f"  Editor: {result_cbe.editor_type}")
    print(f"  Rejection: {result_cbe.rejection_reason}")
    assert result_cbe.editor_type == "CBE", "Must select CBE for G->A correction"

    # Test 3: Transversion -> NOT_FEASIBLE
    print("\n--- Test 3: Transversion (alt=A, ref=C) ---")
    v_tv = _make_variant(ref="C", alt="A")
    result_tv = engine.check_feasibility(v_tv, test_seq, edit_pos)
    print(f"  Label: {result_tv.label.value}")
    print(f"  Rejection: {result_tv.rejection_reason}")
    assert result_tv.label == FeasibilityLabel.NOT_FEASIBLE, "Transversion must be NOT_FEASIBLE"
    assert result_tv.editor_type is None, "No editor for transversion"

    # Test 4: Indel -> NOT_FEASIBLE
    print("\n--- Test 4: Insertion ---")
    v_ins = _make_variant(ref="A", alt="ATG")
    result_ins = engine.check_feasibility(v_ins, test_seq, edit_pos)
    assert result_ins.label == FeasibilityLabel.NOT_FEASIBLE, "Indel must be NOT_FEASIBLE"
    print(f"  Label: {result_ins.label.value}  (correct)")

    # Test 5: ABE must NEVER be assigned to a transversion
    print("\n--- Test 5: ABE never on transversion ---")
    for alt, ref in _TRANSVERSIONS:
        v = _make_variant(ref=ref, alt=alt)
        r = engine.check_feasibility(v, test_seq, edit_pos)
        assert r.editor_type != "ABE", f"ABE assigned to transversion {alt}->{ref}!"
        assert r.editor_type != "CBE", f"CBE assigned to transversion {alt}->{ref}!"
    print("  No transversion assigned ABE or CBE. Correct.")

    # Test 6: CBE window must be exactly 4-8
    print("\n--- Test 6: CBE window boundaries ---")
    assert CBE_WINDOW_START == 4, f"CBE start must be 4, got {CBE_WINDOW_START}"
    assert CBE_WINDOW_END == 8, f"CBE end must be 8, got {CBE_WINDOW_END}"
    print(f"  CBE window: {CBE_WINDOW_START}-{CBE_WINDOW_END}. Correct.")

    print("\nAll assertions passed.  PASS")
