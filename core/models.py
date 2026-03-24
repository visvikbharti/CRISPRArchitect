"""
CRISPRArchitect v2 — Central Data Models
==========================================

All shared dataclasses and enums used across the v2 pipeline.

Design principles:
- Every field has a type annotation and a sensible default
- Enums encode biological categories; strings encode free-form metadata
- NormalizedVariant bridges to v1's Mutation dataclass for backward compat
- All genomic coordinates are 1-based (Ensembl convention)
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Any, Dict, List, Optional, Tuple, TYPE_CHECKING

if TYPE_CHECKING:
    pass  # future forward refs


# ═══════════════════════════════════════════════════════════════════════════
# Enums
# ═══════════════════════════════════════════════════════════════════════════

class ConsequenceType(Enum):
    """Functional consequence of a coding variant.

    Classification follows Ensembl Variation consequence hierarchy.
    Splice site distances follow ACMG/AMP standards:
      - splice_donor / splice_acceptor: within 2 bp of exon boundary
      - splice_region: within 3-8 bp of exon boundary

    References
    ----------
    McLaren et al., Genome Biology, 2016 (Ensembl VEP consequence types)
    Richards et al., Genetics in Medicine, 2015 (ACMG splice definitions)
    """
    SYNONYMOUS = "synonymous"
    MISSENSE = "missense"
    NONSENSE = "nonsense"
    SPLICE_DONOR = "splice_donor"
    SPLICE_ACCEPTOR = "splice_acceptor"
    SPLICE_REGION = "splice_region"
    FRAMESHIFT = "frameshift"
    INFRAME_INSERTION = "inframe_insertion"
    INFRAME_DELETION = "inframe_deletion"
    INTRONIC = "intronic"
    UTR_5 = "5_prime_UTR"
    UTR_3 = "3_prime_UTR"
    NON_CODING = "non_coding"
    UNKNOWN = "unknown"


class EditModality(Enum):
    """Genome editing modality.

    Each value represents a specific editor + donor combination.
    The hierarchy reflects increasing complexity and DSB burden:
    ABE/CBE (no DSB) > PE (nick only) > HDR variants (DSB required).

    References
    ----------
    Komor et al., Nature, 2016 (CBE)
    Gaudelli et al., Nature, 2017 (ABE)
    Anzalone et al., Nature, 2019 (PE)
    """
    ABE = "ABE"
    CBE = "CBE"
    PE = "PE"
    HDR_SSODN = "HDR_ssODN"
    HDR_CSSDNA = "HDR_cssDNA"
    HDR_LSSDNA = "HDR_lssDNA"
    HDR_DSDNA = "HDR_dsDNA"
    EXON_DELETION = "exon_deletion"


class FeasibilityLabel(Enum):
    """Hard feasibility verdict for a modality on a specific variant."""
    FEASIBLE = "feasible"
    MARGINAL = "marginal"
    NOT_FEASIBLE = "not_feasible"


class EvidenceTier(Enum):
    """Confidence tier for a strategy recommendation.

    Tier A: all components PAM-verified, window-verified
    Tier B: feasible but with caveats (bystanders, large donor, etc.)
    Tier C: theoretical only (no PAM found, extrapolated efficiency)
    """
    A = "A"
    B = "B"
    C = "C"


class RiskLevel(Enum):
    """Genomic rearrangement risk level."""
    LOW = "low"
    MODERATE = "moderate"
    HIGH = "high"
    VERY_HIGH = "very_high"


# ═══════════════════════════════════════════════════════════════════════════
# Transcript / Coordinate Models
# ═══════════════════════════════════════════════════════════════════════════

@dataclass
class ExonRecord:
    """One exon in a transcript.

    Attributes
    ----------
    exon_id : str
        Ensembl exon stable ID (e.g., "ENSE00003659301").
    exon_number : int
        1-based exon number within the transcript.
    start : int
        Genomic start coordinate (1-based, inclusive).
    end : int
        Genomic end coordinate (1-based, inclusive).
    strand : int
        1 for forward, -1 for reverse.
    chromosome : str
        Chromosome name (e.g., "17").
    """
    exon_id: str
    exon_number: int
    start: int
    end: int
    strand: int
    chromosome: str

    @property
    def length(self) -> int:
        return self.end - self.start + 1


@dataclass
class TranscriptInfo:
    """Transcript record fetched from Ensembl with CDS-level detail.

    All coordinates are 1-based genomic (Ensembl convention).
    For reverse-strand transcripts, start < end still holds
    (Ensembl always reports start < end).
    """
    transcript_id: str
    gene_symbol: str
    gene_id: str
    chromosome: str
    start: int
    end: int
    strand: int                        # 1 or -1
    biotype: str                       # "protein_coding", etc.
    is_canonical: bool
    exons: List[ExonRecord]

    @property
    def n_exons(self) -> int:
        return len(self.exons)

    @property
    def span_bp(self) -> int:
        return self.end - self.start + 1


@dataclass
class TranscriptCoordinate:
    """Result of mapping a genomic position to transcript context.

    All positions are 1-based. For non-coding positions, cds_position
    is set to -1 and in_cds is False.
    """
    genomic_position: int
    exon_number: int                   # 1-based; -1 if intronic
    transcript_position: int           # 1-based position in spliced exonic seq
    cds_position: int                  # 1-based position in CDS; -1 if non-CDS
    codon_index: int                   # 1-based codon number; -1 if non-CDS
    codon_position: int                # 1, 2, or 3 within codon; -1 if non-CDS
    reference_codon: str               # 3-letter codon (e.g., "ATG"); "" if N/A
    reference_aa: str                  # single-letter AA; "" if N/A
    distance_to_exon_start: int        # bp to 5' end of exon (transcript sense)
    distance_to_exon_end: int          # bp to 3' end of exon (transcript sense)
    in_cds: bool = True
    sequence_context: str = ""         # local sequence around position
    context_center_index: int = -1     # 0-based index of target base in context


# ═══════════════════════════════════════════════════════════════════════════
# Variant Models
# ═══════════════════════════════════════════════════════════════════════════

@dataclass
class GenomicVariantInput:
    """User-supplied variant before normalization.

    This is the raw input. The pipeline normalizes it into a
    NormalizedVariant with full transcript and consequence annotation.

    Attributes
    ----------
    chromosome : str
        Chromosome name (e.g., "17", "X"). No "chr" prefix.
    position : int
        1-based genomic position.
    ref_allele : str
        Reference allele. Use "-" for pure insertion.
    alt_allele : str
        Alternate (patient) allele. Use "-" for pure deletion.
    gene_symbol : str or None
        HGNC gene symbol. Required for transcript lookup.
    transcript_id : str or None
        Ensembl transcript ID. If None, canonical transcript is used.
    name : str or None
        Human-readable label (e.g., ClinVar ID, HGVS notation).
    """
    chromosome: str
    position: int
    ref_allele: str
    alt_allele: str
    gene_symbol: Optional[str] = None
    transcript_id: Optional[str] = None
    name: Optional[str] = None
    species: str = "homo_sapiens"


@dataclass
class CodingAnnotation:
    """Coding consequence annotation for a variant."""
    consequence: ConsequenceType
    hgvs_c: str = ""                   # e.g., "c.910C>T"
    hgvs_p: str = ""                   # e.g., "p.Arg304Ter"
    reference_codon: str = ""
    alternate_codon: str = ""
    reference_aa: str = ""
    alternate_aa: str = ""
    codon_index: int = -1
    codon_position: int = -1
    splice_proximity: Optional[str] = None  # "near_exon_start", "near_exon_end"
    message: str = ""                  # human-readable summary


@dataclass
class ReferenceValidation:
    """Result of validating a ref allele against the genome."""
    is_valid: bool
    expected_ref_genomic: str          # what Ensembl says
    provided_ref_genomic: str          # what the user gave
    expected_ref_transcript: str       # on transcript strand
    provided_ref_transcript: str       # user's ref on transcript strand
    message: str = ""


@dataclass
class NormalizedVariant:
    """Fully annotated variant after normalization + annotation.

    This is the central v2 data structure. It contains everything needed
    to assess editing feasibility: the raw input, transcript coordinate,
    coding consequence, reference validation, and a bridge to v1's
    Mutation dataclass.

    The v1_mutation field allows direct use with v1's StrategyEnumerator
    and StrategyScorer.
    """
    input: GenomicVariantInput
    transcript: TranscriptInfo
    transcript_coord: TranscriptCoordinate
    coding: CodingAnnotation
    ref_validation: ReferenceValidation
    v1_mutation: Any = None            # mosaic.mutation_classifier.Mutation
    local_sequence: str = ""           # genomic sequence +-flank around variant
    local_seq_edit_index: int = -1     # 0-based index of edit in local_sequence


# ═══════════════════════════════════════════════════════════════════════════
# Feasibility Models
# ═══════════════════════════════════════════════════════════════════════════

@dataclass
class GuideCandidate:
    """A candidate sgRNA identified by PAM scanning.

    Attributes
    ----------
    sequence_20mer : str
        20-nt protospacer sequence (5' to 3', PAM-distal to PAM-proximal).
    pam_sequence : str
        The PAM sequence found (e.g., "AGG").
    strand : str
        "+" for sense strand guide, "-" for antisense.
    cut_position : int
        1-based genomic coordinate of predicted cut site.
        For SpCas9: 3 bp upstream of PAM on protospacer strand.
    distance_to_edit : int
        Signed distance from cut to target edit (bp).
        Negative means cut is upstream of edit.
    gc_content : float
        GC fraction of the 20-mer protospacer.
    position_in_window : int
        Position of the target base within the protospacer (1-20).
        Only meaningful for base editing; -1 otherwise.
    has_poly_t : bool
        True if protospacer contains >=4 consecutive T's (Pol III
        terminator signal, reduces expression from U6 promoter).
    score : float
        Composite ranking score (higher = better).
    """
    sequence_20mer: str = ""
    pam_sequence: str = ""
    strand: str = "+"
    cut_position: int = 0
    distance_to_edit: int = 0
    gc_content: float = 0.0
    position_in_window: int = -1
    has_poly_t: bool = False
    score: float = 0.0


@dataclass
class BaseEditingFeasibility:
    """Result of base editing feasibility assessment for one variant.

    Biology
    -------
    ABE converts A>G (correction of G>A patient mutations).
    CBE converts C>T (correction of T>C patient mutations).
    Editing windows (1-indexed in 20-mer, PAM-distal = 1):
      ABE: positions 4-7 (Gaudelli et al., 2017)
      CBE: positions 4-8 (Komor et al., 2016)

    Bystanders are same-type bases within the editing window that may
    be edited unintentionally along with the target.
    """
    label: FeasibilityLabel = FeasibilityLabel.NOT_FEASIBLE
    editor_type: Optional[str] = None  # "ABE" or "CBE"
    best_guide: Optional[GuideCandidate] = None
    target_position_in_window: int = -1
    bystander_count: int = 0
    bystander_positions: List[int] = field(default_factory=list)
    bystander_consequences: List[ConsequenceType] = field(default_factory=list)
    compatible_nucleases: List[str] = field(default_factory=list)
    score: float = 0.0
    warnings: List[str] = field(default_factory=list)
    rejection_reason: str = ""
    metadata: Dict[str, Any] = field(default_factory=dict)


@dataclass
class PrimeEditingFeasibility:
    """Result of prime editing feasibility assessment.

    Biology
    -------
    PE uses a pegRNA with:
    - Spacer (20 nt) for Cas9 targeting
    - PBS (primer binding site, 10-17 nt) for RT priming
    - RT template (10-30 nt) encoding the desired edit

    PE3 adds a second nicking guide 40-100 bp away to improve efficiency.

    References
    ----------
    Anzalone et al., Nature, 2019 (PE design rules)
    Nelson et al., Nature Biotechnology, 2022 (PrimeDesign)
    """
    label: FeasibilityLabel = FeasibilityLabel.NOT_FEASIBLE
    best_guide: Optional[GuideCandidate] = None
    pbs_length: int = 13
    rt_template_length: int = 15
    rt_template_sequence: str = ""
    pe3_nick_guide: Optional[GuideCandidate] = None
    pe3_nick_distance: int = 0
    edit_type: str = ""                # "substitution", "insertion", "deletion"
    score: float = 0.0
    warnings: List[str] = field(default_factory=list)
    rejection_reason: str = ""
    metadata: Dict[str, Any] = field(default_factory=dict)


@dataclass
class HDRFeasibility:
    """Result of HDR feasibility assessment.

    Biology
    -------
    HDR requires a DSB near the mutation and a donor template with
    homology arms flanking the desired correction. Key factors:
    - Cut-to-edit distance: <10 bp optimal, <100 bp feasible
    - Donor type: ssODN for close edits, cssDNA for moderate, lssDNA for far
    - Gene conversion tract must reach from cut to edit site

    References
    ----------
    Paquet et al., Nature, 2016 (cut-to-edit distance)
    Iyer et al., CRISPR Journal, 2022 (cssDNA optimization)
    """
    label: FeasibilityLabel = FeasibilityLabel.NOT_FEASIBLE
    best_guide: Optional[GuideCandidate] = None
    cut_to_edit_distance: int = 0
    recommended_donor_type: str = "cssDNA"
    donor_length_estimate: int = 0
    homology_arm_length: int = 300
    conversion_probability: float = 0.0
    pam_disruption_possible: bool = False
    score: float = 0.0
    warnings: List[str] = field(default_factory=list)
    rejection_reason: str = ""
    metadata: Dict[str, Any] = field(default_factory=dict)


@dataclass
class FeasibilityBundle:
    """All feasibility results for one variant, across all modalities.

    This is what the strategy generator consumes: for each variant in the
    patient, which modalities are feasible and how do they compare?
    """
    variant: NormalizedVariant
    mutation_index: int                # index in the input variant list
    base_editing_results: List[BaseEditingFeasibility] = field(
        default_factory=list
    )
    prime_editing_result: Optional[PrimeEditingFeasibility] = None
    hdr_result: Optional[HDRFeasibility] = None

    def best_base_editing_result(self) -> Optional[BaseEditingFeasibility]:
        """Return the highest-scoring base editing result, or None."""
        feasible = [
            r for r in self.base_editing_results
            if r.label != FeasibilityLabel.NOT_FEASIBLE
        ]
        if not feasible:
            return None
        return max(feasible, key=lambda r: r.score)

    def best_modality(self) -> Optional[EditModality]:
        """Return the recommended modality for this variant."""
        candidates: List[Tuple[float, EditModality]] = []

        be = self.best_base_editing_result()
        if be and be.label != FeasibilityLabel.NOT_FEASIBLE:
            mod = EditModality.ABE if be.editor_type == "ABE" else EditModality.CBE
            candidates.append((be.score, mod))

        pe = self.prime_editing_result
        if pe and pe.label != FeasibilityLabel.NOT_FEASIBLE:
            candidates.append((pe.score, EditModality.PE))

        hdr = self.hdr_result
        if hdr and hdr.label != FeasibilityLabel.NOT_FEASIBLE:
            mod = EditModality.HDR_CSSDNA  # default; refined by donor_type
            if hdr.recommended_donor_type == "ssODN":
                mod = EditModality.HDR_SSODN
            elif hdr.recommended_donor_type == "lssDNA":
                mod = EditModality.HDR_LSSDNA
            candidates.append((hdr.score, mod))

        if not candidates:
            return None
        candidates.sort(reverse=True, key=lambda x: x[0])
        return candidates[0][1]


# ═══════════════════════════════════════════════════════════════════════════
# Strategy Models (v2 enhanced, linking to v1)
# ═══════════════════════════════════════════════════════════════════════════

@dataclass
class StrategyStep:
    """One editing step within a multi-step strategy."""
    modality: EditModality
    target_mutation_index: int         # index into the variant list
    donor_required: bool = False
    editor_name: str = ""              # e.g., "ABE8e", "PE2", "SpCas9"


@dataclass
class Strategy:
    """A complete editing strategy with v2 annotation.

    Links to v1's EditingStrategy via the v1_strategy field, while adding
    consequence-aware scoring and step-level feasibility information.
    """
    name: str
    steps: List[StrategyStep] = field(default_factory=list)
    num_dsbs: int = 0
    simultaneous_dsbs: bool = False
    num_rounds: int = 1
    num_distinct_guides: int = 1
    num_distinct_proteins: int = 1
    num_donors: int = 0
    requires_selection: bool = False
    screening_clones: int = 12
    p53_active: bool = True
    rearrangement_risk: RiskLevel = RiskLevel.LOW
    evidence_tier: EvidenceTier = EvidenceTier.B
    feasibility_results: List[Any] = field(default_factory=list)
    bystander_severity: float = 0.0    # 0-1, 0 = no bystander risk
    donor_feasibility_score: float = 1.0
    modality_prior_score: float = 0.5
    estimated_duration_weeks: int = 4
    included_reasons: List[str] = field(default_factory=list)
    penalties: List[str] = field(default_factory=list)
    rejection_reasons: List[str] = field(default_factory=list)

    @property
    def is_rejected(self) -> bool:
        return len(self.rejection_reasons) > 0


@dataclass
class ScoredStrategy:
    """A strategy with multi-objective scores and final ranking."""
    strategy: Strategy
    safety_score: float = 0.0
    feasibility_score: float = 0.0
    complexity_score: float = 0.0
    risk_score: float = 0.0
    confidence_score: float = 0.0
    consequence_penalty: float = 0.0
    consequence_bonus: float = 0.0
    overall_score: float = 0.0
    rank: int = 0
    annotation_notes: List[str] = field(default_factory=list)

    @property
    def strategy_name(self) -> str:
        return self.strategy.name

    @property
    def confidence(self) -> str:
        return self.strategy.evidence_tier.value


# ═══════════════════════════════════════════════════════════════════════════
# Pipeline Result
# ═══════════════════════════════════════════════════════════════════════════

@dataclass
class PipelineResult:
    """Full output of the v2 strategy pipeline."""
    transcript: TranscriptInfo
    variants: List[NormalizedVariant] = field(default_factory=list)
    bundles: List[FeasibilityBundle] = field(default_factory=list)
    strategies: List[ScoredStrategy] = field(default_factory=list)
    rejected_strategies: List[Strategy] = field(default_factory=list)
    metadata: Dict[str, Any] = field(default_factory=dict)
    warnings: List[str] = field(default_factory=list)

    @property
    def top_strategy(self) -> Optional[ScoredStrategy]:
        if self.strategies:
            return self.strategies[0]
        return None


# ═══════════════════════════════════════════════════════════════════════════
# Benchmark Models
# ═══════════════════════════════════════════════════════════════════════════

@dataclass
class BenchmarkCase:
    """One curated case for benchmarking.

    Each case has:
    - Real genomic coordinates from ClinVar / published studies
    - Expected modality determined by biological reasoning
    - Tiered truth labels (preferred / acceptable / reject)
    """
    case_id: str
    gene_symbol: str
    variants: List[GenomicVariantInput]
    cell_type: str = "iPSC"
    nuclease: str = "SpCas9"
    category: str = ""                 # e.g., "clean_base_editable"
    truth_label: Dict[str, List[str]] = field(default_factory=dict)
    disease_context: str = ""
    source_pmid: str = ""
    rationale: List[str] = field(default_factory=list)
    notes: str = ""


@dataclass
class BenchmarkResult:
    """Result of running one benchmark case through the pipeline."""
    case: BenchmarkCase
    pipeline_result: Optional[PipelineResult] = None
    top_strategy: str = ""
    top3_strategies: List[str] = field(default_factory=list)
    top1_correct: bool = False
    top3_correct: bool = False
    rejected_correctly: bool = False
    pipeline_error: str = ""


@dataclass
class BenchmarkSummary:
    """Aggregate metrics across all benchmark cases."""
    n_cases: int = 0
    top1_accuracy: float = 0.0
    top3_accuracy: float = 0.0
    rejection_accuracy: float = 0.0
    consequence_shift_fraction: float = 0.0
    case_results: List[BenchmarkResult] = field(default_factory=list)

    def to_dict(self) -> Dict[str, Any]:
        return {
            "n_cases": self.n_cases,
            "top1_accuracy": self.top1_accuracy,
            "top3_accuracy": self.top3_accuracy,
            "rejection_accuracy": self.rejection_accuracy,
            "consequence_shift_fraction": self.consequence_shift_fraction,
            "case_results": [
                {
                    "case_id": r.case.case_id,
                    "top_strategy": r.top_strategy,
                    "top3_strategies": r.top3_strategies,
                    "top1_correct": r.top1_correct,
                    "top3_correct": r.top3_correct,
                    "error": r.pipeline_error,
                }
                for r in self.case_results
            ],
        }


# ═══════════════════════════════════════════════════════════════════════════
# Self-test
# ═══════════════════════════════════════════════════════════════════════════

if __name__ == "__main__":
    # Quick smoke test: all dataclasses instantiate with defaults
    print("Testing core/models.py instantiation...")

    t = TranscriptInfo(
        transcript_id="ENST00000358273",
        gene_symbol="NF1",
        gene_id="ENSG00000196712",
        chromosome="17",
        start=31094927,
        end=31377677,
        strand=-1,
        biotype="protein_coding",
        is_canonical=True,
        exons=[
            ExonRecord("ENSE001", 1, 31094927, 31095068, -1, "17"),
        ],
    )
    assert t.n_exons == 1
    assert t.span_bp == 282751

    v = GenomicVariantInput(
        chromosome="17",
        position=31232193,
        ref_allele="C",
        alt_allele="T",
        gene_symbol="NF1",
        name="c.910C>T",
    )
    assert v.ref_allele == "C"

    gc = GuideCandidate(sequence_20mer="ATCGATCGATCGATCGATCG")
    assert len(gc.sequence_20mer) == 20

    be = BaseEditingFeasibility()
    assert be.label == FeasibilityLabel.NOT_FEASIBLE

    pe = PrimeEditingFeasibility()
    hdr = HDRFeasibility()

    s = Strategy(name="test")
    assert not s.is_rejected

    bs = BenchmarkSummary(n_cases=10, top1_accuracy=0.8)
    d = bs.to_dict()
    assert d["n_cases"] == 10

    print("All models instantiate correctly.")
    print(f"ConsequenceType values: {[c.value for c in ConsequenceType]}")
    print(f"EditModality values: {[m.value for m in EditModality]}")
    print("PASS")
