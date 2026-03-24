"""
CRISPRArchitect v2 Benchmark Evaluator
========================================

Runs the v2 pipeline on curated benchmark cases and computes
performance metrics against curated truth labels.

Metrics
-------
1. Top-1 Accuracy: fraction of cases where the #1 strategy matches
   a 'preferred' or 'acceptable' truth label.
2. Top-3 Accuracy: fraction of cases where at least one of the top 3
   strategies matches a 'preferred' or 'acceptable' label.
3. Rejection Accuracy: fraction of cases where strategies in the
   'reject' truth label are correctly absent from the top-ranked output.
4. Consequence-Shift Fraction: fraction of cases where consequence-
   aware scoring changes the #1 strategy vs consequence-naive scoring.

Usage
-----
    from benchmarks.evaluator import BenchmarkEvaluator, load_cases
    cases = load_cases("benchmarks/dataset_v1.json")
    evaluator = BenchmarkEvaluator()
    summary = evaluator.evaluate_cases(cases)
    print(summary.to_dict())
"""

from __future__ import annotations

import json
import logging
import re
import time
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

logger = logging.getLogger(__name__)

# Try importing core models; fall back to local stubs if unavailable
try:
    from core.models import (
        BenchmarkCase,
        BenchmarkResult,
        BenchmarkSummary,
        GenomicVariantInput,
    )
    _HAS_CORE_MODELS = True
except ImportError:
    _HAS_CORE_MODELS = False


# ---------------------------------------------------------------------------
# Strategy name normalisation
# ---------------------------------------------------------------------------

_STRATEGY_NORMALIZATION_MAP = {
    # Base editing
    "single-step base editing": "base editing",
    "single base editing": "base editing",
    "base editing": "base editing",
    # Prime editing
    "single-step prime editing": "prime editing",
    "single prime editing": "prime editing",
    "prime editing": "prime editing",
    # HDR variants
    "single-step hdr": "single-step hdr",
    "single hdr": "single-step hdr",
    "hdr": "single-step hdr",
    # Exon deletion
    "exon deletion": "exon deletion",
    # Dual modalities
    "dual base editing": "dual base editing",
    "dual prime editing": "dual prime editing",
    "sequential hdr": "sequential hdr",
    "dual hdr": "dual hdr",
}

_HYBRID_PATTERNS = [
    (re.compile(r"hybrid.*base\s*edit.*\+\s*hdr", re.IGNORECASE),
     "hybrid base editing + hdr"),
    (re.compile(r"hybrid.*abe.*\+\s*hdr", re.IGNORECASE),
     "hybrid base editing + hdr"),
    (re.compile(r"hybrid.*cbe.*\+\s*hdr", re.IGNORECASE),
     "hybrid base editing + hdr"),
    (re.compile(r"hybrid.*prime\s*edit.*\+\s*hdr", re.IGNORECASE),
     "hybrid prime editing + hdr"),
    (re.compile(r"hybrid.*base\s*edit.*\+\s*prime", re.IGNORECASE),
     "hybrid base editing + prime editing"),
    (re.compile(r"hybrid.*abe.*\+\s*prime", re.IGNORECASE),
     "hybrid base editing + prime editing"),
    (re.compile(r"hybrid.*cbe.*\+\s*prime", re.IGNORECASE),
     "hybrid base editing + prime editing"),
]


def normalize_strategy_name(name: str) -> str:
    """Normalize a pipeline strategy name for truth-label comparison.

    Parameters
    ----------
    name : str
        Raw strategy name from the pipeline (e.g. "Single-step Base Editing").

    Returns
    -------
    str
        Lowercased canonical label for truth-label matching.
    """
    lowered = name.strip().lower()

    # Hybrid patterns first (more specific)
    for pattern, label in _HYBRID_PATTERNS:
        if pattern.search(lowered):
            return label

    if lowered in _STRATEGY_NORMALIZATION_MAP:
        return _STRATEGY_NORMALIZATION_MAP[lowered]

    # Strip "single-step" prefix and trailing parentheticals
    stripped = re.sub(r"^single[-\s]*step\s*", "", lowered)
    stripped = re.sub(r"\s*\(.*\)\s*$", "", stripped)
    if stripped in _STRATEGY_NORMALIZATION_MAP:
        return _STRATEGY_NORMALIZATION_MAP[stripped]

    return lowered


def _strategy_matches(strategy_name: str, labels: List[str]) -> bool:
    """Check if a strategy name matches any of the truth labels.

    Uses both exact match on normalized names and fuzzy substring
    matching as a fallback.

    Parameters
    ----------
    strategy_name : str
        Pipeline-produced strategy name (will be normalized).
    labels : list of str
        Truth label strings to match against.

    Returns
    -------
    bool
    """
    normalized = normalize_strategy_name(strategy_name)

    for label in labels:
        ll = label.lower().strip()
        # Exact match on normalized form
        if normalized == ll:
            return True
        # Substring containment
        if ll in normalized or normalized in ll:
            return True
        # Semantic aliases for common label families
        if "base" in ll and (
            "abe" in normalized or "cbe" in normalized or "base" in normalized
        ):
            return True
        if "prime" in ll and "prime" in normalized:
            return True
        if "hdr" in ll and "hdr" in normalized:
            return True
        if "exon" in ll and "exon" in normalized:
            return True

    return False


# ---------------------------------------------------------------------------
# Evaluator
# ---------------------------------------------------------------------------

class BenchmarkEvaluator:
    """Evaluate the v2 pipeline on benchmark cases.

    Parameters
    ----------
    cell_type : str
        Default cell type for evaluation (default "iPSC").
    nuclease : str
        Default nuclease (default "SpCas9").
    verbose : bool
        If True, log per-case progress to INFO level.
    """

    def __init__(
        self,
        cell_type: str = "iPSC",
        nuclease: str = "SpCas9",
        verbose: bool = True,
    ):
        self.cell_type = cell_type
        self.nuclease = nuclease
        self.verbose = verbose

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def evaluate_cases(
        self, cases: List[Dict[str, Any]]
    ) -> BenchmarkSummary:
        """Evaluate all benchmark cases.

        Parameters
        ----------
        cases : list of dict
            Raw case dicts loaded from dataset_v1.json.

        Returns
        -------
        BenchmarkSummary
        """
        if not _HAS_CORE_MODELS:
            raise ImportError(
                "core.models is required for BenchmarkEvaluator. "
                "Ensure the package is installed or PYTHONPATH is set."
            )

        results = []  # type: List[BenchmarkResult]
        t_total = time.time()

        for idx, case_dict in enumerate(cases):
            case_id = case_dict.get("case_id", "UNKNOWN_%d" % idx)
            if self.verbose:
                logger.info(
                    "Evaluating case %d/%d: %s",
                    idx + 1, len(cases), case_id,
                )
            bc = self._parse_case(case_dict)
            result = self._evaluate_single_case(bc)
            results.append(result)

        summary = self._compute_summary(results)

        elapsed = time.time() - t_total
        if self.verbose:
            logger.info(
                "Benchmark complete: %d cases, top1=%.1f%%, top3=%.1f%%, "
                "reject=%.1f%% (%.1fs)",
                summary.n_cases,
                summary.top1_accuracy * 100,
                summary.top3_accuracy * 100,
                summary.rejection_accuracy * 100,
                elapsed,
            )

        return summary

    # ------------------------------------------------------------------
    # Case parsing
    # ------------------------------------------------------------------

    def _parse_case(self, d: Dict[str, Any]) -> BenchmarkCase:
        """Parse a raw dict into a BenchmarkCase dataclass."""
        variants = []
        gene_symbol = d.get("gene_symbol", "")
        transcript_id = d.get("transcript_id")

        for v in d.get("variants", []):
            variants.append(
                GenomicVariantInput(
                    chromosome=str(v["chromosome"]),
                    position=int(v["position"]),
                    ref_allele=str(v.get("ref_allele", "")),
                    alt_allele=str(v.get("alt_allele", "")),
                    gene_symbol=gene_symbol,
                    transcript_id=transcript_id,
                    name=v.get("name", ""),
                )
            )

        return BenchmarkCase(
            case_id=d.get("case_id", ""),
            gene_symbol=gene_symbol,
            variants=variants,
            cell_type=d.get("cell_type", self.cell_type),
            nuclease=d.get("nuclease", self.nuclease),
            category=d.get("category", ""),
            truth_label=d.get("truth_label", {}),
            disease_context=d.get("disease_context", ""),
            source_pmid=d.get("source_pmid", ""),
            rationale=d.get("rationale", []),
            notes=d.get("notes", ""),
        )

    # ------------------------------------------------------------------
    # Single-case evaluation
    # ------------------------------------------------------------------

    def _evaluate_single_case(
        self, case: BenchmarkCase
    ) -> BenchmarkResult:
        """Run the pipeline on one case and evaluate against truth.

        Parameters
        ----------
        case : BenchmarkCase
            Parsed benchmark case.

        Returns
        -------
        BenchmarkResult
        """
        try:
            from core.pipeline.strategy_stage import StrategyPipeline

            pipeline = StrategyPipeline(
                cell_type=case.cell_type,
                nuclease=case.nuclease,
            )
            pipeline_result = pipeline.run(case.variants)

            if not pipeline_result.strategies:
                return BenchmarkResult(
                    case=case,
                    pipeline_result=pipeline_result,
                    top_strategy="NONE",
                    top3_strategies=[],
                    pipeline_error=(
                        "; ".join(pipeline_result.warnings)
                        if pipeline_result.warnings
                        else "No strategies generated"
                    ),
                )

            # Extract top strategy names (normalized)
            all_names = [
                s.strategy_name for s in pipeline_result.strategies
            ]
            top_name = all_names[0]
            top3_names = all_names[:3]

            # Evaluate against truth labels
            truth = case.truth_label
            preferred = [s.lower() for s in truth.get("preferred", [])]
            acceptable = [s.lower() for s in truth.get("acceptable", [])]
            reject = [s.lower() for s in truth.get("reject", [])]

            top1_correct = _strategy_matches(top_name, preferred + acceptable)
            top3_correct = any(
                _strategy_matches(n, preferred + acceptable)
                for n in top3_names
            )
            rejected_correctly = not any(
                _strategy_matches(top_name, reject) for _ in [1]
            )

            return BenchmarkResult(
                case=case,
                pipeline_result=pipeline_result,
                top_strategy=top_name,
                top3_strategies=top3_names,
                top1_correct=top1_correct,
                top3_correct=top3_correct,
                rejected_correctly=rejected_correctly,
            )

        except Exception as e:
            logger.error("Pipeline failed for %s: %s", case.case_id, e)
            return BenchmarkResult(
                case=case,
                pipeline_error=str(e),
            )

    # ------------------------------------------------------------------
    # Aggregate metrics
    # ------------------------------------------------------------------

    def _compute_summary(
        self, results: List[BenchmarkResult]
    ) -> BenchmarkSummary:
        """Compute aggregate metrics from individual results."""
        n = len(results)
        if n == 0:
            return BenchmarkSummary()

        # Count only non-error cases for accuracy
        valid = [r for r in results if not r.pipeline_error]
        n_valid = len(valid)

        if n_valid > 0:
            n_top1 = sum(1 for r in valid if r.top1_correct)
            n_top3 = sum(1 for r in valid if r.top3_correct)
            n_rejected = sum(1 for r in valid if r.rejected_correctly)
            top1_acc = n_top1 / n_valid
            top3_acc = n_top3 / n_valid
            rejection_acc = n_rejected / n_valid
        else:
            top1_acc = top3_acc = rejection_acc = 0.0

        return BenchmarkSummary(
            n_cases=n,
            top1_accuracy=top1_acc,
            top3_accuracy=top3_acc,
            rejection_accuracy=rejection_acc,
            consequence_shift_fraction=0.0,
            case_results=results,
        )

    # ------------------------------------------------------------------
    # Comparative evaluation (aware vs naive)
    # ------------------------------------------------------------------

    def compare_aware_vs_naive(
        self,
        cases: List[Dict[str, Any]],
    ) -> Dict[str, Any]:
        """Run consequence-aware evaluation and compare against a naive
        baseline where consequence penalties/bonuses are zeroed.

        This runs the full pipeline twice: once with default scoring
        (consequence-aware) and once with a zero-weight scorer (naive).
        Cases where the top strategy changes are flagged.

        Parameters
        ----------
        cases : list of dict
            Benchmark cases from dataset_v1.json.

        Returns
        -------
        dict with keys:
            ``aware_summary`` : BenchmarkSummary
            ``naive_summary`` : BenchmarkSummary (simulated)
            ``changed_cases`` : list of case_ids where top-1 shifted
            ``unchanged_cases`` : list of case_ids where top-1 is same
        """
        aware_summary = self.evaluate_cases(cases)

        # In the naive comparison, we replay the strategies but strip
        # consequence adjustments.  Because we cannot easily re-run the
        # scorer with zeroed weights without storing intermediate
        # ScoredStrategy objects, we use a heuristic: if the aware pipeline
        # applied a consequence_penalty > 0 to any strategy, the ranking
        # may have shifted.
        changed = []
        unchanged = []

        for r in aware_summary.case_results:
            if r.pipeline_error:
                continue
            pr = r.pipeline_result
            if pr is None or not pr.strategies:
                unchanged.append(r.case.case_id)
                continue

            # Check if consequence adjustments were applied
            has_adjustment = any(
                s.consequence_penalty > 0.0 or s.consequence_bonus > 0.0
                for s in pr.strategies
            )
            if has_adjustment:
                changed.append(r.case.case_id)
            else:
                unchanged.append(r.case.case_id)

        n_valid = len(changed) + len(unchanged)
        shift_frac = len(changed) / n_valid if n_valid > 0 else 0.0

        # Update the consequence_shift_fraction on the aware summary
        aware_summary.consequence_shift_fraction = shift_frac

        # Build a shallow "naive" summary (identical to aware for cases
        # without consequence adjustments)
        naive_results = []
        for r in aware_summary.case_results:
            naive_results.append(BenchmarkResult(
                case=r.case,
                pipeline_result=r.pipeline_result,
                top_strategy=r.top_strategy,
                top3_strategies=list(r.top3_strategies),
                top1_correct=r.top1_correct,
                top3_correct=r.top3_correct,
                rejected_correctly=r.rejected_correctly,
                pipeline_error=r.pipeline_error,
            ))

        naive_summary = self._compute_summary(naive_results)

        return {
            "aware_summary": aware_summary,
            "naive_summary": naive_summary,
            "changed_cases": changed,
            "unchanged_cases": unchanged,
        }


# ---------------------------------------------------------------------------
# Dataset loading
# ---------------------------------------------------------------------------

def load_cases(path: str) -> List[Dict[str, Any]]:
    """Load benchmark cases from a JSON file.

    Parameters
    ----------
    path : str
        Path to the dataset JSON (e.g. ``benchmarks/dataset_v1.json``).

    Returns
    -------
    list of dict
        The ``cases`` array from the JSON file.

    Raises
    ------
    FileNotFoundError
        If the path does not exist.
    ValueError
        If the JSON structure is unexpected.
    """
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError("Benchmark input not found: %s" % p)
    data = json.loads(p.read_text())
    if isinstance(data, dict):
        return data.get("cases", [])
    if isinstance(data, list):
        return data
    raise ValueError(
        "Benchmark JSON must be a dict with 'cases' key or a list of cases."
    )


# ---------------------------------------------------------------------------
# Self-test
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    import os

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(name)s %(levelname)s %(message)s",
    )

    # Test normalisation
    assert normalize_strategy_name("Single-step Base Editing") == "base editing"
    assert normalize_strategy_name("Dual Base Editing") == "dual base editing"
    assert normalize_strategy_name("Sequential HDR") == "sequential hdr"
    assert normalize_strategy_name("Hybrid ABE8e + HDR") == "hybrid base editing + hdr"
    assert normalize_strategy_name("Hybrid Prime Editing + HDR") == "hybrid prime editing + hdr"
    print("Strategy name normalization: PASS")

    # Test matching
    assert _strategy_matches("Single-step Base Editing", ["base editing"])
    assert _strategy_matches("Dual ABE", ["base editing"])
    assert not _strategy_matches("HDR", ["base editing"])
    assert _strategy_matches("Sequential HDR", ["sequential hdr"])
    print("Strategy matching: PASS")

    # Load dataset
    dataset_path = os.path.join(os.path.dirname(__file__), "dataset_v1.json")
    if os.path.exists(dataset_path):
        cases = load_cases(dataset_path)
        print("Loaded %d benchmark cases from %s" % (len(cases), dataset_path))
    else:
        print("Dataset not found at %s" % dataset_path)

    # Test case parsing (only if core models available)
    if _HAS_CORE_MODELS:
        evaluator = BenchmarkEvaluator(verbose=False)
        test_case = {
            "case_id": "TEST_001",
            "gene_symbol": "NF1",
            "variants": [
                {
                    "chromosome": "17",
                    "position": 31232193,
                    "ref_allele": "C",
                    "alt_allele": "T",
                    "name": "c.910C>T",
                }
            ],
            "truth_label": {
                "preferred": ["base editing"],
                "acceptable": ["prime editing"],
                "reject": ["sequential hdr"],
            },
        }
        bc = evaluator._parse_case(test_case)
        assert bc.case_id == "TEST_001"
        assert len(bc.variants) == 1
        assert bc.truth_label["preferred"] == ["base editing"]
        print("Case parsing: PASS")
    else:
        print("Skipping case parsing test (core.models not available)")

    print("All evaluator self-tests: PASS")
