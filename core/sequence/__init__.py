"""
Transcript-aware sequence mapping and variant annotation.

This sub-package provides the v2 sequence layer:

- TranscriptFetcher: Ensembl REST API client for transcript structures
- TranscriptMapper: genomic-to-transcript coordinate mapping
- ReferenceValidator: ref allele validation against the genome
- CodingAnnotator: consequence annotation and HGVS notation
- VariantNormalizer: end-to-end orchestrator (the main entry point)
"""

from core.sequence.fetcher import TranscriptFetcher
from core.sequence.transcript_mapper import TranscriptMapper
from core.sequence.reference_validator import ReferenceValidator
from core.sequence.coding_annotation import CodingAnnotator
from core.sequence.variant_normalizer import VariantNormalizer

__all__ = [
    "TranscriptFetcher",
    "TranscriptMapper",
    "ReferenceValidator",
    "CodingAnnotator",
    "VariantNormalizer",
]
