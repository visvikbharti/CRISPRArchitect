"""
HGVS Notation Parser for CRISPRArchitect v3
=============================================

Parses clinical variant notation (HGVS coding DNA format) into
GenomicVariantInput objects that the CRISPRArchitect pipeline can process.

Supported formats
-----------------
- NM_000267.3:c.910C>T      (coding DNA substitution)
- NM_000267.3:c.1234del      (single-base deletion)
- NM_000267.3:c.1234_1236del (multi-base deletion)
- NM_000267.3:c.1234insATG   (insertion)
- NM_000267.3:c.1234delinsATG (delins / indel)
- Gene-prefixed: NF1:c.910C>T (resolved via Ensembl)

The parser converts transcript-relative coordinates (c. notation) to
genomic coordinates via Ensembl REST API.

Limitations
-----------
- Does not handle g. (genomic) or p. (protein) notation directly
- Does not handle complex HGVS (duplications, inversions, etc.)
- Requires network access for Ensembl coordinate mapping

References
----------
den Dunnen et al., Human Mutation, 2016 (HGVS nomenclature v15.11)

Python 3.9 compatible.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

from core.models import GenomicVariantInput

try:
    from core.sequence.fetcher import EnsemblClient
except ImportError:
    EnsemblClient = None  # type: ignore


# ── HGVS parsing regex patterns ──────────────────────────────────────

# NM_000267.3:c.910C>T
_HGVS_SUBSTITUTION = re.compile(
    r'^(?:(?P<gene>[A-Z][A-Z0-9_-]*):)?'   # optional gene symbol prefix
    r'(?P<transcript>[NX][MR]_\d+(?:\.\d+)?)?'  # transcript accession
    r':?c\.(?P<position>-?\d+(?:[+-]\d+)?)'  # c. position (may include intron offset)
    r'(?P<ref>[ACGT])>(?P<alt>[ACGT])$'      # substitution
)

# NM_000267.3:c.1234del or c.1234_1236del
_HGVS_DELETION = re.compile(
    r'^(?:(?P<gene>[A-Z][A-Z0-9_-]*):)?'
    r'(?P<transcript>[NX][MR]_\d+(?:\.\d+)?)?'
    r':?c\.(?P<start>-?\d+)(?:_(?P<end>-?\d+))?'
    r'del(?P<deleted>[ACGT]*)$'
)

# NM_000267.3:c.1234insATG
_HGVS_INSERTION = re.compile(
    r'^(?:(?P<gene>[A-Z][A-Z0-9_-]*):)?'
    r'(?P<transcript>[NX][MR]_\d+(?:\.\d+)?)?'
    r':?c\.(?P<start>-?\d+)_(?P<end>-?\d+)'
    r'ins(?P<inserted>[ACGT]+)$'
)

# NM_000267.3:c.1234delinsATG
_HGVS_DELINS = re.compile(
    r'^(?:(?P<gene>[A-Z][A-Z0-9_-]*):)?'
    r'(?P<transcript>[NX][MR]_\d+(?:\.\d+)?)?'
    r':?c\.(?P<start>-?\d+)(?:_(?P<end>-?\d+))?'
    r'delins(?P<inserted>[ACGT]+)$'
)


@dataclass
class ParsedHGVS:
    """Parsed HGVS notation components."""
    gene_symbol: Optional[str] = None
    transcript_id: Optional[str] = None
    cds_position: int = 0
    cds_end_position: Optional[int] = None
    ref_allele: str = ""
    alt_allele: str = ""
    variant_type: str = ""  # "substitution", "deletion", "insertion", "delins"
    raw_notation: str = ""
    intron_offset: int = 0  # for intronic variants (e.g., c.1234+5)


class HGVSParser:
    """Parse HGVS c. notation into CRISPRArchitect variant inputs.

    Usage
    -----
    >>> parser = HGVSParser()
    >>> parsed = parser.parse("NM_000267.3:c.910C>T")
    >>> parsed.variant_type
    'substitution'
    >>> parsed.cds_position
    910

    For genomic coordinate resolution (requires Ensembl):
    >>> variant = parser.to_genomic_variant("NF1", "NM_000267.3:c.910C>T")
    """

    def parse(self, notation: str) -> ParsedHGVS:
        """Parse an HGVS c. notation string.

        Parameters
        ----------
        notation : str
            HGVS coding DNA notation (e.g., "NM_000267.3:c.910C>T").

        Returns
        -------
        ParsedHGVS
            Parsed components.

        Raises
        ------
        ValueError
            If the notation cannot be parsed.
        """
        notation = notation.strip()

        # Try substitution
        m = _HGVS_SUBSTITUTION.match(notation)
        if m:
            pos_str = m.group("position")
            pos, intron_offset = self._parse_position(pos_str)
            return ParsedHGVS(
                gene_symbol=m.group("gene"),
                transcript_id=m.group("transcript"),
                cds_position=pos,
                ref_allele=m.group("ref"),
                alt_allele=m.group("alt"),
                variant_type="substitution",
                raw_notation=notation,
                intron_offset=intron_offset,
            )

        # Try deletion
        m = _HGVS_DELETION.match(notation)
        if m:
            start, _ = self._parse_position(m.group("start"))
            end_str = m.group("end")
            end = self._parse_position(end_str)[0] if end_str else start
            deleted = m.group("deleted") or ""
            ref = deleted if deleted else "-"
            return ParsedHGVS(
                gene_symbol=m.group("gene"),
                transcript_id=m.group("transcript"),
                cds_position=start,
                cds_end_position=end if end != start else None,
                ref_allele=ref,
                alt_allele="-",
                variant_type="deletion",
                raw_notation=notation,
            )

        # Try insertion
        m = _HGVS_INSERTION.match(notation)
        if m:
            start, _ = self._parse_position(m.group("start"))
            return ParsedHGVS(
                gene_symbol=m.group("gene"),
                transcript_id=m.group("transcript"),
                cds_position=start,
                cds_end_position=self._parse_position(m.group("end"))[0],
                ref_allele="-",
                alt_allele=m.group("inserted"),
                variant_type="insertion",
                raw_notation=notation,
            )

        # Try delins
        m = _HGVS_DELINS.match(notation)
        if m:
            start, _ = self._parse_position(m.group("start"))
            end_str = m.group("end")
            end = self._parse_position(end_str)[0] if end_str else start
            return ParsedHGVS(
                gene_symbol=m.group("gene"),
                transcript_id=m.group("transcript"),
                cds_position=start,
                cds_end_position=end if end != start else None,
                ref_allele="",  # needs sequence lookup
                alt_allele=m.group("inserted"),
                variant_type="delins",
                raw_notation=notation,
            )

        raise ValueError(
            f"Cannot parse HGVS notation: '{notation}'. "
            "Supported formats: NM_xxx:c.NNNRef>Alt, c.NNNdel, c.NNN_NNNinsXXX"
        )

    def to_genomic_variant(
        self,
        gene_symbol: str,
        notation: str,
        transcript_id: Optional[str] = None,
    ) -> GenomicVariantInput:
        """Parse HGVS and resolve to genomic coordinates via Ensembl.

        Parameters
        ----------
        gene_symbol : str
            HGNC gene symbol (e.g., "NF1").
        notation : str
            HGVS c. notation.
        transcript_id : str, optional
            Override transcript ID. If None, uses the one from notation
            or the canonical transcript.

        Returns
        -------
        GenomicVariantInput
            Ready for the CRISPRArchitect pipeline.
        """
        parsed = self.parse(notation)

        # Use transcript from notation or parameter
        tx_id = transcript_id or parsed.transcript_id
        gene = gene_symbol or parsed.gene_symbol

        if not gene:
            raise ValueError("Gene symbol required for genomic coordinate resolution.")

        if EnsemblClient is None:
            raise ImportError("EnsemblClient not available for coordinate resolution.")

        client = EnsemblClient()

        # Resolve transcript
        if tx_id:
            tx = client.fetch_transcript_record(tx_id)
        else:
            tx = client.choose_transcript(gene)

        # Map CDS position to genomic coordinate
        genomic_pos = self._cds_to_genomic(
            tx, parsed.cds_position, parsed.intron_offset
        )

        # Determine chromosome
        chromosome = tx.chromosome

        return GenomicVariantInput(
            chromosome=chromosome,
            position=genomic_pos,
            ref_allele=parsed.ref_allele,
            alt_allele=parsed.alt_allele,
            gene_symbol=gene,
            transcript_id=tx.transcript_id,
            name=parsed.raw_notation,
        )

    def parse_batch(self, notations: List[str], gene_symbol: str) -> List[ParsedHGVS]:
        """Parse a batch of HGVS notations.

        Returns successfully parsed results; logs failures.
        """
        results = []
        for notation in notations:
            try:
                parsed = self.parse(notation)
                if not parsed.gene_symbol:
                    parsed.gene_symbol = gene_symbol
                results.append(parsed)
            except ValueError:
                pass  # skip unparseable
        return results

    @staticmethod
    def _parse_position(pos_str: str) -> Tuple[int, int]:
        """Parse a CDS position string, handling intron offsets.

        Examples: "910" -> (910, 0), "1234+5" -> (1234, 5), "1234-3" -> (1234, -3)
        """
        if "+" in pos_str:
            parts = pos_str.split("+")
            return int(parts[0]), int(parts[1])
        elif "-" in pos_str and pos_str.index("-") > 0:
            # Negative offset, not negative CDS position
            parts = pos_str.rsplit("-", 1)
            if parts[0]:
                return int(parts[0]), -int(parts[1])
        return int(pos_str), 0

    @staticmethod
    def _cds_to_genomic(tx, cds_position: int, intron_offset: int = 0) -> int:
        """Map a CDS position to a genomic coordinate using transcript exons.

        Accounts for UTR regions by clipping exons to CDS boundaries
        (Translation start/end from Ensembl). For intronic offsets (e.g.,
        c.92+5), the offset is applied from the nearest exon boundary
        in the correct genomic direction for the strand.

        For forward-strand transcripts: walk exons 5'->3'.
        For reverse-strand: walk exons 3'->5' (reversed order).
        """
        exons = sorted(tx.exons, key=lambda e: e.start)
        if tx.strand == -1:
            exons = list(reversed(exons))

        cds_start = getattr(tx, 'cds_start', None)
        cds_end = getattr(tx, 'cds_end', None)

        cds_bp_remaining = cds_position
        for exon in exons:
            # Clip exon to CDS boundaries (exclude UTR)
            ex_start = exon.start
            ex_end = exon.end
            if cds_start is not None:
                ex_start = max(ex_start, cds_start)
            if cds_end is not None:
                ex_end = min(ex_end, cds_end)

            cds_length = ex_end - ex_start + 1
            if cds_length <= 0:
                continue  # UTR-only exon, skip

            if cds_bp_remaining <= cds_length:
                if tx.strand == 1:
                    genomic = ex_start + cds_bp_remaining - 1
                else:
                    genomic = ex_end - cds_bp_remaining + 1

                if intron_offset != 0:
                    # Intronic offset: apply from the exon boundary
                    # in the direction away from the exon
                    if intron_offset > 0:
                        # Past the 3' end of this exon (in transcript sense)
                        if tx.strand == 1:
                            genomic = ex_end + intron_offset
                        else:
                            genomic = ex_start - intron_offset
                    else:
                        # Before the 5' start of this exon (in transcript sense)
                        if tx.strand == 1:
                            genomic = ex_start + intron_offset
                        else:
                            genomic = ex_end - intron_offset

                return genomic
            cds_bp_remaining -= cds_length

        raise ValueError(
            f"CDS position {cds_position} exceeds transcript length "
            f"({tx.transcript_id})."
        )


# ═══════════════════════════════════════════════════════════════════════
# ClinVar Batch Processor
# ═══════════════════════════════════════════════════════════════════════

@dataclass
class ClinVarVariant:
    """A variant from ClinVar."""
    clinvar_id: str = ""
    gene_symbol: str = ""
    chromosome: str = ""
    position: int = 0
    ref_allele: str = ""
    alt_allele: str = ""
    clinical_significance: str = ""
    hgvs_c: str = ""
    hgvs_p: str = ""
    review_status: str = ""
    condition: str = ""


class ClinVarProcessor:
    """Process ClinVar variants for batch analysis.

    Reads ClinVar TSV (variant_summary.txt) or VCF format and converts
    to GenomicVariantInput objects for CRISPRArchitect pipeline.

    Usage
    -----
    >>> processor = ClinVarProcessor()
    >>> variants = processor.load_tsv("clinvar_variant_summary.txt",
    ...     gene_filter=["NF1", "BRCA2", "DMD"],
    ...     significance_filter=["Pathogenic", "Likely pathogenic"],
    ... )
    >>> inputs = processor.to_pipeline_inputs(variants)
    """

    def load_tsv(
        self,
        filepath: str,
        gene_filter: Optional[List[str]] = None,
        significance_filter: Optional[List[str]] = None,
        assembly: str = "GRCh38",
        max_variants: int = 10000,
    ) -> List[ClinVarVariant]:
        """Load variants from ClinVar variant_summary.txt.

        Parameters
        ----------
        filepath : str
            Path to the ClinVar TSV file.
        gene_filter : list of str, optional
            Only include variants in these genes.
        significance_filter : list of str, optional
            Only include variants with these clinical significance values.
            Default: ["Pathogenic", "Likely pathogenic"]
        assembly : str
            Genome assembly (default "GRCh38").
        max_variants : int
            Maximum number of variants to return.

        Returns
        -------
        List[ClinVarVariant]
            Parsed variants.
        """
        if significance_filter is None:
            significance_filter = ["Pathogenic", "Likely pathogenic"]

        variants = []
        with open(filepath, "r") as f:
            header = None
            for line in f:
                if line.startswith("#"):
                    continue
                fields = line.strip().split("\t")
                if header is None:
                    header = {name: idx for idx, name in enumerate(fields)}
                    continue

                try:
                    gene = fields[header.get("GeneSymbol", header.get("Gene", 0))]
                    if gene_filter and gene not in gene_filter:
                        continue

                    sig = fields[header.get("ClinicalSignificance", header.get("ClinSig", 0))]
                    if not any(s.lower() in sig.lower() for s in significance_filter):
                        continue

                    asm = fields[header.get("Assembly", 0)]
                    if asm != assembly:
                        continue

                    chrom = fields[header.get("Chromosome", header.get("Chr", 0))]
                    pos_str = fields[header.get("PositionVCF", header.get("Start", 0))]
                    ref = fields[header.get("ReferenceAlleleVCF", header.get("ReferenceAllele", 0))]
                    alt = fields[header.get("AlternateAlleleVCF", header.get("AlternateAllele", 0))]

                    if not pos_str or pos_str == "-1" or not ref or not alt:
                        continue

                    cv = ClinVarVariant(
                        clinvar_id=fields[header.get("VariationID", header.get("#AlleleID", 0))],
                        gene_symbol=gene,
                        chromosome=chrom,
                        position=int(pos_str),
                        ref_allele=ref,
                        alt_allele=alt,
                        clinical_significance=sig,
                        hgvs_c=fields[header.get("HGVS(c)", "")] if "HGVS(c)" in header else "",
                        review_status=fields[header.get("ReviewStatus", 0)] if "ReviewStatus" in header else "",
                        condition=fields[header.get("PhenotypeIDS", header.get("Condition", 0))] if "PhenotypeIDS" in header else "",
                    )
                    variants.append(cv)

                    if len(variants) >= max_variants:
                        break

                except (IndexError, ValueError, KeyError):
                    continue

        return variants

    def load_vcf(
        self,
        filepath: str,
        gene_filter: Optional[List[str]] = None,
        max_variants: int = 10000,
    ) -> List[ClinVarVariant]:
        """Load variants from ClinVar VCF file.

        Parameters
        ----------
        filepath : str
            Path to ClinVar VCF file.
        gene_filter : list of str, optional
            Only include variants in these genes.
        max_variants : int
            Maximum variants to return.
        """
        variants = []
        with open(filepath, "r") as f:
            for line in f:
                if line.startswith("#"):
                    continue
                fields = line.strip().split("\t")
                if len(fields) < 8:
                    continue

                chrom = fields[0].replace("chr", "")
                pos = int(fields[1])
                ref = fields[3]
                alt = fields[4]
                info = fields[7]

                # Parse INFO field
                info_dict = {}
                for item in info.split(";"):
                    if "=" in item:
                        k, v = item.split("=", 1)
                        info_dict[k] = v

                gene = info_dict.get("GENEINFO", "").split(":")[0] if "GENEINFO" in info_dict else ""
                if gene_filter and gene not in gene_filter:
                    continue

                sig = info_dict.get("CLNSIG", "")
                if "pathogenic" not in sig.lower():
                    continue

                # Only SNVs and small indels
                if len(ref) > 50 or len(alt) > 50:
                    continue

                cv = ClinVarVariant(
                    clinvar_id=info_dict.get("ALLELEID", fields[2]),
                    gene_symbol=gene,
                    chromosome=chrom,
                    position=pos,
                    ref_allele=ref,
                    alt_allele=alt,
                    clinical_significance=sig,
                )
                variants.append(cv)

                if len(variants) >= max_variants:
                    break

        return variants

    def to_pipeline_inputs(
        self,
        clinvar_variants: List[ClinVarVariant],
    ) -> List[GenomicVariantInput]:
        """Convert ClinVar variants to GenomicVariantInput objects."""
        inputs = []
        for cv in clinvar_variants:
            inputs.append(GenomicVariantInput(
                chromosome=cv.chromosome,
                position=cv.position,
                ref_allele=cv.ref_allele,
                alt_allele=cv.alt_allele,
                gene_symbol=cv.gene_symbol,
                name=cv.hgvs_c or f"ClinVar:{cv.clinvar_id}",
            ))
        return inputs

    def filter_snvs(
        self,
        variants: List[ClinVarVariant],
    ) -> List[ClinVarVariant]:
        """Filter to only single-nucleotide variants (most BE/PE-amenable)."""
        return [
            v for v in variants
            if len(v.ref_allele) == 1 and len(v.alt_allele) == 1
            and v.ref_allele in "ACGT" and v.alt_allele in "ACGT"
        ]

    def classify_variants(
        self,
        variants: List[ClinVarVariant],
    ) -> Dict[str, List[ClinVarVariant]]:
        """Classify variants by editability potential.

        Returns
        -------
        Dict with keys: "abe_amenable", "cbe_amenable", "transversion",
        "indel", "complex"
        """
        transitions_abe = {("A", "G"), ("G", "A"), ("T", "C"), ("C", "T")}
        transitions_cbe = {("C", "T"), ("T", "C"), ("G", "A"), ("A", "G")}

        result = {
            "abe_amenable": [],
            "cbe_amenable": [],
            "transversion": [],
            "indel": [],
            "complex": [],
        }

        for v in variants:
            if len(v.ref_allele) != 1 or len(v.alt_allele) != 1:
                if len(v.ref_allele) > 50 or len(v.alt_allele) > 50:
                    result["complex"].append(v)
                else:
                    result["indel"].append(v)
                continue

            pair = (v.alt_allele, v.ref_allele)  # correction direction
            if pair in {("A", "G"), ("T", "C")}:
                result["abe_amenable"].append(v)
            elif pair in {("C", "T"), ("G", "A")}:
                result["cbe_amenable"].append(v)
            else:
                result["transversion"].append(v)

        return result
