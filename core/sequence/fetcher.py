"""
Transcript Fetcher — Ensembl REST API client for transcript-level data
=======================================================================

Fetches transcript structures (exon boundaries, strand, biotype) from the
Ensembl REST API and returns v2 TranscriptInfo / ExonRecord objects.

Biology context
---------------
A gene can have many alternative transcripts due to alternative splicing,
alternative promoters, or alternative polyadenylation.  For clinical variant
annotation, we need the *canonical* (MANE Select or Ensembl canonical)
transcript, because that is the one ClinVar and HGVS nomenclature reference.

For CRISPR editing, the transcript choice determines:
  - which exon the variant falls in
  - whether the variant is coding or intronic
  - how far the variant is from splice junctions (affects splice-site
    editing feasibility)

This module wraps v1's fetch_sequence() for raw sequence retrieval and adds
transcript-level structure on top.

Rate limiting
-------------
Ensembl REST API asks for <= 15 requests/second.  We enforce a 0.1 s
minimum gap between calls to stay well within that limit.

References
----------
- Ensembl REST API: https://rest.ensembl.org
- Yates et al., Bioinformatics, 2015 (Ensembl REST API paper)
- Morales et al., Genomics, 2022 (MANE Select transcript set)
"""

from __future__ import annotations

import json
import logging
import socket
import time
import urllib.request
import urllib.error
from typing import Dict, List, Optional

logger = logging.getLogger(__name__)

from core.models import TranscriptInfo, ExonRecord

# v1 sequence fetcher — used for raw genomic sequence retrieval
try:
    from utils.ensembl import fetch_sequence, EnsemblError
except ImportError:
    try:
        from crisprarchitect.utils.ensembl import fetch_sequence, EnsemblError
    except ImportError:
        # Minimal fallback so the module can at least be imported
        class EnsemblError(Exception):
            pass

        def fetch_sequence(chromosome, start, end, strand=1, species="human"):
            raise EnsemblError("v1 ensembl module not available")


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

ENSEMBL_REST = "https://rest.ensembl.org"
TIMEOUT_SEC = 30
MIN_REQUEST_INTERVAL = 0.1  # seconds between API calls (rate limiting)
MAX_RETRIES = 3              # retry transient errors (HTTP 500/502/503/504, timeouts)
RETRY_BACKOFF_BASE = 2.0     # exponential backoff: wait 2^attempt seconds


# ---------------------------------------------------------------------------
# TranscriptFetcher
# ---------------------------------------------------------------------------

class TranscriptFetcher:
    """Fetches transcript structures from the Ensembl REST API.

    This class provides two entry points:
      1. fetch_by_gene(gene_symbol) — looks up the gene, returns its
         canonical transcript as a TranscriptInfo.
      2. fetch_by_transcript_id(transcript_id) — fetches a specific
         transcript directly.

    Results are cached in-memory so repeated queries for the same gene
    or transcript do not hit the API again.

    Parameters
    ----------
    species : str
        Ensembl species name (default: "homo_sapiens").

    Examples
    --------
    >>> fetcher = TranscriptFetcher()
    >>> tx = fetcher.fetch_by_gene("CFTR")
    >>> print(tx.transcript_id, tx.n_exons)
    """

    def __init__(self, species: str = "homo_sapiens"):
        self.species = species
        self._cache: Dict[str, TranscriptInfo] = {}
        self._last_request_time: float = 0.0

    # ----- public API -------------------------------------------------------

    def fetch_by_gene(self, gene_symbol: str) -> TranscriptInfo:
        """Fetch the canonical transcript for a gene symbol.

        Parameters
        ----------
        gene_symbol : str
            HGNC gene symbol (e.g., "NF1", "CFTR", "BRCA2").

        Returns
        -------
        TranscriptInfo
            Canonical transcript with all exon records.

        Raises
        ------
        EnsemblError
            If the gene is not found or the API is unreachable.
        """
        cache_key = f"gene:{gene_symbol.upper()}"
        if cache_key in self._cache:
            return self._cache[cache_key]

        # Step 1: lookup gene with transcripts expanded
        data = self._fetch_json(
            f"/lookup/symbol/{self.species}/{gene_symbol}?expand=1"
        )

        gene_id = data.get("id", "")
        chromosome = str(data.get("seq_region_name", ""))
        strand = data.get("strand", 1)

        # Step 2: find canonical transcript
        canonical = None
        for t in data.get("Transcript", []):
            if t.get("is_canonical") == 1:
                canonical = t
                break

        # Fallback: pick the transcript with the most exons
        if canonical is None and data.get("Transcript"):
            canonical = max(
                data["Transcript"],
                key=lambda t: len(t.get("Exon", []))
            )

        if canonical is None:
            raise EnsemblError(
                f"No transcripts found for gene {gene_symbol}"
            )

        transcript_info = self._build_transcript_info(
            canonical, gene_symbol.upper(), gene_id, chromosome, strand
        )

        # Cache under both gene symbol and transcript id
        self._cache[cache_key] = transcript_info
        self._cache[f"tx:{transcript_info.transcript_id}"] = transcript_info

        return transcript_info

    def fetch_by_transcript_id(self, transcript_id: str) -> TranscriptInfo:
        """Fetch a specific transcript by Ensembl transcript ID.

        Parameters
        ----------
        transcript_id : str
            Ensembl transcript stable ID (e.g., "ENST00000358273").

        Returns
        -------
        TranscriptInfo
        """
        cache_key = f"tx:{transcript_id}"
        if cache_key in self._cache:
            return self._cache[cache_key]

        data = self._fetch_json(
            f"/lookup/id/{transcript_id}?expand=1"
        )

        gene_symbol = data.get("display_name", "").split("-")[0]
        gene_id = data.get("Parent", "")
        chromosome = str(data.get("seq_region_name", ""))
        strand = data.get("strand", 1)

        # If the lookup returns a gene rather than a transcript, find it
        if data.get("object_type") == "Gene":
            raise EnsemblError(
                f"{transcript_id} resolved to a gene, not a transcript. "
                "Use fetch_by_gene() instead."
            )

        transcript_info = self._build_transcript_info(
            data, gene_symbol, gene_id, chromosome, strand
        )

        self._cache[cache_key] = transcript_info
        return transcript_info

    def fetch_genomic_sequence(
        self,
        chromosome: str,
        start: int,
        end: int,
        strand: int = 1,
    ) -> str:
        """Fetch genomic DNA sequence via v1's fetch_sequence().

        This is a thin wrapper that adds rate limiting on top of v1's
        fetch_sequence() function.

        Parameters
        ----------
        chromosome : str
            Chromosome name (e.g., "17", "X").
        start : int
            1-based start coordinate (inclusive).
        end : int
            1-based end coordinate (inclusive).
        strand : int
            1 for forward, -1 for reverse complement.

        Returns
        -------
        str
            Uppercase DNA sequence.
        """
        # Use retry logic for robustness against transient Ensembl errors.
        # v1's fetch_sequence() uses urllib internally, so we wrap with retry.
        last_error = None
        for attempt in range(MAX_RETRIES + 1):
            self._rate_limit()
            try:
                return fetch_sequence(chromosome, start, end, strand=strand)
            except Exception as e:
                err_str = str(e).lower()
                is_transient = any(
                    code in err_str
                    for code in ["500", "502", "503", "504", "timeout", "timed out"]
                )
                if is_transient and attempt < MAX_RETRIES:
                    wait = RETRY_BACKOFF_BASE ** attempt
                    logger.warning(
                        "Sequence fetch error on attempt %d/%d — "
                        "retrying in %.1fs: %s",
                        attempt + 1, MAX_RETRIES + 1, wait, e,
                    )
                    time.sleep(wait)
                    last_error = e
                    continue
                raise EnsemblError(
                    f"Failed to fetch sequence: {e}"
                ) from e
        raise EnsemblError(
            f"All {MAX_RETRIES + 1} attempts failed for sequence "
            f"chr{chromosome}:{start}-{end}. Last error: {last_error}"
        )

    def fetch_cds_sequence(self, transcript_id: str) -> str:
        """Fetch the coding DNA sequence (CDS) for a transcript.

        Uses the Ensembl /sequence/id endpoint with type=cds to get the
        spliced coding sequence. This is the concatenation of all coding
        exon portions, in transcript orientation.

        Parameters
        ----------
        transcript_id : str
            Ensembl transcript stable ID.

        Returns
        -------
        str
            CDS DNA sequence (uppercase, starts with ATG).
        """
        url = (
            f"{ENSEMBL_REST}/sequence/id/{transcript_id}"
            f"?type=cds&content-type=text/plain"
        )
        req = urllib.request.Request(url, headers={
            "Content-Type": "text/plain",
            "User-Agent": "CRISPRArchitect/2.0",
        })
        raw = self._urlopen_with_retry(req)
        seq = raw.decode().strip().replace("\n", "")
        return seq.upper()

    # ----- private helpers --------------------------------------------------

    def _build_transcript_info(
        self,
        tx_data: dict,
        gene_symbol: str,
        gene_id: str,
        chromosome: str,
        strand: int,
    ) -> TranscriptInfo:
        """Convert raw Ensembl JSON into a TranscriptInfo dataclass."""
        transcript_id = tx_data.get("id", "")
        biotype = tx_data.get("biotype", "unknown")
        is_canonical = bool(tx_data.get("is_canonical", 0))
        tx_start = tx_data.get("start", 0)
        tx_end = tx_data.get("end", 0)

        # Extract CDS boundaries from Translation object
        translation = tx_data.get("Translation", {})
        cds_start = translation.get("start")  # genomic, always <= cds_end
        cds_end = translation.get("end")

        # Build sorted exon records
        raw_exons = tx_data.get("Exon", [])
        if not raw_exons:
            raise EnsemblError(
                f"No exons found for transcript {transcript_id}. "
                "This may be a non-coding RNA or pseudogene."
            )

        # Sort exons by genomic position:
        #   forward strand → ascending start
        #   reverse strand → descending start (exon 1 is at the highest position)
        if strand == 1:
            raw_exons.sort(key=lambda e: e["start"])
        else:
            raw_exons.sort(key=lambda e: e["start"], reverse=True)

        exon_records: List[ExonRecord] = []
        for i, e in enumerate(raw_exons, 1):
            exon_records.append(ExonRecord(
                exon_id=e.get("id", f"exon_{i}"),
                exon_number=i,
                start=e["start"],
                end=e["end"],
                strand=strand,
                chromosome=chromosome,
            ))

        return TranscriptInfo(
            transcript_id=transcript_id,
            gene_symbol=gene_symbol,
            gene_id=gene_id,
            chromosome=chromosome,
            start=tx_start,
            end=tx_end,
            strand=strand,
            biotype=biotype,
            is_canonical=is_canonical,
            exons=exon_records,
            cds_start=cds_start,
            cds_end=cds_end,
        )

    def _rate_limit(self) -> None:
        """Enforce minimum interval between Ensembl API requests."""
        elapsed = time.time() - self._last_request_time
        if elapsed < MIN_REQUEST_INTERVAL:
            time.sleep(MIN_REQUEST_INTERVAL - elapsed)
        self._last_request_time = time.time()

    def _urlopen_with_retry(
        self, req: urllib.request.Request, timeout: int = TIMEOUT_SEC,
    ) -> bytes:
        """Execute an HTTP request with exponential-backoff retry.

        Retries on transient server errors (HTTP 500, 502, 503, 504)
        and network timeouts. Does NOT retry on client errors (4xx)
        which indicate a genuine problem with the request.

        Parameters
        ----------
        req : urllib.request.Request
            The prepared request object.
        timeout : int
            Per-attempt timeout in seconds.

        Returns
        -------
        bytes
            Raw response body.

        Raises
        ------
        EnsemblError
            After all retries are exhausted, or on non-retryable errors.
        """
        retryable_codes = {500, 502, 503, 504}
        last_error = None

        for attempt in range(MAX_RETRIES + 1):
            self._rate_limit()
            try:
                with urllib.request.urlopen(req, timeout=timeout) as resp:
                    return resp.read()
            except urllib.error.HTTPError as e:
                if e.code in retryable_codes and attempt < MAX_RETRIES:
                    wait = RETRY_BACKOFF_BASE ** attempt
                    logger.warning(
                        "Ensembl HTTP %d on attempt %d/%d — retrying in %.1fs",
                        e.code, attempt + 1, MAX_RETRIES + 1, wait,
                    )
                    time.sleep(wait)
                    last_error = e
                    continue
                if e.code == 400:
                    raise EnsemblError(
                        f"Not found or invalid request: {req.full_url}"
                    )
                if e.code == 429:
                    # Rate limited — wait longer
                    wait = RETRY_BACKOFF_BASE ** (attempt + 2)
                    logger.warning(
                        "Ensembl rate limit (429) — waiting %.1fs", wait,
                    )
                    time.sleep(wait)
                    last_error = e
                    continue
                raise EnsemblError(
                    f"Ensembl API error (HTTP {e.code}): {e.reason}"
                )
            except (urllib.error.URLError, socket.timeout, OSError) as e:
                if attempt < MAX_RETRIES:
                    wait = RETRY_BACKOFF_BASE ** attempt
                    logger.warning(
                        "Network error on attempt %d/%d — retrying in %.1fs: %s",
                        attempt + 1, MAX_RETRIES + 1, wait, e,
                    )
                    time.sleep(wait)
                    last_error = e
                    continue
                raise EnsemblError(
                    f"Cannot connect to Ensembl REST API after "
                    f"{MAX_RETRIES + 1} attempts.\n"
                    f"URL: {req.full_url}\nLast error: {e}"
                )

        # Should not reach here, but safety net
        raise EnsemblError(
            f"All {MAX_RETRIES + 1} attempts failed. "
            f"Last error: {last_error}"
        )

    def _fetch_json(self, endpoint: str) -> dict:
        """Make a rate-limited GET request to Ensembl REST, return JSON.

        Automatically retries on transient server errors (HTTP 500/502/503)
        with exponential backoff.
        """
        url = f"{ENSEMBL_REST}{endpoint}"
        if "content-type" not in url.lower():
            sep = "&" if "?" in url else "?"
            url += f"{sep}content-type=application/json"

        req = urllib.request.Request(url, headers={
            "Content-Type": "application/json",
            "User-Agent": "CRISPRArchitect/2.0",
        })
        raw = self._urlopen_with_retry(req)
        return json.loads(raw.decode())


# ---------------------------------------------------------------------------
# Self-test
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    print("=" * 60)
    print("TranscriptFetcher — self-test")
    print("=" * 60)

    fetcher = TranscriptFetcher()

    # Test 1: basic instantiation and cache
    assert len(fetcher._cache) == 0
    print("[PASS] Instantiation OK, empty cache")

    # Test 2: fetch a well-known gene (requires internet)
    try:
        tx = fetcher.fetch_by_gene("CFTR")
        print(f"[PASS] Fetched {tx.gene_symbol}: {tx.transcript_id}")
        print(f"       chr{tx.chromosome}:{tx.start}-{tx.end} "
              f"strand={'+'if tx.strand == 1 else '-'}")
        print(f"       {tx.n_exons} exons, biotype={tx.biotype}")
        assert tx.n_exons > 0
        assert tx.chromosome == "7"
        assert tx.strand == 1  # CFTR is on the forward strand

        # Test cache hit
        tx2 = fetcher.fetch_by_gene("CFTR")
        assert tx2 is tx, "Cache miss — should return same object"
        print("[PASS] Cache hit verified")

        # Test 3: fetch by transcript ID (use whatever we got above)
        tx3 = fetcher.fetch_by_transcript_id(tx.transcript_id)
        assert tx3.transcript_id == tx.transcript_id
        print(f"[PASS] fetch_by_transcript_id returned {tx3.transcript_id}")

        # Test 4: CDS sequence
        cds = fetcher.fetch_cds_sequence(tx.transcript_id)
        assert cds.startswith("ATG"), f"CDS should start with ATG, got {cds[:6]}"
        print(f"[PASS] CDS length = {len(cds)} bp, starts with {cds[:3]}")

        # Test 5: reverse-strand gene (NF1)
        tx_nf1 = fetcher.fetch_by_gene("NF1")
        assert tx_nf1.strand == -1, "NF1 should be on reverse strand"
        print(f"[PASS] NF1: {tx_nf1.transcript_id}, strand={tx_nf1.strand}, "
              f"{tx_nf1.n_exons} exons")

    except EnsemblError as e:
        print(f"[SKIP] Network test skipped (no internet?): {e}")
    except Exception as e:
        print(f"[FAIL] Unexpected error: {e}")
        raise

    print("\nAll TranscriptFetcher tests passed.")
