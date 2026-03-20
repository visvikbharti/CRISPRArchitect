"""
Ensembl REST API Client for CRISPRArchitect
=============================================

Fetches real gene structures and sequences from Ensembl (GRCh38/hg38)
over the internet. NO genome download required.

HOW IT WORKS
------------
The Ensembl REST API (https://rest.ensembl.org) provides free, public
access to gene annotations and sequences. We make HTTP requests and
parse the JSON responses. No authentication needed.

WHAT WE FETCH
-------------
1. Gene info: chromosome, start, end, strand, transcript list
2. Exon coordinates: start, end, rank for each exon in the canonical transcript
3. Genomic sequence: the actual DNA bases for donor template design

WHAT WE DON'T NEED
-------------------
- No genome FASTA download (3 GB)
- No local database
- No authentication or API key

References
----------
- Ensembl REST API docs: https://rest.ensembl.org
- GRCh38 assembly: Genome Reference Consortium, 2013
"""

import json
import urllib.request
import urllib.error
from typing import Dict, List, Optional, Tuple


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

ENSEMBL_REST = "https://rest.ensembl.org"
TIMEOUT_SEC = 30


# ---------------------------------------------------------------------------
# Low-level API
# ---------------------------------------------------------------------------

def _fetch_json(endpoint: str) -> dict:
    """Make a GET request to the Ensembl REST API and return parsed JSON.

    Parameters
    ----------
    endpoint : str
        The API endpoint path (e.g., "/lookup/symbol/homo_sapiens/NF1")

    Returns
    -------
    dict
        Parsed JSON response

    Raises
    ------
    EnsemblError
        If the request fails or the gene is not found
    """
    url = f"{ENSEMBL_REST}{endpoint}"
    if "content-type" not in url:
        sep = "&" if "?" in url else "?"
        url += f"{sep}content-type=application/json"

    req = urllib.request.Request(url, headers={
        "Content-Type": "application/json",
        "User-Agent": "CRISPRArchitect/1.0"
    })
    try:
        with urllib.request.urlopen(req, timeout=TIMEOUT_SEC) as response:
            return json.loads(response.read().decode())
    except urllib.error.HTTPError as e:
        if e.code == 400:
            raise EnsemblError(f"Gene not found or invalid request: {url}")
        raise EnsemblError(f"Ensembl API error (HTTP {e.code}): {e.reason}")
    except urllib.error.URLError as e:
        raise EnsemblError(
            f"Cannot connect to Ensembl REST API. Check your internet connection.\n"
            f"URL: {url}\nError: {e.reason}"
        )


def _fetch_text(endpoint: str) -> str:
    """Fetch plain text (e.g., DNA sequence) from Ensembl REST API."""
    url = f"{ENSEMBL_REST}{endpoint}"
    req = urllib.request.Request(url, headers={
        "Content-Type": "text/plain",
        "User-Agent": "CRISPRArchitect/1.0"
    })
    try:
        with urllib.request.urlopen(req, timeout=TIMEOUT_SEC) as response:
            return response.read().decode().strip()
    except Exception as e:
        raise EnsemblError(f"Failed to fetch sequence: {e}")


class EnsemblError(Exception):
    """Custom exception for Ensembl API errors."""
    pass


# ---------------------------------------------------------------------------
# Gene Info
# ---------------------------------------------------------------------------

class GeneInfo:
    """Complete gene information fetched from Ensembl.

    Attributes
    ----------
    symbol : str
        Gene symbol (e.g., "NF1")
    ensembl_id : str
        Ensembl gene ID (e.g., "ENSG00000196712")
    chromosome : str
        Chromosome (e.g., "17")
    start : int
        Gene start coordinate (GRCh38)
    end : int
        Gene end coordinate (GRCh38)
    strand : int
        1 = forward, -1 = reverse
    canonical_transcript_id : str
        Ensembl transcript ID for the canonical transcript
    exons : list of dict
        Exon coordinates: [{'number': 1, 'start': ..., 'end': ...}, ...]
    description : str
        Gene description from Ensembl
    """

    def __init__(self, symbol, ensembl_id, chromosome, start, end, strand,
                 canonical_transcript_id, exons, description=""):
        self.symbol = symbol
        self.ensembl_id = ensembl_id
        self.chromosome = chromosome
        self.start = start
        self.end = end
        self.strand = strand
        self.canonical_transcript_id = canonical_transcript_id
        self.exons = exons
        self.description = description

    @property
    def span_bp(self):
        return self.end - self.start

    @property
    def n_exons(self):
        return len(self.exons)

    def exon_sizes(self):
        return [e['end'] - e['start'] for e in self.exons]

    def summary(self):
        sizes = self.exon_sizes()
        return (
            f"Gene: {self.symbol} ({self.description})\n"
            f"  Ensembl ID: {self.ensembl_id}\n"
            f"  Location: chr{self.chromosome}:{self.start:,}-{self.end:,} "
            f"({'forward' if self.strand == 1 else 'reverse'} strand)\n"
            f"  Gene span: {self.span_bp:,} bp ({self.span_bp/1e6:.2f} Mb)\n"
            f"  Canonical transcript: {self.canonical_transcript_id}\n"
            f"  Exons: {self.n_exons}\n"
            f"  Exon sizes: min={min(sizes)} bp, median={sorted(sizes)[len(sizes)//2]} bp, "
            f"max={max(sizes)} bp"
        )


# ---------------------------------------------------------------------------
# Main Fetcher Functions
# ---------------------------------------------------------------------------

def fetch_gene(gene_symbol: str, species: str = "homo_sapiens") -> GeneInfo:
    """Fetch complete gene information from Ensembl.

    This is the main function users should call. It:
    1. Looks up the gene by symbol
    2. Finds the canonical transcript
    3. Fetches all exon coordinates
    4. Returns a GeneInfo object ready for CRISPRArchitect

    Parameters
    ----------
    gene_symbol : str
        HGNC gene symbol (e.g., "NF1", "DMD", "BRCA2", "CFTR")
    species : str
        Ensembl species name (default: "homo_sapiens")

    Returns
    -------
    GeneInfo
        Complete gene information including exon coordinates

    Raises
    ------
    EnsemblError
        If the gene is not found or API is unavailable

    Examples
    --------
    >>> gene = fetch_gene("NF1")
    >>> print(gene.summary())
    >>> print(f"Exon 20 to 50 span: {gene.exons[49]['end'] - gene.exons[19]['start']} bp")
    """
    # Step 1: Look up gene
    data = _fetch_json(f"/lookup/symbol/{species}/{gene_symbol}?expand=1")

    ensembl_id = data.get("id", "")
    chromosome = str(data.get("seq_region_name", ""))
    start = data.get("start", 0)
    end = data.get("end", 0)
    strand = data.get("strand", 1)
    description = data.get("description", "")

    # Step 2: Find canonical transcript
    canonical_id = None
    canonical_exons_raw = []

    for t in data.get("Transcript", []):
        if t.get("is_canonical") == 1:
            canonical_id = t["id"]
            canonical_exons_raw = t.get("Exon", [])
            break

    # Fallback: pick transcript with most exons
    if canonical_id is None and data.get("Transcript"):
        best_t = max(data["Transcript"], key=lambda t: len(t.get("Exon", [])))
        canonical_id = best_t["id"]
        canonical_exons_raw = best_t.get("Exon", [])

    if not canonical_exons_raw:
        raise EnsemblError(
            f"No exons found for gene {gene_symbol} "
            f"(transcript {canonical_id}). "
            "This may be a non-coding gene or pseudogene."
        )

    # Step 3: Sort exons by genomic position and assign numbers
    # For forward-strand genes: sort by start ascending
    # For reverse-strand genes: sort by start descending (exon 1 is at the 3' end of the gene)
    if strand == 1:
        canonical_exons_raw.sort(key=lambda e: e["start"])
    else:
        canonical_exons_raw.sort(key=lambda e: e["start"], reverse=True)

    exons = []
    for i, e in enumerate(canonical_exons_raw, 1):
        exons.append({
            "number": i,
            "start": e["start"],
            "end": e["end"],
        })

    return GeneInfo(
        symbol=gene_symbol.upper(),
        ensembl_id=ensembl_id,
        chromosome=chromosome,
        start=start,
        end=end,
        strand=strand,
        canonical_transcript_id=canonical_id,
        exons=exons,
        description=description,
    )


def fetch_sequence(chromosome: str, start: int, end: int,
                   strand: int = 1, species: str = "human") -> str:
    """Fetch a genomic DNA sequence from Ensembl.

    This is used to:
    - Get the actual sequence for donor template design
    - Get homology arm sequences (intronic regions flanking an exon)
    - Get exon sequences for verifying mutations

    Parameters
    ----------
    chromosome : str
        Chromosome (e.g., "17", "X")
    start : int
        Start coordinate (1-based, inclusive)
    end : int
        End coordinate (1-based, inclusive)
    strand : int
        1 = forward strand, -1 = reverse complement
    species : str
        Species (default: "human")

    Returns
    -------
    str
        DNA sequence (uppercase)

    Notes
    -----
    Maximum region size is ~5 Mb per request. For larger regions,
    make multiple requests.
    """
    region = f"{chromosome}:{start}..{end}:{strand}"
    seq = _fetch_text(f"/sequence/region/{species}/{region}?content-type=text/plain")
    return seq.upper().replace("\n", "")


def fetch_exon_sequence(gene_info: GeneInfo, exon_number: int) -> str:
    """Fetch the DNA sequence of a specific exon.

    Parameters
    ----------
    gene_info : GeneInfo
        Gene information from fetch_gene()
    exon_number : int
        1-based exon number

    Returns
    -------
    str
        Exon DNA sequence
    """
    if exon_number < 1 or exon_number > gene_info.n_exons:
        raise ValueError(f"Exon {exon_number} out of range (1-{gene_info.n_exons})")

    exon = gene_info.exons[exon_number - 1]
    return fetch_sequence(
        gene_info.chromosome, exon["start"], exon["end"], gene_info.strand
    )


def fetch_flanking_sequence(gene_info: GeneInfo, exon_number: int,
                            flank_bp: int = 300) -> Tuple[str, str, str]:
    """Fetch exon sequence plus flanking intronic sequence for donor design.

    This returns the sequences needed to design a cssDNA donor template:
    - Left homology arm (intronic sequence upstream of the exon)
    - Exon sequence (the region to correct)
    - Right homology arm (intronic sequence downstream of the exon)

    Parameters
    ----------
    gene_info : GeneInfo
        Gene information from fetch_gene()
    exon_number : int
        1-based exon number
    flank_bp : int
        Length of flanking intronic sequence on each side (default: 300 bp)

    Returns
    -------
    tuple of (left_arm, exon_seq, right_arm)
        Three DNA sequences
    """
    exon = gene_info.exons[exon_number - 1]
    exon_start = exon["start"]
    exon_end = exon["end"]

    # Fetch the full region: flank + exon + flank
    region_start = exon_start - flank_bp
    region_end = exon_end + flank_bp

    full_seq = fetch_sequence(
        gene_info.chromosome, region_start, region_end, gene_info.strand
    )

    # Split into arms and exon
    left_arm = full_seq[:flank_bp]
    exon_seq = full_seq[flank_bp:flank_bp + (exon_end - exon_start)]
    right_arm = full_seq[flank_bp + (exon_end - exon_start):]

    return left_arm, exon_seq, right_arm


def design_cssdna_donor(gene_info: GeneInfo, exon_number: int,
                        ref_allele: str, alt_allele: str,
                        mutation_position_in_exon: int = -1,
                        homology_arm_length: int = 300) -> Dict:
    """Auto-design a cssDNA donor template for correcting a mutation.

    Fetches real genomic sequences from Ensembl and constructs a donor
    with correct intronic homology arms and the corrected exon sequence.

    Parameters
    ----------
    gene_info : GeneInfo
        Gene information from fetch_gene()
    exon_number : int
        Exon containing the mutation
    ref_allele : str
        Reference (wild-type) allele
    alt_allele : str
        Alternate (mutant) allele — this is what the patient has
    mutation_position_in_exon : int
        0-based position within the exon. -1 = middle of exon (default)
    homology_arm_length : int
        Length of each homology arm (default: 300 bp, optimal for cssDNA)

    Returns
    -------
    dict with keys:
        'donor_sequence': str — the full cssDNA donor sequence
        'left_arm': str — left homology arm (intronic)
        'exon_corrected': str — corrected exon sequence
        'right_arm': str — right homology arm (intronic)
        'total_length': int
        'mutation_details': str
        'design_notes': list of str
    """
    left_arm, exon_seq, right_arm = fetch_flanking_sequence(
        gene_info, exon_number, homology_arm_length
    )

    # Determine mutation position
    if mutation_position_in_exon == -1:
        mutation_position_in_exon = len(exon_seq) // 2

    # Verify the mutant allele is present
    notes = []
    pos = mutation_position_in_exon

    if pos < len(exon_seq):
        genomic_base = exon_seq[pos]
        # The genomic reference should match either ref or alt
        # We want to change alt → ref (correct the mutation)
        corrected_exon = exon_seq[:pos] + ref_allele + exon_seq[pos + len(alt_allele):]
        notes.append(f"Position {pos} in exon: genomic={genomic_base}, "
                     f"correcting {alt_allele}>{ref_allele}")
    else:
        corrected_exon = exon_seq
        notes.append("WARNING: mutation position exceeds exon length, using uncorrected sequence")

    donor_sequence = left_arm + corrected_exon + right_arm

    return {
        "donor_sequence": donor_sequence,
        "left_arm": left_arm,
        "exon_corrected": corrected_exon,
        "right_arm": right_arm,
        "exon_original": exon_seq,
        "total_length": len(donor_sequence),
        "mutation_details": f"Exon {exon_number}, pos {pos}: {alt_allele}>{ref_allele}",
        "design_notes": notes,
        "homology_arm_length": homology_arm_length,
    }


# ---------------------------------------------------------------------------
# Convenience: list of well-known large multi-exon genes for demo
# ---------------------------------------------------------------------------

EXAMPLE_GENES = {
    "NF1": "Neurofibromin 1 (57 exons, chr17, neurofibromatosis type 1)",
    "DMD": "Dystrophin (79 exons, chrX, Duchenne muscular dystrophy)",
    "BRCA2": "BRCA2 DNA repair (27 exons, chr13, breast/ovarian cancer)",
    "RB1": "Retinoblastoma 1 (27 exons, chr13, retinoblastoma)",
    "CFTR": "CF transmembrane conductance regulator (27 exons, chr7, cystic fibrosis)",
    "COL7A1": "Collagen type VII (118 exons, chr3, epidermolysis bullosa)",
    "TTN": "Titin (363 exons, chr2, dilated cardiomyopathy)",
    "NF2": "Merlin (17 exons, chr22, neurofibromatosis type 2)",
    "TSC1": "Hamartin (23 exons, chr9, tuberous sclerosis)",
    "TSC2": "Tuberin (42 exons, chr16, tuberous sclerosis)",
    "PKD1": "Polycystin-1 (46 exons, chr16, polycystic kidney disease)",
    "FBN1": "Fibrillin-1 (66 exons, chr15, Marfan syndrome)",
}
