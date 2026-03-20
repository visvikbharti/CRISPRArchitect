"""
Sequence Manipulation Utilities
===============================

Helper functions for DNA sequence operations used across CRISPRArchitect modules.
"""

from typing import Dict, List, Tuple


# Standard genetic code: codon -> amino acid
GENETIC_CODE: Dict[str, str] = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}

# Reverse lookup: amino acid -> list of codons
REVERSE_GENETIC_CODE: Dict[str, List[str]] = {}
for codon, aa in GENETIC_CODE.items():
    if aa not in REVERSE_GENETIC_CODE:
        REVERSE_GENETIC_CODE[aa] = []
    REVERSE_GENETIC_CODE[aa].append(codon)


def complement(base: str) -> str:
    """Return the Watson-Crick complement of a DNA base.

    In DNA, bases pair specifically:
    - Adenine (A) pairs with Thymine (T) via 2 hydrogen bonds
    - Guanine (G) pairs with Cytosine (C) via 3 hydrogen bonds

    Parameters
    ----------
    base : str
        Single DNA base character (A, T, G, or C)

    Returns
    -------
    str
        The complementary base
    """
    comp_map = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
                'a': 't', 't': 'a', 'g': 'c', 'c': 'g'}
    return comp_map.get(base, 'N')


def reverse_complement(sequence: str) -> str:
    """Return the reverse complement of a DNA sequence.

    In double-stranded DNA, the two strands run antiparallel (5'->3' and 3'->5').
    The reverse complement of a sequence represents the opposing strand read 5'->3'.

    This is important for donor design because:
    - The sgRNA targets one strand (target strand)
    - The other strand (non-target strand) is displaced by Cas9
    - ssDNA donors complementary to the non-target strand may have higher HDR
      efficiency (Richardson et al., Nature Biotechnology, 2016)

    Parameters
    ----------
    sequence : str
        DNA sequence (5' to 3')

    Returns
    -------
    str
        Reverse complement sequence (5' to 3' of the opposite strand)
    """
    return ''.join(complement(b) for b in reversed(sequence))


def gc_content(sequence: str) -> float:
    """Calculate the GC content of a DNA sequence.

    GC content affects:
    - DNA melting temperature (higher GC = higher Tm)
    - Secondary structure stability (GC base pairs have 3 H-bonds vs 2 for AT)
    - Donor template stability and HDR efficiency
    - G-quadruplex formation potential (high G content)

    Parameters
    ----------
    sequence : str
        DNA sequence

    Returns
    -------
    float
        Fraction of bases that are G or C (0.0 to 1.0)
    """
    seq_upper = sequence.upper()
    gc_count = seq_upper.count('G') + seq_upper.count('C')
    total = len(seq_upper)
    if total == 0:
        return 0.0
    return gc_count / total


def find_pam_sites(sequence: str, pam: str = "NGG", strand: str = "both") -> List[Dict]:
    """Find all PAM (Protospacer Adjacent Motif) sites in a sequence.

    The PAM is a short DNA sequence (2-6 bp) that must be present adjacent to
    the target site for Cas9 to bind and cut. Different Cas proteins recognize
    different PAMs:

    - SpCas9: 5'-NGG-3' (N = any nucleotide)
    - enFnCas9: 5'-NRG-3' or 5'-NGR-3' (R = A or G, broadened PAM)
    - Cas12a: 5'-TTTV-3' (V = A, C, or G)

    Parameters
    ----------
    sequence : str
        DNA sequence to search
    pam : str
        PAM motif using IUPAC codes (N=any, R=A/G, Y=C/T, V=A/C/G)
    strand : str
        "forward", "reverse", or "both"

    Returns
    -------
    list of dict
        Each dict contains: position, strand, pam_sequence, protospacer_end
    """
    import re

    # IUPAC ambiguity codes to regex
    iupac_to_regex = {
        'N': '[ATCG]', 'R': '[AG]', 'Y': '[CT]',
        'S': '[GC]', 'W': '[AT]', 'K': '[GT]',
        'M': '[AC]', 'B': '[CGT]', 'D': '[AGT]',
        'H': '[ACT]', 'V': '[ACG]',
        'A': 'A', 'T': 'T', 'G': 'G', 'C': 'C',
    }
    pattern = ''.join(iupac_to_regex.get(b, b) for b in pam.upper())

    results = []
    seq_upper = sequence.upper()

    if strand in ("forward", "both"):
        for match in re.finditer(f'(?={pattern})', seq_upper):
            results.append({
                'position': match.start(),
                'strand': '+',
                'pam_sequence': seq_upper[match.start():match.start() + len(pam)],
            })

    if strand in ("reverse", "both"):
        rc_seq = reverse_complement(seq_upper)
        for match in re.finditer(f'(?={pattern})', rc_seq):
            # Convert position back to forward strand coordinates
            fwd_pos = len(seq_upper) - match.start() - len(pam)
            results.append({
                'position': fwd_pos,
                'strand': '-',
                'pam_sequence': rc_seq[match.start():match.start() + len(pam)],
            })

    return sorted(results, key=lambda x: x['position'])


def get_synonymous_codons(codon: str) -> List[str]:
    """Get all synonymous codons (encoding the same amino acid).

    Synonymous codons encode the same amino acid but use different nucleotides.
    This is useful for donor template optimization: we can change the DNA sequence
    to disrupt problematic secondary structures (hairpins, G-quadruplexes) without
    changing the encoded protein.

    Parameters
    ----------
    codon : str
        Three-letter DNA codon

    Returns
    -------
    list of str
        All codons encoding the same amino acid (excluding the input codon)
    """
    codon_upper = codon.upper()
    if codon_upper not in GENETIC_CODE:
        return []
    aa = GENETIC_CODE[codon_upper]
    return [c for c in REVERSE_GENETIC_CODE[aa] if c != codon_upper]


def is_transition(ref: str, alt: str) -> bool:
    """Check if a single nucleotide change is a transition.

    Transitions are purine<->purine (A<->G) or pyrimidine<->pyrimidine (C<->T).
    Transversions are purine<->pyrimidine changes.

    This matters for editing strategy because:
    - Cytosine Base Editors (CBE) can correct C>T transitions (or G>A on opposite strand)
    - Adenine Base Editors (ABE) can correct A>G transitions (or T>C on opposite strand)
    - Transversions CANNOT be corrected by standard base editors
    - Transversions require HDR or prime editing

    Parameters
    ----------
    ref : str
        Reference allele (single nucleotide)
    alt : str
        Alternate allele (single nucleotide)

    Returns
    -------
    bool
        True if the change is a transition
    """
    transitions = {('A', 'G'), ('G', 'A'), ('C', 'T'), ('T', 'C')}
    return (ref.upper(), alt.upper()) in transitions


def classify_snv(ref: str, alt: str) -> str:
    """Classify a single nucleotide variant.

    Returns
    -------
    str
        One of: "transition_AG", "transition_CT", "transversion"
    """
    ref, alt = ref.upper(), alt.upper()
    if (ref, alt) in {('A', 'G'), ('G', 'A'), ('T', 'C'), ('C', 'T')}:
        if ref in ('A', 'T'):
            return "transition_AG"  # ABE-amenable (on appropriate strand)
        else:
            return "transition_CT"  # CBE-amenable (on appropriate strand)
    return "transversion"
