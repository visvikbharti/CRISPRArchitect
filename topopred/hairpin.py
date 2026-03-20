"""
Hairpin / Stem-Loop Predictor for cssDNA Donor Templates
=========================================================

BIOLOGICAL BACKGROUND — WHY HAIRPINS MATTER FOR HDR
----------------------------------------------------

Single-stranded DNA (ssDNA) is not truly "single-stranded" in solution. Because
DNA bases can form Watson-Crick pairs (A-T and G-C), any ssDNA molecule with
internal complementary regions will fold back on itself to form HAIRPIN structures
(also called stem-loops). A hairpin consists of:

    5' ------====LLLL====------ 3'
              STEM  LOOP  STEM
              (dsDNA)     (complement of left stem)

  - STEM: A double-stranded region formed by intramolecular base pairing.
    G-C pairs contribute 3 hydrogen bonds (~-2.0 kcal/mol stacking energy),
    while A-T pairs contribute 2 hydrogen bonds (~-1.0 kcal/mol stacking energy).

  - LOOP: The unpaired single-stranded region at the tip of the hairpin.
    The loop must be at least 3 nucleotides — smaller loops are sterically
    impossible because the DNA backbone cannot make a sharp enough turn.
    Loop formation imposes an ENTROPIC PENALTY that opposes hairpin formation.

THERMODYNAMICS OF HAIRPIN FORMATION:
  The free energy (dG) of a hairpin is the balance between:
    (1) FAVORABLE: base-pair stacking in the stem (negative dG, stabilizing)
    (2) UNFAVORABLE: loop closure entropy (positive dG, destabilizing)

  dG_hairpin = dG_stem + dG_loop
  dG_stem = sum of nearest-neighbor stacking energies for each bp in the stem
  dG_loop = +4.0 + 1.4 * ln(loop_length) kcal/mol  (Jacobson-Stockmayer)

  A hairpin is considered "stable" if dG < -3.0 kcal/mol at 37 degrees C.
  This means it persists long enough to interfere with biological processes.

WHY HAIRPINS ARE CRITICAL FOR cssDNA-MEDIATED HDR:
---------------------------------------------------

1. RAD51 REQUIRES SINGLE-STRANDED DNA:
   The RAD51 recombinase — the engine of homology-directed repair — binds
   cooperatively along ssDNA to form a nucleoprotein filament. RAD51 CANNOT
   bind to double-stranded DNA (dsDNA). Therefore, any region of the donor
   that folds into a hairpin stem becomes INVISIBLE to RAD51.

2. HOMOLOGY ARM SEQUESTRATION:
   For HDR to work, the donor's homology arms must be single-stranded and
   available for homology search. If a hairpin forms WITHIN a homology arm,
   the base-paired nucleotides are sequestered in a dsDNA-like structure
   and cannot participate in strand invasion. This directly reduces HDR
   efficiency.

   Example: A 300-nt homology arm with a 20-bp hairpin stem loses 40 nt
   (both sides of the stem) from productive homology search — a 13% reduction
   in effective arm length.

3. RAD51 FILAMENT DISCONTINUITY:
   Even if the hairpin is small, it creates a GAP in the RAD51 filament.
   RAD51 filaments require ~15 bp of continuous ssDNA for stable nucleation
   (see constants.py: RAD51_NUCLEATION_MIN_MONOMERS * RAD51_FOOTPRINT_NT).
   A hairpin that fragments the ssDNA into segments shorter than this minimum
   can completely prevent filament formation.

4. cssDNA IS ESPECIALLY VULNERABLE:
   Circular ssDNA (cssDNA) donors are entirely single-stranded, so they are
   free to fold into complex secondary structures. Unlike linear ssDNA donors
   that may have RPA-coated regions, cssDNA donors delivered to the cell
   initially lack protective ssDNA-binding proteins and are free to fold.

SIMPLIFIED ENERGY MODEL (NO EXTERNAL DEPENDENCIES):
  This module implements a simplified nearest-neighbor model:
    - GC stacking: -2.0 kcal/mol per base pair
    - AT stacking: -1.0 kcal/mol per base pair
    - Loop penalty: +4.0 + 1.4 * ln(loop_length) kcal/mol

  Full nearest-neighbor parameters (SantaLucia, 1998) use 10 unique dinucleotide
  stacking parameters and salt correction terms. Our simplified model captures
  the essential physics: GC-rich stems are more stable than AT-rich stems, and
  larger loops destabilize the hairpin. For donor design purposes, this level
  of accuracy is sufficient to identify problematic structures.

References:
  - SantaLucia, PNAS, 1998 (nearest-neighbor thermodynamic parameters)
  - Zuker, Nucleic Acids Res, 2003 (mfold algorithm)
  - Jacobson & Stockmayer, J Chem Phys, 1950 (loop entropy)
  - Conway & Bhaskaran, Mol Biotechnol, 2019 (ssDNA donor design)
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import List, Optional, Tuple


# ---------------------------------------------------------------------------
# Thermodynamic Constants (Simplified Nearest-Neighbor Model)
# ---------------------------------------------------------------------------

# Stacking free energy per base pair (kcal/mol, at 37 degrees C, 1 M NaCl).
#
# Full nearest-neighbor models use 10 unique dinucleotide stacking parameters,
# but for our purposes, we simplify to two classes:
#
# GC pairs: 3 hydrogen bonds, strong stacking → -2.0 kcal/mol per bp
#   (actual range: -1.84 to -3.42 kcal/mol depending on neighbors)
#
# AT pairs: 2 hydrogen bonds, weaker stacking → -1.0 kcal/mol per bp
#   (actual range: -0.88 to -1.45 kcal/mol depending on neighbors)
#
# These simplified values capture the ~2-fold difference in stability between
# GC and AT pairs that is the primary determinant of hairpin stability.
GC_STACK_ENERGY = -2.0  # kcal/mol per GC base pair
AT_STACK_ENERGY = -1.0  # kcal/mol per AT base pair

# Loop closure free energy (Jacobson-Stockmayer approximation):
#   dG_loop = LOOP_INITIATION + LOOP_LENGTH_FACTOR * ln(loop_length)
#
# The initiation penalty (+4.0 kcal/mol) accounts for the loss of
# conformational entropy when the backbone makes a U-turn at the loop.
# The logarithmic term (+1.4 * ln(N)) reflects the Jacobson-Stockmayer
# theory: the number of conformations accessible to a polymer loop
# scales as N^(-3/2), giving a free energy penalty of ~1.5*R*T*ln(N)
# at 37 degrees C (~1.4 kcal/mol * ln(N)).
LOOP_INITIATION_PENALTY = 4.0  # kcal/mol (entropy cost of forming any loop)
LOOP_LENGTH_FACTOR = 1.4       # kcal/mol, multiplied by ln(loop_length)

# Stability threshold: hairpins with dG below this are considered "stable"
# at physiological temperature (37 degrees C). A dG of -3.0 kcal/mol
# corresponds to a folded:unfolded ratio of approximately 150:1 at 37 C
# (using dG = -RT*ln(K), where R = 0.001987 kcal/(mol*K), T = 310 K).
# This means the hairpin is folded >99% of the time.
STABILITY_THRESHOLD_KCAL = -3.0


# ---------------------------------------------------------------------------
# Watson-Crick Complementarity
# ---------------------------------------------------------------------------

# DNA base pairing rules. In DNA, guanine pairs with cytosine (3 H-bonds)
# and adenine pairs with thymine (2 H-bonds). We also accept G-T wobble
# pairs as a minor contributor (they occur in biological structures but
# are much weaker than canonical pairs). For this simplified model, we
# only count Watson-Crick pairs.
_WC_COMPLEMENT = {
    'A': 'T', 'T': 'A',
    'G': 'C', 'C': 'G',
}


def _is_complement(base1: str, base2: str) -> bool:
    """
    Check if two DNA bases form a Watson-Crick pair.

    Watson-Crick base pairs are the foundation of DNA double helix stability:
      A pairs with T via 2 hydrogen bonds
      G pairs with C via 3 hydrogen bonds

    Parameters
    ----------
    base1, base2 : str
        Single-character DNA bases (uppercase).

    Returns
    -------
    bool
        True if the bases form a Watson-Crick pair.
    """
    return _WC_COMPLEMENT.get(base1) == base2


def _pair_energy(base1: str, base2: str) -> float:
    """
    Return the approximate stacking free energy for a Watson-Crick base pair.

    This is a simplified model where:
      G-C pairs contribute -2.0 kcal/mol (3 hydrogen bonds + strong stacking)
      A-T pairs contribute -1.0 kcal/mol (2 hydrogen bonds + weaker stacking)

    In reality, the energy depends on the neighboring base pairs (nearest-
    neighbor model), but this simplification captures the dominant effect:
    GC-rich stems are approximately twice as stable as AT-rich stems.

    Parameters
    ----------
    base1, base2 : str
        A Watson-Crick base pair (must already be verified as complementary).

    Returns
    -------
    float
        Stacking energy in kcal/mol (negative = stabilizing).
    """
    pair = frozenset((base1, base2))
    if pair == frozenset(('G', 'C')):
        return GC_STACK_ENERGY
    elif pair == frozenset(('A', 'T')):
        return AT_STACK_ENERGY
    else:
        # Should not reach here for canonical Watson-Crick pairs.
        return 0.0


# ---------------------------------------------------------------------------
# Data Classes
# ---------------------------------------------------------------------------

@dataclass
class Hairpin:
    """
    Represents a predicted hairpin (stem-loop) structure in a DNA sequence.

    A hairpin forms when a segment of ssDNA folds back on itself, creating a
    double-stranded stem closed by a single-stranded loop:

        5' ...NNNNSTEM_5'---LOOP---STEM_3'NNNN... 3'
                  ||||               ||||
              base pairs         base pairs
              (the stem 3' is reverse complement of stem 5')

    The stem_start and stem_end define the outer boundaries of the hairpin.
    Within those boundaries:
      - [stem_start, loop_start) is the 5' half of the stem
      - [loop_start, loop_end) is the loop
      - [loop_end, stem_end) is the 3' half of the stem

    Attributes
    ----------
    stem_start : int
        0-based index of the first nucleotide in the 5' stem half.
    stem_end : int
        0-based index one past the last nucleotide in the 3' stem half.
    loop_start : int
        0-based index of the first nucleotide in the loop.
    loop_end : int
        0-based index one past the last nucleotide in the loop.
    stem_sequence_5prime : str
        Sequence of the 5' half of the stem (read 5'->3').
    stem_sequence_3prime : str
        Sequence of the 3' half of the stem (read 5'->3').
        This is the reverse complement of stem_sequence_5prime.
    loop_sequence : str
        Sequence of the unpaired loop region.
    stem_length : int
        Number of base pairs in the stem.
    loop_length : int
        Number of nucleotides in the loop.
    free_energy_kcal_mol : float
        Estimated Gibbs free energy of hairpin formation (kcal/mol).
        Negative values indicate thermodynamically favorable (stable) hairpins.
    is_stable : bool
        True if free_energy_kcal_mol < STABILITY_THRESHOLD_KCAL (-3.0 kcal/mol).
        Stable hairpins are predicted to persist at 37 degrees C and impair HDR.
    gc_content : float
        Fraction of GC pairs in the stem (0.0 to 1.0).
        Higher GC content = more stable stem.
    """
    stem_start: int = 0
    stem_end: int = 0
    loop_start: int = 0
    loop_end: int = 0
    stem_sequence_5prime: str = ""
    stem_sequence_3prime: str = ""
    loop_sequence: str = ""
    stem_length: int = 0
    loop_length: int = 0
    free_energy_kcal_mol: float = 0.0
    is_stable: bool = False
    gc_content: float = 0.0

    @property
    def total_length(self) -> int:
        """Total number of nucleotides involved in this hairpin."""
        return self.stem_end - self.stem_start

    @property
    def structured_positions(self) -> Tuple[range, range]:
        """
        Return the ranges of positions that are BASE-PAIRED (in the stem).

        These are the positions that are NOT accessible to RAD51 — they are
        locked in dsDNA-like structure and cannot participate in homology
        search during strand invasion.

        Returns
        -------
        tuple of (range, range)
            (5_prime_stem_positions, 3_prime_stem_positions)
        """
        return (
            range(self.stem_start, self.loop_start),
            range(self.loop_end, self.stem_end),
        )

    def __repr__(self) -> str:
        return (
            f"Hairpin(pos={self.stem_start}-{self.stem_end}, "
            f"stem={self.stem_length}bp, loop={self.loop_length}nt, "
            f"dG={self.free_energy_kcal_mol:.1f} kcal/mol, "
            f"stable={self.is_stable})"
        )


# ---------------------------------------------------------------------------
# Hairpin Predictor
# ---------------------------------------------------------------------------

class HairpinPredictor:
    """
    Predicts hairpin (stem-loop) structures in ssDNA sequences.

    This predictor uses a sliding-window approach to identify complementary
    regions within a single-stranded DNA sequence that can fold into hairpin
    structures. For each potential hairpin, it computes a simplified free
    energy using nearest-neighbor-inspired parameters and the Jacobson-
    Stockmayer loop entropy model.

    ALGORITHM OVERVIEW:
    ------------------
    For each possible loop position (i, j) where j - i >= min_loop:
      1. Check if bases flanking the loop are complementary (i-1 pairs with j)
      2. Extend the stem outward while bases remain complementary
      3. If stem length >= min_stem, calculate free energy
      4. Record the hairpin if it meets stability criteria

    This is a simplified approach compared to full dynamic programming
    algorithms (like Nussinov or Zuker), but it is:
      - Fast (O(n^2) for typical parameters)
      - Dependency-free (no ViennaRNA needed)
      - Sufficient for identifying problematic structures in donor design

    LIMITATIONS:
      - Does not predict pseudoknots
      - Does not account for multibranch loops or internal loops
      - Uses simplified energy parameters (not full nearest-neighbor)
      - Does not model dangling ends or coaxial stacking

    For high-accuracy structure prediction, use ViennaRNA (RNAfold) or mfold.
    This module is designed for RAPID SCREENING of donor templates.

    Parameters
    ----------
    min_stem : int
        Minimum number of base pairs in the stem (default: 4).
        Stems shorter than 4 bp are generally unstable at 37 degrees C
        because the stacking energy cannot overcome the loop entropy.
    min_loop : int
        Minimum loop size in nucleotides (default: 3).
        Loops of 1-2 nt are sterically impossible — the DNA backbone
        physically cannot make a tight enough turn. The minimum is 3 nt,
        known as a "triloop", which is the tightest possible hairpin.
    max_loop : int
        Maximum loop size in nucleotides (default: 30).
        Very large loops have such high entropic penalty that the hairpin
        is unlikely to form. 30 nt is a generous upper bound; most stable
        hairpins have loops of 3-10 nt.
    """

    def __init__(
        self,
        min_stem: int = 4,
        min_loop: int = 3,
        max_loop: int = 30,
    ):
        self.min_stem = min_stem
        self.min_loop = min_loop
        self.max_loop = max_loop

    # -----------------------------------------------------------------
    # Free energy calculation
    # -----------------------------------------------------------------

    def calculate_free_energy(self, hairpin: Hairpin) -> float:
        """
        Calculate the approximate Gibbs free energy of hairpin formation.

        The free energy determines whether the hairpin will form spontaneously
        at physiological temperature (37 degrees C). The calculation follows
        a simplified nearest-neighbor approach:

        dG_total = dG_stem + dG_loop

        where:
          dG_stem = sum of base-pair stacking energies
                  = (number of GC pairs) * (-2.0) + (number of AT pairs) * (-1.0)
                  Units: kcal/mol

          dG_loop = loop initiation + loop length penalty
                  = 4.0 + 1.4 * ln(loop_length)
                  Units: kcal/mol

        The loop penalty follows the Jacobson-Stockmayer theory for polymer
        ring closure: the probability of a polymer chain of N segments forming
        a ring scales as N^(-3/2), which gives a free energy penalty of
        approximately (3/2) * kT * ln(N). At 37 C, this works out to
        approximately 1.4 * ln(N) kcal/mol.

        Parameters
        ----------
        hairpin : Hairpin
            The hairpin structure to evaluate.

        Returns
        -------
        float
            Free energy in kcal/mol. Negative values are thermodynamically
            favorable (the hairpin wants to form). More negative = more stable.
        """
        # --- Stem energy ---
        # Sum up the stacking energy for each base pair in the stem.
        # The 5' stem half pairs with the 3' stem half in antiparallel
        # orientation. Position i of the 5' half pairs with position
        # (stem_length - 1 - i) of the 3' half, but since we store
        # stem_sequence_3prime in 5'->3' direction, we need to reverse it
        # for the pairing.
        stem_energy = 0.0
        stem_5 = hairpin.stem_sequence_5prime.upper()
        # The 3' arm of the stem pairs in antiparallel: the first base of
        # the 5' arm pairs with the LAST base of the 3' arm.
        stem_3_reversed = hairpin.stem_sequence_3prime.upper()[::-1]

        for b1, b2 in zip(stem_5, stem_3_reversed):
            if _is_complement(b1, b2):
                stem_energy += _pair_energy(b1, b2)
            # Non-complementary positions (mismatches) contribute 0 energy
            # in our simplified model. In reality, mismatches destabilize
            # the stem, but they are typically excluded during hairpin
            # finding (we only extend stems through complementary bases).

        # --- Loop energy ---
        # The entropic cost of forming a loop of N nucleotides.
        # This is always positive (destabilizing) because constraining
        # the backbone into a loop reduces conformational freedom.
        loop_len = hairpin.loop_length
        if loop_len < 1:
            loop_len = 1  # Safety floor (should not happen with min_loop >= 3)
        loop_energy = LOOP_INITIATION_PENALTY + LOOP_LENGTH_FACTOR * math.log(loop_len)

        # --- Total free energy ---
        total_dg = stem_energy + loop_energy
        return total_dg

    # -----------------------------------------------------------------
    # Hairpin finding
    # -----------------------------------------------------------------

    def find_hairpins(
        self,
        sequence: str,
        min_stem: Optional[int] = None,
        min_loop: Optional[int] = None,
        max_loop: Optional[int] = None,
    ) -> List[Hairpin]:
        """
        Find all possible hairpin structures in a DNA sequence.

        This method scans the sequence for all positions where the strand
        can fold back on itself to form a stem-loop. The algorithm:

        1. For each possible loop starting position (i):
           For each possible loop ending position (j = i + min_loop to i + max_loop):
             a. Check if the bases flanking the loop are complementary
             b. Extend the stem outward while complementarity holds
             c. If stem reaches min_stem length, record the hairpin

        This finds ALL possible hairpins, including overlapping ones. In
        practice, overlapping hairpins compete for the same bases, and only
        the most stable one(s) will form. The predict_structure() method
        handles this by selecting non-overlapping hairpins.

        Parameters
        ----------
        sequence : str
            DNA sequence (A, T, C, G only). Case-insensitive.
        min_stem : int, optional
            Override the minimum stem length for this call.
        min_loop : int, optional
            Override the minimum loop length for this call.
        max_loop : int, optional
            Override the maximum loop length for this call.

        Returns
        -------
        list of Hairpin
            All detected hairpins, including potentially overlapping ones.
            Each hairpin has its free energy and stability flag calculated.
        """
        sequence = sequence.upper().strip()
        seq_len = len(sequence)

        # Use instance defaults if not overridden.
        _min_stem = min_stem if min_stem is not None else self.min_stem
        _min_loop = min_loop if min_loop is not None else self.min_loop
        _max_loop = max_loop if max_loop is not None else self.max_loop

        hairpins: List[Hairpin] = []
        # Track found hairpins by (stem_start, stem_end) to avoid duplicates.
        seen: set = set()

        # Iterate over all possible loop regions.
        # The loop occupies positions [loop_start, loop_end).
        # The stem extends outward from the loop.
        for loop_start in range(0, seq_len):
            for loop_len in range(_min_loop, _max_loop + 1):
                loop_end = loop_start + loop_len
                if loop_end >= seq_len:
                    break  # Loop extends beyond the sequence.

                # Try to build a stem extending outward from the loop.
                # The 5' stem extends leftward from loop_start.
                # The 3' stem extends rightward from loop_end.
                stem_len = 0
                left_idx = loop_start - 1  # 5' side, moving leftward
                right_idx = loop_end       # 3' side, moving rightward

                while left_idx >= 0 and right_idx < seq_len:
                    if _is_complement(sequence[left_idx], sequence[right_idx]):
                        stem_len += 1
                        left_idx -= 1
                        right_idx += 1
                    else:
                        break  # Stem extension stops at first mismatch.

                # Check if stem meets minimum length requirement.
                if stem_len < _min_stem:
                    continue

                # Define the hairpin boundaries.
                # stem_start = first nucleotide of the 5' stem half
                # stem_end = one past last nucleotide of the 3' stem half
                stem_start = left_idx + 1
                stem_end = right_idx

                key = (stem_start, stem_end)
                if key in seen:
                    continue
                seen.add(key)

                # Extract sequences.
                stem_5prime = sequence[stem_start:loop_start]
                stem_3prime = sequence[loop_end:stem_end]
                loop_seq = sequence[loop_start:loop_end]

                # Calculate GC content of the stem.
                # This is the fraction of base pairs that are G-C (vs. A-T).
                gc_pairs = sum(
                    1 for b in stem_5prime if b in ('G', 'C')
                )
                gc_content = gc_pairs / max(stem_len, 1)

                # Build the Hairpin object.
                hp = Hairpin(
                    stem_start=stem_start,
                    stem_end=stem_end,
                    loop_start=loop_start,
                    loop_end=loop_end,
                    stem_sequence_5prime=stem_5prime,
                    stem_sequence_3prime=stem_3prime,
                    loop_sequence=loop_seq,
                    stem_length=stem_len,
                    loop_length=loop_len,
                    gc_content=gc_content,
                )

                # Calculate free energy and stability.
                hp.free_energy_kcal_mol = self.calculate_free_energy(hp)
                hp.is_stable = hp.free_energy_kcal_mol < STABILITY_THRESHOLD_KCAL

                hairpins.append(hp)

        return hairpins

    # -----------------------------------------------------------------
    # Structure prediction (non-overlapping, most stable hairpins)
    # -----------------------------------------------------------------

    def predict_structure(self, sequence: str) -> List[Hairpin]:
        """
        Predict the most likely set of hairpin structures in a DNA sequence.

        Unlike find_hairpins() which returns ALL possible hairpins (including
        overlapping ones that cannot co-exist), this method selects a
        non-overlapping set of hairpins that maximizes the total energy
        gained. This approximates the minimum free energy (MFE) structure.

        ALGORITHM:
        ----------
        1. Find all possible hairpins.
        2. Sort by free energy (most stable first, i.e., most negative dG).
        3. Greedily select hairpins, skipping any that overlap with an
           already-selected hairpin.

        This greedy approach does not guarantee the global MFE, but it is
        fast and gives a good approximation for donor design screening.
        For exact MFE prediction, use ViennaRNA (RNAfold).

        Parameters
        ----------
        sequence : str
            DNA sequence to analyze.

        Returns
        -------
        list of Hairpin
            Non-overlapping hairpins sorted by free energy (most stable first).
            Only stable hairpins (dG < -3.0 kcal/mol) are included, as
            unstable hairpins are transient and unlikely to interfere with HDR.
        """
        # Find all possible hairpins.
        all_hairpins = self.find_hairpins(sequence)

        # Filter to only stable hairpins (those likely to persist and
        # interfere with RAD51 loading).
        stable = [hp for hp in all_hairpins if hp.is_stable]

        # Sort by free energy: most stable (most negative dG) first.
        # These are the hairpins most likely to form and most damaging to HDR.
        stable.sort(key=lambda hp: hp.free_energy_kcal_mol)

        # Greedy non-overlapping selection.
        # We pick the most stable hairpin first, then skip any that overlap.
        # Two hairpins overlap if their [stem_start, stem_end) intervals
        # intersect, because the same nucleotide cannot participate in two
        # different hairpin stems simultaneously.
        selected: List[Hairpin] = []
        occupied: set = set()  # Positions already claimed by a hairpin.

        for hp in stable:
            hp_positions = set(range(hp.stem_start, hp.stem_end))
            if hp_positions & occupied:
                # This hairpin overlaps with an already-selected hairpin.
                # Skip it — the overlapping hairpin was more stable.
                continue
            selected.append(hp)
            occupied.update(hp_positions)

        return selected

    # -----------------------------------------------------------------
    # Convenience: get structured positions
    # -----------------------------------------------------------------

    def get_structured_positions(self, sequence: str) -> set:
        """
        Return the set of nucleotide positions involved in stable hairpin stems.

        These are the positions that are BASE-PAIRED (double-stranded) and
        therefore INACCESSIBLE to RAD51 for filament formation. They represent
        the "structured" fraction of the sequence.

        This is useful for the AccessibilityScorer to determine what fraction
        of a homology arm is available for strand invasion.

        Parameters
        ----------
        sequence : str
            DNA sequence to analyze.

        Returns
        -------
        set of int
            0-based positions that are in hairpin stems (base-paired).
        """
        hairpins = self.predict_structure(sequence)
        structured = set()
        for hp in hairpins:
            # Add positions from both halves of the stem.
            # These are the nucleotides locked in base pairs.
            stem_5_range, stem_3_range = hp.structured_positions
            structured.update(stem_5_range)
            structured.update(stem_3_range)
        return structured
