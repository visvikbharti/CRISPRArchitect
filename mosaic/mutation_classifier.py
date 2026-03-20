"""
Mutation Classification for MOSAIC
====================================

This module classifies individual pathogenic mutations to determine which
genome-editing technologies can correct them.  The classification directly
controls which strategies the ``StrategyEnumerator`` will consider.

Background: why mutation type determines editing strategy
---------------------------------------------------------

Not all genome-editing tools are created equal.  The three main classes of
editors — **base editors**, **prime editors**, and **nuclease + HDR** —
have very different capabilities, efficiencies, and safety profiles:

1. **Cytosine base editors (CBEs)**
   - Convert C*G to T*A (or equivalently, G on the target strand to A).
   - Mechanism: a catalytically impaired Cas9 (nickase or dead) fused to
     a cytidine deaminase (e.g., APOBEC1) and a uracil glycosylase
     inhibitor (UGI).  The deaminase converts C to U; the cell reads U
     as T during replication.
   - **No double-strand break (DSB) is introduced** — only a nick on
     the non-edited strand to bias mismatch repair.
   - Editing window: typically positions 4-8 in the protospacer
     (counting the PAM-distal end as position 1).
   - Efficiency: 20-60% in most cell types.
   - Source: Komor et al., Nature, 2016.

2. **Adenine base editors (ABEs)**
   - Convert A*T to G*C (or equivalently, T on the target strand to C).
   - Mechanism: Cas9 nickase fused to an evolved tRNA adenosine
     deaminase (TadA).  The deaminase converts A to inosine (I); the
     cell reads I as G.
   - **No DSB.**
   - Efficiency: 30-70% in most cell types; generally higher purity
     than CBEs.
   - Source: Gaudelli et al., Nature, 2017.

3. **Prime editors (PEs)**
   - Can install any small edit: substitutions, insertions (up to ~40 bp),
     deletions (up to ~80 bp), and combinations thereof.
   - Mechanism: Cas9 H840A nickase fused to an engineered reverse
     transcriptase (RT).  A prime editing guide RNA (pegRNA) contains
     both the spacer (for target binding) and a 3' extension that
     encodes the desired edit plus a primer binding site (PBS).
   - **No DSB** — only a nick on the PAM-containing strand.
   - Efficiency: 5-50%, highly variable by locus and edit type.
   - PE2, PE3, PE4, PE5 variants trade off efficiency vs. indel
     byproducts.
   - Source: Anzalone et al., Nature, 2019; Chen et al., Cell, 2021.

4. **Nuclease + HDR (homology-directed repair)**
   - Can correct *any* mutation, including large insertions, deletions,
     and complex rearrangements.
   - Mechanism: a nuclease (SpCas9, Cas12a, vCas9, etc.) creates a DSB
     at or near the mutation site.  The cell's HDR pathway uses a
     provided donor template (ssODN, cssDNA, plasmid, AAV) to
     incorporate the desired correction.
   - **Requires a DSB**, which triggers:
       * p53-dependent DNA-damage checkpoint (major concern in iPSCs;
         Ihry et al., Nature Medicine, 2018)
       * Risk of unintended NHEJ/MMEJ outcomes (~40-70% of alleles)
       * Risk of large deletions at the cut site (Kosicki et al.,
         Nature Biotechnology, 2018)
       * Translocation risk when two DSBs are introduced
   - Efficiency: highly variable (2-50%), depends on cell type, donor
     topology, homology arm design, and cell cycle phase.

**Key take-away for MOSAIC**: Whenever possible, prefer DSB-free
approaches (base editing, prime editing) because they avoid p53
selection and translocation risk.  Fall back to HDR only when the
mutation type demands it (e.g., transversions that cannot be corrected
by base editors, or large structural changes).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional


# ---------------------------------------------------------------------------
# Mutation dataclass
# ---------------------------------------------------------------------------

@dataclass
class Mutation:
    """A single pathogenic variant in the patient's genome.

    Attributes
    ----------
    exon_number : int
        1-based exon number where this mutation resides.
    position : int
        Genomic coordinate of the mutation (0-based).  For multi-base
        changes, this is the position of the first affected base.
    ref_allele : str
        Reference (wild-type) allele.  For a point mutation this is a
        single nucleotide; for an indel it can be longer.  Use ``"-"``
        for a pure insertion (nothing in the reference).
    alt_allele : str
        Alternate (mutant) allele carried by the patient.  Use ``"-"``
        for a pure deletion.
    mutation_type : str
        Auto-classified mutation type.  Set by ``MutationClassifier``
        after construction.  One of: ``"transition_CT"``,
        ``"transition_AG"``, ``"transversion"``, ``"insertion"``,
        ``"deletion"``, ``"complex"``.
    name : str
        Optional human-readable label (e.g., ``"c.6527insC"``).

    Notes
    -----
    Following HGVS nomenclature, the ``ref_allele`` and ``alt_allele``
    are given in the sense (coding) strand orientation.  However, the
    classifier accounts for complementary base editing on the opposite
    strand — for example, a G>A mutation on the sense strand is a C>T
    change on the antisense strand, which is the natural substrate for
    a CBE if the guide targets the antisense strand.
    """

    exon_number: int
    position: int
    ref_allele: str
    alt_allele: str
    mutation_type: str = ""
    name: str = ""

    def __post_init__(self) -> None:
        """Normalise alleles to uppercase and auto-classify if needed."""
        self.ref_allele = self.ref_allele.upper().strip()
        self.alt_allele = self.alt_allele.upper().strip()
        if not self.mutation_type:
            classifier = MutationClassifier()
            self.mutation_type = classifier.classify(
                self.ref_allele, self.alt_allele
            )
        if not self.name:
            self.name = (
                f"exon{self.exon_number}:{self.ref_allele}>"
                f"{self.alt_allele}"
            )


# ---------------------------------------------------------------------------
# MutationClassifier
# ---------------------------------------------------------------------------

class MutationClassifier:
    """Classify mutations to determine editing-technology compatibility.

    The classifier examines the reference and alternate alleles and
    assigns each mutation to one of six categories.  These categories
    map directly to the editing technologies that can correct the
    mutation:

    +------------------+-------------------+------------------------------+
    | Category         | Example           | Compatible editors           |
    +==================+===================+==============================+
    | transition_CT    | C>T or G>A        | CBE (+ PE, HDR)              |
    +------------------+-------------------+------------------------------+
    | transition_AG    | A>G or T>C        | ABE (+ PE, HDR)              |
    +------------------+-------------------+------------------------------+
    | transversion     | A>C, G>T, etc.    | PE, HDR only                 |
    +------------------+-------------------+------------------------------+
    | insertion        | ->A, ->ACGT       | PE (<~40 bp), HDR            |
    +------------------+-------------------+------------------------------+
    | deletion         | ACG->-            | PE (<~80 bp), HDR            |
    +------------------+-------------------+------------------------------+
    | complex          | AC>GT             | PE (if small), HDR           |
    +------------------+-------------------+------------------------------+

    Biology note — why "correcting" vs. "installing"
    --------------------------------------------------
    The patient carries the *alt_allele* and we want to revert to the
    *ref_allele*.  So the "correction" edit is alt -> ref.  When we ask
    "is this amenable to CBE?", we are asking: can a CBE convert the
    patient's alt allele back to the ref allele?  A CBE converts C>T
    (on the edited strand), so it can revert a T>C mutation (ref=C,
    alt=T: we need to convert the patient's T back to C ... wait, CBE
    goes C->T, not T->C).

    **Important**: we must think about corrections in *both* strand
    orientations.  If the mutation is G>A on the sense strand, then on
    the antisense strand it is C>T.  A CBE targeting the antisense
    strand would convert the C (antisense) to T, installing G>A — that
    is the *mutant* direction.  To *correct* a G>A mutation, we need to
    go A>G, which is what an **ABE** does.

    The classifier therefore labels mutations by the *correction* needed:
      - ``"transition_CT"``: the correction requires C>T (or G>A), i.e.,
        a CBE can do the correction.  This means ref has C (or G) and
        alt has T (or A) — the patient has a C>T (or G>A) mutation, and
        we need to reverse it... Actually, let's be very precise:

        If ref=C, alt=T: patient has T, we want C.  Correction = T>C = ABE.
        If ref=G, alt=A: patient has A, we want G.  Correction = A>G = ABE.

    **Revised logic** (this is critical to get right):

    The *correction* edit is:  alt_allele -> ref_allele.

      * alt=T, ref=C  ->  correction is T>C  ->  ABE can do A>G on the
        opposite strand, which is equivalent to T>C.  **ABE-amenable.**
      * alt=A, ref=G  ->  correction is A>G  ->  ABE directly.  **ABE.**
      * alt=C, ref=T  ->  correction is C>T  ->  CBE directly.  **CBE.**
      * alt=G, ref=A  ->  correction is G>A  ->  CBE on opposite strand
        (C>T on antisense = G>A on sense).  **CBE.**
      * alt=A, ref=C  ->  correction is A>C  ->  transversion.
      * etc.

    So we classify by the *correction direction* (alt -> ref):

      * transition_CT:  correction is C>T or G>A  ->  CBE-amenable
        (ref in {T, A}, alt in {C, G}, and they form a transition pair)
      * transition_AG:  correction is A>G or T>C  ->  ABE-amenable
        (ref in {G, C}, alt in {A, T}, and they form a transition pair)
    """

    # Transition pairs: these are the four possible single-nucleotide
    # transitions (purine<->purine or pyrimidine<->pyrimidine).
    # Complementary pairs for strand-aware analysis:
    #   A<->G (purine transition) and C<->T (pyrimidine transition)
    _TRANSITIONS = {
        frozenset({"A", "G"}),
        frozenset({"C", "T"}),
    }

    def classify(self, ref: str, alt: str) -> str:
        """Classify a mutation based on reference and alternate alleles.

        Parameters
        ----------
        ref : str
            Reference (wild-type) allele.  Single nucleotide for SNVs;
            multiple nucleotides or ``"-"`` for indels.
        alt : str
            Alternate (patient) allele.

        Returns
        -------
        str
            One of ``"transition_CT"``, ``"transition_AG"``,
            ``"transversion"``, ``"insertion"``, ``"deletion"``,
            ``"complex"``.

        Examples
        --------
        >>> mc = MutationClassifier()
        >>> mc.classify("C", "T")   # Patient has T, want C back -> ABE
        'transition_AG'
        >>> mc.classify("A", "G")   # Patient has G, want A back -> CBE
        'transition_CT'
        """
        ref = ref.upper().strip()
        alt = alt.upper().strip()

        # --- Handle indels ---
        if ref == "-" or len(ref) == 0:
            # Pure insertion (nothing in reference, something in alt)
            return "insertion"
        if alt == "-" or len(alt) == 0:
            # Pure deletion (something in reference, nothing in alt)
            return "deletion"
        if len(ref) != len(alt):
            # Length-changing: could be a complex indel
            if len(ref) < len(alt):
                return "insertion"
            else:
                return "deletion"

        # --- Handle multi-nucleotide variants (MNVs) ---
        if len(ref) > 1:
            return "complex"

        # --- Single nucleotide variant (SNV) ---
        # Check if this is a transition (purine<->purine or pyrimidine<->pyrimidine).
        pair = frozenset({ref, alt})
        if pair in self._TRANSITIONS:
            # It is a transition.  Now classify by the *correction* direction.
            # Correction = alt -> ref (we want to revert the patient's
            # allele back to reference).

            # Determine what the correction edit is:
            correction_from = alt  # What the patient has
            correction_to = ref    # What we want

            # CBE performs: C>T (on the edited strand)
            # Equivalently on the opposite strand: G>A
            # So CBE can accomplish corrections where the net change is
            # C>T or G>A.
            #
            # ABE performs: A>G (on the edited strand)
            # Equivalently on the opposite strand: T>C
            # So ABE can accomplish corrections where the net change is
            # A>G or T>C.

            if (correction_from == "C" and correction_to == "T") or \
               (correction_from == "G" and correction_to == "A"):
                # Correction is C>T or G>A -> CBE can do this
                return "transition_CT"
            elif (correction_from == "A" and correction_to == "G") or \
                 (correction_from == "T" and correction_to == "C"):
                # Correction is A>G or T>C -> ABE can do this
                return "transition_AG"
            else:
                # Should not reach here for transitions, but be safe
                return "transversion"
        else:
            # Transversion (purine <-> pyrimidine)
            return "transversion"

    def base_editing_amenable(self, mutation: Mutation) -> bool:
        """Can this mutation be corrected by a base editor (CBE or ABE)?

        Base editing is the safest class of genome editing because it
        introduces **no double-strand break**.  Instead, a chemically
        modified base (deaminated C or A) is processed by the cell's
        mismatch repair machinery to install the desired change.

        A mutation is base-editing-amenable if and only if the
        *correction* (alt -> ref) is a transition:
          - C>T or G>A correction -> CBE
          - A>G or T>C correction -> ABE

        Transversions, insertions, deletions, and complex changes
        **cannot** be corrected by current base editors.

        Parameters
        ----------
        mutation : Mutation
            The mutation to evaluate.

        Returns
        -------
        bool
            ``True`` if a CBE or ABE can correct this mutation.

        References
        ----------
        - Komor et al., "Programmable editing of a target base in genomic
          DNA without double-stranded DNA cleavage," Nature, 2016.
        - Gaudelli et al., "Programmable base editing of A*T to G*C in
          genomic DNA without DNA cleavage," Nature, 2017.
        """
        return mutation.mutation_type in ("transition_CT", "transition_AG")

    def prime_editing_amenable(self, mutation: Mutation) -> bool:
        """Can this mutation be corrected by prime editing?

        Prime editing uses a Cas9 nickase fused to a reverse
        transcriptase (RT).  The prime editing guide RNA (pegRNA)
        encodes the desired edit in a 3' extension.  After nicking the
        target strand, the RT uses the pegRNA extension as a template
        to synthesize a DNA flap containing the edit, which is then
        incorporated by the cell's DNA repair machinery.

        Prime editing can install:
          - Any single-nucleotide substitution (transitions AND
            transversions)
          - Small insertions (up to ~44 bp demonstrated reliably)
          - Small deletions (up to ~80 bp demonstrated reliably)
          - Combinations of the above

        The practical size limit for the edit is roughly **50 bp** of
        net sequence change — beyond this, pegRNA stability and RT
        processivity become limiting.

        Prime editing does NOT introduce a DSB, making it much safer
        than nuclease + HDR, especially in iPSCs.

        Parameters
        ----------
        mutation : Mutation
            The mutation to evaluate.

        Returns
        -------
        bool
            ``True`` if prime editing can correct this mutation (any
            small change <50 bp).

        References
        ----------
        - Anzalone et al., "Search-and-replace genome editing without
          double-strand breaks or donor DNA," Nature, 2019.
        - Chen et al., "Enhanced prime editing systems by manipulating
          cellular determinants of editing outcomes," Cell, 2021.

        Notes
        -----
        Even though prime editing *can* correct transitions that base
        editors also handle, base editing is generally preferred for
        transitions because:
          1. Base editing efficiency is often higher (30-70% vs 5-50%).
          2. Base editing produces fewer indel byproducts.
          3. The technology is more mature and widely validated.
        """
        # Prime editing can handle any small mutation.
        # We define "small" as total edit size < 50 bp.
        # For SNVs: always amenable.
        # For indels: check size.

        if mutation.mutation_type in (
            "transition_CT", "transition_AG", "transversion"
        ):
            # Any single-nucleotide substitution: prime-editable
            return True

        # For insertions/deletions/complex: check the size of the change.
        # The "size" of the edit is max(len(ref), len(alt)) — the larger
        # of the two alleles determines the length of the pegRNA extension
        # needed.
        ref_len = len(mutation.ref_allele.replace("-", ""))
        alt_len = len(mutation.alt_allele.replace("-", ""))
        edit_size = max(ref_len, alt_len)

        # Practical limit for prime editing: ~50 bp of edit
        # Source: Anzalone et al., Nature, 2019; empirical observations
        # from multiple labs suggest reliable editing up to ~40-50 bp.
        MAX_PRIME_EDIT_SIZE_BP = 50
        return edit_size <= MAX_PRIME_EDIT_SIZE_BP

    def hdr_required(self, mutation: Mutation) -> bool:
        """Does this mutation *require* HDR-based correction?

        A mutation requires HDR if it cannot be corrected by any DSB-free
        approach (base editing or prime editing).  In practice, this
        means:
          - Transversions that are too large for prime editing
          - Insertions > ~50 bp
          - Deletions > ~80 bp
          - Complex structural variants

        When HDR is required, MOSAIC must plan for:
          - A DSB at or near the mutation site
          - A donor template (cssDNA, ssODN, plasmid, or AAV)
          - Higher screening burden (due to competing NHEJ/MMEJ)
          - p53-related toxicity in iPSCs

        Parameters
        ----------
        mutation : Mutation
            The mutation to evaluate.

        Returns
        -------
        bool
            ``True`` if neither base editing nor prime editing can
            correct this mutation, forcing an HDR-based approach.

        Biology
        -------
        HDR is the "last resort" correction strategy because it
        introduces a DSB, which:
          1. Activates p53, causing apoptosis or G1 arrest in iPSCs
             (Ihry et al., Nature Medicine, 2018).
          2. Can lead to large deletions (>100 bp, sometimes >10 kb)
             at the cut site (Kosicki et al., Nature Biotechnology, 2018).
          3. Produces NHEJ/MMEJ outcomes in 40-70% of alleles, meaning
             only a fraction of surviving cells will have the desired
             HDR correction.
          4. When two DSBs are made simultaneously, chromosomal
             translocations can occur.
        """
        return (
            not self.base_editing_amenable(mutation)
            and not self.prime_editing_amenable(mutation)
        )

    def describe_mutation(self, mutation: Mutation) -> str:
        """Return a human-readable description of the mutation and its
        correction options.

        Parameters
        ----------
        mutation : Mutation
            The mutation to describe.

        Returns
        -------
        str
            Multi-line description string.
        """
        lines = [
            f"Mutation: {mutation.name}",
            f"  Location:        Exon {mutation.exon_number}, "
            f"position {mutation.position:,}",
            f"  Change:          {mutation.ref_allele} > {mutation.alt_allele}",
            f"  Classification:  {mutation.mutation_type}",
            f"  Base-editable:   {'Yes' if self.base_editing_amenable(mutation) else 'No'}",
        ]
        if self.base_editing_amenable(mutation):
            if mutation.mutation_type == "transition_CT":
                lines.append(
                    f"    -> CBE can correct (C>T on edited strand)"
                )
            elif mutation.mutation_type == "transition_AG":
                lines.append(
                    f"    -> ABE can correct (A>G on edited strand)"
                )
        lines.append(
            f"  Prime-editable:  "
            f"{'Yes' if self.prime_editing_amenable(mutation) else 'No'}"
        )
        lines.append(
            f"  Requires HDR:    "
            f"{'Yes' if self.hdr_required(mutation) else 'No'}"
        )
        return "\n".join(lines)


# =========================================================================
# Quick self-test
# =========================================================================

if __name__ == "__main__":
    mc = MutationClassifier()

    # Example: compound heterozygous mutations in a large gene
    # Mutation 1: G>A in exon 40 (transition — ABE can correct A>G)
    m1 = Mutation(
        exon_number=40,
        position=5_000_000,
        ref_allele="G",
        alt_allele="A",
        name="c.5000G>A (exon 40)",
    )

    # Mutation 2: 2-bp deletion in exon 60
    m2 = Mutation(
        exon_number=60,
        position=8_000_000,
        ref_allele="CT",
        alt_allele="-",
        name="c.8000delCT (exon 60)",
    )

    for m in [m1, m2]:
        print(mc.describe_mutation(m))
        print()
