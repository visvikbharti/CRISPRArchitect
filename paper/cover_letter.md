Date: [Date of submission]

Dr. [Editor Name]
Editor-in-Chief
PLOS Computational Biology

Dear Editor,

We are pleased to submit our manuscript entitled **"A Monte Carlo framework for predicting HDR gene conversion outcomes and optimizing multi-site genome editing strategies"** for consideration as a Research Article in *PLOS Computational Biology*.

**The problem.** Correcting multiple pathogenic mutations in large multi-exon genes using CRISPR technologies requires choosing among several editing modalities — homology-directed repair (HDR), base editing, prime editing, and hybrid approaches. This decision is particularly consequential in human induced pluripotent stem cells (iPSCs), where double-strand breaks activate p53-mediated apoptosis. Despite the importance of this decision, no computational tool currently helps researchers navigate this multi-modality, multi-site strategy space.

**What we did.** We developed CRISPRArchitect, a computational framework with two core components. ConversionSim is a Monte Carlo simulator that models the sequential steps of HDR — end resection, RAD51 filament formation, strand invasion, and synthesis-dependent strand annealing — to predict gene conversion tract length distributions. MOSAIC is a strategy optimizer that enumerates all feasible editing strategies for a given multi-site correction problem and ranks them by efficiency, safety, time, and cost.

**Key findings.** We validated ConversionSim against four independent published datasets. The model correctly predicts the approximately two-fold enhancement of circular over linear single-stranded DNA donors (predicted 2.07-fold versus observed 1.5–2.1-fold; Iyer et al., *CRISPR Journal*, 2022) and the enhancement from staggered nuclease cuts (predicted 1.82-fold versus observed 1.4–2.8-fold; Chauhan et al., *PNAS*, 2023). The validation also revealed a systematic discrepancy for short oligonucleotide donors, attributable to the single-strand template repair (SSTR) pathway — a finding that delineates the model's applicability and highlights the mechanistic distinction between long-donor SDSA and short-donor SSTR pathways. We benchmarked MOSAIC against strategies employed in 14 published genome editing studies, achieving 71% concordance with expert decisions and 100% concordance for base editing applications.

**Why PLOS Computational Biology.** We believe this work is well suited to *PLOS Computational Biology* for several reasons. First, ConversionSim represents a new class of Monte Carlo simulation for genome editing biology, combining known biochemical parameters into a predictive framework that produces testable quantitative predictions. Second, the work directly addresses a practical gap — no existing tool compares editing modalities for multi-site correction — with clear implications for the gene therapy community. Third, we provide complete open-source software (https://github.com/visvikbharti/CRISPRArchitect), a web interface, Docker container, and full validation scripts, aligning with the journal's commitment to reproducibility and open science.

**Competing interests.** The authors declare no competing interests.

**Ethical statement.** This study is entirely computational and did not involve human subjects, animal experiments, or biological materials.

We confirm that this manuscript has not been published elsewhere and is not under consideration by another journal. All authors have approved the manuscript for submission.

We suggest the following potential reviewers with relevant expertise:
1. [Reviewer suggestion 1 — to be filled by authors]
2. [Reviewer suggestion 2 — to be filled by authors]
3. [Reviewer suggestion 3 — to be filled by authors]

Thank you for considering our manuscript.

Sincerely,

Vishal Bharti and Debojyoti Chakraborty

CSIR-Institute of Genomics and Integrative Biology (IGIB)
New Delhi, 110025, India

Correspondence: vishalvikashbharti@gmail.com; debojyoti.chakraborty@igib.in
