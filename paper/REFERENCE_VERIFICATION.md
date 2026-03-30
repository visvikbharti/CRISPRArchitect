# Reference Verification Report — v3 Manuscript

**Manuscript:** CRISPRArchitect v3: multi-nuclease decision support for genome editing strategy design with TOPSIS ranking and sensitivity analysis

**Authors:** Vishal Bharti and Debojyoti Chakraborty

**Verification date:** 2026-03-30

**Method:** Each reference verified against PubMed (PMID lookup), DOI resolution, and publisher databases.

---

## Summary

- **Total references:** 20
- **VERIFIED:** 17
- **CORRECTED (this session):** 3
- **Previously known errors (v1 manuscript):** All 5 resolved (4 removed, 1 reformatted)

---

## Corrections Applied (2026-03-30)

### Ref 9 — Arbab et al. (2020): WRONG JOURNAL
- **Was:** *Nature* **584**, 268–276 (2020) — **FABRICATED journal/volume/pages**
- **Corrected to:** *Cell* **182**, 463–480.e30 (2020)
- **PMID:** 32533916 | **DOI:** 10.1016/j.cell.2020.05.037

### Ref 14 — enFnCas9 paper: WRONG TITLE, WRONG FIRST AUTHOR, WRONG FORMAT
- **Was:** Chakraborty, D. *et al.* (Hirano, S. *et al.*) "enFnCas9: an engineered *Francisella novicida* Cas9..." *Nat. Commun.* **15**, 1–14 (2024) — **FABRICATED title**
- **Corrected to:** Acharya, S. *et al.* "PAM-flexible Engineered FnCas9 variants for robust and ultra-precise genome editing and diagnostics." *Nat. Commun.* **15**, 5471 (2024)
- **PMID:** 38942756 | **DOI:** 10.1038/s41467-024-49233-w
- **Notes:** First author is Acharya S (not Chakraborty D, who is corresponding/last author). Article number 5471, not pages 1–14.

### Ref 13 — Walton et al. (2020): STYLE FIX
- **Was:** *Science* **368**, eaba8853 (2020)
- **Corrected to:** *Science* **368**, 290–296 (2020)
- **PMID:** 32217751 | **DOI:** 10.1126/science.aba8853
- **Notes:** eLocator technically valid but inconsistent with other refs using page numbers.

---

## Full Verification Table

| Ref | Citation | PMID | Status |
|-----|----------|------|--------|
| 1 | Komor AC et al. Nature 533, 420-424 (2016) | 27096365 | VERIFIED |
| 2 | Gaudelli NM et al. Nature 551, 464-471 (2017) | 29160308 | VERIFIED |
| 3 | Anzalone AV et al. Nature 576, 149-157 (2019) | 31634902 | VERIFIED |
| 4 | Paquet D et al. Nature 533, 125-129 (2016) | 27120160 | VERIFIED |
| 5 | Ihry RJ et al. Nat Med 24, 939-946 (2018) | 29892062 | VERIFIED |
| 6 | Hwang C-L & Yoon K. Springer-Verlag (1981) | Book | VERIFIED |
| 7 | Iyer S et al. CRISPR J 5, 685-701 (2022) | 36070530 | VERIFIED |
| 8 | Richards S et al. Genet Med 17, 405-424 (2015) | 25741868 | VERIFIED |
| 9 | Arbab M et al. Cell 182, 463-480 (2020) | 32533916 | **CORRECTED** |
| 10 | Rees HA & Liu DR. Nat Rev Genet 19, 770-788 (2018) | 30323312 | VERIFIED |
| 11 | Richter MF et al. Nat Biotechnol 38, 883-891 (2020) | 32433547 | VERIFIED |
| 12 | Nishimasu H et al. Science 361, 1259-1262 (2018) | 30166441 | VERIFIED |
| 13 | Walton RT et al. Science 368, 290-296 (2020) | 32217751 | **CORRECTED** |
| 14 | Acharya S et al. Nat Commun 15, 5471 (2024) | 38942756 | **CORRECTED** |
| 15 | Doench JG et al. Nat Biotechnol 34, 184-191 (2016) | 26780180 | VERIFIED |
| 16 | Hsu PD et al. Nat Biotechnol 31, 827-832 (2013) | 23873081 | VERIFIED |
| 17 | Koblan LW et al. Nat Biotechnol 36, 843-846 (2018) | 29813047 | VERIFIED |
| 18 | den Dunnen JT et al. Hum Mutat 37, 564-569 (2016) | 26931183 | VERIFIED |
| 19 | Elliott B et al. Mol Cell Biol 18, 93-101 (1998) | 9418857 | VERIFIED |
| 20 | Concordet J-P & Haeussler M. NAR 46, W242-W245 (2018) | 29762716 | VERIFIED |

---

## Resolution of Previously Flagged Errors (v1 Manuscript)

The v1 REFERENCE_VERIFICATION.md (2026-03-20) identified 7 errors in 22 references, including 3 fabricated author names. Status in v3:

| v1 Error | Resolution |
|----------|------------|
| Ref 13 (Paquet): "Bhatt S" wrong author | Resolved — v3 uses "Paquet D *et al.*" |
| Ref 14 (Symington): wrong year/volume | Resolved — reference removed from v3 |
| Ref 15 (Cejka): wrong journal | Resolved — reference removed from v3 |
| Ref 19 (Kan): "Taber S" wrong author | Resolved — reference removed from v3 |
| Ref 22 (Aymard): "Clapier P" wrong author | Resolved — reference removed from v3 |

---

## Code-Level Reference Updates (2026-03-30)

The enFnCas9 reference was also corrected in `utils/constants.py` across all occurrences:
- "Chakraborty et al., Nat Commun, 2024" → "Acharya et al., Nat Commun, 2024 (15:5471)"
- Applied to NUCLEASE_PARAMS["enFnCas9"], BASE_EDITOR_PROFILES["ABE8e-enFnCas9"], and BASE_EDITOR_PROFILES["BE4max-enFnCas9"]
