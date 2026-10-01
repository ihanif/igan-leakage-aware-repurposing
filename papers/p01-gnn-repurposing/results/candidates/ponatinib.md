# Candidate dossier: ponatinib

**Generated:** 2026-06-05  |  **Pipeline:** IgA Nephropathy GNN drug repurposing

---

## Pipeline scores

| Metric | Value |
|--------|-------|
| Final rank | 12 |
| Priority score | 5 / 9 |
| Network proximity rank | 59 / 1,989 (z = -1.588) |
| Open Targets genetic rank | 9 / 3,005 (score 0.5052) |
| GNN link prediction rank (full R-GCN) | 641 / 7,957 (score 5.22e-07) |
| Hit alignment | Hit 2 — autoantibody/B-cell (LYN/BAFF/APRIL axis) |
| GWAS target | LYN |
| Clinical phase | Launched |
| Mechanism of action | Bcr-Abl kinase inhibitor | FLT3 inhibitor | PDGFR tyrosine kinase receptor inhibitor |

---

## GWAS locus context

LYN encodes a Src-family tyrosine kinase at the apex of the B-cell receptor signalling cascade. Kiryluk et al. 2023 identified LYN as a genome-wide significant IgAN susceptibility locus (rank 22 of 3,005 by Open Targets genetic association score, score 0.505). Downstream cascade: LYN -> SYK -> BTK -> NF-kB -> BAFF/APRIL. Fostamatinib (SYK inhibitor) is in Phase 2 for IgAN, confirming pathway validity. Drugs targeting LYN directly have not been tested in IgAN.

---

## CKD safety profile

| Parameter | Assessment |
|-----------|------------|
| eGFR adjustment | monitor |
| Nephrotoxicity risk | moderate |
| Hypertension risk | high |

Arterial occlusive events (14% of patients), hypertension (67%), hepatotoxicity. Significant cardiovascular burden inappropriate for IgAN without careful risk-benefit. Flag for specialist cardiovascular review before IgAN use.


---

## PubMed evidence (drug + IgA Nephropathy/kidney fibrosis)

**Cross-Domain Text Mining to Predict Adverse Events from Tyrosine Kinase Inhibitors for Chronic Myeloid Leukemia.**
*Cancers* (2022) | PMID: 36230609
Tyrosine kinase inhibitors (TKIs) are prescribed for chronic myeloid leukemia (CML) and some other cancers. The objective was to predict and rank TKI-related adverse events (AEs), including under-reported or preclinical AEs, using novel text mining. First, k-means clustering of 2575 clinical CML TKI...

**Identification of biomarkers and drug repurposing candidates based on an immune-, inflammation- and membranous glomerulonephritis-associated triplets network for membranous glomerulonephritis.**
*BMC medical genomics* (2020) | PMID: 31910852
Membranous glomerulonephritis (MGN) is a common kidney disease. Despite many evidences support that many immune and inflammation-related genes could serve as effective biomarkers and treatment targets for MGN patients, the potential associations among MGN-, immune- and inflammation-related genes hav...

---

## Open Targets known indications

| Disease | Phase |
|---------|-------|
| small cell lung carcinoma | PHASE_2 |
| medullary thyroid gland carcinoma | PHASE_2 |
| non-small cell lung carcinoma | PHASE_2_3 |
| leukemia | PHASE_2 |
| acute lymphoblastic leukemia | PHASE_3 |
| endometrial neoplasm | UNKNOWN |
| lung adenocarcinoma | PHASE_2 |
| glioblastoma multiforme | PHASE_2 |
| acute myeloid leukemia | PHASE_2 |
| chronic myelogenous leukemia | PHASE_3 |

IgA Nephropathy not listed as a known indication — confirms novelty of this repurposing prediction.

---

## ClinicalTrials.gov active and completed trials

| NCT ID | Title | Status | Phase | Condition | Start |
|--------|-------|--------|-------|-----------|-------|
| [NCT06175702](https://clinicaltrials.gov/study/NCT06175702) | Treatment Protocol for Newky Diagnosed Adult Ph Positive ALL | NOT_YET_RECRUITING | N/A | Lymphoblastic Leukemia; Philadelphia-Positive ALL; Adult ALL | 2023-12-25 |
| [NCT04530565](https://clinicaltrials.gov/study/NCT04530565) | Testing the Use of Steroids and Tyrosine Kinase Inhibitors With Blinatumomab or Chemotherapy for Newly Diagnosed BCR-ABL | ACTIVE_NOT_RECRUITING | PHASE3 | B Acute Lymphoblastic Leukemia With t(9;22)(q34.1;q11.2); BCR-ABL1 | 2021-01-25 |
| [NCT01641107](https://clinicaltrials.gov/study/NCT01641107) | Phase II Front-line Ponatinib in Adult Philadelphia+/BCR-ABL+ Acute Lymphoblastic Leukemia. | COMPLETED | PHASE2 | Philadelphia Positive; BCR-ABL Positive; Acute Lymphoblastic Leukemia | 2014-12-04 |
| [NCT07188428](https://clinicaltrials.gov/study/NCT07188428) | Treatment Patterns and Outcomes Among Patients With Chronic Myeloid Leukemia (CML) in All Lines of Treatment | COMPLETED | N/A | Leukemia, Chronic Myeloid | 2024-02-27 |
| [NCT02265341](https://clinicaltrials.gov/study/NCT02265341) | Ponatinib Hydrochloride in Treating Patients With Advanced Biliary Cancer With FGFR2 Fusions | COMPLETED | PHASE2 | Malignant Hepatobiliary Neoplasm | 2014-12 |
| [NCT02829840](https://clinicaltrials.gov/study/NCT02829840) | Dose-Escalation Study of Ponatinib, a FLT3 Inhibitor, With and Without Combination of 5-Azacytidine, in Patients With FL | WITHDRAWN | PHASE1 | Leukemia; FLT3-Mutated Acute Myeloid Leukemia; FLT3-Mutated High-Risk Myelodyspl | 2016-09 |
| [NCT03589326](https://clinicaltrials.gov/study/NCT03589326) | A Study of Ponatinib Versus Imatinib in Adults With Acute Lymphoblastic Leukemia | ACTIVE_NOT_RECRUITING | PHASE3 | Philadelphia Chromosome Positive Acute Lymphoblastic Leukemia (Ph+ALL) | 2018-10-04 |
| [NCT02981784](https://clinicaltrials.gov/study/NCT02981784) | Comparative Evaluation of Results of Allogeneic Hematopoietic Stem Cells Versus Ponatinib in CML Patients Carrying a Mut | COMPLETED | N/A | Leukemia | 2000-01 |

No kidney or IgA Nephropathy trials found for this drug — confirms this is a novel repurposing prediction.

---

## Next steps for validation

1. **In silico:** protein-ligand co-folding against LYN (Boltz-2 results for the five top-ranked candidates are in results/boltz/)
2. **In vitro:** Test in primary human mesangial cells, or the mouse mesangial cell line MES-13
3. **In vivo:** ddY mouse model (spontaneous IgAN) — measure proteinuria, UPCR, Gd-IgA1
4. **Clinical:** IgAN-specific Phase 2 trial (primary endpoint: 40% UPCR reduction at 12 months)