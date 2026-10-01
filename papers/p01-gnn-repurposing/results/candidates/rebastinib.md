# Candidate dossier: rebastinib

**Generated:** 2026-06-05  |  **Pipeline:** IgA Nephropathy GNN drug repurposing

---

## Pipeline scores

| Metric | Value |
|--------|-------|
| Final rank | 10 |
| Priority score | 6 / 9 |
| Network proximity rank | 55 / 1,989 (z = -2.091) |
| Open Targets genetic rank | 10 / 3,005 (score 0.5052) |
| GNN link prediction rank (full R-GCN) | 1,362 / 7,957 (score 2.38e-08) |
| Hit alignment | Hit 2 — autoantibody/B-cell (LYN/BAFF/APRIL axis) |
| GWAS target | LYN |
| Clinical phase | Phase 1/Phase 2 |
| Mechanism of action | Bcr-Abl kinase inhibitor | TIE tyrosine kinase inhibitor | VEGFR inhibitor |

---

## GWAS locus context

LYN encodes a Src-family tyrosine kinase at the apex of the B-cell receptor signalling cascade. Kiryluk et al. 2023 identified LYN as a genome-wide significant IgAN susceptibility locus (rank 22 of 3,005 by Open Targets genetic association score, score 0.505). Downstream cascade: LYN -> SYK -> BTK -> NF-kB -> BAFF/APRIL. Fostamatinib (SYK inhibitor) is in Phase 2 for IgAN, confirming pathway validity. Drugs targeting LYN directly have not been tested in IgAN.

---

## CKD safety profile

| Parameter | Assessment |
|-----------|------------|
| eGFR adjustment | monitor |
| Nephrotoxicity risk | low |
| Hypertension risk | low |

Phase 1/2 only; limited safety data. No renal contraindication identified. Author should verify most recent trial safety reporting.


---

## PubMed evidence (drug + IgA Nephropathy/kidney fibrosis)

No PubMed results found for this drug + IgA Nephropathy/kidney fibrosis query.

---

## Open Targets known indications

| Disease | Phase |
|---------|-------|
| breast cancer | PHASE_1 |
| chronic myelogenous leukemia | PHASE_1_2 |
| neoplasm | PHASE_1_2 |

IgA Nephropathy not listed as a known indication — confirms novelty of this repurposing prediction.

---

## ClinicalTrials.gov active and completed trials

| NCT ID | Title | Status | Phase | Condition | Start |
|--------|-------|--------|-------|-----------|-------|
| [NCT00827138](https://clinicaltrials.gov/study/NCT00827138) | Study Safety and Preliminary Efficacy of DCC-2036 in Patients With Leukemias (Ph+ CML With T315I Mutation) | COMPLETED | PHASE1 | Chronic Myeloid Leukemia | 2009-03 |
| [NCT03717415](https://clinicaltrials.gov/study/NCT03717415) | A Study of Rebastinib (DCC-2036) in Combination With Carboplatin in Patients With Advanced or Metastatic Solid Tumors | TERMINATED | PHASE1 | Locally Advanced or Metastatic Solid Tumor | 2019-01-02 |
| [NCT03601897](https://clinicaltrials.gov/study/NCT03601897) | A Phase 1b/2 Study of Rebastinib (DCC-2036) in Combination With Paclitaxel in Patients With Advanced or Metastatic Solid | TERMINATED | PHASE1 | Locally Advanced or Metastatic Solid Tumor | 2018-10-25 |
| [NCT02824575](https://clinicaltrials.gov/study/NCT02824575) | Rebastinib Plus Antitubulin Therapy With Paclitaxel or Eribulin in Metastatic Breast Cancer | TERMINATED | PHASE1 | Breast Cancer; Breast Adenocarcinoma; Human Epidermal Growth Factor 2 Negative C | 2016-07 |

No kidney or IgA Nephropathy trials found for this drug — confirms this is a novel repurposing prediction.

---

## Next steps for validation

1. **In silico:** protein-ligand co-folding against LYN (Boltz-2 results for the five top-ranked candidates are in results/boltz/)
2. **In vitro:** Test in primary human mesangial cells, or the mouse mesangial cell line MES-13
3. **In vivo:** ddY mouse model (spontaneous IgAN) — measure proteinuria, UPCR, Gd-IgA1
4. **Clinical:** IgAN-specific Phase 2 trial (primary endpoint: 40% UPCR reduction at 12 months)