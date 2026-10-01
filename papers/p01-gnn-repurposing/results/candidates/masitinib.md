# Candidate dossier: masitinib

**Generated:** 2026-06-05  |  **Pipeline:** IgA Nephropathy GNN drug repurposing

---

## Pipeline scores

| Metric | Value |
|--------|-------|
| Final rank | 3 |
| Priority score | 8 / 9 |
| Network proximity rank | 53 / 1,989 (z = -2.430) |
| Open Targets genetic rank | 14 / 3,005 (score 0.5052) |
| GNN link prediction rank (full R-GCN) | 2,175 / 7,957 (score 1.59e-08) |
| Hit alignment | Hit 2 — autoantibody/B-cell (LYN/BAFF/APRIL axis) |
| GWAS target | LYN |
| Clinical phase | Launched |
| Mechanism of action | KIT inhibitor | PDGFR tyrosine kinase receptor inhibitor | SRC inhibitor |

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

No CKD-specific dose adjustment published. Tested in ALS, MS, pancreatic cancer. Mild edema in ~15% of patients. Renal safety data from non-oncology trials favourable.


---

## PubMed evidence (drug + IgA Nephropathy/kidney fibrosis)

No PubMed results found for this drug + IgA Nephropathy/kidney fibrosis query.

---

## Open Targets known indications

| Disease | Phase |
|---------|-------|
| Malignant Pancreatic Neoplasm | PHASE_3 |
| pancreatic carcinoma | PHASE_3 |
| amyotrophic lateral sclerosis | PHASE_3 |
| metastatic melanoma | PHASE_3 |
| asthma | PHASE_3 |
| neoplasm | APPROVAL |
| Alzheimer disease | PHASE_3 |
| multiple sclerosis | PHASE_3 |
| rheumatoid arthritis | PHASE_2_3 |
| pancreatic neoplasm | PHASE_3 |

IgA Nephropathy not listed as a known indication — confirms novelty of this repurposing prediction.

---

## ClinicalTrials.gov active and completed trials

| NCT ID | Title | Status | Phase | Condition | Start |
|--------|-------|--------|-------|-----------|-------|
| [NCT00814073](https://clinicaltrials.gov/study/NCT00814073) | Masitinib in Severe Indolent or Smoldering Systemic Mastocytosis | COMPLETED | PHASE3 | Indolent Systemic Mastocytosis | 2008-12 |
| [NCT03556956](https://clinicaltrials.gov/study/NCT03556956) | Masitinib in Combination With FOLFIRI in Third or Fourth Line of Treatment of Patients With Metastatic Colorectal Cancer | COMPLETED | PHASE2 | Colorectal Cancer Metastatic | 2015-07 |
| [NCT01266369](https://clinicaltrials.gov/study/NCT01266369) | Masitinib in Patients With Mastocytosis With Handicap and Bearing the D816V Mutation | COMPLETED | PHASE2 | Mastocytosis | 2007-02 |
| [NCT04333108](https://clinicaltrials.gov/study/NCT04333108) | Masitinib in Severe Indolent or Smoldering Systemic Mastocytosis Unresponsive to Optimal Symptomatic Treatment | UNKNOWN | PHASE3 | Indolent Systemic Mastocytosis | 2020-07-01 |
| [NCT01449162](https://clinicaltrials.gov/study/NCT01449162) | Masitinib in Treatment of Patients With Severe Persistent Asthma Treated With Oral Corticosteroids | COMPLETED | PHASE3 | Asthma | 2011-01 |
| [NCT03127267](https://clinicaltrials.gov/study/NCT03127267) | Efficacy and Safety of Masitinib Versus Placebo in the Treatment of ALS Patients | RECRUITING | PHASE3 | Amyotrophic Lateral Sclerosis | 2021-02-02 |
| [NCT02588677](https://clinicaltrials.gov/study/NCT02588677) | Masitinib in Combination With Riluzole for the Treatment of Patients Suffering From Amyotrophic Lateral Sclerosis (ALS) | COMPLETED | PHASE2 | Amyotrophic Lateral Sclerosis (ALS) | 2013-04 |
| [NCT01280565](https://clinicaltrials.gov/study/NCT01280565) | Masitinib in Non-Resectable or Metastatic Stage 3/4 Melanoma Carrying a Mutation in the Juxta Membrane Domain of c-Kit | TERMINATED | PHASE3 | Metastatic Melanoma | 2011-01 |

No kidney or IgA Nephropathy trials found for this drug — confirms this is a novel repurposing prediction.

---

## Next steps for validation

1. **In silico:** protein-ligand co-folding against LYN (Boltz-2 results for the five top-ranked candidates are in results/boltz/)
2. **In vitro:** Test in primary human mesangial cells, or the mouse mesangial cell line MES-13
3. **In vivo:** ddY mouse model (spontaneous IgAN) — measure proteinuria, UPCR, Gd-IgA1
4. **Clinical:** IgAN-specific Phase 2 trial (primary endpoint: 40% UPCR reduction at 12 months)