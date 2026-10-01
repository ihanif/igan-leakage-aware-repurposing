# Candidate dossier: bosutinib

**Generated:** 2026-06-05  |  **Pipeline:** IgA Nephropathy GNN drug repurposing

---

## Pipeline scores

| Metric | Value |
|--------|-------|
| Final rank | 2 |
| Priority score | 8 / 9 |
| Network proximity rank | 56 / 1,989 (z = -2.061) |
| Open Targets genetic rank | 7 / 3,005 (score 0.5052) |
| GNN link prediction rank (full R-GCN) | 253 / 7,957 (score 7.11e-06) |
| Hit alignment | Hit 2 — autoantibody/B-cell (LYN/BAFF/APRIL axis) |
| GWAS target | LYN |
| Clinical phase | Launched |
| Mechanism of action | Abl kinase inhibitor | Bcr-Abl kinase inhibitor | SRC inhibitor |

---

## GWAS locus context

LYN encodes a Src-family tyrosine kinase at the apex of the B-cell receptor signalling cascade. Kiryluk et al. 2023 identified LYN as a genome-wide significant IgAN susceptibility locus (rank 22 of 3,005 by Open Targets genetic association score, score 0.505). Downstream cascade: LYN -> SYK -> BTK -> NF-kB -> BAFF/APRIL. Fostamatinib (SYK inhibitor) is in Phase 2 for IgAN, confirming pathway validity. Drugs targeting LYN directly have not been tested in IgAN.

---

## CKD safety profile

| Parameter | Assessment |
|-----------|------------|
| eGFR adjustment | dose_reduce_egfr<30 |
| Nephrotoxicity risk | low |
| Hypertension risk | low |

Reduce dose to 200 mg/day for eGFR < 30 (FDA label). Hepatic metabolism; renal excretion minimal. No CKD contraindication above eGFR 30.


---

## PubMed evidence (drug + IgA Nephropathy/kidney fibrosis)

**Cross-Domain Text Mining to Predict Adverse Events from Tyrosine Kinase Inhibitors for Chronic Myeloid Leukemia.**
*Cancers* (2022) | PMID: 36230609
Tyrosine kinase inhibitors (TKIs) are prescribed for chronic myeloid leukemia (CML) and some other cancers. The objective was to predict and rank TKI-related adverse events (AEs), including under-reported or preclinical AEs, using novel text mining. First, k-means clustering of 2575 clinical CML TKI...

**Identification of biomarkers and drug repurposing candidates based on an immune-, inflammation- and membranous glomerulonephritis-associated triplets network for membranous glomerulonephritis.**
*BMC medical genomics* (2020) | PMID: 31910852
Membranous glomerulonephritis (MGN) is a common kidney disease. Despite many evidences support that many immune and inflammation-related genes could serve as effective biomarkers and treatment targets for MGN patients, the potential associations among MGN-, immune- and inflammation-related genes hav...

**Neutrophil Fc&#x3b3;RIIA promotes IgG-mediated glomerular neutrophil capture via Abl/Src kinases.**
*The Journal of clinical investigation* (2017) | PMID: 28891817
The kidney glomerular capillaries are frequent sites of immune complex deposition and subsequent neutrophil accumulation in post-infectious and rapidly progressive glomerulonephritis. However, the mechanisms of neutrophil recruitment remain enigmatic, and there is no targeted therapeutic to avert th...

---

## Open Targets known indications

| Disease | Phase |
|---------|-------|
| pancreatic carcinoma | PHASE_1 |
| leukemia | PHASE_2 |
| chronic myelogenous leukemia | PHASE_3 |
| Lewy body dementia | PHASE_2 |
| breast cancer | APPROVAL |
| glioblastoma multiforme | PHASE_2 |
| Autosomal dominant polycystic kidney disease | PHASE_2 |
| myeloid leukemia | APPROVAL |
| breast carcinoma | PHASE_1 |
| amyotrophic lateral sclerosis | PHASE_1_2 |

---

## ClinicalTrials.gov active and completed trials

| NCT ID | Title | Status | Phase | Condition | Start |
|--------|-------|--------|-------|-----------|-------|
| [NCT01374139](https://clinicaltrials.gov/study/NCT01374139) | Bioequivalence And Food Effect Study Of Bosutinib In Healthy Subjects | COMPLETED | PHASE1 | Philadelphia Chromosome Positive (Ph+) Chronic Myeloid Leukemia (CML) | 2011-08 |
| [NCT07188428](https://clinicaltrials.gov/study/NCT07188428) | Treatment Patterns and Outcomes Among Patients With Chronic Myeloid Leukemia (CML) in All Lines of Treatment | COMPLETED | N/A | Leukemia, Chronic Myeloid | 2024-02-27 |
| [NCT03023319](https://clinicaltrials.gov/study/NCT03023319) | Bosutinib in Combination With Pemetrexed in Patients With Selected Metastatic Solid Tumors | COMPLETED | PHASE1 | Carcinoma, Non-Small-Cell Lung; Mesothelioma; Bladder Cancer | 2019-12-10 |
| [NCT00914121](https://clinicaltrials.gov/study/NCT00914121) | Study Evaluating The Effect of Bosutinib (SKI-606) On Cardiac Repolarization (Rhythms Of The Heart) | COMPLETED | PHASE1 | Healthy Subjects | 2009-06 |
| [NCT03106779](https://clinicaltrials.gov/study/NCT03106779) | Study of Efficacy of CML-CP Patients Treated With ABL001 Versus Bosutinib, Previously Treated With 2 or More TKIs | COMPLETED | PHASE3 | Chronic Myelogenous Leukemia | 2017-10-26 |
| [NCT07091019](https://clinicaltrials.gov/study/NCT07091019) | A Real-world Chronic Myelogenous Leukemia (CML) Patient Disease Registry to Describe Patient Experience and Clinical Out | NOT_YET_RECRUITING | N/A | Chronic Myelogenous Leukemia - Chronic Phase | 2025-11-15 |
| [NCT00574873](https://clinicaltrials.gov/study/NCT00574873) | Compare Bosutinib To Imatinib In Subjects With Newly Diagnosed Chronic Phase Philadelphia Chromosome Positive CML | COMPLETED | PHASE3 | Chronic Myeloid Leukemia | 2008-02-05 |
| [NCT06423911](https://clinicaltrials.gov/study/NCT06423911) | Study of Olverembatinib (HQP1351) in Patients With CP-CML | RECRUITING | PHASE3 | Chronic Myeloid Leukemia; CML; CML, Chronic Phase | 2024-02-05 |

No kidney or IgA Nephropathy trials found for this drug — confirms this is a novel repurposing prediction.

---

## Next steps for validation

1. **In silico:** protein-ligand co-folding against LYN (Boltz-2 results for the five top-ranked candidates are in results/boltz/)
2. **In vitro:** Test in primary human mesangial cells, or the mouse mesangial cell line MES-13
3. **In vivo:** ddY mouse model (spontaneous IgAN) — measure proteinuria, UPCR, Gd-IgA1
4. **Clinical:** IgAN-specific Phase 2 trial (primary endpoint: 40% UPCR reduction at 12 months)