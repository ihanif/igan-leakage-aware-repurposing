# Candidate dossier: dasatinib

**Generated:** 2026-06-05  |  **Pipeline:** IgA Nephropathy GNN drug repurposing

---

## Pipeline scores

| Metric | Value |
|--------|-------|
| Final rank | 6 |
| Priority score | 7 / 9 |
| Network proximity rank | 61 / 1,989 (z = -0.913) |
| Open Targets genetic rank | 15 / 3,005 (score 0.5052) |
| GNN link prediction rank (full R-GCN) | 644 / 7,957 (score 5.2e-07) |
| Hit alignment | Hit 2 — autoantibody/B-cell (LYN/BAFF/APRIL axis) |
| GWAS target | LYN |
| Clinical phase | Launched |
| Mechanism of action | Bcr-Abl kinase inhibitor | ephrin inhibitor | KIT inhibitor | PDGFR tyrosine kinase receptor inhibitor | SRC inhibitor | tyrosine kinase inhibitor |

---

## GWAS locus context

LYN encodes a Src-family tyrosine kinase at the apex of the B-cell receptor signalling cascade. Kiryluk et al. 2023 identified LYN as a genome-wide significant IgAN susceptibility locus (rank 22 of 3,005 by Open Targets genetic association score, score 0.505). Downstream cascade: LYN -> SYK -> BTK -> NF-kB -> BAFF/APRIL. Fostamatinib (SYK inhibitor) is in Phase 2 for IgAN, confirming pathway validity. Drugs targeting LYN directly have not been tested in IgAN.

---

## CKD safety profile

| Parameter | Assessment |
|-----------|------------|
| eGFR adjustment | monitor |
| Nephrotoxicity risk | low |
| Hypertension risk | moderate |

No renal dose adjustment required (hepatic metabolism). Hypertension in 11% of patients; monitor BP closely in IgAN. Pleural effusion at oncology doses (100 mg); lower doses (20–50 mg) under investigation for autoimmune indications.


---

## PubMed evidence (drug + IgA Nephropathy/kidney fibrosis)

**Cross-Domain Text Mining to Predict Adverse Events from Tyrosine Kinase Inhibitors for Chronic Myeloid Leukemia.**
*Cancers* (2022) | PMID: 36230609
Tyrosine kinase inhibitors (TKIs) are prescribed for chronic myeloid leukemia (CML) and some other cancers. The objective was to predict and rank TKI-related adverse events (AEs), including under-reported or preclinical AEs, using novel text mining. First, k-means clustering of 2575 clinical CML TKI...

**A Human Conditionally Immortalized Proximal Tubule Epithelial Cell Line as a Novel Model for Studying Senescence and Response to Senolytics.**
*Frontiers in pharmacology* (2022) | PMID: 35422705
Accumulating evidence suggests that senescence of kidney tubule epithelial cells leads to fibrosis. These cells secrete senescence-associated secretory phenotype (SASP) factors that are involved in diverse signaling pathways, influencing kidney fibrosis. Here, we investigated whether our previously ...

**Dasatinib-induced nephrotic syndrome in a patient with chronic myelogenous leukemia: a case report.**
*BMC nephrology* (2019) | PMID: 30845905
Dasatinib is a second-generation tyrosine kinase inhibitor that is indicated for the treatment of patients with chronic myeloid leukemia. Here, we report the case of a man with nephrotic syndrome that was caused by dasatinib.

---

## Open Targets known indications

| Disease | Phase |
|---------|-------|
| liver disease | UNKNOWN |
| gastrointestinal stromal tumor | UNKNOWN |
| chronic myelogenous leukemia | APPROVAL |
| acute myeloid leukemia | UNKNOWN |
| myelodysplastic syndrome | UNKNOWN |
| multiple myeloma | PHASE_2 |
| bladder transitional cell carcinoma | UNKNOWN |
| leukemia | UNKNOWN |
| acute lymphoblastic leukemia | APPROVAL |

IgA Nephropathy not listed as a known indication — confirms novelty of this repurposing prediction.

---

## ClinicalTrials.gov active and completed trials

| NCT ID | Title | Status | Phase | Condition | Start |
|--------|-------|--------|-------|-----------|-------|
| [NCT01660971](https://clinicaltrials.gov/study/NCT01660971) | Gemcitabine Hydrochloride, Dasatinib, and Erlotinib Hydrochloride in Treating Patients With Pancreatic Cancer That Is Me | ACTIVE_NOT_RECRUITING | PHASE1 | Metastatic Pancreatic Adenocarcinoma; Recurrent Pancreatic Carcinoma; Stage III  | 2012-07-30 |
| [NCT02222272](https://clinicaltrials.gov/study/NCT02222272) | Effect of 2nd Gen TKI in CML | COMPLETED | N/A | Myeloid Leukemia, Chronic | 2010-01-01 |
| [NCT06175702](https://clinicaltrials.gov/study/NCT06175702) | Treatment Protocol for Newky Diagnosed Adult Ph Positive ALL | NOT_YET_RECRUITING | N/A | Lymphoblastic Leukemia; Philadelphia-Positive ALL; Adult ALL | 2023-12-25 |
| [NCT04530565](https://clinicaltrials.gov/study/NCT04530565) | Testing the Use of Steroids and Tyrosine Kinase Inhibitors With Blinatumomab or Chemotherapy for Newly Diagnosed BCR-ABL | ACTIVE_NOT_RECRUITING | PHASE3 | B Acute Lymphoblastic Leukemia With t(9;22)(q34.1;q11.2); BCR-ABL1 | 2021-01-25 |
| [NCT04603872](https://clinicaltrials.gov/study/NCT04603872) | CAR-T Cells Combined With Dasatinib for Patients With Relapsed and/or Refractory B-cell Hematological Malignancies | RECRUITING | EARLY_PHASE1 | Multiple Myeloma in Relapse; Multiple Myeloma, Refractory; Acute Lymphoblastic L | 2020-11-01 |
| [NCT06124157](https://clinicaltrials.gov/study/NCT06124157) | A Study Testing the Combination of Dasatinib or Imatinib to Chemotherapy Treatment With Blinatumomab for Children, Adole | RECRUITING | PHASE2 | B Acute Lymphoblastic Leukemia | 2025-05-30 |
| [NCT02011945](https://clinicaltrials.gov/study/NCT02011945) | A Phase 1B Study to Investigate the Safety and Preliminary Efficacy for the Combination of Dasatinib Plus Nivolumab in P | COMPLETED | PHASE1 | Chronic Myeloid Leukemia | 2014-02-07 |
| [NCT06645886](https://clinicaltrials.gov/study/NCT06645886) | A Study to Investigate the Safety and Efficacy of KQB198 as Monotherapy and in Combination in Participants With Advanced | ACTIVE_NOT_RECRUITING | PHASE1 | Hematologic Malignancies; Adult | 2024-12-09 |

No kidney or IgA Nephropathy trials found for this drug — confirms this is a novel repurposing prediction.

---

## Next steps for validation

1. **In silico:** protein-ligand co-folding against LYN (Boltz-2 results for the five top-ranked candidates are in results/boltz/)
2. **In vitro:** Test in primary human mesangial cells, or the mouse mesangial cell line MES-13
3. **In vivo:** ddY mouse model (spontaneous IgAN) — measure proteinuria, UPCR, Gd-IgA1
4. **Clinical:** IgAN-specific Phase 2 trial (primary endpoint: 40% UPCR reduction at 12 months)