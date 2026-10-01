# Candidate dossier: bafetinib

**Generated:** 2026-06-05  |  **Pipeline:** IgA Nephropathy GNN drug repurposing

---

## Pipeline scores

| Metric | Value |
|--------|-------|
| Final rank | 8 |
| Priority score | 6 / 9 |
| Network proximity rank | 35 / 1,989 (z = -3.924) |
| Open Targets genetic rank | 12 / 3,005 (score 0.5052) |
| GNN link prediction rank (full R-GCN) | — (drug not in PrimeKG drug node set) |
| Hit alignment | Hit 2 — autoantibody/B-cell (LYN/BAFF/APRIL axis) |
| GWAS target | LYN |
| Clinical phase | Phase 2 |
| Mechanism of action | Bcr-Abl kinase inhibitor | LYN tyrosine kinase inhibitor |

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

Phase 2 safety data limited; hepatic metabolism; no renal contraindication identified. Most LYN-selective agent in shortlist. Author should verify Phase 2 safety reports.


---

## PubMed evidence (drug + IgA Nephropathy/kidney fibrosis)

No PubMed results found for this drug + IgA Nephropathy/kidney fibrosis query.

---

## Open Targets known indications

| Disease | Phase |
|---------|-------|
| prostate cancer | PHASE_2 |
| chronic lymphocytic leukemia | PHASE_2 |
| bone disease | PHASE_2 |

IgA Nephropathy not listed as a known indication — confirms novelty of this repurposing prediction.

---

## ClinicalTrials.gov active and completed trials

| NCT ID | Title | Status | Phase | Condition | Start |
|--------|-------|--------|-------|-----------|-------|
| [NCT01215799](https://clinicaltrials.gov/study/NCT01215799) | Study of Bafetinib (INNO-406) as Treatment for Patients With Hormone-Refractory Prostate Cancer | COMPLETED | PHASE2 | Hormone Refractory Prostate Cancer | 2010-08 |
| [NCT00352677](https://clinicaltrials.gov/study/NCT00352677) | Safety and Tolerability Study of INNO-406 to Treat Chronic Myeloid Leukemia or Acute Lymphocytic Leukemia | COMPLETED | PHASE1 | Chronic Myeloid Leukemia; Acute Lymphocytic Leukemia | 2006-07 |
| [NCT01234740](https://clinicaltrials.gov/study/NCT01234740) | Bafetinib in Treating Patients With Recurrent High-Grade Glioma or Brain Metastases | COMPLETED | PHASE1 | Adult Anaplastic Astrocytoma; Adult Anaplastic Ependymoma; Adult Anaplastic Olig | 2010-12 |
| [NCT01144260](https://clinicaltrials.gov/study/NCT01144260) | Study of Bafetinib as Treatment for Relapsed or Refractory B-Cell Chronic Lymphocytic Leukemia (B-CLL) | COMPLETED | PHASE2 | B-Cell Chronic Lymphocytic Leukemia | 2010-06 |

No kidney or IgA Nephropathy trials found for this drug — confirms this is a novel repurposing prediction.

---

## Next steps for validation

1. **In silico:** protein-ligand co-folding against LYN (Boltz-2 results for the five top-ranked candidates are in results/boltz/)
2. **In vitro:** Test in primary human mesangial cells, or the mouse mesangial cell line MES-13
3. **In vivo:** ddY mouse model (spontaneous IgAN) — measure proteinuria, UPCR, Gd-IgA1
4. **Clinical:** IgAN-specific Phase 2 trial (primary endpoint: 40% UPCR reduction at 12 months)