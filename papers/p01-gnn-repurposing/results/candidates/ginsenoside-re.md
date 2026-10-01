# Candidate dossier: ginsenoside-re

**Generated:** 2026-06-05  |  **Pipeline:** IgA Nephropathy GNN drug repurposing

---

## Pipeline scores

| Metric | Value |
|--------|-------|
| Final rank | 11 |
| Priority score | 5 / 9 |
| Network proximity rank | 30 / 1,989 (z = -4.265) |
| Open Targets genetic rank | 16 / 3,005 (score 0.4961) |
| GNN link prediction rank (full R-GCN) | — (drug not in PrimeKG drug node set) |
| Hit alignment | Hit 4 — fibrosis and mesangial damage (TGF-beta/NF-kB) |
| GWAS target | REL |
| Clinical phase | Phase 1 |
| Mechanism of action | anti-inflammatory agent |

---

## GWAS locus context

REL encodes NF-kB proto-oncogene c-Rel, a transcription factor subunit controlling inflammatory cytokine production and B-cell survival. Kiryluk et al. 2023 identified REL as a genome-wide significant IgAN locus (rank 23 of 3,005, score 0.496). NF-kB drives BAFF and APRIL expression, connecting REL to the autoantibody axis (Hit 2) and fibrotic gene transcription (Hit 4).

---

## CKD safety profile

| Parameter | Assessment |
|-----------|------------|
| eGFR adjustment | none |
| Nephrotoxicity risk | none |
| Hypertension risk | none |

Phase 1 anti-inflammatory. Natural compound; low clinical tractability for IgAN repurposing. Include as mechanistic validation of REL/NF-κB axis.


---

## PubMed evidence (drug + IgA Nephropathy/kidney fibrosis)

**Integrative informatics analysis identifies that ginsenoside Re improves renal fibrosis through regulation of autophagy.**
*Journal of natural medicines* (2024) | PMID: 38683298
We previously demonstrated that ginsenoside Re (G-Re) has protective effects on acute kidney injury. However, the underlying mechanism is still unclear. In this study, we conducted a meta-analysis and pathway enrichment analysis of all published transcriptome data to identify differentially expresse...

---

## Open Targets known indications

No Open Targets indication data retrieved.

---

## ClinicalTrials.gov active and completed trials

| NCT ID | Title | Status | Phase | Condition | Start |
|--------|-------|--------|-------|-----------|-------|
| [NCT04069715](https://clinicaltrials.gov/study/NCT04069715) | The Effect of Farlong® NotoGinseng™ (Ginseng Plus®) on Cholesterol and Blood Pressure | COMPLETED | PHASE2 | Hyperlipidemias; Hypertension | 2016-07-20 |
| [NCT02204826](https://clinicaltrials.gov/study/NCT02204826) | Effects of Korean Red Ginseng on Semen Parameters in Male Infertility Patients: a Randomized, Placebo-controlled, Double | COMPLETED | PHASE4 | Male Infertility | 2011-04 |
| [NCT04142151](https://clinicaltrials.gov/study/NCT04142151) | Double Randomized and Placebo Controlled Trail of Sanchitongshu Combined Antiplatelet Drug to Prevent Ischemic Stroke | UNKNOWN | NA | Ischemic Stroke; Antiplatelet Effect | 2019-09-01 |
| [NCT00781534](https://clinicaltrials.gov/study/NCT00781534) | A Clinical Trial of Ginseng in Diabetes | COMPLETED | EARLY_PHASE1 | Diabetes | 2003-09 |
| [NCT02975076](https://clinicaltrials.gov/study/NCT02975076) | Sanchitongtshu Plus Asprine for Minor Ischemic Stroke or Transient Ischemic Attack: A Randomized Double-blind Study | UNKNOWN | NA | Stroke, Acute; Ischemic Attack, Transient | 2016-12 |
| [NCT02413099](https://clinicaltrials.gov/study/NCT02413099) | The Efficacy and Safety of New Herbal Formula (KBMSI-2) in the Treatment of Erectile Dysfunction | COMPLETED | PHASE4 | Erectile Dysfunction | 2012-02 |

No kidney or IgA Nephropathy trials found for this drug — confirms this is a novel repurposing prediction.

---

## Next steps for validation

1. **In silico:** protein-ligand co-folding against REL (Boltz-2 results for the five top-ranked candidates are in results/boltz/)
2. **In vitro:** Test in primary human mesangial cells, or the mouse mesangial cell line MES-13
3. **In vivo:** ddY mouse model (spontaneous IgAN) — measure proteinuria, UPCR, Gd-IgA1
4. **Clinical:** IgAN-specific Phase 2 trial (primary endpoint: 40% UPCR reduction at 12 months)