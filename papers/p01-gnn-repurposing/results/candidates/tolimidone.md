# Candidate dossier: tolimidone

**Generated:** 2026-06-05  |  **Pipeline:** IgA Nephropathy GNN drug repurposing

---

## Pipeline scores

| Metric | Value |
|--------|-------|
| Final rank | 9 |
| Priority score | 6 / 9 |
| Network proximity rank | 40 / 1,989 (z = -3.793) |
| Open Targets genetic rank | 11 / 3,005 (score 0.5052) |
| GNN link prediction rank (full R-GCN) | — (drug not in PrimeKG drug node set) |
| Hit alignment | Hit 2 — autoantibody/B-cell (LYN/BAFF/APRIL axis) |
| GWAS target | LYN |
| Clinical phase | Phase 2 |
| Mechanism of action | SRC activator |

---

## GWAS locus context

LYN encodes a Src-family tyrosine kinase at the apex of the B-cell receptor signalling cascade. Kiryluk et al. 2023 identified LYN as a genome-wide significant IgAN susceptibility locus (rank 22 of 3,005 by Open Targets genetic association score, score 0.505). Downstream cascade: LYN -> SYK -> BTK -> NF-kB -> BAFF/APRIL. Fostamatinib (SYK inhibitor) is in Phase 2 for IgAN, confirming pathway validity. Drugs targeting LYN directly have not been tested in IgAN.

---

## CKD safety profile

| Parameter | Assessment |
|-----------|------------|
| eGFR adjustment | none |
| Nephrotoxicity risk | none |
| Hypertension risk | none |

Phase 2; development discontinued (antidiabetic). SRC activator — unusual mechanism for IgAN. Mechanistic rationale requires further scrutiny. Flag for review.


---

## PubMed evidence (drug + IgA Nephropathy/kidney fibrosis)

No PubMed results found for this drug + IgA Nephropathy/kidney fibrosis query.

---

## Open Targets known indications

| Disease | Phase |
|---------|-------|
| type 2 diabetes mellitus | PHASE_2 |
| obesity | PHASE_2 |

IgA Nephropathy not listed as a known indication — confirms novelty of this repurposing prediction.

---

## ClinicalTrials.gov active and completed trials

| NCT ID | Title | Status | Phase | Condition | Start |
|--------|-------|--------|-------|-----------|-------|
| [NCT06474598](https://clinicaltrials.gov/study/NCT06474598) | An Adaptive Design Study of MTX228 | RECRUITING | PHASE2 | Type 1 Diabetes Mellitus | 2024-11-14 |

No kidney or IgA Nephropathy trials found for this drug — confirms this is a novel repurposing prediction.

---

## Next steps for validation

1. **In silico:** protein-ligand co-folding against LYN (Boltz-2 results for the five top-ranked candidates are in results/boltz/)
2. **In vitro:** Test in primary human mesangial cells, or the mouse mesangial cell line MES-13
3. **In vivo:** ddY mouse model (spontaneous IgAN) — measure proteinuria, UPCR, Gd-IgA1
4. **Clinical:** IgAN-specific Phase 2 trial (primary endpoint: 40% UPCR reduction at 12 months)