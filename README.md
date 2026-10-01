# Leakage-aware GWAS-informed drug repurposing for IgA nephropathy

Code and results for the manuscript "Leakage-aware retrospective validation of a GWAS-informed drug
repurposing pipeline for IgA nephropathy" (Hanif Rahman, 2026).

The pipeline builds an IgA nephropathy disease module from 30 genome-wide association loci, scores
drugs by network proximity, Open Targets genetic association, a relational graph convolutional network
(R-GCN) on PrimeKG and LINCS L1000 transcriptomic reversal, and checks the leading candidates with
Boltz-2 protein-ligand co-folding. Retrospective recall is credited only when a drug's scoring targets
are annotated in a dated (March 2020) release of the Drug Repurposing Hub, and GNN recall is measured
by leave-one-out.

## Setup

Requires [uv](https://docs.astral.sh/uv/). Dependency versions are pinned in `uv.lock`.

```bash
uv sync                 # core pipeline
uv sync --extra gnn     # adds torch and torch-geometric for scripts 06 and 06b
```

## Data

All primary sources are open access. Raw downloads are not redistributed here; place them in
`papers/p01-gnn-repurposing/data/raw/` and verify each against the matching `.sha256` file.

| Source | Location in `data/raw/` | Notes |
|---|---|---|
| PrimeKG (Harvard Dataverse) | `primekg/nodes.tab`, `edges.csv`, `drug_features.tab` | Knowledge graph for the GNN |
| Drug Repurposing Hub, file date 18 Aug 2025 | `drug_repurposing_hub/repurposing_drugs.tsv` | Drug-target annotations used for scoring |
| Drug Repurposing Hub, release 24 Mar 2020 | `drug_repurposing_hub/repurposing_drugs_20200324.txt` | Annotation-vintage check only; [download](https://s3.amazonaws.com/data.clue.io/repurposing/downloads/repurposing_drugs_20200324.txt) |
| STRING v12 | queried via REST API by script 01 | Confidence threshold 700 |
| Open Targets Platform | queried via GraphQL API by script 05 | EFO_0004194 |
| GEO GSE93798 (Liu et al. 2017) | downloaded by `cmap_analysis.py` | Glomerular biopsies, 20 IgAN vs 22 controls |
| LINCS L1000 | queried via Enrichr by `cmap_analysis.py` | |

Scripts 01 and 05 and `cmap_analysis.py` call live APIs, so a rerun reflects the current state of
STRING, Open Targets and Enrichr. The outputs reported in the manuscript are committed under
`results/` and `data/processed/`.

## Pipeline

Run from the repository root:

```bash
S=papers/p01-gnn-repurposing/scripts
uv run python $S/01_build_gwas_module.py          # disease module (953 proteins, 1,388 interactions)
uv run python $S/04_network_proximity.py          # network proximity, 1,989 drugs
uv run python $S/05_opentargets_query.py          # Open Targets genetic association
uv run python $S/06b_gnn_loo_validation.py        # R-GCN: leave-one-out recall + full-model candidate ranks
uv run python $S/06b_gnn_loo_validation.py --encoder graphsage   # architecture ablation
uv run python $S/07_consensus_shortlist.py        # network proximity ∩ Open Targets top-100
uv run python $S/08_clinical_filter.py            # clinical feasibility filter, 15 final candidates
uv run python $S/cmap_analysis.py                 # transcriptomic reversal
uv run python $S/09_candidate_dossier.py          # per-candidate evidence dossiers
uv run python $S/11_boltz_predict.py              # Boltz-2 co-folding (runs boltz via uvx)
uv run python $S/12_boltz_analysis.py
uv run python $S/13_check_annotation_vintage.py   # annotation-vintage check for the recall claim
```

`06_gnn_train.py` is the original single 80/20 split model. It ranks known drugs that were in its
training split, which is why the reported GNN results come from `06b_gnn_loo_validation.py`. It is
kept so that comparison can be reproduced; its output (`results/gnn_ranked.tsv`) is not used for any
reported number.

Figures:

```bash
F=papers/p01-gnn-repurposing/manuscript/figures
uv run python $F/generate_figures.py              # fig1_gwas_coverage.png, fig2_pipeline.png
uv run python $F/fig03_method_convergence.py      # fig3_method_convergence.png
uv run python $F/fig04_loo_ablation.py            # fig4_loo_ablation.png
```

## Where each reported result comes from

Paths are relative to `papers/p01-gnn-repurposing/`.

| Manuscript item | File | Produced by |
|---|---|---|
| Figure 1, pipeline | `manuscript/figures/fig2_pipeline.png` | `generate_figures.py` |
| Figure 2, GWAS coverage map | `manuscript/figures/fig1_gwas_coverage.png` | `generate_figures.py` |
| Figure 3, leave-one-out recall | `results/kg_loo_recall_rgcn.tsv`, `kg_loo_recall_graphsage.tsv` | `06b` |
| Figure 4, method convergence | `manuscript/figures/fig3_method_convergence.png` | `fig03_method_convergence.py` |
| Table 1 and retrospective ranks | `results/network_proximity_ranked.tsv` | `04` |
| Table 2, final candidates | `results/final_candidates.tsv` | `07`, `08` |
| Table 3, transcriptomic reversal | `results/cmap_analysis/p01_candidates_cmap.tsv` | `cmap_analysis.py` |
| Supplementary Table 1, co-folding | `results/boltz/confidence_summary.tsv`, `contact_residues.tsv` | `11`, `12` |
| GNN ranks for candidates | `results/kg_full_candidate_ranked_rgcn.tsv` | `06b` |
| Genetic association ranks | `results/opentargets_genes.tsv`, `opentargets_ranked.tsv` | `05` |
| Annotation-vintage protocol | console output | `13` |
| Candidate dossiers | `results/candidates/*.md` | `09` |

The dossiers' PubMed, Open Targets and ClinicalTrials.gov sections are snapshots from the date shown
in each file.

## License

Code: MIT (see `LICENSE`). Result tables, figures and dossiers in this repository: CC BY 4.0.
Third-party data remain under their providers' terms.
