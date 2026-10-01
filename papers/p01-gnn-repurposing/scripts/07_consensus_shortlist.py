"""
Step 7: Build consensus shortlist from three ranked drug lists.

Takes the top-N candidates from each of three independent ranking methods:
  - results/network_proximity_ranked.tsv  (Guney 2016 z-scores)
  - results/opentargets_ranked.tsv        (Open Targets genetic association × DRH)
  - results/gnn_ranked.tsv                (PrimeKG GNN link prediction score)

Candidates appearing in all three top-N lists form the consensus shortlist.
The intersection requirement reduces false-positives vs. any single method.

Usage:
    # IgAN (default):
    uv run python scripts/07_consensus_shortlist.py [--top N]

    # Any disease:
    uv run python scripts/07_consensus_shortlist.py --config shared/configs/fsgs.yaml

    Default: --top 100 (intersection of top 100 from each method)

Output (in config repurposing.output_dir/results/):
    consensus_shortlist.tsv

Validation checkpoint:
    Shortlist should contain at least 10 novel candidates. If fewer, increase --top.
"""

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "shared/scripts"))
from pipeline_utils import find_repo_root, load_config

import pandas as pd

DEFAULT_CONFIG = "shared/configs/igan.yaml"


def load_ranked(
    path: Path,
    drug_col: str,
    score_col: str,
    source_label: str,
    ascending: bool = False,
    dedup_by: str | None = None,
    dedup_agg: str = "max",
) -> pd.DataFrame | None:
    """
    Load a ranked drug list from a TSV.

    ascending=False → higher score is better (OT genetic score, GNN score).
    ascending=True  → lower score is better (network proximity z-score: most
                       negative = most proximal = rank 1).
    dedup_by        → column to group by before keeping best score per drug
                       (used for OT file which has one row per drug-gene pair).
    """
    if not path.exists():
        print(f"  {source_label}: file not found ({path})", file=sys.stderr)
        return None
    df = pd.read_csv(path, sep="\t")

    if dedup_by and dedup_by in df.columns:
        # Deduplicate: keep best score per drug across multiple gene rows
        agg_fn = "min" if ascending else "max"
        df = df.groupby(dedup_by, as_index=False)[score_col].agg(agg_fn)
        drug_col = dedup_by

    df = df[[drug_col, score_col]].copy()
    df.columns = ["drug", "score"]
    df["drug"] = df["drug"].str.lower().str.strip()
    df["source"] = source_label
    # Rank: ascending=True means lower score → rank 1
    df = df.sort_values("score", ascending=ascending).reset_index(drop=True)
    df["rank"] = df.index + 1
    return df


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default=DEFAULT_CONFIG,
                        help="Path to disease config YAML (relative to repo root or absolute)")
    parser.add_argument("--top", type=int, default=100,
                        help="Top-N from each method for intersection")
    parser.add_argument("--np-only", action="store_true",
                        help="Use network proximity alone (skip OT intersection). "
                             "Appropriate for Mendelian diseases with weak GWAS signal.")
    args = parser.parse_args()

    cfg, repo_root = load_config(args.config)
    rep = cfg["repurposing"]
    disease_id   = cfg["disease"]["id"]
    disease_name = cfg["disease"]["name"]

    output_dir  = repo_root / rep["output_dir"]
    results_dir = output_dir / "results"
    results_dir.mkdir(parents=True, exist_ok=True)

    np_file  = results_dir / "network_proximity_ranked.tsv"
    ot_file  = results_dir / "opentargets_ranked.tsv"
    gnn_file = results_dir / "gnn_ranked.tsv"
    out_file = results_dir / "consensus_shortlist.tsv"

    known_drugs = {d.lower() for d in rep.get("known_trial_drugs", [])}
    top_n = args.top

    print(f"Disease: {disease_name} ({disease_id})")
    print(f"Output:  {output_dir}")
    mode = "NP-only" if args.np_only else f"NP ∩ OT (top {top_n})"
    print(f"Mode:    {mode}\n")

    # Network proximity: lower z-score = more proximal = better (ascending=True)
    np_df = load_ranked(np_file, "drug", "z_score", "network_proximity", ascending=True)
    # Open Targets: higher genetic score = better; deduplicate per drug
    ot_df = load_ranked(ot_file, "pert_iname", "ot_genetic_score", "opentargets",
                        ascending=False, dedup_by="pert_iname")
    # GNN: annotation only — DrugBank names differ from DRH names used by NP/OT
    gnn_df = load_ranked(gnn_file, "drug", "gnn_score", "gnn", ascending=False)

    if np_df is None:
        print("network_proximity file missing — run script 04 first.", file=sys.stderr)
        sys.exit(1)

    if args.np_only:
        # NP-only mode: take all drugs with negative z-score (proximal to disease module)
        # Appropriate for Mendelian diseases where OT GWAS signal is insufficient.
        proximal = np_df[np_df["score"] < 0]
        consensus_drugs = set(proximal["drug"])
        print(f"NP-only: {len(consensus_drugs)} drugs with z < 0 (proximal to module)")
        if not consensus_drugs:
            print("No proximal candidates (z < 0). Check disease module and DRH.", file=sys.stderr)
            sys.exit(1)
    else:
        if ot_df is None:
            print("opentargets file missing — run script 05 first, or use --np-only.", file=sys.stderr)
            sys.exit(1)
        np_top = set(np_df[np_df["rank"] <= top_n]["drug"])
        ot_top = set(ot_df[ot_df["rank"] <= top_n]["drug"])
        consensus_drugs = np_top & ot_top
        print(f"Candidates in top {top_n}:")
        print(f"  network_proximity: {len(np_top)}")
        print(f"  opentargets:       {len(ot_top)}")
        if gnn_df is not None:
            gnn_top = set(gnn_df[gnn_df["rank"] <= top_n]["drug"])
            print(f"  gnn (annotation):  {len(gnn_top)} (DrugBank names — not used in intersection)")
        print(f"NP ∩ OT consensus: {len(consensus_drugs)} candidates")
        if not consensus_drugs:
            print("\nNo consensus candidates. Try --np-only for Mendelian diseases.", file=sys.stderr)
            sys.exit(1)

    np_rank_map  = dict(zip(np_df["drug"], np_df["rank"]))
    np_score_map = dict(zip(np_df["drug"], np_df["score"]))
    ot_rank_map  = dict(zip(ot_df["drug"], ot_df["rank"])) if ot_df is not None else {}
    ot_score_map = dict(zip(ot_df["drug"], ot_df["score"])) if ot_df is not None else {}

    # Load OT annotations for hit/moa/phase/gene metadata (optional in np-only mode)
    ot_meta = pd.DataFrame()
    if ot_file.exists():
        ot_full = pd.read_csv(ot_file, sep="\t")
        ot_full["drug_l"] = ot_full["pert_iname"].str.lower().str.strip()
        ot_best = ot_full.sort_values("ot_genetic_score", ascending=False).drop_duplicates("drug_l")
        ot_meta = ot_best.set_index("drug_l")[["hit_alignment", "moa", "clinical_phase", "gene_symbol"]]

    # DRH fallback: provides clinical_phase and moa for drugs missing from OT metadata
    # This is critical in --np-only mode where OT GWAS signal may be absent.
    drh_meta: dict[str, dict] = {}
    shared_data_dir = repo_root / rep["shared_data_dir"]
    drh_path = shared_data_dir / "drug_repurposing_hub/repurposing_drugs.tsv"
    if drh_path.exists():
        drh_raw = pd.read_csv(drh_path, sep="\t", comment="!")
        drh_raw.columns = [c.strip().lower().replace(" ", "_") for c in drh_raw.columns]
        drh_raw["drug_l"] = drh_raw["pert_iname"].str.lower().str.strip()
        phase_col = "clinical_phase" if "clinical_phase" in drh_raw.columns else None
        moa_col   = "moa"            if "moa"            in drh_raw.columns else None
        target_col = "target"        if "target"         in drh_raw.columns else None
        for _, r in drh_raw.drop_duplicates("drug_l").iterrows():
            drh_meta[r["drug_l"]] = {
                "clinical_phase": str(r[phase_col]) if phase_col else "",
                "moa":            str(r[moa_col])   if moa_col   else "",
                "target":         str(r[target_col]) if target_col else "",
            }

    rows = []
    for drug in consensus_drugs:
        row = {
            "drug": drug,
            "np_rank": int(np_rank_map.get(drug, 9999)),
            "np_z_score": round(np_score_map.get(drug, 0.0), 4),
            "ot_rank": int(ot_rank_map.get(drug, 9999)),
            "ot_genetic_score": round(ot_score_map.get(drug, 0.0), 4),
            "mean_rank": round(
                (np_rank_map.get(drug, 9999) + ot_rank_map.get(drug, 9999)) / 2, 1
            ),
        }
        ot = ot_meta.loc[drug] if (not ot_meta.empty and drug in ot_meta.index) else {}
        drh = drh_meta.get(drug, {})
        row["hit_alignment"]  = int(ot.get("hit_alignment", 0)) if hasattr(ot, "get") else 0
        # Target gene: prefer OT gene_symbol, fall back to first DRH target
        ot_gene = str(ot.get("gene_symbol", "")) if hasattr(ot, "get") else ""
        if ot_gene not in ("", "nan"):
            row["target_gene"] = ot_gene
        else:
            drh_target = drh.get("target", "")
            row["target_gene"] = drh_target.split("|")[0].strip() if drh_target else ""
        # Prefer OT phase/moa; fall back to DRH when OT metadata is absent
        ot_phase = str(ot.get("clinical_phase", "")) if hasattr(ot, "get") else ""
        ot_moa   = str(ot.get("moa", ""))            if hasattr(ot, "get") else ""
        row["moa"]            = ot_moa   if ot_moa   not in ("", "nan") else drh.get("moa", "")
        row["clinical_phase"] = ot_phase if ot_phase not in ("", "nan") else drh.get("clinical_phase", "")
        row["novel"] = drug not in known_drugs
        rows.append(row)

    shortlist = pd.DataFrame(rows).sort_values("mean_rank")
    shortlist.to_csv(out_file, sep="\t", index=False)
    print(f"\nConsensus shortlist saved to {out_file}")

    novel_count = int(shortlist["novel"].sum())
    print(f"\nNovel candidates (no known {disease_id.upper()} trial): {novel_count}/{len(shortlist)}")
    print(f"\n{mode} shortlist (sorted by mean NP rank):")
    display_cols = ["drug", "np_rank", "np_z_score", "ot_rank", "ot_genetic_score",
                    "hit_alignment", "target_gene", "clinical_phase", "novel"]
    print(shortlist[display_cols].to_string(index=False))

    if novel_count < 10:
        print(f"\nWARNING: Only {novel_count} novel candidates (target: >= 10).", file=sys.stderr)
        if not args.np_only:
            print("Try --np-only for Mendelian diseases with weak OT GWAS signal.", file=sys.stderr)
        else:
            print("Consider lowering STRING_CONFIDENCE_CUTOFF in script 01.", file=sys.stderr)


if __name__ == "__main__":
    main()
