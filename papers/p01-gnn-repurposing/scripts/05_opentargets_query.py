"""
Step 5: Open Targets quick-win query.

Queries the Open Targets Platform GraphQL API for all gene associations with the
target disease, then cross-joins with the Drug Repurposing Hub to produce a
preliminary ranked drug list.

This is the fastest path to a preliminary ranked list — no model training required.
Serves as the naive baseline that GNN and network-proximity methods improve on.

Usage:
    # IgAN (default):
    uv run python scripts/05_opentargets_query.py

    # Any disease:
    uv run python scripts/05_opentargets_query.py --config shared/configs/fsgs.yaml

Output (in config repurposing.output_dir/results/):
    opentargets_genes.tsv
    opentargets_ranked.tsv

Validation checkpoint:
    Known approved drugs (from config repurposing.known_approved_drugs) should
    appear in the top 50. If not, verify the EFO ID against Open Targets.

Data sources (no download required — all API):
    Open Targets Platform API: https://platform.opentargets.org/api
    Drug Repurposing Hub: download to <output_dir>/data/raw/drug_repurposing_hub/
"""

import argparse
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "shared/scripts"))
from pipeline_utils import find_repo_root, load_config

import pandas as pd
import requests

OT_API = "https://api.platform.opentargets.org/api/v4/graphql"

DEFAULT_CONFIG = "shared/configs/igan.yaml"

OT_QUERY = """
query DiseaseAssociations($efoId: String!, $size: Int!, $index: Int!) {
  disease(efoId: $efoId) {
    id
    name
    associatedTargets(page: { size: $size, index: $index }) {
      count
      rows {
        target {
          id
          approvedSymbol
          approvedName
          biotype
        }
        score
        datatypeScores {
          id
          score
        }
      }
    }
  }
}
"""


def query_ot_associations(disease_id: str, disease_name: str,
                          page_size: int = 100) -> list[dict]:
    """Fetch all gene associations for a disease from Open Targets."""
    rows = []
    index = 0
    total = None

    while True:
        resp = requests.post(
            OT_API,
            json={"query": OT_QUERY, "variables": {"efoId": disease_id,
                  "size": page_size, "index": index}},
            headers={"Content-Type": "application/json"},
            timeout=30,
        )
        resp.raise_for_status()
        data = resp.json()

        if "errors" in data:
            print(f"GraphQL errors: {data['errors']}", file=sys.stderr)
            break

        disease_data = data["data"]["disease"]
        if disease_data is None:
            print(f"Disease {disease_id} not found in Open Targets.", file=sys.stderr)
            break

        assoc = disease_data["associatedTargets"]
        if total is None:
            total = assoc["count"]
            print(f"Total gene associations for {disease_name}: {total}")

        batch = assoc["rows"]
        rows.extend(batch)
        print(f"  Fetched {len(rows)}/{total} associations...")

        if len(rows) >= total:
            break

        index += 1
        time.sleep(0.2)

    return rows


def parse_associations(rows: list[dict], hit_map: dict[str, int]) -> pd.DataFrame:
    records = []
    for row in rows:
        target = row["target"]
        gene = target["approvedSymbol"]
        genetic_score = next(
            (d["score"] for d in row["datatypeScores"] if d["id"] == "genetic_association"),
            0.0,
        )
        records.append({
            "ensembl_id": target["id"],
            "gene_symbol": gene,
            "gene_name": target["approvedName"],
            "biotype": target["biotype"],
            "ot_overall_score": row["score"],
            "ot_genetic_score": genetic_score,
            "hit_alignment": hit_map.get(gene, 0),
        })
    return pd.DataFrame(records).sort_values("ot_genetic_score", ascending=False)


def load_drh(path: Path) -> pd.DataFrame | None:
    """Load Drug Repurposing Hub. Returns None if not yet downloaded."""
    if not path.exists():
        print(
            f"\nDrug Repurposing Hub not found at {path}.\n"
            "Download from https://clue.io/repurposing and save as:\n"
            f"  {path}\n"
            "Continuing with gene-only output.",
            file=sys.stderr,
        )
        return None

    df = pd.read_csv(path, sep="\t", comment="!")
    # Normalise column names (DRH has inconsistent headers across versions)
    df.columns = [c.strip().lower().replace(" ", "_") for c in df.columns]
    return df


def cross_join_drh(gene_df: pd.DataFrame, drh: pd.DataFrame) -> pd.DataFrame:
    """
    Join gene associations with Drug Repurposing Hub on target gene symbol.
    DRH 'target' column contains semicolon-separated gene symbols.
    """
    # Explode multi-target entries
    drh = drh.copy()
    if "target" in drh.columns:
        # DRH drug annotation file uses " | " as multi-target separator
        drh["gene_symbol"] = drh["target"].str.split(r"\s*\|\s*")
        drh = drh.explode("gene_symbol")
        drh["gene_symbol"] = drh["gene_symbol"].str.strip()
    elif "gene_symbol" in drh.columns:
        pass
    else:
        # Try to find a column that looks like gene symbols
        candidates = [c for c in drh.columns if "gene" in c or "target" in c or "symbol" in c]
        if candidates:
            drh["gene_symbol"] = drh[candidates[0]]
        else:
            print("Cannot find target gene column in DRH.", file=sys.stderr)
            return gene_df

    merged = gene_df.merge(drh, on="gene_symbol", how="inner")
    merged = merged.sort_values(["ot_genetic_score", "ot_overall_score"], ascending=False)
    return merged


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default=DEFAULT_CONFIG,
                        help="Path to disease config YAML (relative to repo root or absolute)")
    args = parser.parse_args()

    cfg, repo_root = load_config(args.config)
    rep = cfg["repurposing"]
    disease_id   = cfg["disease"]["id"]
    disease_name = cfg["disease"]["name"]

    output_dir      = repo_root / rep["output_dir"]
    shared_data_dir = repo_root / rep["shared_data_dir"]
    results_dir     = output_dir / "results"
    drh_path        = shared_data_dir / "drug_repurposing_hub/repurposing_drugs.tsv"
    results_dir.mkdir(parents=True, exist_ok=True)

    efo_id  = rep["efo_id"]
    hit_map = {locus["gene"]: locus["hit"] for locus in rep["gwas_loci"]}
    known_drugs = {d.lower() for d in rep.get("known_approved_drugs", [])}

    print(f"Disease: {disease_name} ({disease_id})")
    print(f"EFO ID:  {efo_id}")
    print(f"Output:  {output_dir}\n")
    print(f"Querying Open Targets for {disease_name} gene associations...")

    rows = query_ot_associations(efo_id, disease_name)
    if not rows:
        print("No associations returned. Check EFO ID.", file=sys.stderr)
        sys.exit(1)

    gene_df = parse_associations(rows, hit_map)
    print(f"\nTop 20 genes by genetic association score:")
    print(gene_df[["gene_symbol", "ot_genetic_score", "ot_overall_score",
                   "hit_alignment"]].head(20).to_string(index=False))

    # Save gene-level output
    gene_out = results_dir / "opentargets_genes.tsv"
    gene_df.to_csv(gene_out, sep="\t", index=False)
    print(f"\nGene associations saved to {gene_out}")

    # Cross-join with Drug Repurposing Hub
    drh = load_drh(drh_path)
    if drh is not None:
        merged = cross_join_drh(gene_df, drh)
        out_path = results_dir / "opentargets_ranked.tsv"
        merged.to_csv(out_path, sep="\t", index=False)
        print(f"Drug-level ranked list saved to {out_path} ({len(merged)} rows)")

        # Validation checkpoint — known approved drugs should appear near top
        print(f"\nValidation checkpoint ({disease_name} known drugs):")
        gene_list = gene_df["gene_symbol"].tolist()
        for locus in rep["gwas_loci"]:
            gene = locus["gene"]
            if gene not in gene_list:
                continue
            gene_rank_idx = gene_list.index(gene) + 1
            score = gene_df[gene_df["gene_symbol"] == gene]["ot_genetic_score"].values[0]
            # Check if any known drug targets this gene
            gene_drugs = merged[merged["gene_symbol"] == gene]["pert_iname"].str.lower() \
                if "pert_iname" in merged.columns else pd.Series([], dtype=str)
            known_hits = [d for d in gene_drugs if d in known_drugs]
            if known_hits:
                print(f"  {gene} (rank {gene_rank_idx}, score {score:.4f}): "
                      f"known drugs = {known_hits}")
    else:
        print("\nSkipped drug cross-join (DRH not downloaded).")
        print(f"Download DRH to {drh_path} and re-run for full output.")


if __name__ == "__main__":
    main()
