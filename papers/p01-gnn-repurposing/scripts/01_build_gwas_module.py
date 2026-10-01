"""
Step 1: Build the disease module from GWAS loci.

Maps GWAS seed genes onto STRING v12 protein-protein interactions to construct a
disease module — a NetworkX subgraph of seed proteins and their first-degree
interaction partners (confidence >= 700).

This module is the foundation of network proximity scoring (script 04) and the
disease enrichment layer in GNN training (script 06).

Usage:
    # IgAN (default):
    uv run python scripts/01_build_gwas_module.py

    # Any disease with a config file:
    uv run python scripts/01_build_gwas_module.py --config shared/configs/fsgs.yaml

Output (paths from config repurposing.output_dir):
    data/processed/gwas_loci_proteins.tsv      -- seed loci table
    data/processed/<disease_id>_disease_module.json -- NetworkX node-link JSON

Validation checkpoint:
    Module should contain 500-2000 proteins. Proteins in config
    repurposing.validation_proteins must all be present. If module is too small,
    lower STRING_CONFIDENCE_CUTOFF to 500.

Data sources:
    STRING v12 human PPI (download manually):
        https://string-db.org/cgi/download?sessionId=&species_text=Homo+sapiens
        File: 9606.protein.links.detailed.v12.0.txt.gz
        Save to: <output_dir>/data/raw/string_v12/

    If STRING files are absent, falls back to the STRING REST API.
"""

import argparse
import gzip
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "shared/scripts"))
from pipeline_utils import find_repo_root, load_config

import networkx as nx
import pandas as pd
import requests

STRING_API = "https://string-db.org/api/json/network"
STRING_CONFIDENCE_CUTOFF = 700
STRING_SPECIES = 9606  # Homo sapiens

DEFAULT_CONFIG = "shared/configs/igan.yaml"


def fetch_string_api(gene_list: list[str], caller_id: str) -> pd.DataFrame:
    """
    Fetch STRING interactions via REST API with first-degree neighbourhood expansion.

    Two-step:
    1. Get all interaction partners for seed genes (interaction_partners endpoint)
    2. Get seed-seed interactions (network endpoint) — catches missed intra-seed edges
    Merges and deduplicates both result sets.
    """
    print("Using STRING REST API (slower). Download STRING files for faster runs.")

    identifiers = "%0d".join(gene_list)
    all_rows: list[dict] = []

    # Step 1: expand to first-degree neighbours
    partners_url = "https://string-db.org/api/json/interaction_partners"
    params = {
        "identifiers": identifiers,
        "species": STRING_SPECIES,
        "required_score": STRING_CONFIDENCE_CUTOFF,
        "caller_identity": caller_id,
    }
    print("  Querying STRING interaction_partners (full first-degree expansion)...")
    resp = requests.post(partners_url, data=params, timeout=90)
    resp.raise_for_status()
    partners = resp.json()
    print(f"  interaction_partners returned {len(partners)} edges")
    all_rows.extend(partners)

    # Step 2: seed-seed network
    network_params = {
        "identifiers": identifiers,
        "species": STRING_SPECIES,
        "required_score": STRING_CONFIDENCE_CUTOFF,
        "caller_identity": caller_id,
    }
    print("  Querying STRING network (seed-seed edges)...")
    resp2 = requests.post(STRING_API, data=network_params, timeout=60)
    resp2.raise_for_status()
    seed_seed = resp2.json()
    print(f"  network returned {len(seed_seed)} seed-seed edges")
    all_rows.extend(seed_seed)

    if not all_rows:
        return pd.DataFrame(columns=["preferredName_A", "preferredName_B", "score"])

    df = pd.DataFrame(all_rows)[["preferredName_A", "preferredName_B", "score"]]
    # Deduplicate (same edge, different order)
    df["pair"] = df.apply(lambda r: tuple(sorted([r["preferredName_A"], r["preferredName_B"]])), axis=1)
    df = df.sort_values("score", ascending=False).drop_duplicates(subset="pair").drop(columns="pair")
    return df


def load_string_local(seed_genes: set[str], data_raw: Path) -> pd.DataFrame:
    """Load STRING interactions from local files, filtered to seed genes."""
    string_links = data_raw / "9606.protein.links.detailed.v12.0.txt.gz"
    string_info  = data_raw / "9606.protein.info.v12.0.txt.gz"
    if not string_links.exists() or not string_info.exists():
        return pd.DataFrame()  # Signal to fall back to API

    print(f"Loading STRING v12 from {string_links}...")

    # Load protein info for gene-name mapping
    with gzip.open(string_info, "rt") as f:
        info = pd.read_csv(f, sep="\t")
    # Columns: #string_protein_id, preferred_name, protein_size, annotation
    string_id_to_gene = dict(zip(info["#string_protein_id"], info["preferred_name"]))

    seed_ids = {v for k, v in {g: k for k, g in string_id_to_gene.items()}.items() if v in seed_genes}

    print(f"Loading PPI links (confidence >= {STRING_CONFIDENCE_CUTOFF})...")
    with gzip.open(string_links, "rt") as f:
        links = pd.read_csv(f, sep=" ")

    # Filter: at least one partner must be a seed gene
    links = links[links["combined_score"] >= STRING_CONFIDENCE_CUTOFF]
    seed_links = links[links["protein1"].isin(seed_ids) | links["protein2"].isin(seed_ids)]

    # Map STRING IDs back to gene symbols
    seed_links = seed_links.copy()
    seed_links["preferredName_A"] = seed_links["protein1"].map(string_id_to_gene)
    seed_links["preferredName_B"] = seed_links["protein2"].map(string_id_to_gene)
    seed_links["score"] = seed_links["combined_score"]

    return seed_links[["preferredName_A", "preferredName_B", "score"]].dropna()


def build_module(loci: list[dict], interactions: pd.DataFrame) -> nx.Graph:
    """Construct disease module graph from GWAS seeds and interaction data."""
    G = nx.Graph()

    # Add seed nodes
    for locus in loci:
        G.add_node(locus["gene"], hit=locus["hit"], seed=True, note=locus["note"])

    # Add interaction edges (seed + first-degree neighbours)
    seed_genes = {l["gene"] for l in loci}
    for _, row in interactions.iterrows():
        a, b = row["preferredName_A"], row["preferredName_B"]
        if a in seed_genes or b in seed_genes:
            score = float(row.get("score", 0))
            if not G.has_node(a):
                G.add_node(a, hit=0, seed=False, note="STRING neighbour")
            if not G.has_node(b):
                G.add_node(b, hit=0, seed=False, note="STRING neighbour")
            G.add_edge(a, b, weight=score / 1000.0)

    return G


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default=DEFAULT_CONFIG,
                        help="Path to disease config YAML (relative to repo root or absolute)")
    args = parser.parse_args()

    cfg, repo_root = load_config(args.config)
    rep = cfg["repurposing"]
    disease_id   = cfg["disease"]["id"]
    disease_name = cfg["disease"]["name"]
    caller_id    = f"{disease_id}_repurposing_pipeline"

    output_dir = repo_root / rep["output_dir"]
    data_raw   = output_dir / "data/raw/string_v12"
    data_proc  = output_dir / "data/processed"
    data_proc.mkdir(parents=True, exist_ok=True)

    loci = rep["gwas_loci"]
    seed_genes = [l["gene"] for l in loci]
    validation_proteins = set(rep.get("validation_proteins", []))

    print(f"Disease: {disease_name} ({disease_id})")
    print(f"Output:  {output_dir}")
    print(f"Loci:    {len(loci)} GWAS seed genes\n")

    # Save loci table
    loci_df = pd.DataFrame(loci)
    loci_path = data_proc / "gwas_loci_proteins.tsv"
    loci_df.to_csv(loci_path, sep="\t", index=False)
    print(f"Saved {len(loci_df)} GWAS loci to {loci_path}")

    # Load STRING interactions
    interactions = load_string_local(set(seed_genes), data_raw)
    if interactions.empty:
        print("Local STRING files not found — falling back to API...")
        interactions = fetch_string_api(seed_genes, caller_id)

    print(f"Loaded {len(interactions)} STRING interactions involving seed genes")

    # Build NetworkX graph
    G = build_module(loci, interactions)
    n_nodes = G.number_of_nodes()
    n_edges = G.number_of_edges()
    print(f"\nDisease module: {n_nodes} proteins, {n_edges} interactions")

    # Validation checkpoint
    present = validation_proteins & set(G.nodes)
    missing = validation_proteins - present
    if missing:
        print(f"WARNING: Required proteins missing: {missing}", file=sys.stderr)
        print("Consider lowering STRING_CONFIDENCE_CUTOFF to 500.", file=sys.stderr)
    else:
        print(f"Checkpoint passed: {validation_proteins} all present in module")

    if not (500 <= n_nodes <= 2000):
        print(f"WARNING: Module size {n_nodes} outside expected range 500-2000",
              file=sys.stderr)

    # Save module
    out_path = data_proc / f"{disease_id}_disease_module.json"
    module_data = nx.node_link_data(G)
    with open(out_path, "w") as f:
        json.dump(module_data, f, indent=2)
    print(f"\nDisease module saved to {out_path}")

    # Print hit distribution
    hit_labels = rep.get("hit_labels", {})
    hit_counts: dict[int, int] = {}
    for _, attr in G.nodes(data=True):
        h = attr.get("hit", 0)
        hit_counts[h] = hit_counts.get(h, 0) + 1
    print("\nNode distribution by hit:")
    for h in sorted(hit_counts):
        label = hit_labels.get(h, hit_labels.get(str(h), f"Hit {h}"))
        print(f"  Hit {h} ({label}): {hit_counts[h]} proteins")


if __name__ == "__main__":
    main()
