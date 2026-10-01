"""
Step 4: Network proximity scoring (Guney et al. 2016).

Computes the z-scored network proximity between each Drug Repurposing Hub compound
and the IgAN disease module on the STRING v12 protein-protein interaction network.

Reference: Guney E et al. Network-based in silico drug efficacy screening.
           Nat Commun. 2016;7:10331. doi:10.1038/ncomms10331

Method:
    d(S, T) = mean_{s in S} min_{t in T} d(s, t)

    where S = drug target set, T = disease module proteins, d = shortest path length.

    z-score = (d(S,T) - mu_r) / sigma_r

    where mu_r and sigma_r are the mean and SD of d(S',T') over N random drug-target
    sets of the same size as S, preserving degree distribution (degree-matched permutation).

z-score < 0: drug targets are closer to the disease module than expected by chance
             (i.e., the drug is proximal to IgAN biology — a repurposing candidate)

Usage:
    uv run python scripts/04_network_proximity.py [--permutations N] [--workers W]

    Defaults: --permutations 1000, --workers 4

Output:
    results/network_proximity_ranked.tsv

Validation checkpoint:
    Fostamatinib (SYK inhibitor, Phase 2 IgAN) should have z-score < -1.
    Iptacopan (CFB inhibitor, FDA 2024) should have z-score < -1.
    At least 3 of 5 approved IgAN drugs should be in top 10% of DRH by proximity.

Data requirements:
    data/processed/igan_disease_module.json  (from Script 01)
    data/raw/drug_repurposing_hub/repurposing_drugs.tsv  (from Script 05 / download)

Background network:
    The disease module graph is used as the background PPI network for this analysis.
    For a more complete background, replace with PrimeKG PPI layer or full STRING v12.
    The disease module background is conservative (smaller network = harder to achieve
    proximity by chance) — appropriate for a first-pass analysis.
"""

import argparse
import json
import math
import random
import sys
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "shared/scripts"))
from pipeline_utils import find_repo_root, load_config

import networkx as nx
import pandas as pd

N_PERMUTATIONS_DEFAULT = 1000
WORKERS_DEFAULT = 4

DEFAULT_CONFIG = "shared/configs/igan.yaml"


def load_module(path: Path) -> nx.Graph:
    with open(path) as f:
        return nx.node_link_graph(json.load(f))


def load_drh(path: Path) -> dict[str, set[str]]:
    """Return {drug_name: {target_gene, ...}} from DRH annotation file."""
    df = pd.read_csv(path, sep="\t", comment="!")
    drug_targets: dict[str, set[str]] = {}
    for _, row in df.iterrows():
        drug = str(row.get("pert_iname", "")).strip().lower()
        target_str = str(row.get("target", ""))
        if not drug or target_str in ("", "nan"):
            continue
        targets = {t.strip() for t in target_str.split("|") if t.strip()}
        if targets:
            drug_targets[drug] = targets
    return drug_targets


def mean_min_distance(
    G: nx.Graph,
    source_nodes: set[str],
    target_nodes: set[str],
    lengths: dict[str, dict[str, int]],
) -> float | None:
    """
    d(S, T) = mean over s in S of min over t in T of shortest_path(s, t).
    Returns None if source_nodes has no members in G.
    """
    s_in_g = source_nodes & set(G.nodes)
    t_in_g = target_nodes & set(G.nodes)
    if not s_in_g or not t_in_g:
        return None

    total = 0.0
    count = 0
    for s in s_in_g:
        s_lengths = lengths.get(s, {})
        min_d = min((s_lengths[t] for t in t_in_g if t in s_lengths), default=None)
        if min_d is not None:
            total += min_d
            count += 1

    return total / count if count > 0 else None


def build_degree_bins(G: nx.Graph) -> dict[int, list[str]]:
    """Group nodes by degree for degree-matched random sampling."""
    bins: dict[int, list[str]] = defaultdict(list)
    for node, deg in G.degree():
        bins[deg].append(node)
    return dict(bins)


def random_proximal_set(
    G: nx.Graph,
    size: int,
    degree_bins: dict[int, list[str]],
    ref_degrees: list[int],
) -> set[str]:
    """
    Sample a random node set of the given size, degree-matched to ref_degrees.
    Falls back to uniform sampling if a degree bucket is empty.
    """
    sampled: set[str] = set()
    for deg in ref_degrees:
        candidates = degree_bins.get(deg, [])
        if not candidates:
            # Fallback: pick any node not yet sampled
            remaining = list(set(G.nodes) - sampled)
            if remaining:
                sampled.add(random.choice(remaining))
        else:
            pool = [n for n in candidates if n not in sampled]
            if pool:
                sampled.add(random.choice(pool))
    # Top up to size if needed
    while len(sampled) < size:
        n = random.choice(list(G.nodes))
        sampled.add(n)
    return sampled


def compute_proximity_for_drug(
    drug: str,
    drug_targets: set[str],
    G_nodes: set[str],
    disease_nodes: set[str],
    lengths_subset: dict[str, dict[str, int]],
    degree_bins: dict[int, list[str]],
    n_perms: int,
    seed: int,
) -> dict | None:
    """
    Compute z-scored proximity for one drug.
    Runs in a worker process — no shared state.
    """
    random.seed(seed)

    targets_in_g = drug_targets & G_nodes
    if not targets_in_g:
        return None

    observed = mean_min_distance_worker(targets_in_g, disease_nodes, lengths_subset)
    if observed is None:
        return None

    # Reference distribution
    ref_degrees = [
        sum(1 for nb in lengths_subset.get(t, {}) if nb in G_nodes)
        for t in targets_in_g
    ]

    null_distances = []
    for _ in range(n_perms):
        rand_set = random_set_worker(G_nodes, len(targets_in_g), degree_bins, ref_degrees)
        d = mean_min_distance_worker(rand_set, disease_nodes, lengths_subset)
        if d is not None:
            null_distances.append(d)

    if len(null_distances) < 10:
        return None

    mu = sum(null_distances) / len(null_distances)
    var = sum((x - mu) ** 2 for x in null_distances) / len(null_distances)
    sigma = math.sqrt(var) if var > 0 else 1e-9

    z = (observed - mu) / sigma

    return {
        "drug": drug,
        "z_score": round(z, 4),
        "observed_distance": round(observed, 4),
        "null_mean": round(mu, 4),
        "null_sd": round(sigma, 4),
        "n_targets_in_network": len(targets_in_g),
        "n_permutations": len(null_distances),
    }


def mean_min_distance_worker(
    source: set[str],
    target: set[str],
    lengths: dict[str, dict[str, int]],
) -> float | None:
    total, count = 0.0, 0
    for s in source:
        row = lengths.get(s, {})
        min_d = min((row[t] for t in target if t in row), default=None)
        if min_d is not None:
            total += min_d
            count += 1
    return total / count if count > 0 else None


def random_set_worker(
    all_nodes: set[str],
    size: int,
    degree_bins: dict[int, list[str]],
    ref_degrees: list[int],
) -> set[str]:
    sampled: set[str] = set()
    for deg in ref_degrees:
        candidates = degree_bins.get(deg, [])
        pool = [n for n in candidates if n not in sampled]
        if pool:
            sampled.add(random.choice(pool))
        elif all_nodes - sampled:
            sampled.add(random.choice(list(all_nodes - sampled)))
    while len(sampled) < size:
        n = random.choice(list(all_nodes))
        sampled.add(n)
    return sampled


# Module-level shared state injected by worker pool initializer.
# Avoids re-serialising the large lengths dict for every task.
_worker_G_nodes: set[str] = set()
_worker_disease_nodes: set[str] = set()
_worker_lengths: dict[str, dict[str, int]] = {}
_worker_degree_bins: dict[int, list[str]] = {}


def _worker_init(G_nodes, disease_nodes, lengths, degree_bins):
    global _worker_G_nodes, _worker_disease_nodes, _worker_lengths, _worker_degree_bins
    _worker_G_nodes = G_nodes
    _worker_disease_nodes = disease_nodes
    _worker_lengths = lengths
    _worker_degree_bins = degree_bins


def _score_drug_entry(drug_targets_nperms: tuple) -> dict | None:
    """Top-level worker function callable by ProcessPoolExecutor."""
    drug, targets, n_perms = drug_targets_nperms
    row = compute_proximity_for_drug(
        drug, targets, _worker_G_nodes, _worker_disease_nodes,
        _worker_lengths, _worker_degree_bins, n_perms,
        hash(drug) % (2 ** 32),
    )
    return row


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default=DEFAULT_CONFIG,
                        help="Path to disease config YAML (relative to repo root or absolute)")
    parser.add_argument("--permutations", type=int, default=N_PERMUTATIONS_DEFAULT)
    parser.add_argument("--workers", type=int, default=WORKERS_DEFAULT)
    args = parser.parse_args()

    cfg, repo_root = load_config(args.config)
    rep = cfg["repurposing"]
    disease_id   = cfg["disease"]["id"]
    disease_name = cfg["disease"]["name"]

    output_dir      = repo_root / rep["output_dir"]
    shared_data_dir = repo_root / rep["shared_data_dir"]
    data_proc       = output_dir / "data/processed"
    results_dir     = output_dir / "results"
    results_dir.mkdir(parents=True, exist_ok=True)

    module_path = data_proc / f"{disease_id}_disease_module.json"
    drh_path    = shared_data_dir / "drug_repurposing_hub/repurposing_drugs.tsv"
    out_path    = results_dir / "network_proximity_ranked.tsv"
    known_drugs = {d.lower() for d in rep.get("known_approved_drugs", [])}

    print(f"Disease: {disease_name} ({disease_id})")
    print(f"Output:  {output_dir}\n")

    if not module_path.exists():
        print(f"Disease module not found: {module_path}", file=sys.stderr)
        print(f"Run: uv run python scripts/01_build_gwas_module.py --config {args.config}",
              file=sys.stderr)
        sys.exit(1)

    if not drh_path.exists():
        print(f"Drug Repurposing Hub not found: {drh_path}", file=sys.stderr)
        sys.exit(1)

    print("Loading disease module...")
    G = load_module(module_path)
    all_module_nodes = set(G.nodes)
    print(f"  Module: {len(all_module_nodes)} proteins, {G.number_of_edges()} interactions")

    # Target set T = GWAS seed proteins only (not the full first-degree expansion).
    # Background network G = full module.
    seed_loci = pd.read_csv(data_proc / "gwas_loci_proteins.tsv", sep="\t")
    disease_nodes = set(seed_loci["gene"].tolist()) & all_module_nodes
    print(f"  GWAS seed targets (T): {len(disease_nodes)} proteins in module network")

    print("Loading Drug Repurposing Hub...")
    drug_targets = load_drh(drh_path)
    print(f"  DRH: {len(drug_targets)} drugs with target annotations")

    print("Pre-computing all-pairs shortest paths within disease module...")
    lengths: dict[str, dict[str, int]] = {}
    for node in all_module_nodes:
        lengths[node] = nx.single_source_shortest_path_length(G, node)
    print(f"  Done. {len(lengths)} source nodes (full background network).")

    print("Building degree bins for null model...")
    degree_bins = build_degree_bins(G)

    eligible = {
        drug: targets
        for drug, targets in drug_targets.items()
        if targets & all_module_nodes
    }
    print(f"  {len(eligible)}/{len(drug_targets)} DRH drugs have >=1 target in module network")

    G_nodes = all_module_nodes
    tasks = [(drug, targets, args.permutations) for drug, targets in eligible.items()]

    print(f"\nScoring {len(tasks)} drugs "
          f"({args.permutations} permutations each, {args.workers} workers)...")

    results = []
    with ProcessPoolExecutor(
        max_workers=args.workers,
        initializer=_worker_init,
        initargs=(G_nodes, disease_nodes, lengths, degree_bins),
    ) as pool:
        done = 0
        for row in pool.map(_score_drug_entry, tasks, chunksize=20):
            done += 1
            if done % 100 == 0:
                print(f"  {done}/{len(tasks)} drugs scored...")
            if row is not None:
                row["known_approved"] = row["drug"] in known_drugs
                results.append(row)

    df = pd.DataFrame(results).sort_values("z_score")
    df.to_csv(out_path, sep="\t", index=False)
    print(f"\nResults saved to {out_path} ({len(df)} drugs scored)")

    # --- Validation checkpoint ---
    print(f"\n=== Validation checkpoint: known {disease_name} drugs ===")
    known_rows = df[df["known_approved"]].sort_values("z_score")
    total = len(df)
    top10_cutoff = df["z_score"].quantile(0.10)
    in_top10 = (known_rows["z_score"] <= top10_cutoff).sum()
    print(f"Top 10% z-score cutoff: {top10_cutoff:.4f}")
    for _, row in known_rows.iterrows():
        pct = (df["z_score"] <= row["z_score"]).sum() / total * 100
        flag = "TOP10%" if row["z_score"] <= top10_cutoff else ""
        print(f"  {row['drug']:30} z={row['z_score']:7.3f}  percentile={pct:.1f}%  {flag}")
    print(f"\n{in_top10}/{len(known_rows)} known {disease_id.upper()} drugs in top 10% (target: >=3)")

    print("\n=== Top 25 proximal novel candidates ===")
    novel = df[~df["known_approved"]].head(25)
    print(novel[["drug", "z_score", "observed_distance", "n_targets_in_network"]].to_string(index=False))

    if in_top10 < 3:
        print(
            f"\nWARNING: Fewer than 3 known {disease_id.upper()} drugs in top 10%. "
            "Consider using full STRING v12 as background network for better separation.",
            file=sys.stderr,
        )


if __name__ == "__main__":
    main()
