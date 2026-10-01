"""
Step 6b: Leave-one-out cross-validation and R-GCN vs GraphSAGE architecture
ablation for the GNN link predictor.

06_gnn_train.py trains a single R-GCN with one random 80/20 split of the known
positive drugs and reports rank for ALL drugs, including the 80% the model was
fit on. That conflates a drug's training-set membership with genuine recall for
any drug in the 80% split. This script implements the protocol described in
the manuscript's Methods (GNN link prediction, leave-one-out cross-validation):
for each of the known IgAN drugs in shared/configs/igan.yaml's
known_drugs_drugbank list, train a fresh model with that ONE drug excluded from
the positive supervision set (all known drugs are excluded from the negative
pool in every fold, so a held-out positive is never mistakenly sampled as a
negative label either), then record that drug's rank among all drug nodes.
Every known drug gets exactly one genuinely blind evaluation. After the folds,
one additional "full" model is trained on all known positives (no holdout);
that model is what scores the NOVEL candidates (nintedanib, bosutinib, F351,
etc.), which were never positive labels in the first place and so were never
subject to this circularity. This separates two different questions cleanly:

    1. kg_loo_recall_<encoder>.tsv          - did the model ever genuinely
                                               predict a known drug it wasn't
                                               trained on?
    2. kg_full_candidate_ranked_<encoder>.tsv - where do this paper's 15 final
                                               candidates rank under a
                                               properly-trained model?

Also implements the R-GCN vs GraphSAGE architecture ablation reported in
Results (relation-aware vs relation-agnostic message passing on identical
data and identical LOO protocol): pass --encoder graphsage to reproduce it
with a relation-agnostic SAGEConv encoder that ignores edge relation type.

Usage:
    uv run python scripts/06b_gnn_loo_validation.py [--encoder {rgcn,graphsage}] [--epochs N --hidden D --lr LR --seed S]

    Default encoder is rgcn (matches 06_gnn_train.py). Pass --encoder graphsage
    to run the architecture ablation with a relation-agnostic SAGEConv encoder
    on identical data and identical LOO protocol.

Output:
    results/kg_loo_recall_<encoder>.tsv
    results/kg_full_candidate_ranked_<encoder>.tsv

Data requirements: same as 06_gnn_train.py (data/raw/primekg/nodes.tab,
data/raw/primekg/edges.csv).
"""

import argparse
import random
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parent.parent
REPO_ROOT = ROOT.parent.parent
RESULTS = ROOT / "results"

PRIMEKG_DIR = ROOT / "data" / "raw" / "primekg"
IGAN_CONFIG = REPO_ROOT / "shared" / "configs" / "igan.yaml"

RELEVANT_RELATIONS = {
    "drug_protein", "protein_drug", "protein_protein", "protein_disease",
    "disease_protein", "drug_effect",
}
RELATION_TO_IDX = {rel: i for i, rel in enumerate(sorted(RELEVANT_RELATIONS))}
NUM_RELATIONS = len(RELATION_TO_IDX)
EMB_DIM = 32  # matches 06_gnn_train.py


def check_dependencies() -> bool:
    try:
        import torch  # noqa: F401
        import torch_geometric  # noqa: F401
        return True
    except ImportError:
        return False


def load_igan_config() -> tuple[int, list[str]]:
    with open(IGAN_CONFIG) as f:
        cfg = yaml.safe_load(f)
    rep = cfg["repurposing"]
    return rep["primekg_node_index"], rep["known_drugs_drugbank"]


def load_primekg(nodes_path: Path, edges_path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    print("Loading PrimeKG nodes...")
    nodes = pd.read_csv(nodes_path, sep="\t")
    nodes["node_index"] = nodes["node_index"].astype(int)
    print(f"  {len(nodes):,} nodes")

    print("Loading PrimeKG edges (may take ~30s for 4M rows)...")
    edges = pd.read_csv(edges_path)
    edges["x_index"] = edges["x_index"].astype(int)
    edges["y_index"] = edges["y_index"].astype(int)
    print(f"  {len(edges):,} edges")
    return nodes, edges


def build_graph(nodes: pd.DataFrame, edges: pd.DataFrame, disease_idx: int):
    """Identical construction to 06_gnn_train.py's PrimeKG branch."""
    import torch
    from torch_geometric.data import Data

    rel_edges = edges[edges["relation"].isin(RELEVANT_RELATIONS)].copy()
    drug_nodes = nodes[nodes["node_type"] == "drug"]["node_index"].tolist()
    gene_nodes = nodes[nodes["node_type"] == "gene/protein"]["node_index"].tolist()
    active_nodes = set(drug_nodes) | set(gene_nodes) | {disease_idx}

    mask = rel_edges["x_index"].isin(active_nodes) & rel_edges["y_index"].isin(active_nodes)
    rel_edges = rel_edges[mask]
    print(f"  Relevant edges after node filtering: {len(rel_edges):,}")

    all_active = sorted(active_nodes)
    idx_map = {orig: new for new, orig in enumerate(all_active)}
    n_nodes = len(all_active)

    src = [idx_map[x] for x in rel_edges["x_index"]]
    dst = [idx_map[y] for y in rel_edges["y_index"]]
    rel_type = [RELATION_TO_IDX[r] for r in rel_edges["relation"]]
    edge_index = torch.tensor([src + dst, dst + src], dtype=torch.long)
    edge_type = torch.tensor(rel_type + rel_type, dtype=torch.long)

    x = torch.zeros(n_nodes, 3)  # one-hot: [gene, drug, disease]
    drug_set = set(drug_nodes)
    for orig in all_active:
        i = idx_map[orig]
        if orig in drug_set:
            x[i, 1] = 1.0
        elif orig == disease_idx:
            x[i, 2] = 1.0
        else:
            x[i, 0] = 1.0

    data = Data(x=x, edge_index=edge_index, edge_type=edge_type)
    data.node_ids = torch.arange(n_nodes, dtype=torch.long)
    data.n_nodes = n_nodes
    data.idx_map = idx_map
    data.rev_map = {v: k for k, v in idx_map.items()}
    data.drug_nodes_remapped = [idx_map[n] for n in drug_nodes if n in idx_map]
    data.disease_remapped = idx_map[disease_idx]
    return data


def make_model(data, hidden, encoder="rgcn"):
    """Build a link predictor with either an R-GCN or GraphSAGE encoder.

    Both share the same decoder and embedding layer so the comparison is
    controlled: the only variable is whether the message-passing step is
    relation-aware (R-GCN) or relation-agnostic (GraphSAGE). Both accept
    the same encode(x, edge_index, edge_type, node_ids) call signature so
    train_once() is unchanged; GraphSAGE simply ignores edge_type.
    """
    import torch
    import torch.nn as nn
    import torch.nn.functional as F

    if encoder == "rgcn":
        from torch_geometric.nn import RGCNConv

        class LinkPredictor(nn.Module):
            def __init__(self, n_nodes, feat_channels, hidden_channels):
                super().__init__()
                self.node_emb = nn.Embedding(n_nodes, EMB_DIM)
                nn.init.xavier_uniform_(self.node_emb.weight)
                in_ch = EMB_DIM + feat_channels
                self.conv1 = RGCNConv(in_ch, hidden_channels, num_relations=NUM_RELATIONS)
                self.conv2 = RGCNConv(hidden_channels, hidden_channels, num_relations=NUM_RELATIONS)
                self.decoder = nn.Sequential(
                    nn.Linear(hidden_channels * 2, hidden_channels),
                    nn.ReLU(),
                    nn.Linear(hidden_channels, 1),
                )

            def encode(self, x, edge_index, edge_type, node_ids):
                emb = self.node_emb(node_ids)
                x_in = torch.cat([emb, x], dim=-1)
                h = F.relu(self.conv1(x_in, edge_index, edge_type))
                h = F.dropout(h, p=0.3, training=self.training)
                return self.conv2(h, edge_index, edge_type)

            def decode(self, z, src_idx, tgt_indices):
                src_emb = z[src_idx].unsqueeze(0).expand(len(tgt_indices), -1)
                return self.decoder(torch.cat([src_emb, z[tgt_indices]], dim=-1)).squeeze(-1)

    elif encoder == "graphsage":
        from torch_geometric.nn import SAGEConv

        class LinkPredictor(nn.Module):
            def __init__(self, n_nodes, feat_channels, hidden_channels):
                super().__init__()
                self.node_emb = nn.Embedding(n_nodes, EMB_DIM)
                nn.init.xavier_uniform_(self.node_emb.weight)
                in_ch = EMB_DIM + feat_channels
                self.conv1 = SAGEConv(in_ch, hidden_channels)
                self.conv2 = SAGEConv(hidden_channels, hidden_channels)
                self.decoder = nn.Sequential(
                    nn.Linear(hidden_channels * 2, hidden_channels),
                    nn.ReLU(),
                    nn.Linear(hidden_channels, 1),
                )

            def encode(self, x, edge_index, edge_type, node_ids):  # edge_type unused
                emb = self.node_emb(node_ids)
                x_in = torch.cat([emb, x], dim=-1)
                h = F.relu(self.conv1(x_in, edge_index))
                h = F.dropout(h, p=0.3, training=self.training)
                return self.conv2(h, edge_index)

            def decode(self, z, src_idx, tgt_indices):
                src_emb = z[src_idx].unsqueeze(0).expand(len(tgt_indices), -1)
                return self.decoder(torch.cat([src_emb, z[tgt_indices]], dim=-1)).squeeze(-1)

    else:
        raise ValueError(f"Unknown encoder: {encoder!r}. Choose 'rgcn' or 'graphsage'.")

    return LinkPredictor(data.n_nodes, data.x.size(1), hidden)


def train_once(data, disease_new: int, train_pos: list[int], neg_pool: list[int],
                n_epochs: int, hidden: int, lr: float, encoder: str = "rgcn"):
    """Train one GNN model (R-GCN or GraphSAGE). train_pos is the positive
    supervision set for THIS run only — callers control what's excluded."""
    import torch
    import torch.nn.functional as F

    model = make_model(data, hidden, encoder=encoder)
    optimizer = torch.optim.Adam(model.parameters(), lr=lr)

    for epoch in range(1, n_epochs + 1):
        model.train()
        optimizer.zero_grad()
        z = model.encode(data.x, data.edge_index, data.edge_type, data.node_ids)
        neg_sample = random.sample(neg_pool, min(len(train_pos) * 3, len(neg_pool)))
        pos_scores = model.decode(z, disease_new, train_pos)
        neg_scores = model.decode(z, disease_new, neg_sample)
        labels = torch.cat([torch.ones(len(train_pos)), torch.zeros(len(neg_sample))])
        scores = torch.cat([pos_scores, neg_scores])
        loss = F.binary_cross_entropy_with_logits(scores, labels)
        loss.backward()
        optimizer.step()

    model.eval()
    with torch.no_grad():
        z = model.encode(data.x, data.edge_index, data.edge_type, data.node_ids)
        drug_indices = data.drug_nodes_remapped
        scores = torch.sigmoid(model.decode(z, disease_new, drug_indices)).numpy()
    return drug_indices, scores


def rank_table(data, nodes: pd.DataFrame, drug_indices, scores) -> pd.DataFrame:
    orig_indices = [data.rev_map[i] for i in drug_indices]
    name_map = nodes.set_index("node_index")["node_name"].to_dict()
    drug_names = [name_map.get(orig, str(orig)) for orig in orig_indices]
    df = pd.DataFrame({"drug": drug_names, "node_index": orig_indices, "score": scores})
    df = df.sort_values("score", ascending=False).reset_index(drop=True)
    df["rank"] = range(1, len(df) + 1)
    df["percentile"] = df["rank"] / len(df) * 100
    return df


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--encoder", choices=["rgcn", "graphsage"], default="rgcn",
                        help="GNN encoder architecture (default: rgcn)")
    parser.add_argument("--epochs", type=int, default=200)
    parser.add_argument("--hidden", type=int, default=128)
    parser.add_argument("--lr", type=float, default=0.001)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    if not check_dependencies():
        print("ERROR: torch / torch_geometric not installed.", file=sys.stderr)
        sys.exit(1)
    import torch

    print(f"Encoder: {args.encoder}")
    random.seed(args.seed)
    np.random.seed(args.seed)
    torch.manual_seed(args.seed)

    disease_node_idx, known_drugs = load_igan_config()
    print(f"PrimeKG disease node index: {disease_node_idx}")
    print(f"Known IgAN drugs (positive label set, from shared/configs/igan.yaml): "
          f"{len(known_drugs)}\n")

    for path, label in [(PRIMEKG_DIR / "nodes.tab", "nodes.tab"), (PRIMEKG_DIR / "edges.csv", "edges.csv")]:
        if not path.exists():
            print(f"ERROR: PrimeKG {label} not found at {path}", file=sys.stderr)
            sys.exit(1)

    nodes, edges = load_primekg(PRIMEKG_DIR / "nodes.tab", PRIMEKG_DIR / "edges.csv")
    print(f"\nBuilding drug-protein-disease subgraph (disease node {disease_node_idx})...")
    data = build_graph(nodes, edges, disease_node_idx)
    print(f"  Graph: {data.n_nodes:,} nodes, {data.edge_index.size(1)//2:,} edges (undirected)\n")

    drug_node_set = set(data.drug_nodes_remapped)
    pos_orig = nodes[(nodes["node_type"] == "drug") & nodes["node_name"].isin(known_drugs)]
    pos_name_by_new = {}
    pos_new_all = []
    for _, row in pos_orig.iterrows():
        orig = row["node_index"]
        if orig in data.idx_map:
            new = data.idx_map[orig]
            pos_new_all.append(new)
            pos_name_by_new[new] = row["node_name"]
    missing = sorted(set(known_drugs) - set(pos_name_by_new.values()))
    if missing:
        print(f"WARNING: {len(missing)} known drugs not found as PrimeKG drug nodes: {missing}", file=sys.stderr)
    print(f"Known drugs resolved in current PrimeKG: {len(pos_new_all)}/{len(known_drugs)}\n")

    neg_pool_base = [n for n in drug_node_set if n not in set(pos_new_all)]
    disease_new = data.disease_remapped

    # ---- Phase 1: leave-one-out over the known positives ----
    print("=" * 70)
    print(f"Phase 1: leave-one-out recall test ({len(pos_new_all)} folds, encoder={args.encoder})")
    print("=" * 70)
    loo_rows = []
    t0 = time.time()
    for fold_i, held_out_new in enumerate(pos_new_all, start=1):
        held_out_name = pos_name_by_new[held_out_new]
        train_pos = [n for n in pos_new_all if n != held_out_new]
        t_fold = time.time()
        drug_indices, scores = train_once(
            data, disease_new, train_pos, neg_pool_base, args.epochs, args.hidden, args.lr,
            encoder=args.encoder)
        ranked = rank_table(data, nodes, drug_indices, scores)
        row = ranked[ranked["node_index"] == held_out_new]
        if row.empty:
            print(f"  [{fold_i}/{len(pos_new_all)}] {held_out_name:25} ERROR: held-out node missing from score output")
            continue
        r = row.iloc[0]
        elapsed = time.time() - t_fold
        loo_rows.append({
            "drug": held_out_name,
            "node_index": held_out_new,
            "held_out_score": r["score"],
            "held_out_rank": int(r["rank"]),
            "total_drugs": len(ranked),
            "held_out_percentile": r["percentile"],
        })
        print(f"  [{fold_i:2}/{len(pos_new_all)}] {held_out_name:25} rank={int(r['rank']):5}/{len(ranked)} "
              f"({r['percentile']:.2f}th pctile)  score={r['score']:.4f}  ({elapsed:.0f}s)")

    loo_df = pd.DataFrame(loo_rows).sort_values("held_out_percentile")
    loo_path = RESULTS / f"kg_loo_recall_{args.encoder}.tsv"
    loo_df.to_csv(loo_path, sep="\t", index=False)
    print(f"\nLOO recall results saved to {loo_path.relative_to(REPO_ROOT)}")
    in_top10 = (loo_df["held_out_percentile"] <= 10.0).sum()
    print(f"{in_top10}/{len(loo_df)} known drugs land in top 10% of all drugs when genuinely held out")
    print(f"Phase 1 total time: {time.time() - t0:.0f}s\n")

    # ---- Phase 2: one full model (all positives, no holdout) for candidate scoring ----
    print("=" * 70)
    print(f"Phase 2: full model (all known positives trained) for novel-candidate ranking, encoder={args.encoder}")
    print("=" * 70)
    print("(Novel candidates were never positive labels, so scoring them with a model")
    print(" trained on all known positives is not circular - only re-scoring a KNOWN")
    print(" drug with a model trained on itself would be.)\n")
    drug_indices, scores = train_once(
        data, disease_new, pos_new_all, neg_pool_base, args.epochs, args.hidden, args.lr,
        encoder=args.encoder)
    full_ranked = rank_table(data, nodes, drug_indices, scores)
    full_ranked["is_known_positive"] = full_ranked["drug"].isin(known_drugs)
    full_path = RESULTS / f"kg_full_candidate_ranked_{args.encoder}.tsv"
    full_ranked.to_csv(full_path, sep="\t", index=False)
    print(f"Full-model candidate rankings saved to {full_path.relative_to(REPO_ROOT)} ({len(full_ranked)} drugs)\n")

    print("Final candidates (results/final_candidates.tsv) - rank under this properly-trained full model:")
    final_candidates_path = ROOT / "results" / "final_candidates.tsv"
    if final_candidates_path.exists():
        p01_candidates = pd.read_csv(final_candidates_path, sep="\t")["drug"].tolist()
    else:
        p01_candidates = [
            "nintedanib", "bosutinib", "masitinib", "tranilast", "adaprev", "dasatinib",
            "f351", "bafetinib", "tolimidone", "rebastinib", "ginsenoside-re", "ponatinib",
            "d-4476", "ft011", "on123300",
        ]
    fr_lower = full_ranked.copy()
    fr_lower["drug_lower"] = fr_lower["drug"].str.lower()
    for cand in p01_candidates:
        match = fr_lower[fr_lower["drug_lower"] == cand.lower()]
        if match.empty:
            print(f"  {cand:18} ABSENT from current PrimeKG's drug nodes")
        else:
            r = match.iloc[0]
            print(f"  {cand:18} rank={int(r['rank']):5}/{len(full_ranked)} "
                  f"({r['percentile']:.2f}th pctile)  score={r['score']:.4f}")

    print("\nComparison against 06_gnn_train.py's single-split sparsentan/atrasentan result:")
    for cand in ["sparsentan", "atrasentan"]:
        loo_row = loo_df[loo_df["drug"].str.lower() == cand]
        if not loo_row.empty:
            r = loo_row.iloc[0]
            print(f"  {cand:18} genuine LOO-held-out rank={int(r['held_out_rank'])}/{int(r['total_drugs'])} "
                  f"({r['held_out_percentile']:.2f}th pctile)")
        else:
            print(f"  {cand:18} not found in LOO results")


if __name__ == "__main__":
    main()
