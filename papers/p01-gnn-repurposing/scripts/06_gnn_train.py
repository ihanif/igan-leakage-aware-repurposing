"""
Step 6: GNN link prediction on PrimeKG.

Trains a relational GCN (R-GCN) link predictor on the PrimeKG knowledge graph
to predict drug-disease associations for IgAN. Uses a retrospective validation
design: known IgAN drug-disease edges are masked during training and used as a
held-out test set.

Architecture:
    - 2-layer RGCNConv encoder — one weight matrix per relation type, preserving
      the semantic distinction between drug_protein, protein_disease, etc.
    - Dot-product link prediction decoder
    - Binary cross-entropy loss (positive = known IgAN drug, negative = random drug)
    - Trained on drug-protein and protein-disease subgraph of PrimeKG

Usage:
    uv run python scripts/06_gnn_train.py [--epochs N] [--hidden D] [--lr LR]

    Defaults: --epochs 200, --hidden 128, --lr 0.001

Output:
    results/gnn_ranked.tsv

Note: this script's single random 80/20 split reports rank for ALL known drugs,
including the 80% the model was trained on. That is not a recall result for
those drugs. The manuscript's leave-one-out cross-validation and R-GCN vs
GraphSAGE architecture ablation (Methods, GNN link prediction; Results,
"13 of 15 known drugs in the top 10%... versus 10 of 15") are produced by
06b_gnn_loo_validation.py instead, which trains a fresh model per held-out
drug so no drug is ever ranked by a model that saw its own label.

Validation checkpoint:
    AUC on held-out IgAN-drug links >= 0.70.
    At least 3 of 5 approved IgAN drugs in top 10% of all scored drugs.

Data requirements:
    data/raw/primekg/nodes.tab  (from Harvard Dataverse, file ID 6180617)
    data/raw/primekg/edges.csv  (from Harvard Dataverse, file ID 6180616, ~370 MB)

    Download:
        curl -L "https://dataverse.harvard.edu/api/access/datafile/6180617" -o data/raw/primekg/nodes.tab
        curl -L "https://dataverse.harvard.edu/api/access/datafile/6180616" -o data/raw/primekg/edges.csv
"""

import argparse
import random
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "shared/scripts"))
from pipeline_utils import find_repo_root, load_config

import numpy as np
import pandas as pd

# PrimeKG edge types to keep (covers drug-target and gene-disease associations)
RELEVANT_RELATIONS = {
    "drug_protein",        # drug → protein (target interaction)
    "protein_drug",        # protein → drug (reverse)
    "protein_protein",     # PPI
    "protein_disease",     # gene → disease association
    "disease_protein",     # disease → gene (reverse)
    "drug_effect",         # drug → side effect / phenotype
}

# Stable integer index for each relation type (used by RGCNConv)
RELATION_TO_IDX: dict[str, int] = {rel: i for i, rel in enumerate(sorted(RELEVANT_RELATIONS))}
NUM_RELATIONS = len(RELATION_TO_IDX)

DEFAULT_CONFIG = "shared/configs/igan.yaml"


def check_dependencies() -> bool:
    try:
        import torch
        import torch_geometric
        return True
    except ImportError:
        return False


def install_dependencies():
    import subprocess
    print("Installing PyTorch and PyTorch Geometric...")
    subprocess.run(["uv", "pip", "install", "torch", "--quiet"], check=True)
    subprocess.run(
        ["uv", "pip", "install", "torch-geometric", "--quiet"],
        check=True,
    )


def load_primekg(nodes_path: Path, edges_path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    print("Loading PrimeKG nodes...")
    nodes = pd.read_csv(nodes_path, sep="\t")
    nodes["node_index"] = nodes["node_index"].astype(int)
    print(f"  {len(nodes):,} nodes ({nodes['node_type'].value_counts().to_dict()})")

    print("Loading PrimeKG edges (may take ~30s for 4M rows)...")
    edges = pd.read_csv(edges_path)
    edges["x_index"] = edges["x_index"].astype(int)
    edges["y_index"] = edges["y_index"].astype(int)
    print(f"  {len(edges):,} edges ({edges['relation'].value_counts().head(8).to_dict()})")

    return nodes, edges


def load_drug_molecular_features(drug_feat_path: Path, drug_node_indices: list[int]) -> dict[int, list[float]]:
    """
    Load PrimeKG drug_features.tab and extract numerical molecular descriptors
    (molecular_weight, tpsa, clogp). Returns {node_index: [mw, tpsa, clogp]}.
    Values are z-score normalised; missing filled with 0.
    """
    if not drug_feat_path.exists():
        return {}
    df = pd.read_csv(drug_feat_path, sep="\t", usecols=["node_index", "molecular_weight", "tpsa", "clogp"])
    df["node_index"] = df["node_index"].astype(int)
    df = df[df["node_index"].isin(set(drug_node_indices))]

    for col in ["molecular_weight", "tpsa", "clogp"]:
        df[col] = pd.to_numeric(df[col], errors="coerce").fillna(0.0)
        mu, sd = df[col].mean(), df[col].std()
        if sd > 0:
            df[col] = (df[col] - mu) / sd

    result = {}
    for _, row in df.iterrows():
        result[int(row["node_index"])] = [row["molecular_weight"], row["tpsa"], row["clogp"]]
    return result


def build_drug_graph(nodes: pd.DataFrame, edges: pd.DataFrame, disease_idx: int,
                     drug_feat_path: Path):
    """
    Build graph for drug-disease link prediction.

    Node features: learnable embedding index (nn.Embedding) + 3 molecular descriptors.
    Edges: drug-protein, protein-protein, protein-disease (filtered to active nodes).
    Task: predict drug → disease node.
    """
    import torch
    from torch_geometric.data import Data

    rel_edges = edges[edges["relation"].isin(RELEVANT_RELATIONS)].copy()
    print(f"  Relevant edges: {len(rel_edges):,}")

    drug_nodes = nodes[nodes["node_type"] == "drug"]["node_index"].tolist()
    gene_nodes = nodes[nodes["node_type"] == "gene/protein"]["node_index"].tolist()
    active_nodes = set(drug_nodes) | set(gene_nodes) | {disease_idx}

    mask = rel_edges["x_index"].isin(active_nodes) & rel_edges["y_index"].isin(active_nodes)
    rel_edges = rel_edges[mask]
    print(f"  Edges after node filtering: {len(rel_edges):,}")

    all_active = sorted(active_nodes)
    idx_map = {orig: new for new, orig in enumerate(all_active)}
    n_nodes = len(all_active)

    src = [idx_map[x] for x in rel_edges["x_index"]]
    dst = [idx_map[y] for y in rel_edges["y_index"]]
    rel_type = [RELATION_TO_IDX[r] for r in rel_edges["relation"]]
    # Undirected: each edge appears in both directions; relation type is preserved per direction
    edge_index = torch.tensor([src + dst, dst + src], dtype=torch.long)
    edge_type = torch.tensor(rel_type + rel_type, dtype=torch.long)

    mol_feats = load_drug_molecular_features(drug_feat_path, drug_nodes)
    print(f"  Molecular features loaded for {len(mol_feats)}/{len(drug_nodes)} drug nodes")

    feat_dim_total = 3 + 3  # type one-hot [gene, drug, disease] + mol features
    x = torch.zeros(n_nodes, feat_dim_total)
    drug_set = set(drug_nodes)
    for orig in all_active:
        i = idx_map[orig]
        if orig in drug_set:
            x[i, 1] = 1.0
            mol = mol_feats.get(orig, [0.0, 0.0, 0.0])
            x[i, 3], x[i, 4], x[i, 5] = mol[0], mol[1], mol[2]
        elif orig == disease_idx:
            x[i, 2] = 1.0
        else:
            x[i, 0] = 1.0

    node_ids = torch.arange(n_nodes, dtype=torch.long)
    data = Data(x=x, edge_index=edge_index, edge_type=edge_type)
    data.node_ids = node_ids
    data.n_nodes = n_nodes
    data.idx_map = idx_map
    data.rev_map = {v: k for k, v in idx_map.items()}
    data.drug_nodes_remapped = [idx_map[n] for n in drug_nodes if n in idx_map]
    data.disease_remapped = idx_map[disease_idx]
    return data


def find_known_disease_edges(nodes: pd.DataFrame, known_names: set[str]) -> list[int]:
    """Return PrimeKG node indices of known drugs (positive training examples)."""
    drug_nodes = nodes[
        (nodes["node_type"] == "drug") &
        nodes["node_name"].isin(known_names)
    ]
    return drug_nodes["node_index"].tolist()


def train_gnn(data, known_pos_indices: list[int], n_epochs: int, hidden: int, lr: float):
    """Train R-GCN link predictor."""
    import torch
    import torch.nn as nn
    import torch.nn.functional as F
    from torch_geometric.nn import RGCNConv

    igan_new = data.disease_remapped
    drug_new = set(data.drug_nodes_remapped)

    # Positive examples: known IgAN drugs (remapped)
    pos_drugs = [data.idx_map[orig] for orig in known_pos_indices if orig in data.idx_map]
    neg_pool = [n for n in drug_new if n not in set(pos_drugs)]

    if not pos_drugs:
        print("WARNING: No known IgAN drugs found in graph. Check drug name matching.", file=sys.stderr)
        return None

    print(f"  Positive examples: {len(pos_drugs)} known IgAN drugs")
    print(f"  Negative pool: {len(neg_pool)} other drugs")

    EMB_DIM = 32  # learnable embedding dimension per node

    class LinkPredictor(nn.Module):
        def __init__(self, n_nodes: int, feat_channels: int, hidden_channels: int):
            super().__init__()
            # Learnable node embeddings — each node gets its own vector
            self.node_emb = nn.Embedding(n_nodes, EMB_DIM)
            nn.init.xavier_uniform_(self.node_emb.weight)
            # R-GCN operates on concatenation of embedding + structural features
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
            h = self.conv2(h, edge_index, edge_type)
            return h

        def decode(self, z, src_idx: int, tgt_indices: list[int]):
            src_emb = z[src_idx].unsqueeze(0).expand(len(tgt_indices), -1)
            tgt_emb = z[tgt_indices]
            return self.decoder(torch.cat([src_emb, tgt_emb], dim=-1)).squeeze(-1)

    model = LinkPredictor(data.n_nodes, data.x.size(1), hidden)
    optimizer = torch.optim.Adam(model.parameters(), lr=lr)

    # 80/20 train/test split on positives
    random.shuffle(pos_drugs)
    split = max(1, int(len(pos_drugs) * 0.8))
    train_pos = pos_drugs[:split]
    test_pos = pos_drugs[split:]

    print(f"  Train positives: {len(train_pos)}, Test positives: {len(test_pos)}")

    best_auc = 0.0
    best_state = None

    for epoch in range(1, n_epochs + 1):
        model.train()
        optimizer.zero_grad()

        z = model.encode(data.x, data.edge_index, data.edge_type, data.node_ids)

        # Sample equal negatives
        neg_sample = random.sample(neg_pool, min(len(train_pos) * 3, len(neg_pool)))

        pos_scores = model.decode(z, igan_new, train_pos)
        neg_scores = model.decode(z, igan_new, neg_sample)

        labels = torch.cat([
            torch.ones(len(train_pos)),
            torch.zeros(len(neg_sample)),
        ])
        scores = torch.cat([pos_scores, neg_scores])
        loss = F.binary_cross_entropy_with_logits(scores, labels)
        loss.backward()
        optimizer.step()

        if epoch % 20 == 0:
            model.eval()
            with torch.no_grad():
                z = model.encode(data.x, data.edge_index, data.edge_type, data.node_ids)
                if test_pos:
                    test_neg = random.sample(neg_pool, min(len(test_pos) * 5, len(neg_pool)))
                    tp_scores = torch.sigmoid(model.decode(z, igan_new, test_pos)).numpy()
                    tn_scores = torch.sigmoid(model.decode(z, igan_new, test_neg)).numpy()
                    all_scores = np.concatenate([tp_scores, tn_scores])
                    all_labels = np.concatenate([np.ones(len(tp_scores)), np.zeros(len(tn_scores))])
                    # Compute AUC manually
                    order = np.argsort(-all_scores)
                    auc = 0.0
                    n_pos = int(all_labels.sum())
                    n_neg = len(all_labels) - n_pos
                    if n_pos > 0 and n_neg > 0:
                        tp_cum = 0
                        for label in all_labels[order]:
                            if label == 1:
                                tp_cum += 1
                            else:
                                auc += tp_cum
                        auc /= (n_pos * n_neg)
                    if auc > best_auc:
                        best_auc = auc
                        best_state = {k: v.clone() for k, v in model.state_dict().items()}
                    print(f"    Epoch {epoch:3d}  loss={loss.item():.4f}  AUC={auc:.4f}  (best={best_auc:.4f})")
                else:
                    print(f"    Epoch {epoch:3d}  loss={loss.item():.4f}")

    if best_state:
        model.load_state_dict(best_state)
    return model, best_auc


def score_all_drugs(model, data, nodes: pd.DataFrame) -> pd.DataFrame:
    """Score all drug nodes for IgAN link prediction."""
    import torch

    model.eval()
    with torch.no_grad():
        z = model.encode(data.x, data.edge_index, data.edge_type, data.node_ids)
        drug_indices = data.drug_nodes_remapped
        scores = torch.sigmoid(model.decode(z, data.disease_remapped, drug_indices)).numpy()

    # Map back to drug names
    orig_indices = [data.rev_map[i] for i in drug_indices]
    name_map = nodes.set_index("node_index")["node_name"].to_dict()
    drug_names = [name_map.get(orig, str(orig)) for orig in orig_indices]

    df = pd.DataFrame({
        "drug": drug_names,
        "node_index": orig_indices,
        "gnn_score": scores,
    }).sort_values("gnn_score", ascending=False).reset_index(drop=True)
    df["rank"] = range(1, len(df) + 1)
    return df


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default=DEFAULT_CONFIG,
                        help="Path to disease config YAML (relative to repo root or absolute)")
    parser.add_argument("--epochs", type=int, default=200)
    parser.add_argument("--hidden", type=int, default=128)
    parser.add_argument("--lr", type=float, default=0.001)
    parser.add_argument("--seed", type=int, default=42,
                        help="Random seed for reproducibility")
    args = parser.parse_args()

    random.seed(args.seed)
    np.random.seed(args.seed)
    try:
        import torch
        torch.manual_seed(args.seed)
    except ImportError:
        pass

    cfg, repo_root = load_config(args.config)
    rep = cfg["repurposing"]
    disease_id   = cfg["disease"]["id"]
    disease_name = cfg["disease"]["name"]

    output_dir      = repo_root / rep["output_dir"]
    shared_data_dir = repo_root / rep["shared_data_dir"]
    results_dir     = output_dir / "results"
    results_dir.mkdir(parents=True, exist_ok=True)

    nodes_path    = shared_data_dir / "primekg/nodes.tab"
    edges_path    = shared_data_dir / "primekg/edges.csv"
    drug_feat_path = shared_data_dir / "primekg/drug_features.tab"
    out_path      = results_dir / "gnn_ranked.tsv"

    disease_node_idx   = rep["primekg_node_index"]
    known_drugs_drugbank = set(rep.get("known_drugs_drugbank", []))

    print(f"Disease: {disease_name} ({disease_id})")
    print(f"PrimeKG node index: {disease_node_idx}")
    print(f"Output:  {output_dir}\n")

    for path, label in [(nodes_path, "nodes.tab"), (edges_path, "edges.csv")]:
        if not path.exists():
            print(f"PrimeKG {label} not found at {path}", file=sys.stderr)
            print("Download:", file=sys.stderr)
            print("  curl -L 'https://dataverse.harvard.edu/api/access/datafile/6180617'"
                  f" -o {shared_data_dir}/primekg/nodes.tab", file=sys.stderr)
            print("  curl -L 'https://dataverse.harvard.edu/api/access/datafile/6180616'"
                  f" -o {shared_data_dir}/primekg/edges.csv", file=sys.stderr)
            sys.exit(1)

    if not check_dependencies():
        install_dependencies()

    nodes, edges = load_primekg(nodes_path, edges_path)

    print(f"\nBuilding drug-protein-disease subgraph for {disease_name}...")
    data = build_drug_graph(nodes, edges, disease_node_idx, drug_feat_path)
    print(f"  Graph: {data.n_nodes:,} nodes, {data.edge_index.size(1)//2:,} edges (undirected)")

    known_pos = find_known_disease_edges(nodes, known_drugs_drugbank)
    print(f"  Known {disease_id.upper()} drugs in PrimeKG: {len(known_pos)}")

    print(f"\nTraining R-GCN ({args.epochs} epochs, hidden={args.hidden}, lr={args.lr}, "
          f"num_relations={NUM_RELATIONS})...")
    result = train_gnn(data, known_pos, args.epochs, args.hidden, args.lr)
    if result is None:
        sys.exit(1)
    model, best_auc = result

    print(f"\nBest test AUC: {best_auc:.4f}")
    if best_auc < 0.70:
        print("WARNING: AUC < 0.70. Consider more epochs or a larger hidden dimension.",
              file=sys.stderr)

    print("\nScoring all drug nodes...")
    ranked = score_all_drugs(model, data, nodes)
    ranked.to_csv(out_path, sep="\t", index=False)
    print(f"GNN rankings saved to {out_path} ({len(ranked)} drugs)")

    print(f"\nValidation checkpoint — known {disease_name} drugs:")
    total = len(ranked)
    top10_cutoff = ranked["gnn_score"].quantile(0.90)
    in_top10 = 0
    for drug in known_drugs_drugbank:
        match = ranked[ranked["drug"].str.lower() == drug.lower()]
        if not match.empty:
            r = match.iloc[0]
            flag = "TOP10%" if r["gnn_score"] >= top10_cutoff else ""
            if flag:
                in_top10 += 1
            print(f"  {drug:30} rank={int(r['rank']):5}/{total}  score={r['gnn_score']:.4f}  {flag}")
    print(f"\n{in_top10} known {disease_id.upper()} drugs in top 10% (target: >=3)")

    print("\nTop 20 novel candidates:")
    known_lower = {d.lower() for d in known_drugs_drugbank}
    novel = ranked[~ranked["drug"].str.lower().isin(known_lower)].head(20)
    print(novel[["rank", "drug", "gnn_score"]].to_string(index=False))


if __name__ == "__main__":
    main()
