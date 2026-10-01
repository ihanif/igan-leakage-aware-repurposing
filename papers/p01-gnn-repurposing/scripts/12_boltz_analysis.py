"""
Step 12: Post-prediction analysis of Boltz-2 protein-ligand structures.

Generates three outputs from the CIF predictions produced by 11_boltz_predict.py:

1. Figure: ligand_ipTM bar chart (GNN-ranked, per-model dots + mean bar)
2. Table: contact residues per drug (protein residues within 4 Å of ligand heavy atoms)
3. Table: inter-model ligand RMSD (binding pose consistency across 3 diffusion samples)

All outputs written to results/boltz/:
  figures/fig_boltz2_iptm.png          — ipTM bar chart (manuscript Figure)
  figures/fig_boltz2_iptm.pdf          — vector version
  contact_residues.tsv                 — per-drug binding-site residue table
  inter_model_rmsd.tsv                 — per-drug RMSD between model pairs

Usage:
    uv run python papers/p01-gnn-repurposing/scripts/12_boltz_analysis.py
"""

import json
import math
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "shared/scripts"))
from pipeline_utils import find_repo_root, load_config

DEFAULT_CONFIG = "shared/configs/igan.yaml"

# Canonical drug order (by GNN/pipeline final_rank, not Boltz rank)
DRUGS = [
    {"drug": "nintedanib", "target": "LYN",   "final_rank": 1, "cluster": "LYN"},
    {"drug": "bosutinib",  "target": "LYN",   "final_rank": 2, "cluster": "LYN"},
    {"drug": "masitinib",  "target": "LYN",   "final_rank": 3, "cluster": "LYN"},
    {"drug": "tranilast",  "target": "TGFB1", "final_rank": 4, "cluster": "TGFB1"},
    {"drug": "dasatinib",  "target": "LYN",   "final_rank": 6, "cluster": "LYN"},
]

CLUSTER_COLOURS = {"LYN": "#2166ac", "TGFB1": "#d73027"}


# ---------------------------------------------------------------------------
# mmCIF atom_site parser (pure Python, no biopython dependency)
# ---------------------------------------------------------------------------

def parse_atom_site(cif_path: Path) -> list[dict]:
    """Parse _atom_site loop from an mmCIF file into a list of dicts."""
    text = cif_path.read_text()
    # Locate the _atom_site loop
    loop_start = text.find("loop_\n_atom_site.")
    if loop_start == -1:
        raise ValueError(f"No _atom_site loop in {cif_path}")
    block = text[loop_start:]
    lines = block.splitlines()

    # Collect column names
    cols = []
    data_start = 0
    for i, line in enumerate(lines):
        line_stripped = line.strip()
        if line_stripped.startswith("_atom_site."):
            cols.append(line_stripped.removeprefix("_atom_site."))
        elif cols and not line_stripped.startswith("_") and line_stripped and line_stripped != "loop_":
            data_start = i
            break

    atoms = []
    for line in lines[data_start:]:
        line_stripped = line.strip()
        if not line_stripped or line_stripped.startswith("#") or line_stripped.startswith("loop_") or line_stripped.startswith("_"):
            break
        parts = line_stripped.split()
        if len(parts) == len(cols):
            atoms.append(dict(zip(cols, parts)))
    return atoms


def ligand_heavy_atoms(atoms: list[dict]) -> list[tuple[float, float, float]]:
    """Return (x, y, z) for all non-hydrogen atoms in the ligand chain (label_entity_id=2)."""
    coords = []
    for a in atoms:
        if a.get("label_entity_id") == "2" and a.get("type_symbol") != "H":
            coords.append((float(a["Cartn_x"]), float(a["Cartn_y"]), float(a["Cartn_z"])))
    return coords


def protein_residues_near_ligand(atoms: list[dict], ligand_coords: list[tuple], cutoff: float = 4.0) -> list[str]:
    """Return sorted list of unique 'RES_SEQID' labels for protein residues within cutoff Å of any ligand heavy atom."""
    contact_set = set()
    lig_arr = np.array(ligand_coords)  # (N_lig, 3)
    for a in atoms:
        if a.get("label_entity_id") != "1":
            continue
        if a.get("type_symbol") == "H":
            continue
        pa = np.array([float(a["Cartn_x"]), float(a["Cartn_y"]), float(a["Cartn_z"])])
        dists = np.linalg.norm(lig_arr - pa, axis=1)
        if dists.min() <= cutoff:
            res_id = f"{a['label_comp_id']}_{a['label_seq_id']}"
            contact_set.add(res_id)
    return sorted(contact_set, key=lambda r: int(r.split("_")[1]))


def ligand_rmsd(coords_a: list[tuple], coords_b: list[tuple]) -> float | None:
    """RMSD between two sets of matched ligand heavy-atom coordinates."""
    if len(coords_a) != len(coords_b) or not coords_a:
        return None
    arr_a = np.array(coords_a)
    arr_b = np.array(coords_b)
    return float(np.sqrt(np.mean(np.sum((arr_a - arr_b) ** 2, axis=1))))


# ---------------------------------------------------------------------------
# Confidence loading
# ---------------------------------------------------------------------------

def load_confidence(pred_dir: Path, drug: str, target: str, model_idx: int) -> dict:
    job = f"{drug}_{target}"
    p = pred_dir / f"boltz_results_{job}" / "predictions" / job / f"confidence_{job}_model_{model_idx}.json"
    if not p.exists():
        return {}
    with open(p) as f:
        return json.load(f)


# ---------------------------------------------------------------------------
# Task 1: ipTM bar chart
# ---------------------------------------------------------------------------

def make_iptm_figure(pred_dir: Path, out_dir: Path) -> None:
    n_models = 3
    drug_labels, means, per_model, colours, ranks = [], [], [], [], []

    for d in DRUGS:
        drug, target = d["drug"], d["target"]
        vals = []
        for m in range(n_models):
            conf = load_confidence(pred_dir, drug, target, m)
            vals.append(conf.get("ligand_iptm", conf.get("iptm", float("nan"))))
        mean_val = float(np.nanmean(vals))
        drug_labels.append(f"{drug.capitalize()}\n(rank {d['final_rank']}, {d['cluster']})")
        means.append(mean_val)
        per_model.append(vals)
        colours.append(CLUSTER_COLOURS[d["cluster"]])
        ranks.append(d["final_rank"])

    # Sort descending by mean ipTM for visual clarity
    order = sorted(range(len(means)), key=lambda i: means[i], reverse=True)

    fig, ax = plt.subplots(figsize=(8, 5))
    y_pos = np.arange(len(order))

    for yi, oi in enumerate(order):
        ax.barh(yi, means[oi], color=colours[oi], alpha=0.85, height=0.55, zorder=2)
        # Per-model dots
        for mv in per_model[oi]:
            if not math.isnan(mv):
                ax.plot(mv, yi, "o", color="white", markeredgecolor=colours[oi],
                        markeredgewidth=1.2, markersize=5, zorder=3)

    ax.set_yticks(y_pos)
    ax.set_yticklabels([drug_labels[i] for i in order], fontsize=10)
    ax.set_xlabel("Boltz-2 ligand ipTM (mean of 3 diffusion samples)", fontsize=10)
    ax.set_xlim(0, 1.05)
    ax.axvline(0.8, color="gray", linestyle="--", linewidth=0.8, label="ipTM = 0.8 (good)")
    ax.axvline(0.5, color="lightgray", linestyle=":", linewidth=0.8, label="ipTM = 0.5 (marginal)")
    ax.set_title("Boltz-2 structural binding confidence\nIgAN top-5 drug-protein pairs", fontsize=11)
    ax.grid(axis="x", alpha=0.3, zorder=1)

    lyn_patch = mpatches.Patch(color=CLUSTER_COLOURS["LYN"], alpha=0.85, label="LYN cluster")
    tgf_patch = mpatches.Patch(color=CLUSTER_COLOURS["TGFB1"], alpha=0.85, label="TGFB1 cluster")
    ax.legend(handles=[lyn_patch, tgf_patch], loc="lower right", fontsize=9)

    plt.tight_layout()
    out_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_dir / "fig_boltz2_iptm.png", dpi=200, bbox_inches="tight")
    fig.savefig(out_dir / "fig_boltz2_iptm.pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out_dir / 'fig_boltz2_iptm.png'}")


# ---------------------------------------------------------------------------
# Task 2: Contact residues
# ---------------------------------------------------------------------------

def make_contact_table(pred_dir: Path, out_dir: Path) -> None:
    rows = [["drug", "target", "model", "n_contacts", "contact_residues"]]
    for d in DRUGS:
        drug, target = d["drug"], d["target"]
        for m in range(3):
            cif_path = (pred_dir / f"boltz_results_{drug}_{target}"
                        / "predictions" / f"{drug}_{target}"
                        / f"{drug}_{target}_model_{m}.cif")
            if not cif_path.exists():
                rows.append([drug, target, str(m), "NA", "NA"])
                continue
            atoms = parse_atom_site(cif_path)
            lig_coords = ligand_heavy_atoms(atoms)
            contacts = protein_residues_near_ligand(atoms, lig_coords, cutoff=4.0)
            rows.append([drug, target, str(m), str(len(contacts)), ";".join(contacts)])
    out_path = out_dir / "contact_residues.tsv"
    with open(out_path, "w") as f:
        for row in rows:
            f.write("\t".join(row) + "\n")
    print(f"  Saved: {out_path}")


# ---------------------------------------------------------------------------
# Task 3: Inter-model RMSD
# ---------------------------------------------------------------------------

def contact_jaccard(contacts_a: list[str], contacts_b: list[str]) -> float:
    """Jaccard similarity between two contact residue lists."""
    a, b = set(contacts_a), set(contacts_b)
    if not a and not b:
        return 1.0
    return len(a & b) / len(a | b)


def make_rmsd_table(pred_dir: Path, out_dir: Path) -> None:
    """
    Raw ligand RMSD across Boltz-2 diffusion samples requires backbone superposition
    to be interpretable (each sample has an independent global coordinate frame).
    Contact residue Jaccard similarity is the appropriate binding-mode consistency
    metric here — it measures whether the same binding pocket is engaged regardless
    of absolute orientation.
    """
    rows = [["drug", "target",
             "contact_jaccard_01", "contact_jaccard_02", "contact_jaccard_12", "mean_jaccard",
             "n_contacts_m0", "n_contacts_m1", "n_contacts_m2",
             "raw_rmsd_01_unaligned", "raw_rmsd_02_unaligned", "raw_rmsd_12_unaligned",
             "note"]]
    for d in DRUGS:
        drug, target = d["drug"], d["target"]
        cif_data = []  # list of (coords, contacts) or None
        for m in range(3):
            cif_path = (pred_dir / f"boltz_results_{drug}_{target}"
                        / "predictions" / f"{drug}_{target}"
                        / f"{drug}_{target}_model_{m}.cif")
            if cif_path.exists():
                atoms = parse_atom_site(cif_path)
                coords = ligand_heavy_atoms(atoms)
                contacts = protein_residues_near_ligand(atoms, coords, cutoff=4.0)
                cif_data.append((coords, contacts))
            else:
                cif_data.append(None)

        # Contact Jaccard
        jaccards = []
        for a, b in [(0, 1), (0, 2), (1, 2)]:
            if cif_data[a] is not None and cif_data[b] is not None:
                j = contact_jaccard(cif_data[a][1], cif_data[b][1])
                jaccards.append(j)
            else:
                jaccards.append(None)
        valid_j = [j for j in jaccards if j is not None]
        mean_j = round(float(np.mean(valid_j)), 3) if valid_j else None

        # Raw unaligned RMSD (for completeness; requires superposition to be meaningful)
        rmsds = []
        for a, b in [(0, 1), (0, 2), (1, 2)]:
            if cif_data[a] is not None and cif_data[b] is not None:
                r = ligand_rmsd(cif_data[a][0], cif_data[b][0])
                rmsds.append(r)
            else:
                rmsds.append(None)

        n_contacts = [str(len(cif_data[m][1])) if cif_data[m] is not None else "NA" for m in range(3)]

        rows.append([
            drug, target,
            f"{jaccards[0]:.3f}" if jaccards[0] is not None else "NA",
            f"{jaccards[1]:.3f}" if jaccards[1] is not None else "NA",
            f"{jaccards[2]:.3f}" if jaccards[2] is not None else "NA",
            f"{mean_j:.3f}" if mean_j is not None else "NA",
        ] + n_contacts + [
            f"{rmsds[0]:.2f}" if rmsds[0] is not None else "NA",
            f"{rmsds[1]:.2f}" if rmsds[1] is not None else "NA",
            f"{rmsds[2]:.2f}" if rmsds[2] is not None else "NA",
            "raw_unaligned_rmsd_requires_superposition",
        ])

        mj = f"{mean_j:.3f}" if mean_j is not None else "NA"
        print(f"  {drug}_{target}: mean contact Jaccard={mj} "
              f"(contacts: {n_contacts[0]}/{n_contacts[1]}/{n_contacts[2]})")

    out_path = out_dir / "inter_model_rmsd.tsv"
    with open(out_path, "w") as f:
        for row in rows:
            f.write("\t".join(row) + "\n")
    print(f"  Saved: {out_path}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default=DEFAULT_CONFIG)
    args = parser.parse_args()

    _, repo_root = load_config(args.config)
    base = repo_root / "papers/p01-gnn-repurposing/results/boltz"
    pred_dir = base / "predictions"
    fig_dir = base / "figures"

    print("=== Boltz-2 post-prediction analysis ===\n")

    print("[1/3] Generating ipTM bar chart...")
    make_iptm_figure(pred_dir, fig_dir)

    print("\n[2/3] Extracting contact residues (4 Å cutoff)...")
    make_contact_table(pred_dir, base)

    print("\n[3/3] Computing inter-model ligand RMSD...")
    make_rmsd_table(pred_dir, base)

    print("\nDone. Outputs:")
    print(f"  {fig_dir}/fig_boltz2_iptm.{{png,pdf}}")
    print(f"  {base}/contact_residues.tsv")
    print(f"  {base}/inter_model_rmsd.tsv")


if __name__ == "__main__":
    main()
