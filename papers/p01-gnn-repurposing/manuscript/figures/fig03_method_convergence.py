"""Generate the method-convergence evidence heatmap. Built entirely from existing result files:
results/final_candidates.tsv, results/kg_full_candidate_ranked_rgcn.tsv (the full R-GCN model
from scripts/06b_gnn_loo_validation.py, trained on all known positives), and
results/cmap_analysis/cmap_scores.tsv.
No numbers are invented; drugs absent from a given scoring universe (e.g. not in LINCS L1000)
are shown as a distinct "not scored" cell rather than assigned a value.
"""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.patches as mpatches
import pandas as pd
import numpy as np
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).parent

final = pd.read_csv(ROOT / "results" / "final_candidates.tsv", sep="\t")
gnn = pd.read_csv(ROOT / "results" / "kg_full_candidate_ranked_rgcn.tsv", sep="\t")
cmap = pd.read_csv(ROOT / "results" / "cmap_analysis" / "cmap_scores.tsv", sep="\t")

GNN_TOTAL = 7957
CMAP_TOTAL = 3344

final = final.sort_values("final_rank").reset_index(drop=True)
gnn_lookup = {d.lower(): r for d, r in zip(gnn["drug"], gnn["rank"])}
cmap_lookup = {d.lower(): r for d, r in zip(cmap["drug"], cmap["rank"])}

rows = []
for _, row in final.iterrows():
    drug = row["drug"]
    key = drug.lower()
    np_pct = row["np_rank"] / 1989 * 100
    ot_pct = row["ot_rank"] / 3005 * 100
    gnn_rank = gnn_lookup.get(key)
    gnn_pct = (gnn_rank / GNN_TOTAL * 100) if gnn_rank is not None else None
    cmap_rank = cmap_lookup.get(key)
    cmap_pct = (cmap_rank / CMAP_TOTAL * 100) if cmap_rank is not None else None
    rows.append({
        "drug": drug,
        "cluster": "LYN" if row["target_gene"] == "LYN" else ("TGF-β" if row["target_gene"] == "TGFB1" else row["target_gene"]),
        "NP": np_pct,
        "OT": ot_pct,
        "GNN": gnn_pct,
        "CMap": cmap_pct,
    })

df = pd.DataFrame(rows)

methods = ["NP", "OT", "GNN", "CMap"]
mat = df[methods].to_numpy(dtype=float)  # NaN where not scored

fig, ax = plt.subplots(figsize=(6.0, 5.2))

cmap_colors = LinearSegmentedColormap.from_list("evidence", ["#08519c", "#c6dbef", "#f7f7f7"])
masked = np.ma.masked_invalid(mat)
im = ax.imshow(masked, cmap=cmap_colors, vmin=0, vmax=30, aspect="auto")

# grey hatch for not-scored cells
for i in range(mat.shape[0]):
    for j in range(mat.shape[1]):
        if np.isnan(mat[i, j]):
            ax.add_patch(mpatches.Rectangle((j - 0.5, i - 0.5), 1, 1, facecolor="#e0e0e0",
                                             edgecolor="white", hatch="////", linewidth=0.5))
        else:
            val = mat[i, j]
            txt_color = "white" if val < 12 else "black"
            ax.text(j, i, f"{val:.1f}", ha="center", va="center", fontsize=8, color=txt_color)

ax.set_xticks(range(len(methods)))
ax.set_xticklabels(["Network\nproximity", "Genetic\nassociation", "GNN link\nprediction", "Transcriptomic\nreversal"], fontsize=8)
ax.set_yticks(range(len(df)))
ax.set_yticklabels([f"{d} ({c})" for d, c in zip(df["drug"], df["cluster"])], fontsize=8)
ax.set_xlabel("Percentile rank in each method's scored universe (lower = stronger)", fontsize=8)
ax.set_title("Convergence of four independent scoring methods\non the 15 final candidates", fontsize=9)

cbar = fig.colorbar(im, ax=ax, shrink=0.7, pad=0.02)
cbar.set_label("Percentile rank (%)", fontsize=8)
cbar.ax.tick_params(labelsize=7)

legend_patch = mpatches.Patch(facecolor="#e0e0e0", edgecolor="white", hatch="////", label="Not in this method's scored universe")
ax.legend(handles=[legend_patch], loc="upper center", bbox_to_anchor=(0.5, -0.18), fontsize=7, frameon=False)

plt.tight_layout()
plt.savefig(OUT / "fig3_method_convergence.png", dpi=300, bbox_inches="tight")
print("Wrote", OUT / "fig3_method_convergence.png")

# print the underlying numbers for manual verification against the manuscript text
print(df.to_string(index=False))
