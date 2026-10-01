"""Generate the leave-one-out validation figure (R-GCN vs GraphSAGE ablation). Built entirely from results/kg_loo_recall_rgcn.tsv and
results/kg_loo_recall_graphsage.tsv, produced by scripts/06b_gnn_loo_validation.py.
No numbers are invented.
"""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).parent

rgcn = pd.read_csv(ROOT / "results" / "kg_loo_recall_rgcn.tsv", sep="\t")
sage = pd.read_csv(ROOT / "results" / "kg_loo_recall_graphsage.tsv", sep="\t")

merged = rgcn.merge(sage, on="drug", suffixes=("_rgcn", "_sage"))
merged = merged.sort_values("held_out_percentile_rgcn", ascending=False).reset_index(drop=True)

fig, ax = plt.subplots(figsize=(6.5, 5.0))

y = range(len(merged))
C_RGCN = "#08519c"
C_SAGE = "#d7191c"

for i, row in merged.iterrows():
    ax.plot([row["held_out_percentile_rgcn"], row["held_out_percentile_sage"]], [i, i],
            color="#bbbbbb", linewidth=1, zorder=1)

ax.scatter(merged["held_out_percentile_rgcn"], y, color=C_RGCN, s=45, zorder=2, label="R-GCN (relation-aware)")
ax.scatter(merged["held_out_percentile_sage"], y, color=C_SAGE, s=45, zorder=2, marker="^", label="GraphSAGE (relation-agnostic)")

ax.axvline(10.0, color="black", linestyle="--", linewidth=1, zorder=0)
ax.text(10.3, len(merged) - 0.3, "top decile", fontsize=7, va="top")

ax.set_yticks(list(y))
ax.set_yticklabels(merged["drug"], fontsize=8)
ax.set_xlabel("Leave-one-out held-out percentile rank (of 7,957 PrimeKG drug nodes)", fontsize=8)
ax.set_xscale("log")
ax.set_xlim(0.3, 100)
ax.set_title("Leave-one-out recall: R-GCN vs GraphSAGE architecture ablation", fontsize=9)
ax.legend(loc="upper left", fontsize=8, frameon=False, bbox_to_anchor=(0.02, 0.99))
ax.spines[["top", "right"]].set_visible(False)
ax.set_ylim(-0.8, len(merged) + 1.6)

n_rgcn_top10 = (merged["held_out_percentile_rgcn"] <= 10).sum()
n_sage_top10 = (merged["held_out_percentile_sage"] <= 10).sum()

plt.tight_layout()
plt.savefig(OUT / "fig4_loo_ablation.png", dpi=300, bbox_inches="tight")
print("Wrote", OUT / "fig4_loo_ablation.png")
print(f"R-GCN top-decile recall: {n_rgcn_top10}/{len(merged)}")
print(f"GraphSAGE top-decile recall: {n_sage_top10}/{len(merged)}")
print(merged[["drug", "held_out_rank_rgcn", "held_out_percentile_rgcn",
              "held_out_rank_sage", "held_out_percentile_sage"]].to_string(index=False))
