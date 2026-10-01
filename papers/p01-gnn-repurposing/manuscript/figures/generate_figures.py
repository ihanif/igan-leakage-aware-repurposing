"""Generate Figure 1 (GWAS coverage map) and Figure 2 (pipeline flowchart) for DRAFT_v1."""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
from pathlib import Path

OUT = Path(__file__).parent

# ── colour palette (colourblind-safe) ────────────────────────────────────────
C_APPROVED  = "#2166ac"   # blue   – FDA-approved drug exists
C_CLINICAL  = "#4dac26"   # green  – Phase 2/3 drug, no approval
C_NOVEL     = "#d7191c"   # red    – GWAS-validated, no IgAN drug
C_LOW       = "#999999"   # grey   – lower OT genetic rank

# ── Figure 1: GWAS coverage map ──────────────────────────────────────────────
loci = [
    # (gene, hit, status, ot_rank, approved_drug)
    ("TNFSF13",   2, "approved",  1,   "sibeprenlimab"),
    ("CFH",       3, "clinical",  2,   ""),
    ("LYN",       2, "novel",     22,  ""),
    ("REL",       4, "novel",     23,  ""),
    ("TGFB1",     4, "novel",     34,  ""),
    ("TNFRSF13B", 2, "approved",  20,  "atacicept"),
    ("TNFSF13B",  2, "clinical",  134, "telitacicept"),
    ("CFB",       3, "approved",  133, "iptacopan"),
    ("EDNRA",     4, "approved",  132, "atrasentan"),
    ("AGTR1",     4, "approved",  131, "sparsentan"),
    ("FCAR",      2, "novel",     53,  ""),
    ("IRF8",      2, "novel",     63,  ""),
    ("C5",        3, "clinical",  139, "ravulizumab"),
    ("MASP2",     3, "clinical",  135, "narsoplimab"),
    ("CFD",       3, "clinical",  133, ""),
    ("CD28",      2, "clinical",  None,"abatacept"),
    ("TNFSF4",    2, "clinical",  None,""),
    ("VAV3",      2, "low",       None,""),
    ("CCR6",      2, "novel",     None,""),
    ("HLA-DQB1",  2, "approved",  None,"budesonide"),
    ("GALNT2",    1, "novel",     1013,""),
    ("C1GALT1",   1, "novel",     None,""),
    ("COSMC",     1, "novel",     None,""),
    ("ITGA1",     4, "low",       None,""),
    ("HORMAD2",   0, "low",       None,""),
    ("DEFA",      3, "low",       None,""),
    ("ST6GALNAC2",1, "low",       None,""),
    ("RELA",      4, "clinical",  None,""),
    ("CFHR1",     3, "novel",     None,""),
    ("CFHR5",     3, "novel",     None,""),
]

hit_labels = {0: "Unknown", 1: "Hit 1\nGlyco-\nsylation", 2: "Hit 2\nAutoAb/\nB-cell",
              3: "Hit 3\nComple-\nment", 4: "Hit 4\nFibrosis/\nNF-κB"}

# sort by hit then status priority
status_order = {"approved": 0, "clinical": 1, "novel": 2, "low": 3}
loci.sort(key=lambda x: (x[1], status_order[x[2]]))

genes    = [l[0] for l in loci]
hits     = [l[1] for l in loci]
statuses = [l[2] for l in loci]
colors   = {"approved": C_APPROVED, "clinical": C_CLINICAL,
            "novel": C_NOVEL, "low": C_LOW}
bar_colors = [colors[s] for s in statuses]

fig, ax = plt.subplots(figsize=(12, 6))
y = np.arange(len(genes))
bars = ax.barh(y, [1]*len(genes), color=bar_colors, edgecolor="white", height=0.7)

ax.set_yticks(y)
ax.set_yticklabels(genes, fontsize=9, fontfamily="sans-serif")
ax.set_xticks([])
ax.set_xlim(0, 1)
ax.invert_yaxis()

# hit alignment dividers
hit_changes = []
prev = None
for i, h in enumerate(hits):
    if h != prev:
        hit_changes.append((i, h))
        prev = h

for idx, h in hit_changes:
    if idx > 0:
        ax.axhline(idx - 0.5, color="#cccccc", linewidth=0.8, linestyle="--")
    lbl = hit_labels.get(h, "")
    ax.text(1.02, idx + (hits.count(h) - 1) / 2, lbl,
            va="center", ha="left", fontsize=7.5, color="#555555",
            fontfamily="sans-serif")

# legend
patches = [
    mpatches.Patch(color=C_APPROVED, label="FDA-approved for IgAN"),
    mpatches.Patch(color=C_CLINICAL, label="Phase 2/3 clinical candidate"),
    mpatches.Patch(color=C_NOVEL,    label="GWAS-validated — no IgAN drug (pipeline target)"),
    mpatches.Patch(color=C_LOW,      label="Locus present — lower OT evidence"),
]
ax.legend(handles=patches, loc="upper center", fontsize=8, framealpha=0.95,
          edgecolor="#cccccc", bbox_to_anchor=(0.5, -0.04), ncol=2)

ax.set_title("GWAS coverage map: therapeutic status of 30 IgAN susceptibility loci\n"
             "(Kiryluk et al. 2023, Nature Genetics)",
             fontsize=10, fontfamily="sans-serif", pad=10)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.spines["bottom"].set_visible(False)

plt.tight_layout()
plt.savefig(OUT / "fig1_gwas_coverage.png", dpi=300, bbox_inches="tight")
plt.close()
print("Figure 1 saved.")

# ── Figure 2: pipeline flowchart ─────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(10, 8))
ax.set_xlim(0, 10)
ax.set_ylim(0, 10)
ax.axis("off")

def box(ax, x, y, w, h, text, facecolor="#e8f4fd", edgecolor="#2166ac", fontsize=9):
    rect = mpatches.FancyBboxPatch((x - w/2, y - h/2), w, h,
                                    boxstyle="round,pad=0.1",
                                    facecolor=facecolor, edgecolor=edgecolor,
                                    linewidth=1.2)
    ax.add_patch(rect)
    ax.text(x, y, text, ha="center", va="center", fontsize=fontsize,
            fontfamily="sans-serif", wrap=True,
            multialignment="center")

def arrow(ax, x1, y1, x2, y2):
    ax.annotate("", xy=(x2, y2), xytext=(x1, y1),
                arrowprops=dict(arrowstyle="->", color="#555555", lw=1.2))

# Input data boxes
box(ax, 2.5, 9.2, 3.5, 0.8,  "Kiryluk 2023 GWAS\n30 IgAN loci",          "#fff3cd", "#856404")
box(ax, 7.5, 9.2, 3.5, 0.8,  "STRING v12 PPI\n(confidence ≥ 700)",        "#d4edda", "#155724")

# Disease module
box(ax, 5.0, 7.8, 5.0, 0.9,
    "IgAN disease module\n953 proteins · 1,388 interactions\n(Script 01)",
    "#d4edda", "#155724")
arrow(ax, 2.5, 8.8, 3.8, 8.25)
arrow(ax, 7.5, 8.8, 6.2, 8.25)

# Three analysis streams
box(ax, 1.5, 6.2, 2.6, 0.85,
    "Network proximity\nGuney 2016 z-scores\n1,989 drugs  (Script 04)",
    "#e8f4fd", "#2166ac", fontsize=8)
box(ax, 5.0, 6.2, 2.6, 0.85,
    "Open Targets\nGenetic assoc. score\n3,236 drug-gene pairs\n(Script 05)",
    "#e8f4fd", "#2166ac", fontsize=8)
box(ax, 8.5, 6.2, 2.6, 0.85,
    "PrimeKG R-GCN\nlink prediction, leave-one-out\n7,957 drugs  (Scripts 06, 06b)",
    "#e8f4fd", "#2166ac", fontsize=8)

arrow(ax, 5.0, 7.35, 1.5, 6.63)
arrow(ax, 5.0, 7.35, 5.0, 6.63)
arrow(ax, 5.0, 7.35, 8.5, 6.63)

# Consensus intersection
box(ax, 5.0, 4.7, 4.0, 0.85,
    "Consensus shortlist\nNP ∩ OT top-100\n17 novel candidates  (Script 07)",
    "#fde8e8", "#9e1c1c")

# Arrows from method boxes to consensus — all three methods
arrow(ax, 1.5, 5.78, 3.2, 5.13)   # NP → consensus (left edge)
arrow(ax, 5.0, 5.78, 5.0, 5.13)   # OT → consensus (centre)
# GNN scores candidates but does not select them (shortlist is NP ∩ OT), so no arrow

# Retrospective validation callout — positioned right of centre, clear of all arrows
ax.text(8.5, 4.7, "Annotation-safe recall (NP):\natrasentan rank 8 · fostamatinib rank 188\nsparsentan excluded (annotation backfill)",
        ha="center", va="center", fontsize=7, color="#444444",
        fontfamily="sans-serif", style="italic",
        bbox=dict(boxstyle="round,pad=0.25", facecolor="#fffef0",
                  edgecolor="#b8860b", linewidth=0.8, alpha=0.9))

# Clinical filter
box(ax, 5.0, 3.35, 4.0, 0.85,
    "Clinical feasibility filter\neGFR · nephrotoxicity · BP risk\n15 final candidates  (Script 08)",
    "#fde8e8", "#9e1c1c")
arrow(ax, 5.0, 4.28, 5.0, 3.78)

# Evidence dossiers
box(ax, 5.0, 2.0, 4.0, 0.85,
    "Candidate dossiers\nPubMed · literature cross-ref · OT\nTop-10 candidates  (Script 09)",
    "#f3e8fd", "#4a0e6b")
arrow(ax, 5.0, 2.93, 5.0, 2.43)

# Output highlight
box(ax, 5.0, 0.75, 6.5, 0.85,
    "15 novel candidates  ·  6 launched drugs  ·  3 mechanism clusters\n"
    "Nintedanib (rank 1)  ·  F351 (rank 7)  ·  Bafetinib (rank 8)",
    "#fff3cd", "#856404", fontsize=8)
arrow(ax, 5.0, 1.57, 5.0, 1.18)

ax.set_title("Computational pipeline for IgAN drug repurposing",
             fontsize=10, fontfamily="sans-serif", pad=6)

plt.tight_layout()
plt.savefig(OUT / "fig2_pipeline.png", dpi=300, bbox_inches="tight")
plt.close()
print("Figure 2 saved.")
