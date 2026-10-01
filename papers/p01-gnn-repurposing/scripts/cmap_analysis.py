"""
CMAP connectivity analysis: nintedanib and p01 candidates vs IgAN disease signature.

Method:
  1. Download GEO dataset GSE93798 (Liu et al. 2017, J Am Soc Nephrol; IgAN glomerular transcriptome)
  2. Identify IgAN vs. healthy samples; run differential expression (Mann-Whitney U)
  3. Extract up-regulated and down-regulated gene signatures
  4. Query Enrichr LINCS_L1000_Chem_Pert libraries for reversal (connectivity) scores
  5. Rank all perturbagens; report p01 candidates' positions

Logic:
  - CMAP reversal: a drug that DOWNregulates IgAN-UPregulated genes (and vice versa)
    is a candidate to reverse the disease transcriptome.
  - Connectivity score = geometric mean of rank percentiles from both Enrichr queries.
    Higher = stronger reversal signal.

Output:
  results/cmap_analysis/igan_signature.tsv       -- top DEGs used as query
  results/cmap_analysis/cmap_scores.tsv          -- all ranked drugs with scores
  results/cmap_analysis/p01_candidates_cmap.tsv  -- p01 candidates' CMAP ranks

Usage:
    uv run python scripts/cmap_analysis.py

Fallback:
    If GEO download fails, a curated IgAN gene signature from published meta-analyses
    is used automatically.
"""

import argparse
import json
import math
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "shared/scripts"))
from pipeline_utils import find_repo_root, load_config

import pandas as pd
import requests
from scipy import stats

ENRICHR_URL = "https://maayanlab.cloud/Enrichr"
DEFAULT_CONFIG = "shared/configs/igan.yaml"

# ── Curated fallback IgAN signature ─────────────────────────────────────────────
# Used only if the GEO download fails. From: Liu et al. 2017 (GSE93798), Berthier et al. 2012, Zhou et al. 2025,
#       IJBS 2025 TGF-beta review. Glomerular compartment, IgAN vs healthy.
CURATED_UP = [
    # NF-κB / inflammatory
    "NFKB1", "RELA", "REL", "NFKBIA", "TNF", "IL6", "IL1B", "CXCL10",
    "CCL2", "CCL5", "ICAM1", "VCAM1",
    # TGF-β / fibrosis
    "TGFB1", "TGFBR1", "TGFBR2", "SMAD2", "SMAD3", "FN1",
    "COL1A1", "COL1A2", "COL3A1", "COL4A1", "ACTA2", "VIM",
    # Complement
    "C3", "CFB", "CFH", "C5",
    # BAFF/APRIL axis
    "TNFSF13", "TNFSF13B", "TNFRSF13B",
    # Podocyte stress markers (upregulated in injured podocytes)
    "NPHS1", "CD2AP",  # paradoxically upregulated in early injury
    # Macrophage / monocyte
    "CD68", "MRC1", "LYN", "SYK",
    # Additional IgAN markers
    "EDIL3", "MALAT1", "GADD45B",
]
CURATED_DOWN = [
    # Podocyte structural integrity (lost in IgAN)
    "SYNPO", "NPHS2", "PODXL", "WT1", "PTPRO",
    # Glomerular basement membrane
    "LAMB2", "LAMA5", "COL4A3", "COL4A4",
    # Metabolic / mitochondrial
    "COQ2", "HMGCR",
    # Anti-inflammatory
    "PPARA", "PPARG", "PPARGC1A",
    # Tight junction / barrier
    "CLDN1", "OCLN", "TJP1",
]


# ── GEO download ─────────────────────────────────────────────────────────────────
def download_geo_signature(gse_id: str = "GSE93798", geo_dir: Path | None = None) -> tuple[list[str], list[str]] | None:
    """
    Download GEO dataset and compute IgAN vs. healthy differential expression.
    Returns (upregulated_genes, downregulated_genes) or None if download fails.
    """
    if geo_dir is None:
        geo_dir = Path(__file__).parent.parent / "data/raw/geo"
    try:
        import GEOparse
        print(f"Downloading {gse_id} from GEO (may take 1-3 minutes)...")
        gse = GEOparse.get_GEO(
            geo=gse_id,
            destdir=str(geo_dir),
            silent=True,
        )

        # Get expression matrix from first GPL platform
        gsms = gse.gsms
        platform = list(gse.gpls.keys())[0]
        gpl = gse.gpls[platform]

        # Build expression table
        expr_rows = {}
        sample_labels = {}
        for gsm_id, gsm in gsms.items():
            title = gsm.metadata.get("title", [""])[0].lower()
            char = " ".join(gsm.metadata.get("characteristics_ch1", []))
            is_igan = "igan" in title or "iga" in title or "igan" in char.lower()
            is_ctrl = "control" in title or "healthy" in title or "normal" in title
            if not (is_igan or is_ctrl):
                continue
            sample_labels[gsm_id] = "IgAN" if is_igan else "Control"
            if gsm.table is not None and not gsm.table.empty:
                tbl = gsm.table
                # Identify ID column and value column robustly
                id_col  = tbl.columns[0]
                val_col = next((c for c in tbl.columns[1:] if tbl[c].dtype in ["float64","object"] or c == "VALUE"), tbl.columns[1])
                try:
                    expr_rows[gsm_id] = tbl.set_index(id_col)[val_col].astype(float)
                except Exception:
                    pass

        if not expr_rows or not sample_labels:
            print("  Could not parse sample labels — using curated signature.")
            return None

        expr = pd.DataFrame(expr_rows)
        igan_cols = [c for c, l in sample_labels.items() if l == "IgAN" and c in expr.columns]
        ctrl_cols = [c for c, l in sample_labels.items() if l == "Control" and c in expr.columns]

        if len(igan_cols) < 3 or len(ctrl_cols) < 3:
            print(f"  Too few samples (IgAN={len(igan_cols)}, Ctrl={len(ctrl_cols)}) — using curated.")
            return None

        print(f"  IgAN samples: {len(igan_cols)}, Control samples: {len(ctrl_cols)}")

        # Map probe IDs to gene symbols via GPL annotation
        # Column name varies by platform: "Gene Symbol", "Symbol", "GENE_SYMBOL", etc.
        probe_to_gene: dict[str, str] = {}
        if gpl.table is not None:
            sym_col = next(
                (c for c in gpl.table.columns
                 if c.lower() in ("symbol", "gene symbol", "gene_symbol",
                                  "gene.symbol", "genesymbol")),
                None,
            )
            id_col = gpl.table.columns[0]  # probe ID column (e.g. "ID")
            if sym_col:
                for _, row in gpl.table.iterrows():
                    genes = str(row.get(sym_col, "")).strip()
                    if genes and genes not in ("nan", "---", ""):
                        # Some platforms use " /// " as multi-gene separator
                        gene = genes.split("///")[0].split(" // ")[0].strip()
                        if gene:
                            probe_to_gene[str(row[id_col])] = gene

        # Differential expression: Mann-Whitney U per probe
        results = []
        for probe in expr.index:
            igan_vals = expr.loc[probe, igan_cols].dropna().values
            ctrl_vals = expr.loc[probe, ctrl_cols].dropna().values
            if len(igan_vals) < 3 or len(ctrl_vals) < 3:
                continue
            try:
                _, pval = stats.mannwhitneyu(igan_vals, ctrl_vals, alternative="two-sided")
            except Exception:
                continue
            fc = float(igan_vals.mean()) - float(ctrl_vals.mean())
            gene = probe_to_gene.get(str(probe), "")
            if gene and len(gene) <= 20 and not gene.startswith("AFFX"):  # skip control probes
                results.append({"probe": str(probe), "gene": gene, "log2fc": fc, "pval": float(pval)})

        if not results:
            print("  No DEA results with gene symbols — using curated signature.")
            return None

        de = pd.DataFrame(results)
        de = de.sort_values("pval").drop_duplicates("gene")

        # Top 200 up and down (z-score < 0.05)
        sig = de[de["pval"] < 0.05]
        up_genes   = sig[sig["log2fc"] > 0].head(200)["gene"].tolist()
        down_genes = sig[sig["log2fc"] < 0].head(200)["gene"].tolist()

        print(f"  Significant DEGs: {len(up_genes)} up, {len(down_genes)} down")
        return up_genes, down_genes

    except Exception as e:
        print(f"  GEO download failed ({e}) — using curated signature.")
        return None


# ── Enrichr CMAP query ───────────────────────────────────────────────────────────
def enrichr_query(gene_list: list[str], lib: str) -> dict[str, dict]:
    """
    Query Enrichr with a gene list against a LINCS L1000 library.
    Returns {term: {pval, rank, combined_score}} dict.
    """
    # Step 1: submit gene list (multipart/form-data required by Enrichr)
    r = requests.post(
        f"{ENRICHR_URL}/addList",
        files={
            "list":        (None, "\n".join(gene_list)),
            "description": (None, f"IgAN_query_{lib}"),
        },
        timeout=30,
    )
    r.raise_for_status()
    user_list_id = r.json()["userListId"]
    time.sleep(0.5)

    # Step 2: get enrichment results
    r2 = requests.get(
        f"{ENRICHR_URL}/enrich",
        params={"userListId": user_list_id, "backgroundType": lib},
        timeout=60,
    )
    r2.raise_for_status()
    data = r2.json().get(lib, [])

    results = {}
    for row in data:
        # row: [rank, term, pval, zscore, combined_score, genes, adj_pval, ...]
        rank      = row[0]
        term      = row[1]
        pval      = row[2]
        combined  = row[4]
        results[term] = {"rank": rank, "pval": pval, "combined_score": combined}
    return results


def parse_drug_name(term: str) -> str:
    """
    Extract drug name from Enrichr LINCS L1000 term.
    Format: "BATCH CELL TIME-DRUG_NAME-DOSE"
    e.g. "LJP006 SKBR3 24H-nintedanib-1.0" → "nintedanib"
         "LJP008 HA1E 24H-GSK-461364-3.33"  → "gsk-461364"
    """
    parts = term.split("-")
    if len(parts) < 2:
        return term.strip().lower()
    # parts[0] = "BATCH CELL TIME"
    # parts[1:] = drug segments + dose
    drug_parts = parts[1:]
    # Drop the last segment if it looks like a concentration (float)
    if drug_parts:
        try:
            float(drug_parts[-1])
            drug_parts = drug_parts[:-1]
        except ValueError:
            pass
    return "-".join(drug_parts).strip().lower()


def compute_connectivity(
    up_genes: list[str],
    down_genes: list[str],
    n_drugs: int = 500,
) -> pd.DataFrame:
    """
    Query both LINCS libraries and compute combined connectivity score per drug.

    Connectivity logic:
      - IgAN-UP genes enriched in a drug's DOWN-signature → drug reverses up-regulation
      - IgAN-DOWN genes enriched in a drug's UP-signature → drug reverses down-regulation
    Combined score = geometric mean of combined_scores from both queries.
    """
    print("\nQuerying Enrichr LINCS_L1000_Chem_Pert_down with IgAN-UP genes...")
    down_results = enrichr_query(up_genes[:150], "LINCS_L1000_Chem_Pert_down")
    print(f"  Got {len(down_results)} terms")
    time.sleep(1)

    print("Querying Enrichr LINCS_L1000_Chem_Pert_up with IgAN-DOWN genes...")
    up_results = enrichr_query(down_genes[:150], "LINCS_L1000_Chem_Pert_up")
    print(f"  Got {len(up_results)} terms")

    # Aggregate by drug name: take max combined_score per drug across cell lines
    drug_down: dict[str, float] = {}  # drug → best score when downregulating IgAN-UP genes
    for term, vals in down_results.items():
        drug = parse_drug_name(term)
        drug_down[drug] = max(drug_down.get(drug, 0), vals["combined_score"])

    drug_up: dict[str, float] = {}    # drug → best score when upregulating IgAN-DOWN genes
    for term, vals in up_results.items():
        drug = parse_drug_name(term)
        drug_up[drug] = max(drug_up.get(drug, 0), vals["combined_score"])

    # Combine: only drugs appearing in BOTH queries get a combined score
    all_drugs = set(drug_down) | set(drug_up)
    rows = []
    for drug in all_drugs:
        s_down = drug_down.get(drug, 0.0)
        s_up   = drug_up.get(drug, 0.0)
        if s_down > 0 and s_up > 0:
            combined = math.sqrt(s_down * s_up)   # geometric mean
        else:
            combined = max(s_down, s_up) * 0.5    # penalty for single direction
        rows.append({
            "drug": drug,
            "score_down_query": round(s_down, 2),
            "score_up_query":   round(s_up, 2),
            "connectivity_score": round(combined, 2),
        })

    df = pd.DataFrame(rows).sort_values("connectivity_score", ascending=False).reset_index(drop=True)
    df["rank"] = df.index + 1
    return df


# ── Main ─────────────────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default=DEFAULT_CONFIG,
                        help="Path to disease config YAML (relative to repo root or absolute)")
    args = parser.parse_args()

    cfg, repo_root = load_config(args.config)
    rep = cfg["repurposing"]
    disease_id   = cfg["disease"]["id"]
    disease_name = cfg["disease"]["name"]

    output_dir = repo_root / rep["output_dir"]
    out_dir = output_dir / "results/cmap_analysis"
    out_dir.mkdir(parents=True, exist_ok=True)
    geo_dir = output_dir / "data/raw/geo"

    # Read candidates from final_candidates.tsv when available; fall back to config list
    final_cands_path = output_dir / "results/final_candidates.tsv"
    if final_cands_path.exists():
        final_df = pd.read_csv(final_cands_path, sep="\t")
        candidates = final_df["drug"].str.lower().str.strip().tolist()
        print(f"Loaded {len(candidates)} candidates from {final_cands_path}")
    else:
        candidates = [c.lower() for c in rep.get("cmap_candidates", [])]
        if not candidates:
            print("WARNING: final_candidates.tsv not found and no cmap_candidates in config.",
                  file=sys.stderr)
            print("Run scripts 07 and 08 first, or add 'cmap_candidates' to your config.",
                  file=sys.stderr)

    known_drugs = [d.lower() for d in rep.get("known_approved_drugs", [])]

    print(f"=== CMAP Connectivity Analysis: {disease_name} vs. pipeline candidates ===\n")

    # Step 1: Get disease gene signature
    geo_result = download_geo_signature("GSE93798", geo_dir=geo_dir)
    if geo_result:
        up_genes, down_genes = geo_result
        sig_source = "GEO:GSE93798 (Liu et al. 2017)"
    else:
        print(f"Using curated {disease_name} gene signature from published meta-analyses.")
        up_genes, down_genes = CURATED_UP, CURATED_DOWN
        sig_source = "Curated (Liu 2017, Berthier 2012, Zhou 2025, IJBS 2025)"

    print(f"\nSignature source: {sig_source}")
    print(f"Up-regulated genes ({len(up_genes)}): {', '.join(up_genes[:10])}{'...' if len(up_genes)>10 else ''}")
    print(f"Down-regulated genes ({len(down_genes)}): {', '.join(down_genes[:10])}{'...' if len(down_genes)>10 else ''}")

    # Save signature
    sig_df = pd.DataFrame({
        "gene":      up_genes + down_genes,
        "direction": ["up"] * len(up_genes) + ["down"] * len(down_genes),
    })
    sig_df.to_csv(out_dir / f"{disease_id}_signature.tsv", sep="\t", index=False)
    print(f"\nSignature saved to {out_dir / f'{disease_id}_signature.tsv'}")

    # Step 2: Query CMAP
    scores = compute_connectivity(up_genes, down_genes)
    scores.to_csv(out_dir / "cmap_scores.tsv", sep="\t", index=False)
    print(f"\nFull CMAP results saved ({len(scores)} drugs): {out_dir / 'cmap_scores.tsv'}")

    # Step 3: Extract pipeline candidates
    n_total = len(scores)
    p01_rows = []
    for cand in candidates:
        # Exact match first, then prefix match (handles salt forms like "bosutinib hcl")
        exact = scores[scores["drug"] == cand.lower()]
        match = exact if not exact.empty else \
                scores[scores["drug"].str.startswith(cand.lower(), na=False)]
        if not match.empty:
            r = match.iloc[0]
            pct = round(r["rank"] / n_total * 100, 1)
            top10 = r["rank"] <= n_total * 0.1
            p01_rows.append({
                "candidate": cand,
                "cmap_rank": int(r["rank"]),
                "total_drugs": n_total,
                "percentile": pct,
                "connectivity_score": r["connectivity_score"],
                "top_10pct": top10,
            })
        else:
            p01_rows.append({
                "candidate": cand,
                "cmap_rank": None,
                "total_drugs": n_total,
                "percentile": None,
                "connectivity_score": 0.0,
                "top_10pct": False,
            })

    candidates_tsv = out_dir / f"{disease_id}_candidates_cmap.tsv"
    p01_df = pd.DataFrame(p01_rows).sort_values("cmap_rank")
    p01_df.to_csv(candidates_tsv, sep="\t", index=False)

    # Step 4: Print summary
    print(f"\n{'='*65}")
    print(f"CMAP CONNECTIVITY RESULTS — {disease_name} candidates")
    print(f"Signature: {sig_source}")
    print(f"Total drugs scored: {n_total}")
    print(f"{'='*65}")
    print(f"\n{'Candidate':<22} {'Rank':>6} {'Percentile':>11} {'Score':>8} {'Top 10%':>8}")
    print("-" * 60)
    for _, row in p01_df.iterrows():
        rank_str  = str(int(row["cmap_rank"])) if pd.notna(row["cmap_rank"]) else "NOT FOUND"
        pct_str   = f"{row['percentile']:.1f}%" if row["percentile"] else "—"
        score_str = f"{row['connectivity_score']:.1f}"
        flag      = "✓" if row["top_10pct"] else ""
        print(f"{row['candidate']:<22} {rank_str:>6} {pct_str:>11} {score_str:>8} {flag:>8}")

    in_top10 = p01_df["top_10pct"].sum()
    found    = p01_df["cmap_rank"].notna().sum()
    print(f"\n{in_top10}/{found} {disease_id.upper()} candidates in top 10% of CMAP connectivity scores")

    # Top 20 overall for context
    print(f"\n{'='*65}")
    print(f"TOP 20 OVERALL CONNECTIVITY (CMAP reversers of {disease_name} signature)")
    print(f"{'='*65}")
    print(scores[["rank", "drug", "connectivity_score", "score_down_query",
                  "score_up_query"]].head(20).to_string(index=False))

    print(f"\n\nAll outputs saved to: {out_dir}/")
    print(f"  {disease_id}_signature.tsv  — {len(sig_df)} genes used as query")
    print(f"  cmap_scores.tsv        — full ranked drug list")
    print(f"  {disease_id}_candidates_cmap.tsv — candidates' CMAP ranks")


if __name__ == "__main__":
    main()
