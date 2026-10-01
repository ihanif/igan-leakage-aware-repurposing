"""
Step 8: Clinical feasibility filter.

Applies three successive filters to the consensus shortlist:
  1. Hard removal  — non-therapeutic compounds (contrast agents, inorganic tools)
  2. Safety flags  — eGFR-based dose requirements, nephrotoxicity, hypertension risk
  3. Priority score — phase, safety profile, mechanistic fit, config-driven bonuses

Safety information is drawn from FDA prescribing information summaries and
published CKD pharmacokinetics literature. All flags require author verification
against the primary source before publication.

The SAFETY lookup is a global pharmacological database keyed by drug name.
Drugs not found in SAFETY receive a default "No curated safety data" note.

Usage:
    # IgAN (default):
    uv run python scripts/08_clinical_filter.py

    # Any disease:
    uv run python scripts/08_clinical_filter.py --config shared/configs/fsgs.yaml

Input (from config repurposing.output_dir/results/):
    consensus_shortlist.tsv

Output:
    final_candidates.tsv

Validation checkpoint:
    Final list should contain >= 8 candidates.
"""

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "shared/scripts"))
from pipeline_utils import find_repo_root, load_config

import pandas as pd

DEFAULT_CONFIG = "shared/configs/igan.yaml"

# ── Hard removal ──────────────────────────────────────────────────────────────
# Compounds that are not therapeutic drug candidates regardless of network score.
HARD_REMOVE = {
    "iodipamide",                          # contrast agent; not a systemic therapeutic
    "barium-6-o-phosphonato-d-glucose",    # inorganic research tool; no human formulation
}

# ── Safety profile lookup ─────────────────────────────────────────────────────
# Manually curated from FDA labels and CKD pharmacokinetics literature.
# Each entry: (egfr_flag, nephrotox_risk, bp_risk, note)
#   egfr_flag:      none | monitor | dose_reduce_egfr<30 | dose_reduce_egfr<60 | ci_egfr<30
#   nephrotox_risk: none | low | moderate | high
#   bp_risk:        none | low | moderate | high (hypertension risk — critical in IgAN)
#   note:           brief clinical note for manuscript
SAFETY = {
    "nintedanib": (
        "monitor",
        "low",
        "none",
        "No formal dose adjustment for eGFR > 30 (hepatic metabolism). Used in IPF patients "
        "with comorbid CKD. Mild GI toxicity manageable. Best CKD safety profile in cluster.",
    ),
    "bosutinib": (
        "dose_reduce_egfr<30",
        "low",
        "low",
        "Reduce dose to 200 mg/day for eGFR < 30 (FDA label). Hepatic metabolism; renal "
        "excretion minimal. No CKD contraindication above eGFR 30.",
    ),
    "dasatinib": (
        "monitor",
        "low",
        "moderate",
        "No renal dose adjustment required (hepatic metabolism). Hypertension in 11% of "
        "patients; monitor BP closely in IgAN. Pleural effusion at oncology doses (100 mg); "
        "lower doses (20–50 mg) under investigation for autoimmune indications.",
    ),
    "masitinib": (
        "monitor",
        "low",
        "low",
        "No CKD-specific dose adjustment published. Tested in ALS, MS, pancreatic cancer. "
        "Mild edema in ~15% of patients. Renal safety data from non-oncology trials favourable.",
    ),
    "bafetinib": (
        "monitor",
        "low",
        "low",
        "Phase 2 safety data limited; hepatic metabolism; no renal contraindication identified. "
        "Most LYN-selective agent in shortlist. Author should verify Phase 2 safety reports.",
    ),
    "ponatinib": (
        "monitor",
        "moderate",
        "high",
        "Arterial occlusive events (14% of patients), hypertension (67%), hepatotoxicity. "
        "Significant cardiovascular burden inappropriate for IgAN without careful risk-benefit. "
        "Flag for specialist cardiovascular review before IgAN use.",
    ),
    "rebastinib": (
        "monitor",
        "low",
        "low",
        "Phase 1/2 only; limited safety data. No renal contraindication identified. "
        "Author should verify most recent trial safety reporting.",
    ),
    "tranilast": (
        "monitor",
        "low",
        "none",
        "Launched antiallergic/anti-fibrotic. Well-tolerated in long-term use (Japan). "
        "Rare hepatic effects. No CKD contraindication; used in Japanese patients with CKD.",
    ),
    "f351": (
        "monitor",
        "low",
        "none",
        "Phase 2 liver fibrosis trial. Anti-fibrotic mechanism. Safety profile emerging; "
        "no renal-specific contraindication identified. Most clinically advanced novel Hit 4 prediction.",
    ),
    "adaprev": (
        "monitor",
        "low",
        "none",
        "Phase 3 (DME indication). TGF-β receptor inhibitor. Systemic safety data limited "
        "for non-ocular route. Author should verify systemic dosing safety.",
    ),
    "d-4476": (
        "none",
        "none",
        "none",
        "Preclinical only. CK1α + TGF-β dual inhibitor. No human safety data. "
        "Include as mechanistic anchor; flag for in vitro IgAN validation.",
    ),
    "ft011": (
        "none",
        "none",
        "none",
        "Preclinical anti-fibrotic. Published preclinical CKD data (diabetic nephropathy model). "
        "No human safety data. Strong mechanistic rationale; needs Phase 1 development.",
    ),
    "on123300": (
        "none",
        "none",
        "none",
        "Preclinical CDK/LYN inhibitor. No human safety data. Include as mechanistic anchor.",
    ),
    "ginsenoside-re": (
        "none",
        "none",
        "none",
        "Phase 1 anti-inflammatory. Natural compound; low clinical tractability for "
        "IgAN repurposing. Include as mechanistic validation of REL/NF-κB axis.",
    ),
    "tolimidone": (
        "none",
        "none",
        "none",
        "Phase 2; development discontinued (antidiabetic). SRC activator — unusual mechanism "
        "for IgAN. Mechanistic rationale requires further scrutiny. Flag for review.",
    ),
}

# ── Priority score components ─────────────────────────────────────────────────
PHASE_SCORE = {
    "Launched": 4,
    "Phase 3": 3,
    "Phase 2": 2,
    "Phase 1/Phase 2": 2,
    "Phase 1": 1,
    "Preclinical": 0,
}
SAFETY_SCORE = {"none": 2, "low": 2, "moderate": 1, "high": 0}
BP_SCORE = {"none": 2, "low": 2, "moderate": 1, "high": 0}


def phase_score(phase: str) -> int:
    for k, v in PHASE_SCORE.items():
        if k.lower() in str(phase).lower():
            return v
    return 0


def priority(row: pd.Series, safety_row: tuple) -> float:
    _, nephrotox, bp_risk, _ = safety_row
    p = phase_score(row["clinical_phase"])
    s = SAFETY_SCORE.get(nephrotox, 1)
    b = BP_SCORE.get(bp_risk, 1)
    return float(p + s + b)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default=DEFAULT_CONFIG,
                        help="Path to disease config YAML (relative to repo root or absolute)")
    args = parser.parse_args()

    cfg, repo_root = load_config(args.config)
    rep = cfg["repurposing"]
    disease_id   = cfg["disease"]["id"]
    disease_name = cfg["disease"]["name"]

    output_dir  = repo_root / rep["output_dir"]
    results_dir = output_dir / "results"
    in_path     = results_dir / "consensus_shortlist.tsv"
    out_path    = results_dir / "final_candidates.tsv"

    priority_bonuses  = rep.get("priority_bonuses", {})
    mechanism_clusters = rep.get("mechanism_clusters", [])

    print(f"Disease: {disease_name} ({disease_id})")
    print(f"Output:  {output_dir}\n")

    df = pd.read_csv(in_path, sep="\t")
    print(f"Loaded {len(df)} candidates from consensus shortlist")

    # Step 1: hard removal
    before = len(df)
    df = df[~df["drug"].isin(HARD_REMOVE)].copy()
    print(f"Hard removal: {before - len(df)} non-therapeutic compounds removed ({len(df)} remain)")

    # Step 2: safety annotation
    rows = []
    for _, row in df.iterrows():
        drug = row["drug"]
        safety_row = SAFETY.get(
            drug,
            ("none", "none", "none", "No curated safety data — author review required."),
        )
        egfr_flag, nephrotox, bp_risk, note = safety_row
        rows.append({
            **row.to_dict(),
            "egfr_flag": egfr_flag,
            "nephrotoxicity_risk": nephrotox,
            "hypertension_risk": bp_risk,
            "safety_note": note,
            "priority_score": priority(row, safety_row),
        })

    annotated = pd.DataFrame(rows)

    # Step 2b: apply config-driven priority bonuses (e.g. dual-hit drugs)
    for drug_name, bonus in priority_bonuses.items():
        mask = annotated["drug"] == drug_name.lower()
        if mask.any():
            annotated.loc[mask, "priority_score"] += float(bonus)

    # Step 3: rank by priority score, then mean_rank
    annotated = annotated.sort_values(
        ["priority_score", "mean_rank"], ascending=[False, True]
    ).reset_index(drop=True)
    annotated["final_rank"] = range(1, len(annotated) + 1)

    annotated.to_csv(out_path, sep="\t", index=False)
    print(f"\nFinal candidates saved to {out_path} ({len(annotated)} drugs)")

    # ── Summary ───────────────────────────────────────────────────────────────
    print("\n=== Final candidate list ===")
    print(f"\n{'#':<3} {'Drug':<33} {'Phase':<18} {'Hit':<4} {'eGFR':<22} {'Neph':<8} {'BP':<8} {'Priority'}")
    print("-" * 110)
    for _, r in annotated.iterrows():
        print(
            f"{int(r['final_rank']):<3} {r['drug']:<33} {str(r['clinical_phase']):<18} "
            f"{int(r['hit_alignment']):<4} {r['egfr_flag']:<22} {r['nephrotoxicity_risk']:<8} "
            f"{r['hypertension_risk']:<8} {r['priority_score']:.0f}"
        )

    launched = annotated[annotated["clinical_phase"].str.lower().str.contains("launch", na=False)]
    high_bp  = annotated[annotated["hypertension_risk"] == "high"]

    print(f"\nSummary:")
    print(f"  Total final candidates: {len(annotated)}")
    print(f"  Launched drugs:         {len(launched)}")
    for cluster in mechanism_clusters:
        n = len(annotated[annotated["target_gene"] == cluster["target_gene"]])
        print(f"  {cluster['name']}: {n}")
    print(f"  High BP risk flag:      {len(high_bp)}"
          + (f" ({', '.join(high_bp['drug'].tolist())})" if len(high_bp) else ""))

    if len(annotated) < 8:
        print(f"\nWARNING: Only {len(annotated)} final candidates (target: >= 8).")


if __name__ == "__main__":
    main()
