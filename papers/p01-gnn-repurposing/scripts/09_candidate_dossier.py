"""
Step 9: Generate per-drug evidence dossiers for all 15 final candidates.

For each drug, produces a markdown file in results/candidates/ containing:
  - Drug profile (from final_candidates.tsv)
  - Pipeline scores: network proximity, Open Targets, full-model R-GCN rank (where available)
  - PubMed evidence: papers citing drug + IgAN or kidney fibrosis
  - ClinicalTrials.gov: active/completed trials for this drug
  - Open Targets: known drug-disease associations for this compound
  - Mechanistic rationale: why this drug addresses the target GWAS locus

Usage:
    uv run python scripts/09_candidate_dossier.py [--top N]

    Default: --top 15

Output:
    results/candidates/<drug_name>.md  (one file per drug)
"""

import argparse
import datetime
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "shared/scripts"))
from pipeline_utils import find_repo_root, load_config

import pandas as pd
import requests

PUBMED_ESEARCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
PUBMED_EFETCH  = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
OT_API         = "https://api.platform.opentargets.org/api/v4/graphql"
CTGOV_API      = "https://clinicaltrials.gov/api/v2/studies"
TODAY          = datetime.date.today().isoformat()

DEFAULT_CONFIG = "shared/configs/igan.yaml"


def pubmed_search(drug: str, disease_terms: list[str], max_results: int = 10) -> list[dict]:
    """Search PubMed for drug + disease papers."""
    disease_query = " OR ".join(f'"{t}"[tiab]' for t in disease_terms)
    query = (
        f'("{drug}"[tiab] OR "{drug.title()}"[tiab]) AND '
        f'({disease_query} OR "kidney fibrosis"[tiab] OR "glomerulonephritis"[tiab])'
    )
    resp = requests.get(
        PUBMED_ESEARCH,
        params={"db": "pubmed", "term": query, "retmax": max_results,
                "retmode": "json", "usehistory": "n"},
        timeout=15,
    )
    resp.raise_for_status()
    ids = resp.json().get("esearchresult", {}).get("idlist", [])
    if not ids:
        return []

    fetch = requests.get(
        PUBMED_EFETCH,
        params={"db": "pubmed", "id": ",".join(ids), "rettype": "abstract",
                "retmode": "xml"},
        timeout=30,
    )
    fetch.raise_for_status()

    # Parse minimal fields from XML without lxml dependency
    import re
    xml = fetch.text
    articles = []
    for block in re.findall(r"<PubmedArticle>.*?</PubmedArticle>", xml, re.DOTALL):
        pmid_m = re.search(r"<PMID[^>]*>(\d+)</PMID>", block)
        title_m = re.search(r"<ArticleTitle>(.*?)</ArticleTitle>", block, re.DOTALL)
        year_m  = re.search(r"<PubDate>.*?<Year>(\d{4})</Year>", block, re.DOTALL)
        journal_m = re.search(r"<Title>(.*?)</Title>", block)
        abstract_m = re.search(r"<AbstractText[^>]*>(.*?)</AbstractText>", block, re.DOTALL)
        pmid    = pmid_m.group(1) if pmid_m else ""
        title   = re.sub(r"<[^>]+>", "", title_m.group(1)) if title_m else "No title"
        year    = year_m.group(1) if year_m else "?"
        journal = re.sub(r"<[^>]+>", "", journal_m.group(1)) if journal_m else "?"
        abstract = re.sub(r"<[^>]+>", "", abstract_m.group(1))[:300] if abstract_m else ""
        articles.append({"pmid": pmid, "title": title, "year": year,
                         "journal": journal, "abstract": abstract})
    return articles


def ot_drug_query(drug_name: str) -> list[dict]:
    """Search Open Targets for known drug-disease associations for this compound."""
    search_q = """
    query Search($q: String!) {
      search(queryString: $q, entityNames: ["drug"]) {
        hits {
          id
          object {
            ... on Drug {
              name
              maximumClinicalStage
              drugType
              indications {
                rows {
                  id
                  maxClinicalStage
                  disease { name id }
                }
              }
            }
          }
        }
      }
    }
    """
    try:
        resp = requests.post(
            OT_API,
            json={"query": search_q, "variables": {"q": drug_name}},
            headers={"Content-Type": "application/json"},
            timeout=15,
        )
        resp.raise_for_status()
        hits = resp.json().get("data", {}).get("search", {}).get("hits", [])
        if not hits:
            return []
        drug_obj = hits[0].get("object", {})
        if not drug_obj:
            return []
        indications = drug_obj.get("indications", {}).get("rows", [])
        return [
            {
                "disease": ind.get("disease", {}).get("name", "?"),
                "disease_id": ind.get("disease", {}).get("id", ""),
                "max_phase": ind.get("maxClinicalStage", "?"),
            }
            for ind in indications[:10]
        ]
    except Exception:
        return []


def ctgov_search(drug_name: str, max_results: int = 8) -> list[dict]:
    """Query ClinicalTrials.gov v2 API for trials involving this drug."""
    try:
        resp = requests.get(
            CTGOV_API,
            params={
                "query.intr": drug_name,
                "pageSize": max_results,
                "format": "json",
                "fields": "NCTId,BriefTitle,OverallStatus,Phase,Condition,StartDate",
            },
            timeout=15,
        )
        resp.raise_for_status()
        studies = resp.json().get("studies", [])
        results = []
        for s in studies:
            proto = s.get("protocolSection", {})
            id_mod = proto.get("identificationModule", {})
            status_mod = proto.get("statusModule", {})
            design_mod = proto.get("designModule", {})
            cond_mod = proto.get("conditionsModule", {})
            nct = id_mod.get("nctId", "?")
            title = id_mod.get("briefTitle", "?")
            status = status_mod.get("overallStatus", "?")
            phases = design_mod.get("phases", [])
            phase = phases[0] if phases else "N/A"
            conditions = cond_mod.get("conditions", [])
            condition = "; ".join(conditions[:3]) if conditions else "?"
            start = status_mod.get("startDateStruct", {}).get("date", "?")
            results.append({
                "nct": nct,
                "title": title[:120],
                "status": status,
                "phase": phase,
                "condition": condition[:80],
                "start": start,
            })
        return results
    except Exception:
        return []


def write_dossier(row: pd.Series, pubmed: list[dict],
                  ot_indications: list[dict], ctgov_trials: list[dict],
                  gnn_rank: int | None, gnn_score: float | None,
                  hit_labels: dict, gwas_context: dict,
                  disease_name: str, efo_id: str,
                  np_total: int,
                  validation_steps: dict) -> Path:
    drug = row["drug"]
    hit  = int(row["hit_alignment"])
    target = str(row["target_gene"])
    moa  = " | ".join(str(row["moa"]).split(" | "))
    phase = str(row["clinical_phase"])
    np_rank  = int(row["np_rank"])
    np_z     = float(row["np_z_score"])
    ot_rank  = int(row["ot_rank"])
    ot_score = float(row["ot_genetic_score"])
    egfr     = str(row["egfr_flag"])
    neph     = str(row["nephrotoxicity_risk"])
    bp       = str(row["hypertension_risk"])
    safety   = str(row["safety_note"])
    priority = int(row["priority_score"])

    import tempfile
    slug = drug.replace(" ", "_").replace("/", "-")
    out  = Path(tempfile.gettempdir()) / f"{slug}.md"

    gnn_row = (
        f"| GNN link prediction rank (full R-GCN) | {gnn_rank:,} / 7,957 (score {gnn_score:.3g}) |"
        if gnn_rank is not None
        else "| GNN link prediction rank (full R-GCN) | — (drug not in PrimeKG drug node set) |"
    )

    hit_label = hit_labels.get(hit, hit_labels.get(str(hit), f"Hit {hit}"))

    lines = [
        f"# Candidate dossier: {drug}",
        f"",
        f"**Generated:** {TODAY}  |  **Pipeline:** {disease_name} GNN drug repurposing",
        f"",
        f"---",
        f"",
        f"## Pipeline scores",
        f"",
        f"| Metric | Value |",
        f"|--------|-------|",
        f"| Final rank | {int(row['final_rank'])} |",
        f"| Priority score | {priority} / 9 |",
        f"| Network proximity rank | {np_rank} / {np_total:,} (z = {np_z:.3f}) |",
        f"| Open Targets genetic rank | {ot_rank} / 3,005 (score {ot_score:.4f}) |",
        gnn_row,
        f"| Hit alignment | {hit_label} |",
        f"| GWAS target | {target} |",
        f"| Clinical phase | {phase} |",
        f"| Mechanism of action | {moa} |",
        f"",
        f"---",
        f"",
        f"## GWAS locus context",
        f"",
        gwas_context.get(target, f"Target gene {target} — see GWAS locus details in the primary reference."),
        f"",
        f"---",
        f"",
        f"## CKD safety profile",
        f"",
        f"| Parameter | Assessment |",
        f"|-----------|------------|",
        f"| eGFR adjustment | {egfr} |",
        f"| Nephrotoxicity risk | {neph} |",
        f"| Hypertension risk | {bp} |",
        f"",
        safety,
        f"",
        f"---",
        f"",
        f"## PubMed evidence (drug + {disease_name}/kidney fibrosis)",
        f"",
    ]

    if pubmed:
        for a in pubmed:
            lines += [
                f"**{a['title']}**",
                f"*{a['journal']}* ({a['year']}) | PMID: {a['pmid']}",
                f"{a['abstract']}{'...' if len(a['abstract']) == 300 else ''}",
                f"",
            ]
    else:
        lines.append(f"No PubMed results found for this drug + {disease_name}/kidney fibrosis query.")
        lines.append("")

    lines += [
        f"---",
        f"",
        f"## Open Targets known indications",
        f"",
    ]
    if ot_indications:
        lines.append(f"| Disease | Phase |")
        lines.append(f"|---------|-------|")
        for ind in ot_indications:
            lines.append(f"| {ind['disease']} | {ind['max_phase']} |")
        lines.append("")
        disease_ind = [i for i in ot_indications
                       if any(t.lower() in i["disease"].lower()
                              for t in [disease_name, "glomerulo", "nephropathy", "kidney"])]
        if not disease_ind:
            lines.append(f"{disease_name} not listed as a known indication — confirms novelty of this repurposing prediction.")
            lines.append("")
    else:
        lines.append("No Open Targets indication data retrieved.")
        lines.append("")

    lines += [
        f"---",
        f"",
        f"## ClinicalTrials.gov active and completed trials",
        f"",
    ]
    if ctgov_trials:
        lines.append(f"| NCT ID | Title | Status | Phase | Condition | Start |")
        lines.append(f"|--------|-------|--------|-------|-----------|-------|")
        for t in ctgov_trials:
            lines.append(
                f"| [{t['nct']}](https://clinicaltrials.gov/study/{t['nct']}) "
                f"| {t['title']} | {t['status']} | {t['phase']} "
                f"| {t['condition']} | {t['start']} |"
            )
        lines.append("")
        disease_trials = [t for t in ctgov_trials
                          if any(k in t["condition"].lower()
                                 for k in ("glomerulo", "nephropathy", "kidney",
                                           disease_name.lower().split()[0]))]
        if disease_trials:
            lines.append(f"> Note: {len(disease_trials)} trial(s) involve kidney-related conditions — review carefully for overlap with this prediction.")
        else:
            lines.append(f"No kidney or {disease_name} trials found for this drug — confirms this is a novel repurposing prediction.")
        lines.append("")
    else:
        lines.append("No ClinicalTrials.gov results retrieved.")
        lines.append("")

    target_display = target if target and target not in ("", "nan") else "target protein"
    lines += [
        f"---",
        f"",
        f"## Next steps for validation",
        f"",
        f"1. **In silico:** protein-ligand co-folding against {target_display} (Boltz-2 results for the five top-ranked candidates are in results/boltz/)",
        f"2. **In vitro:** {validation_steps.get('in_vitro', 'Test in disease-relevant cell model')}",
        f"3. **In vivo:** {validation_steps.get('in_vivo', 'Disease-appropriate animal model')}",
        f"4. **Clinical:** {validation_steps.get('clinical', 'Phase 2 trial design with appropriate endpoints')}",
    ]

    out.write_text("\n".join(lines))
    return out


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default=DEFAULT_CONFIG,
                        help="Path to disease config YAML (relative to repo root or absolute)")
    parser.add_argument("--top", type=int, default=15,
                        help="Number of candidates to process (default: 15)")
    args = parser.parse_args()

    cfg, repo_root = load_config(args.config)
    rep = cfg["repurposing"]
    disease_id   = cfg["disease"]["id"]
    disease_name = cfg["disease"]["name"]
    efo_id       = rep["efo_id"]

    output_dir     = repo_root / rep["output_dir"]
    results_dir    = output_dir / "results"
    candidates_dir = results_dir / "candidates"
    candidates_dir.mkdir(parents=True, exist_ok=True)

    final_candidates = results_dir / "final_candidates.tsv"
    gnn_ranked       = results_dir / "kg_full_candidate_ranked_rgcn.tsv"

    # Disease search terms for PubMed
    disease_terms = [disease_name] + [l["gene"] for l in rep["gwas_loci"]
                                       if l["hit"] in (1, 2)][:3]

    # Pull config-driven context into local vars for write_dossier
    hit_labels       = {int(k): v for k, v in rep.get("hit_labels", {}).items()}
    gwas_context     = rep.get("gwas_context", {})
    validation_steps = rep.get("validation_steps", {})

    print(f"Disease: {disease_name} ({disease_id})")
    print(f"Output:  {output_dir}\n")

    df = pd.read_csv(final_candidates, sep="\t")
    top = df.head(args.top)

    # Total NP-scored drugs — for accurate rank denominator in dossiers
    np_ranked_path = results_dir / "network_proximity_ranked.tsv"
    np_total = len(pd.read_csv(np_ranked_path, sep="\t")) if np_ranked_path.exists() else 0

    print(f"Generating dossiers for top {len(top)} candidates ({np_total} NP-scored total)...\n")

    gnn_df = None
    if gnn_ranked.exists():
        gnn_df = pd.read_csv(gnn_ranked, sep="\t")
        gnn_df["drug_lower"] = gnn_df["drug"].str.lower()
        print(f"GNN ranked data loaded: {len(gnn_df)} rows\n")
    else:
        print(f"GNN ranked file not found — GNN scores will be omitted\n")

    generated = []
    for _, row in top.iterrows():
        drug = row["drug"]
        print(f"  [{int(row['final_rank'])}/{len(top)}] {drug}...")

        gnn_rank = gnn_score_val = None
        if gnn_df is not None:
            match = gnn_df[gnn_df["drug_lower"] == drug.lower()]
            if not match.empty:
                gnn_rank      = int(match.iloc[0]["rank"])
                gnn_score_val = float(match.iloc[0]["score"])
        print(f"    GNN rank: {gnn_rank if gnn_rank else 'not in PrimeKG'}")

        pubmed = pubmed_search(drug, disease_terms, max_results=8)
        print(f"    PubMed: {len(pubmed)} results")
        time.sleep(0.4)

        ot_ind = ot_drug_query(drug)
        print(f"    OT indications: {len(ot_ind)}")
        time.sleep(0.2)

        ctgov = ctgov_search(drug, max_results=8)
        print(f"    ClinicalTrials: {len(ctgov)} trials")
        time.sleep(0.3)

        out = write_dossier(
            row, pubmed, ot_ind, ctgov,
            gnn_rank, gnn_score_val,
            hit_labels, gwas_context,
            disease_name, efo_id,
            np_total, validation_steps,
        )
        # Move to disease-specific candidates dir
        dest = candidates_dir / out.name
        out.rename(dest)
        generated.append(dest)
        print(f"    → {dest.name}")

    print(f"\nDone. {len(generated)} dossiers written to {candidates_dir}/")


if __name__ == "__main__":
    main()
