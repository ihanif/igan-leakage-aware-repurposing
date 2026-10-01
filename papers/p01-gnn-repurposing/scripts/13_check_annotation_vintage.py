"""
Step 13: Annotation-vintage check for the retrospective recall claim.

For each known IgAN drug, lists the Drug Repurposing Hub targets that fall inside
the disease module (the only targets that enter its network proximity score) and
whether each one is already annotated in the release dated 24 March 2020. A known
drug counts toward recall only if all of its in-module targets pass.

Usage:
    uv run python papers/p01-gnn-repurposing/scripts/13_check_annotation_vintage.py

Data requirements:
    data/raw/drug_repurposing_hub/repurposing_drugs.tsv           (file date 18 Aug 2025, used for scoring)
    data/raw/drug_repurposing_hub/repurposing_drugs_20200324.txt  (dated release, vintage check only)
        https://s3.amazonaws.com/data.clue.io/repurposing/downloads/repurposing_drugs_20200324.txt
    data/processed/igan_disease_module.json                       (from Script 01)
Checksums for both Drug Repurposing Hub files are in data/raw/drug_repurposing_hub/*.sha256.
"""

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "shared/scripts"))
from pipeline_utils import load_config

DEFAULT_CONFIG = "shared/configs/igan.yaml"


def read_targets(path: Path) -> dict[str, set[str]]:
    """Map lower-case drug name -> set of target genes from a Drug Repurposing Hub file."""
    targets = {}
    header = None
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith("!") or not line.strip():
            continue
        cols = line.split("\t")
        if header is None:
            header = cols
            continue
        row = dict(zip(header, cols))
        genes = {g.strip() for g in row.get("target", "").split("|") if g.strip()}
        targets[row["pert_iname"].strip().lower()] = genes
    return targets


def main():
    cfg, repo_root = load_config(DEFAULT_CONFIG)
    rep = cfg["repurposing"]
    base = repo_root / rep["output_dir"]
    drh = repo_root / rep["shared_data_dir"] / "drug_repurposing_hub"

    current = read_targets(drh / "repurposing_drugs.tsv")
    dated = read_targets(drh / "repurposing_drugs_20200324.txt")
    module = {n["id"] for n in json.loads((base / "data/processed/igan_disease_module.json").read_text())["nodes"]}

    print(f"{'drug':<24}{'in-module targets (2025 file)':<32}{'annotated Mar 2020':<22}annotation-safe")
    verdict = {}
    for drug in dict.fromkeys(d.lower() for d in rep["known_approved_drugs"]):
        if drug not in current:
            print(f"{drug:<24}{'not in Drug Repurposing Hub':<32}")
            continue
        scored = sorted(current[drug] & module)
        ok_2020 = [t for t in scored if t in dated.get(drug, set())]
        safe = bool(scored) and len(ok_2020) == len(scored)
        verdict[drug] = safe
        print(f"{drug:<24}{', '.join(scored) or '-':<32}{', '.join(ok_2020) or '-':<22}{'yes' if safe else 'no'}")

    # The three cases the manuscript's validation protocol rests on.
    assert verdict["atrasentan"] and verdict["fostamatinib"], "annotation-safe drugs changed"
    assert not verdict["sparsentan"], "sparsentan unexpectedly annotation-safe"
    print("\nManuscript claims hold: atrasentan and fostamatinib annotation-safe; sparsentan excluded.")


if __name__ == "__main__":
    main()
